/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   TCHSStokesFEM-impl.h
 *  Created on November 29, 2017, 2:03 PM
 */

#pragma once

#include "TCHSStokesFEM.h"

#include <boost/math/special_functions/pow.hpp>

#include "MathLib/KelvinVector.h"
#include "NumLib/Fem/CoordinatesMapping/NaturalNodeCoordinates.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"

namespace ProcessLib
{
namespace TCHSStokes
{
template <typename ShapeFunctionVelocity, typename ShapeFunctionPressure,
          typename IntegrationMethod, int VelocityDim>
TCHSStokesLocalAssembler<ShapeFunctionVelocity, ShapeFunctionPressure,
                         IntegrationMethod, VelocityDim>::
    TCHSStokesLocalAssembler(MeshLib::Element const& e,
                             std::size_t const /*local_matrix_size*/,
                             bool const is_axially_symmetric,
                             unsigned const integration_order,
                             TCHSStokesProcessData<VelocityDim>& process_data)
    : _process_data(process_data),
      _integration_method(integration_order),
      _element(e),
      _is_axially_symmetric(is_axially_symmetric)
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    _ip_data.resize(n_integration_points);

    auto const shape_matrices_v =
        initShapeMatrices<ShapeFunctionVelocity, ShapeMatricesTypeVelocity,
                          IntegrationMethod, VelocityDim>(
            e, is_axially_symmetric, _integration_method);

    auto const shape_matrices_p =
        initShapeMatrices<ShapeFunctionPressure, ShapeMatricesTypePressure,
                          IntegrationMethod, VelocityDim>(
            e, is_axially_symmetric, _integration_method);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        // velocity (subscript u)
        auto& ip_data = _ip_data[ip];
        auto const& sm_v = shape_matrices_v[ip];
        _ip_data[ip].integration_weight =
            _integration_method.getWeightedPoint(ip).getWeight() *
            sm_v.integralMeasure * sm_v.detJ;

        ip_data.H = ShapeMatricesTypeVelocity::template MatrixType<
            VelocityDim, velocity_size>::Zero(VelocityDim, velocity_size);
        for (int i = 0; i < VelocityDim; ++i)
            ip_data.H
                .template block<1, velocity_size / VelocityDim>(
                    i, i * velocity_size / VelocityDim)
                .noalias() = sm_v.N;

        ip_data.N_v = sm_v.N;
        ip_data.dNdx_v = sm_v.dNdx;

        ip_data.N_p = shape_matrices_p[ip].N;
        ip_data.dNdx_p = shape_matrices_p[ip].dNdx;
    }
}

template <typename ShapeFunctionVelocity, typename ShapeFunctionPressure,
          typename IntegrationMethod, int VelocityDim>
void TCHSStokesLocalAssembler<
    ShapeFunctionVelocity, ShapeFunctionPressure, IntegrationMethod,
    VelocityDim>::assemble(double const t, std::vector<double> const& local_x,
                           std::vector<double>& /*local_M_data*/,
                           std::vector<double>& local_K_data,
                           std::vector<double>& local_rhs_data)
{
    assert(local_x.size() == pressure_size + velocity_size);

    auto nodal_v =
        Eigen::Map<typename ShapeMatricesTypeVelocity::template VectorType<
            velocity_size> const>(local_x.data() + velocity_index,
                                  velocity_size);

    auto local_K = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeVelocity::template MatrixType<
            velocity_size + pressure_size, velocity_size + pressure_size>>(
        local_K_data, velocity_size + pressure_size,
        velocity_size + pressure_size);

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesTypeVelocity::template VectorType<velocity_size +
                                                                pressure_size>>(
        local_rhs_data, velocity_size + pressure_size);

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;

        auto const& H = _ip_data[ip].H;

        auto const& N_v = _ip_data[ip].N_v;
        auto const& dNdx_v = _ip_data[ip].dNdx_v;

        auto const& N_p = _ip_data[ip].N_p;

        auto const x_coord =
            interpolateXCoordinate<ShapeFunctionVelocity,
                                   ShapeMatricesTypeVelocity>(_element, N_v);
        auto const B =
            LinearBMatrix::computeBMatrix<VelocityDim,
                                          ShapeFunctionVelocity::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_v, N_v, x_coord, _is_axially_symmetric);
        auto const& I =
            MathLib::KelvinVector::Invariants<KelvinVectorSize>::identity2;

        auto const mat_id = _process_data.material_ids[_element.getID()];

        auto const mu = _process_data.fluid_viscosity(t, x_position)[0];
        double porosity;
        Eigen::Matrix<double, VelocityDim, 1> grad_porosity =
            Eigen::Matrix<double, VelocityDim, 1>::Zero();
        double mu_eff;
        double f_1;
        double f_2;

        if (_element.getID() == 0 && ip == 0)
        {
            auto const rho_GR = _process_data.fluid_density(t, x_position)[0];
            auto const Re0 = _process_data.average_darcy_velocity *
                             _process_data.pellet_diameter * rho_GR / mu;

            INFO("Reynolds number is %g.", Re0);
            INFO("Weight is %g.", w);
        }

        if (mat_id == TCHSStokesProcessData<VelocityDim>::MATID_VOID)
        {
            porosity = 1.0;
            mu_eff = mu;
            f_1 = 0.0;
            f_2 = 0.0;
        }
        else if (mat_id == TCHSStokesProcessData<VelocityDim>::MATID_BED)
        {
            auto const r_bed = _process_data.bed_radius;
            auto const d_pel = _process_data.pellet_diameter;
            auto const poro_inf = _process_data.homogeneous_porosity;
            auto const exp_term = std::exp(-5.0 * (r_bed - x_coord) / d_pel);
            porosity = poro_inf + poro_inf * 1.36 * exp_term;

            grad_porosity[0] = poro_inf * 1.36 * 5.0 / d_pel * exp_term;

            auto const poro3 = boost::math::pow<3>(porosity);
            auto const rho_GR = _process_data.fluid_density(t, x_position)[0];
            f_1 = 150.0 * boost::math::pow<2>(1.0 - porosity) / poro3 * mu /
                  d_pel / d_pel;
            f_2 = 1.75 * (1.0 - porosity) / poro3 * rho_GR / d_pel;

            mu_eff = (*_process_data.effective_fluid_viscosity)(t, mu, rho_GR);
        }
        else
            OGS_FATAL("wrong material id: %d", mat_id);

        // p, v_Darcy no div ///////////////////////////////////////////////////

        // K_pv (from mass balance)
        local_K
            .template block<pressure_size, velocity_size>(pressure_index,
                                                          velocity_index)
            .noalias() += N_p.transpose() * I.transpose() * B * w;

        auto const& dNdx_p = _ip_data[ip].dNdx_p;
        // K_vp
        local_K
            .template block<velocity_size, pressure_size>(velocity_index,
                                                          pressure_index)
            .noalias() -= H.transpose() * (porosity * w) * dNdx_p;

        Eigen::Matrix<double, VelocityDim, 1> v = H * nodal_v;
        auto const& P_dev = MathLib::KelvinVector::Invariants<
            KelvinVectorSize>::deviatoric_projection;

        // K_vv
        local_K
            .template block<velocity_size, velocity_size>(velocity_index,
                                                          velocity_index)
            .noalias() -=
            B.transpose() * (2 * mu_eff * w) * P_dev * B +
            H.transpose() * (porosity * (f_1 + f_2 * v.norm()) * w) * H;

        // rhs_v
        auto const two_sym_v_grad_phi = [&]() {
            Eigen::Matrix3d m = Eigen::Matrix3d::Zero();
            m.topLeftCorner<VelocityDim, VelocityDim>().noalias() =
                v * grad_porosity.transpose() / porosity;
            m += m.transpose();
            return MathLib::KelvinVector::tensorToKelvin<VelocityDim>(m);
        };
        local_rhs.template segment<velocity_size>(velocity_index).noalias() -=
            mu_eff * B.transpose() * P_dev * two_sym_v_grad_phi() * w;
    }
}

template <typename ShapeFunctionVelocity, typename ShapeFunctionPressure,
          typename IntegrationMethod, int VelocityDim>
void TCHSStokesLocalAssembler<ShapeFunctionVelocity, ShapeFunctionPressure,
                              IntegrationMethod, VelocityDim>::
    postNonLinearSolverConcrete(std::vector<double> const&
                                /*local_x*/,
                                double const /*t*/,
                                bool const /*use_monolithic_scheme*/)
{
}

template <typename ShapeFunctionVelocity, typename ShapeFunctionPressure,
          typename IntegrationMethod, int VelocityDim>
void TCHSStokesLocalAssembler<
    ShapeFunctionVelocity, ShapeFunctionPressure, IntegrationMethod,
    VelocityDim>::preOutputConcrete(std::vector<double> const& local_x,
                                    const double /*t*/)
{
    auto const p =
        Eigen::Map<typename ShapeMatricesTypeVelocity::template VectorType<
            pressure_size> const>(local_x.data() + pressure_index,
                                  pressure_size);

    using FemType = NumLib::TemplateIsoparametric<ShapeFunctionPressure,
                                                  ShapeMatricesTypePressure>;

    FemType fe(*static_cast<const typename ShapeFunctionPressure::MeshElement*>(
        &_element));
    int const number_base_nodes = _element.getNumberOfBaseNodes();
    int const number_all_nodes = _element.getNumberOfNodes();

    for (int n = 0; n < number_base_nodes; ++n)
    {
        std::size_t const global_index = _element.getNodeIndex(n);
        (*_process_data.mesh_prop_nodal_p)[global_index] = p[n];
    }

    for (int n = number_base_nodes; n < number_all_nodes; ++n)
    {
        // Evaluated at higher order nodes' coordinates.
        typename ShapeMatricesTypePressure::ShapeMatrices shape_matrices_p{
            ShapeFunctionPressure::DIM, VelocityDim,
            ShapeFunctionPressure::NPOINTS};

        fe.computeShapeFunctions(
            NumLib::NaturalCoordinates<
                typename ShapeFunctionVelocity::MeshElement>::coordinates[n]
                .data(),
            shape_matrices_p, VelocityDim, _is_axially_symmetric);

        auto const& N_p = shape_matrices_p.N;

        std::size_t const global_index = _element.getNodeIndex(n);
        (*_process_data.mesh_prop_nodal_p)[global_index] = N_p * p;
    }
}

}  // namespace TCHSStokes
}  // namespace ProcessLib
