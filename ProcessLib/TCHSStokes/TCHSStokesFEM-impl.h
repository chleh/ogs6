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

    auto const shape_matrices_quadratic =
        initShapeMatrices<ShapeFunctionVelocity, ShapeMatricesTypeVelocity,
                          IntegrationMethod, VelocityDim>(
            e, is_axially_symmetric, _integration_method);

    auto const shape_matrices_linear =
        initShapeMatrices<ShapeFunctionPressure, ShapeMatricesTypePressure,
                          IntegrationMethod, VelocityDim>(
            e, is_axially_symmetric, _integration_method);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        // velocity (subscript u)
        auto& ip_data = _ip_data[ip];
        auto const& sm_2 = shape_matrices_quadratic[ip];
        _ip_data[ip].integration_weight =
            _integration_method.getWeightedPoint(ip).getWeight() *
            sm_2.integralMeasure * sm_2.detJ;

        ip_data.H = ShapeMatricesTypeVelocity::template MatrixType<
            VelocityDim, Block::size(Block::V)>::Zero(VelocityDim,
                                                      Block::size(Block::V));
        for (int i = 0; i < VelocityDim; ++i)
            ip_data.H
                .template block<1, Block::size(Block::V) / VelocityDim>(
                    i, i * Block::size(Block::V) / VelocityDim)
                .noalias() = sm_2.N;

        ip_data.N_2 = sm_2.N;
        ip_data.dNdx_2 = sm_2.dNdx;

        ip_data.N_1 = shape_matrices_linear[ip].N;
        ip_data.dNdx_1 = shape_matrices_linear[ip].dNdx;
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
    assert(local_x.size() == Block::index(Block::V) + Block::size(Block::V));

    // make Eigen::Map<>s
    auto map_vector_segment = [](std::vector<double> const& v, auto b) {
        return Eigen::Map<typename ShapeMatricesTypeVelocity::
                              template VectorType<Block::size(b)> const>(
            v.data() + Block::index(b), Block::size(b));
    };
    auto const nodal_p = map_vector_segment(local_x, Block::P);
    auto const nodal_T = map_vector_segment(local_x, Block::T);
    auto const nodal_xmV = map_vector_segment(local_x, Block::X);
    auto const nodal_v = map_vector_segment(local_x, Block::V);

    auto local_K = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeVelocity::template MatrixType<
            Block::index(Block::V) + Block::size(Block::V),
            Block::index(Block::V) + Block::size(Block::V)>>(
        local_K_data, Block::index(Block::V) + Block::size(Block::V),
        Block::index(Block::V) + Block::size(Block::V));

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesTypeVelocity::template VectorType<
            Block::index(Block::V) + Block::size(Block::V)>>(
        local_rhs_data, Block::index(Block::V) + Block::size(Block::V));

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        // shape functions and weights /////////////////////////////////////////
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;

        auto const& H = _ip_data[ip].H;

        auto const& N_2 = _ip_data[ip].N_2;
        auto const& dNdx_2 = _ip_data[ip].dNdx_2;

        auto const& N_1 = _ip_data[ip].N_1;
        auto const& dNdx_p = _ip_data[ip].dNdx_1;

        auto const x_coord =
            interpolateXCoordinate<ShapeFunctionVelocity,
                                   ShapeMatricesTypeVelocity>(_element, N_2);
        auto const B =
            LinearBMatrix::computeBMatrix<VelocityDim,
                                          ShapeFunctionVelocity::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_2, N_2, x_coord, _is_axially_symmetric);

        // special vectors and tensors /////////////////////////////////////////
        auto const& I =
            MathLib::KelvinVector::Invariants<KelvinVectorSize>::identity2;
        auto const& P_dev = MathLib::KelvinVector::Invariants<
            KelvinVectorSize>::deviatoric_projection;

        // interpolate nodal values ////////////////////////////////////////////
        Eigen::Matrix<double, VelocityDim, 1> const v = H * nodal_v;
        double const p = N_1 * nodal_p;
        double const T = N_1 * nodal_T;
        double const xmV = N_1 * nodal_xmV;

        // material parameters /////////////////////////////////////////////////
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

        double const rho_GR = 0;  // TODO

        // assemble local matrices /////////////////////////////////////////////

        // K_p? (total mass balance) ///////////////////////////////////////////
        // K_pp
        Block::block(local_K, Block::P, Block::P).noalias() +=
            N_1.transpose() * v.transpose() * (rho_GR / p * w) * dNdx_p;

        // K_pv
        Block::block(local_K, Block::P, Block::V).noalias() +=
            N_1.transpose() * I.transpose() * B * w;

        // K_v? (gas momentum balance //////////////////////////////////////////
        // K_vp
        Block::block(local_K, Block::V, Block::P).noalias() -=
            H.transpose() * (porosity * w) * dNdx_p;

        // K_vv
        Block::block(local_K, Block::V, Block::V).noalias() -=
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
        Block::segment(local_rhs, Block::V).noalias() -=
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
            Block::size(Block::P)> const>(
            local_x.data() + Block::index(Block::P), Block::size(Block::P));

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

        auto const& N_1 = shape_matrices_p.N;

        std::size_t const global_index = _element.getNodeIndex(n);
        (*_process_data.mesh_prop_nodal_p)[global_index] = N_1 * p;
    }
}

}  // namespace TCHSStokes
}  // namespace ProcessLib
