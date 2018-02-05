/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   IncompressibleStokesBrinkmanFEM-impl.h
 *  Created on November 29, 2017, 2:03 PM
 */

#pragma once

#include "IncompressibleStokesBrinkmanFEM.h"

#include "MaterialLib/SolidModels/KelvinVector.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"

namespace ProcessLib
{
namespace IncompressibleStokesBrinkman
{
template <typename ShapeFunctionVelocity, typename ShapeFunctionPressure,
          typename IntegrationMethod, int VelocityDim>
IncompressibleStokesBrinkmanLocalAssembler<ShapeFunctionVelocity,
                                           ShapeFunctionPressure,
                                           IntegrationMethod, VelocityDim>::
    IncompressibleStokesBrinkmanLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        IncompressibleStokesBrinkmanProcessData<VelocityDim>& process_data)
    : _process_data(process_data),
      _integration_method(integration_order),
      _element(e),
      _is_axially_symmetric(is_axially_symmetric)
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    _ip_data.resize(n_integration_points);

    auto const shape_matrices_u =
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
        auto const& sm_u = shape_matrices_u[ip];
        _ip_data[ip].integration_weight =
            _integration_method.getWeightedPoint(ip).getWeight() *
            sm_u.integralMeasure * sm_u.detJ;

        ip_data.H = ShapeMatricesTypeVelocity::template MatrixType<
            VelocityDim, velocity_size>::Zero(VelocityDim, velocity_size);
        for (int i = 0; i < VelocityDim; ++i)
            ip_data.H
                .template block<1, velocity_size / VelocityDim>(
                    i, i * velocity_size / VelocityDim)
                .noalias() = sm_u.N;

        ip_data.N_v = sm_u.N;
        ip_data.dNdx_v = sm_u.dNdx;

        ip_data.N_p = shape_matrices_p[ip].N;
        ip_data.dNdx_p = shape_matrices_p[ip].dNdx;
    }
}

template <typename ShapeFunctionVelocity, typename ShapeFunctionPressure,
          typename IntegrationMethod, int VelocityDim>
void IncompressibleStokesBrinkmanLocalAssembler<
    ShapeFunctionVelocity, ShapeFunctionPressure, IntegrationMethod,
    VelocityDim>::assemble(double const t, std::vector<double> const& local_x,
                           std::vector<double>& /*local_M_data*/,
                           std::vector<double>& local_K_data,
                           std::vector<double>& local_rhs_data)
{
    assert(local_x.size() == pressure_size + velocity_size);

    auto p = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
        pressure_size> const>(local_x.data() + pressure_index, pressure_size);

    auto v = Eigen::Map<typename ShapeMatricesTypeVelocity::template VectorType<
        velocity_size> const>(local_x.data() + velocity_index, velocity_size);

    auto local_K = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeVelocity::template MatrixType<
            velocity_size + pressure_size, velocity_size + pressure_size>>(
        local_K_data, velocity_size + pressure_size,
        velocity_size + pressure_size);

#if 0
    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesTypeVelocity::template VectorType<velocity_size +
                                                                pressure_size>>(
        local_rhs_data, velocity_size + pressure_size);
#endif

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;

        auto const& H = _ip_data[ip].H;

        auto const& N_u = _ip_data[ip].N_v;
        auto const& dNdx_u = _ip_data[ip].dNdx_v;

        auto const& N_p = _ip_data[ip].N_p;
        // auto const& dNdx_p = _ip_data[ip].dNdx_p;

        auto const x_coord =
            interpolateXCoordinate<ShapeFunctionVelocity,
                                   ShapeMatricesTypeVelocity>(_element, N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<VelocityDim,
                                          ShapeFunctionVelocity::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, _is_axially_symmetric);
        auto const& I =
            MaterialLib::SolidModels::Invariants<KelvinVectorSize>::identity2;

        auto const porosity = _process_data.porosity(t, x_position)[0];
        auto const mu_eff = _process_data.mu_eff(t, x_position)[0];
        auto const lambda_eff = _process_data.lambda_eff(t, x_position)[0];
        auto const f_1 = _process_data.f_1(t, x_position)[0];
        auto const f_2 = _process_data.f_2(t, x_position)[0];

#if 0
        // K_pp
        local_K
            .template block<velocity_size, pressure_size>(pressure_index,
                                                          pressure_index)
            .noalias() += 0;
#endif

        // K_pv
        local_K
            .template block<pressure_size, velocity_size>(pressure_index,
                                                          velocity_index)
            .noalias() += N_p.transpose() * I.transpose() * B * w;

        // K_vp
        local_K
            .template block<velocity_size, pressure_size>(velocity_index,
                                                          pressure_index)
            .noalias() += B.transpose() * I * porosity * N_p * w;

        // K_vv
        local_K
            .template block<velocity_size, velocity_size>(velocity_index,
                                                          velocity_index)
            .noalias() +=
            B.transpose() *
                (2 * mu_eff * B +
                 (lambda_eff - 2.0 * mu_eff / 3.0) * I * I.transpose() * B) *
                w +
            H.transpose() * (porosity * (-f_1 - f_2 * v.norm())) * H * w;
    }
}

template <typename ShapeFunctionVelocity, typename ShapeFunctionPressure,
          typename IntegrationMethod, int VelocityDim>
void IncompressibleStokesBrinkmanLocalAssembler<
    ShapeFunctionVelocity, ShapeFunctionPressure, IntegrationMethod,
    VelocityDim>::postNonLinearSolverConcrete(std::vector<double> const&
                                                  local_x,
                                              double const t,
                                              bool const use_monolithic_scheme)
{
}

}  // namespace IncompressibleStokesBrinkman
}  // namespace ProcessLib
