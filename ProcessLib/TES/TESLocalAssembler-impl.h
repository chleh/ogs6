/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "MaterialLib/Adsorption/Adsorption.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "TESDielectricByKraus.h"
#include "TESLocalAssembler.h"

namespace
{
enum class MatOutType
{
    OGS5,
    PYTHON
};

const MatOutType MATRIX_OUTPUT_FORMAT = MatOutType::PYTHON;

// TODO move to some location in the OGS core.
template <typename Mat>
void ogs5OutMat(const Mat& mat)
{
    for (unsigned r = 0; r < mat.rows(); ++r)
    {
        switch (MATRIX_OUTPUT_FORMAT)
        {
            case MatOutType::OGS5:
                if (r != 0)
                    std::printf("\n");
                std::printf("|");
                break;
            case MatOutType::PYTHON:
                if (r != 0)
                    std::printf(",\n");
                std::printf("[");
                break;
        }

        for (unsigned c = 0; c < mat.cols(); ++c)
        {
            switch (MATRIX_OUTPUT_FORMAT)
            {
                case MatOutType::OGS5:
                    std::printf(" %.16e", mat(r, c));
                    break;
                case MatOutType::PYTHON:
                    if (c != 0)
                        std::printf(",");
                    std::printf(" %23.16g", mat(r, c));
                    break;
            }
        }

        switch (MATRIX_OUTPUT_FORMAT)
        {
            case MatOutType::OGS5:
                std::printf(" | ");
                break;
            case MatOutType::PYTHON:
                std::printf(" ]");
                break;
        }
    }
    std::printf("\n");
}

template <typename Vec>
void ogs5OutVec(const Vec& vec)
{
    for (unsigned r = 0; r < vec.size(); ++r)
    {
        switch (MATRIX_OUTPUT_FORMAT)
        {
            case MatOutType::OGS5:
                if (r != 0)
                    std::printf("\n");
                std::printf("| %.16e | ", vec[r]);
                break;
            case MatOutType::PYTHON:
                if (r != 0)
                    std::printf(",\n");
                std::printf("[ %23.16g ]", vec[r]);
                break;
        }
    }
    std::printf("\n");
}

}  // anonymous namespace

namespace ProcessLib
{
namespace TES
{
template <typename ShapeFunction_, typename IntegrationMethod_,
          unsigned GlobalDim>
TESLocalAssembler<ShapeFunction_, IntegrationMethod_, GlobalDim>::
    TESLocalAssembler(MeshLib::Element const& e,
                      std::size_t const /*local_matrix_size*/,
                      bool is_axially_symmetric,
                      unsigned const integration_order,
                      AssemblyParams const& asm_params)
    : _element(e),
      _integration_method(integration_order),
      _shape_matrices(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                        IntegrationMethod_, GlobalDim>(
          e, is_axially_symmetric, _integration_method)),
      _d(asm_params, e.getID(), _integration_method.getNumberOfPoints(),
         GlobalDim)
{
}

template <typename ShapeFunction_, typename IntegrationMethod_,
          unsigned GlobalDim>
void TESLocalAssembler<ShapeFunction_, IntegrationMethod_, GlobalDim>::assemble(
    double const t, std::vector<double> const& local_x,
    std::vector<double>& local_M_data, std::vector<double>& local_K_data,
    std::vector<double>& local_b_data)
{
    auto const local_matrix_size = local_x.size();
    // This assertion is valid only if all nodal d.o.f. use the same shape matrices.
    assert(local_matrix_size == ShapeFunction::NPOINTS * NODAL_DOF);

    auto local_M = MathLib::createZeroedMatrix<NodalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::createZeroedMatrix<NodalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);
    auto local_b = MathLib::createZeroedVector<NodalVectorType>(local_b_data,
                                                            local_matrix_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& sm = _shape_matrices[ip];
        auto const& wp = _integration_method.getWeightedPoint(ip);
        auto const weight = wp.getWeight();

        _d.assembleIntegrationPoint(ip, t, local_x, sm, weight, local_M,
                                    local_K, local_b);
    }

    if (_d.getAssemblyParameters().output_element_matrices)
    {
        std::puts("### Element: ?");
        std::printf("\n---Mass matrix: \n");
        ogs5OutMat(local_M);
        std::printf("\n");

        std::printf("---Laplacian + Advective + Content matrix: \n");
        ogs5OutMat(local_K);
        std::printf("\n");

        std::printf("---RHS: \n");
        ogs5OutVec(local_b);
        std::printf("\n");
    }
}

template <typename ShapeFunction_, typename IntegrationMethod_,
          unsigned GlobalDim>
std::vector<double> const&
TESLocalAssembler<ShapeFunction_, IntegrationMethod_, GlobalDim>::
    getIntPtSolidDensity(const double /*t*/,
                         GlobalVector const& /*current_solution*/,
                         NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
                         std::vector<double>& cache) const
{
    auto const& d = _d.getData();
    cache.clear();
    cache.reserve(d.reactive_solid_state.size());
    for (auto& s : d.reactive_solid_state)
    {
        cache.emplace_back(d.ap.reactive_solid->getSolidDensity(*s));
    }
    return cache;
}

template <typename ShapeFunction_, typename IntegrationMethod_,
          unsigned GlobalDim>
std::vector<double> const&
TESLocalAssembler<ShapeFunction_, IntegrationMethod_, GlobalDim>::
    getIntPtReactionRate(const double /*t*/,
                         GlobalVector const& /*current_solution*/,
                         NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
                         std::vector<double>& cache) const
{
    auto& d = _d.getData();
    cache.clear();
    cache.reserve(d.reaction_rate.size());

    for (auto& qR : d.reaction_rate)
    {
        cache.emplace_back(d.ap.reactive_solid->getOverallRate(*qR));
    }
    return cache;
}

template <typename ShapeFunction_, typename IntegrationMethod_,
          unsigned GlobalDim>
std::vector<double> const&
TESLocalAssembler<ShapeFunction_, IntegrationMethod_, GlobalDim>::
    getIntPtDarcyVelocity(const double t,
                          GlobalVector const& current_solution,
                          NumLib::LocalToGlobalIndexMap const& dof_table,
                          std::vector<double>& cache) const
{
    auto const n_integration_points = _integration_method.getNumberOfPoints();

    auto const indices = NumLib::getIndices(_element.getID(), dof_table);
    assert(!indices.empty());
    auto const local_x = current_solution.get(indices);
    // local_x is ordered by component, local_x_mat is row major
    // TODO fixed size matrix ShapeFunction_::NPOINTS columns. Conflicts with
    // RowMajor for (presumably) 0D element.
    auto const local_x_mat = MathLib::toMatrix<
        Eigen::Matrix<double, NODAL_DOF, Eigen::Dynamic, Eigen::RowMajor>>(
        local_x, NODAL_DOF, ShapeFunction_::NPOINTS);

    cache.clear();
    auto cache_mat = MathLib::createZeroedMatrix<
        Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, GlobalDim, n_integration_points);

    auto const& perm = *_d.getAssemblyParameters().permeability;

    SpatialPosition pos;
    pos.setElementID(_d.getData().element_id);

    for (unsigned i = 0; i < n_integration_points; ++i)
    {
        double p, T, x;
        NumLib::shapeFunctionInterpolate(local_x, _shape_matrices[i].N, p, T,
                                         x);
        const double eta_GR = fluid_viscosity(p, T, x);

        pos.setIntegrationPoint(i);
        auto const k = perm(t, pos).front();

        cache_mat.col(i).noalias() =
            k * (_shape_matrices[i].dNdx *
                 local_x_mat.row(COMPONENT_ID_PRESSURE).transpose()) /
            -eta_GR;
    }

    return cache;
}

template <typename ShapeFunction_, typename IntegrationMethod_,
          unsigned GlobalDim>
std::vector<double> const& TESLocalAssembler<
    ShapeFunction_, IntegrationMethod_,
    GlobalDim>::getIntPtMassFlux(const double t,
                                 GlobalVector const& current_solution,
                                 NumLib::LocalToGlobalIndexMap const& dof_table,
                                 std::vector<double>& cache) const
{
    // auto const num_nodes = ShapeFunction_::NPOINTS;
    auto const num_intpts = _shape_matrices.size();

    auto const indices = NumLib::getIndices(_d.getData().element_id, dof_table);
    assert(!indices.empty());
    auto const local_x = current_solution.get(indices);
    // local_x is ordered by component, local_x_mat is row major
    // TODO fixed size matrix ShapeFunction_::NPOINTS columns. Conflicts with
    // RowMajor for (presumably) 0D element.
    auto const local_x_mat = MathLib::toMatrix<
        Eigen::Matrix<double, NODAL_DOF, Eigen::Dynamic, Eigen::RowMajor>>(
        local_x, NODAL_DOF, ShapeFunction_::NPOINTS);

    cache.clear();
    auto cache_mat = MathLib::createZeroedMatrix<
        Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, GlobalDim, num_intpts);

    auto const& perm = *_d.getAssemblyParameters().permeability;

    SpatialPosition pos;
    pos.setElementID(_d.getData().element_id);

    for (unsigned i = 0; i < num_intpts; ++i)
    {
        double p, T, x;
        NumLib::shapeFunctionInterpolate(local_x, _shape_matrices[i].N, p, T,
                                         x);
        const double eta_GR = fluid_viscosity(p, T, x);
        const double rho_GR = fluid_density(p, T, x);

        pos.setIntegrationPoint(i);
        auto const k = perm(t, pos).front();

        cache_mat.col(i).noalias() =
            k *
            (_shape_matrices[i].dNdx *
             local_x_mat.row(COMPONENT_ID_PRESSURE).transpose()) *
            (-rho_GR / eta_GR);
    }

    return cache;
}

template <typename ShapeFunction_, typename IntegrationMethod_,
          unsigned GlobalDim>
std::vector<double> const&
TESLocalAssembler<ShapeFunction_, IntegrationMethod_, GlobalDim>::
    getIntPtConductiveHeatFlux(
        const double t,
        GlobalVector const& current_solution,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& cache) const
{
    // auto const dim = _shape_matrices[0].dNdx.rows();
    // auto const num_nodes = ShapeFunction_::NPOINTS;
    auto const num_intpts = _shape_matrices.size();

    auto const indices = NumLib::getIndices(
        _d.getData().element_id, *_d.getAssemblyParameters().dof_table);
    auto const local_x = current_solution.get(indices);
    // local_x is ordered by component, local_x_mat is row major
    auto const local_x_mat =
        MathLib::toMatrix(local_x, NODAL_DOF, ShapeFunction_::NPOINTS);

    cache.clear();
    auto cache_mat = MathLib::createZeroedMatrix<
        Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, GlobalDim, num_intpts);

    SpatialPosition pos;
    pos.setElementID(_d.getData().element_id);

    const double lambda_S = _d.getAssemblyParameters().solid_heat_cond;

    for (unsigned i = 0; i < _shape_matrices.size(); ++i)
    {
        double p, T, x;
        NumLib::shapeFunctionInterpolate(local_x, _shape_matrices[i].N, p, T,
                                         x);

        pos.setIntegrationPoint(i);
        auto const poro = (*_d.getAssemblyParameters().poro)(t, pos)[0];

        const double lambda_F = fluid_heat_conductivity(p, T, x);

        const double lambda = poro * lambda_F + (1.0 - poro) * lambda_S;

        cache_mat.col(i).noalias() =
            -lambda * _shape_matrices[i].dNdx *
            local_x_mat.row(COMPONENT_ID_TEMPERATURE).transpose();
        /*
        cache_vec.row(i).noalias() = -lambda *
                                     local_x_mat.row(COMPONENT_ID_TEMPERATURE) *
                                     _shape_matrices[i].dNdx.transpose();
                                     */
    }

    return cache;
}

template <typename ShapeFunction_, typename IntegrationMethod_,
          unsigned GlobalDim>
const TESLocalAssemblerData&
TESLocalAssembler<ShapeFunction_, IntegrationMethod_,
                  GlobalDim>::getLocalAssemblerData() const
{
    return _d.getData();
}

template <typename ShapeFunction_, typename IntegrationMethod_,
          unsigned GlobalDim>
void TESLocalAssembler<ShapeFunction_, IntegrationMethod_, GlobalDim>::
    initializeSolidDensity(MeshLib::MeshItemType item_type,
                           std::vector<double> const& values)
{
    switch (item_type)
    {
        case MeshLib::MeshItemType::Cell:
        {
            assert(values.size() == 1);

            for (auto& s : _d.getData().reactive_solid_state)
            {
                assert(s->conversion().size() == 1);
                s->conversion()[0] = values.front();
            }
            break;
        }
        case MeshLib::MeshItemType::Node:
        {
            assert((int)values.size() == _shape_matrices[0].N.rows());

            auto& state = _d.getData().reactive_solid_state;

            for (std::size_t i = 0; i < _shape_matrices.size(); ++i)
            {
                auto& s = state[i];
                assert(s->conversion().size() == 1);
                NumLib::shapeFunctionInterpolate(values, _shape_matrices[i].N,
                                                 s->conversion()[0]);
            }
            break;
        }
        default:
            ERR("Unhandled mesh item type for initialization of secondary "
                "variable.");
            std::abort();
    }
}

template <typename ShapeFunction_, typename IntegrationMethod_,
          unsigned GlobalDim>
void TESLocalAssembler<ShapeFunction_, IntegrationMethod_, GlobalDim>::
    preTimestep(std::size_t const /*mesh_item_id*/,
                NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
                GlobalVector const& /*x*/, double const /*t*/,
                double const /*delta_t*/)
{
    _d.preTimestep();
}

}  // namespace TES
}  // namespace ProcessLib
