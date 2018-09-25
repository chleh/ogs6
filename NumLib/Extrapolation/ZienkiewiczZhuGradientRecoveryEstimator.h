/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <map>

#include "BaseLib/DUNEConfig.h"
#include "BaseLib/Error.h"

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"

namespace NumLib
{
static std::size_t ZienkiewiczZhuGradientRecoveryEstimator_DebugOutputCounter =
    0;

/*!
 * Gradient recovery according to the method of Zienkiewicz and Zhu, 1987.
 *
 * Zienkiewicz, O.C., Zhu, J.Z., 1987. A simple error estimator and adaptive
 * procedure for practical engineering analysis. International Journal for
 * Numerical Methods in Engineering 24, 337–357. doi:10.1002/nme.1620240206
 */
class ZienkiewiczZhuGradientRecoveryEstimator
{
public:
    //! Computes the local error estimate and stores the result internally.
    template <typename Grid, typename LocalAssemblers, typename IntPtDataGetter>
    void estimate(Grid const& grid,
                  LocalAssemblers const& local_assemblers,
                  IntPtDataGetter&& int_pt_data_getter);

    //! Returns the local and global error estimates.
    //!
    //! \param global_relative_error The estimated global relative error.
    //!
    //! The returned vector contains the error estimate for each cell.
    //! \pre estimate() must have been called before.
    GlobalVector const& getErrorEstimate(double& global_relative_error) const
    {
        global_relative_error = _global_relative_error;
        return *_error_estimate;
    }

private:
    std::unique_ptr<GlobalVector> _error_estimate;
    double _global_relative_error;

    // Avoids frequent reallocations.
    std::vector<double> _integration_point_values_cache;
};

namespace detail
{
// TODO this is copied code
template <int DisplacementDim>
decltype(auto) makeScalarBasis(
    typename BaseLib::DUNEGridType<DisplacementDim>::LeafGridView const&
        gridView)
{
    namespace DFB = Dune::Functions::BasisBuilder;
    return DFB::makeBasis(gridView, DFB::lagrange<1>());
}
}  // namespace detail

template <typename Grid, typename LocalAssemblers, typename IntPtDataGetter>
void ZienkiewiczZhuGradientRecoveryEstimator::estimate(
    Grid const& grid,
    LocalAssemblers const& local_assemblers,
    IntPtDataGetter&& int_pt_data_getter)
{
    auto const gridView = grid.leafGridView();
    auto const& indexSet = gridView.indexSet();

    auto const basis = detail::makeScalarBasis<Grid::dimension>(gridView);
    auto localView = basis.localView();
    auto localIndexSet = basis.localIndexSet();

    auto constexpr dimension = Grid::dimension;
    auto const num_nodes = gridView.size(dimension);

    Eigen::VectorXd A_lumped = Eigen::VectorXd::Zero(num_nodes);
    Eigen::MatrixXd rhs;

    std::vector<Dune::FieldVector<double, 1>> shape_functions;
    std::vector<double> shape_functions_;
    std::vector<double> A_local;
    std::vector<double> rhs_local;

    Eigen::Map<Eigen::MatrixXd> A_local_mat(nullptr, 1, 1);
    Eigen::Map<Eigen::MatrixXd> rhs_local_mat(nullptr, 1, 1);

    // Step 1: assemble global equation system

    double total_norm_squared = 0;
    std::size_t num_components;
    for (auto const& element : Dune::elements(gridView))
    {
        auto const idx = indexSet.index(element);
        auto const& loc_asm = *local_assemblers[idx];

        auto const int_pt_data =
            (loc_asm.*int_pt_data_getter)(_integration_point_values_cache);
        num_components = int_pt_data.rows();
        auto const num_int_pts = static_cast<std::size_t>(int_pt_data.cols());
        if (rhs.size() == 0)
            rhs = Eigen::MatrixXd::Zero(num_nodes, num_components);

        localView.bind(element);
        localIndexSet.bind(localView);

        const auto& localFiniteElement = localView.tree().finiteElement();

        auto const order =
            2 * (dimension * localFiniteElement.localBasis().order() - 1);
        auto const& quadrature = Dune::QuadratureRules<double, dimension>::rule(
            element.type(), order);
        assert(num_int_pts == quadrature.size());

        auto const geometry = element.geometry();

        for (std::size_t ip = 0; ip < num_int_pts; ++ip)
        {
            auto const quad_pos = quadrature[ip].position();
            auto const integration_element =
                geometry.integrationElement(quad_pos);

            localFiniteElement.localBasis().evaluateFunction(quad_pos,
                                                             shape_functions);
            shape_functions_.clear();
            for (auto v : shape_functions)
                shape_functions_.push_back(v);
            auto const N =
                MathLib::toVector(shape_functions_);  // column vector!
            auto const num_nodes_element = shape_functions.size();

            if (ip == 0)
            {
                A_local.clear();
                A_local.resize(num_nodes_element * num_nodes_element);
                new (&A_local_mat) Eigen::Map<Eigen::MatrixXd>(
                    A_local.data(), num_nodes_element, num_nodes_element);
                rhs_local.clear();
                rhs_local.resize(num_nodes_element * num_components);
                new (&rhs_local_mat) Eigen::Map<Eigen::MatrixXd>(
                    rhs_local.data(), num_nodes_element, num_components);
            }

            auto const weight = (quadrature[ip].weight() * integration_element);
            A_local_mat += N * N.transpose() * weight;
            rhs_local_mat += N * int_pt_data.col(ip).transpose() * weight;

            // Note: here we implicitly use the Kelvin mapping in the case of
            // symmetric tensors, i.e., that in 3D the 6-vector squared is the
            // tensor's Frobenius norm.
            total_norm_squared += int_pt_data.col(ip).squaredNorm() * weight;
        }

        for (std::size_t i = 0; i < localView.size(); ++i)
        {
            assert(localIndexSet.index(i).size() == 1);
            auto const global_index =
                static_cast<GlobalIndexType>(localIndexSet.index(i)[0]);
            A_lumped[global_index] += A_local_mat.row(i).sum();
            rhs.row(global_index).noalias() += rhs_local_mat.row(i);
        }
    }

    // Step 2: solve global equation system
    Eigen::MatrixXd nodal_values =
        (rhs.array().colwise() / A_lumped.array()).matrix();

    // Step 3: interpolate nodal values to integration points to compute element
    // residuals (errors)

    Eigen::VectorXd int_pt_data_interpolated(num_components);
    std::vector<double> nodal_data;
    _error_estimate.reset(new MathLib::EigenVector(gridView.size(0)));
    double total_residual_squared = 0;

    for (auto const& element : Dune::elements(gridView))
    {
        auto const idx = indexSet.index(element);
        auto const& loc_asm = *local_assemblers[idx];

        auto const int_pt_data =
            (loc_asm.*int_pt_data_getter)(_integration_point_values_cache);
        assert(num_components == static_cast<std::size_t>(int_pt_data.rows()));
        auto const num_int_pts = static_cast<std::size_t>(int_pt_data.cols());

        localView.bind(element);
        localIndexSet.bind(localView);

        const auto& localFiniteElement = localView.tree().finiteElement();

        auto const order =
            2 * (dimension * localFiniteElement.localBasis().order() - 1);
        auto const& quadrature = Dune::QuadratureRules<double, dimension>::rule(
            element.type(), order);
        assert(num_int_pts == quadrature.size());

        auto const geometry = element.geometry();

        nodal_data.resize(num_components * localView.size());
        auto nodal_data_mat =
            MathLib::toMatrix(nodal_data, num_components, localView.size());

        for (std::size_t i = 0; i < localView.size(); ++i)
        {
            assert(localIndexSet.index(i).size() == 1);
            auto const global_index =
                static_cast<GlobalIndexType>(localIndexSet.index(i)[0]);

            nodal_data_mat.col(i) = nodal_values.row(global_index).transpose();
        }

        double residual_squared = 0;
        for (std::size_t ip = 0; ip < num_int_pts; ++ip)
        {
            auto const quad_pos = quadrature[ip].position();
            auto const integration_element =
                geometry.integrationElement(quad_pos);

            localFiniteElement.localBasis().evaluateFunction(quad_pos,
                                                             shape_functions);
            shape_functions_.clear();
            for (auto v : shape_functions)
                shape_functions_.push_back(v);
            auto const N =
                MathLib::toVector(shape_functions_);  // column vector!

            int_pt_data_interpolated.noalias() = nodal_data_mat * N;

            // Note: here we implicitly use the Kelvin mapping in the case of
            // symmetric tensors, i.e., that in 3D the 6-vector squared is the
            // tensor's Frobenius norm.
            double const squared_norm_diff =
                (int_pt_data.col(ip) - int_pt_data_interpolated).squaredNorm();

            residual_squared += squared_norm_diff *
                                (quadrature[ip].weight() * integration_element);
        }

        total_residual_squared += residual_squared;

        // Cf.
        // Zienkiewicz, O.C., Zhu, J.Z., 1987. A simple error estimator and
        // adaptive procedure for practical engineering analysis. International
        // Journal for Numerical Methods in Engineering 24, 337–357.
        // doi:10.1002/nme.1620240206
        // ||e||_i / sqrt(||u||^2 + ||e||^2)
        if (residual_squared != 0)
            residual_squared /= (total_norm_squared + residual_squared);

        _error_estimate->set(indexSet.index(element),
                             std::sqrt(residual_squared));
    }

    _global_relative_error = std::sqrt(
        total_residual_squared / (total_residual_squared + total_norm_squared));

    INFO(
        "ZZ gradient recovery estimator: %d elements, %d vertices, total "
        "norm^2: %g, total error^2: %g, global rel error^2: %g, global rel "
        "error: %g.",
        gridView.size(0),
        gridView.size(Grid::dimension),
        total_norm_squared,
        total_residual_squared,
        total_residual_squared / (total_residual_squared + total_norm_squared),
        _global_relative_error);

    // TODO debug output
    Dune::VTKWriter<decltype(gridView)> vtkWriter(gridView);
    vtkWriter.addCellData(_error_estimate->getRawVector(), "error_estimate", 1);
    vtkWriter.addCellData(
        _error_estimate->getRawVector() * std::sqrt(gridView.size(0)),
        "error_estimate_times_sqrt_num_cells",
        1);
    vtkWriter.write(
        "zz_grad_error_estimate_" +
            std::to_string(
                ++ZienkiewiczZhuGradientRecoveryEstimator_DebugOutputCounter),
        Dune::VTK::OutputType::base64);
}

}  // namespace NumLib
