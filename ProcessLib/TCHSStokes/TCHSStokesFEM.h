/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>

#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "LocalAssemblerInterface.h"
#include "TCHSStokesProcessData.h"

namespace ProcessLib
{
namespace TCHSStokes
{
template <typename BMatricesType, typename ShapeMatrixTypeVelocity,
          typename ShapeMatricesTypePressure, int VelocityDim, int NPoints>
struct IntegrationPointData final
{
    typename ShapeMatrixTypeVelocity::template MatrixType<VelocityDim,
                                                          NPoints * VelocityDim>
        H;

    typename ShapeMatrixTypeVelocity::NodalRowVectorType N_v;
    typename ShapeMatrixTypeVelocity::GlobalDimNodalMatrixType dNdx_v;

    typename ShapeMatricesTypePressure::NodalRowVectorType N_p;
    typename ShapeMatricesTypePressure::GlobalDimNodalMatrixType dNdx_p;
    double integration_weight;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

/// Used by for extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>> N_u;
};

template <typename ShapeFunctionVelocity, typename ShapeFunctionPressure,
          typename IntegrationMethod, int VelocityDim>
class TCHSStokesLocalAssembler : public LocalAssemblerInterface
{
public:
    using ShapeMatricesTypeVelocity =
        ShapeMatrixPolicyType<ShapeFunctionVelocity, VelocityDim>;

    // Types for pressure.
    using ShapeMatricesTypePressure =
        ShapeMatrixPolicyType<ShapeFunctionPressure, VelocityDim>;

    static int const KelvinVectorSize =
        MathLib::KelvinVector::KelvinVectorDimensions<VelocityDim>::value;
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;

    TCHSStokesLocalAssembler(TCHSStokesLocalAssembler const&) = delete;
    TCHSStokesLocalAssembler(TCHSStokesLocalAssembler&&) = delete;

    TCHSStokesLocalAssembler(MeshLib::Element const& e,
                             std::size_t const /*local_matrix_size*/,
                             bool const is_axially_symmetric,
                             unsigned const integration_order,
                             TCHSStokesProcessData<VelocityDim>& process_data);

    void assemble(double const /*t*/, std::vector<double> const& /*local_x*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_rhs_data*/) override;

    void assembleWithJacobian(double const /*t*/,
                              std::vector<double> const& /*local_x*/,
                              std::vector<double> const& /*local_xdot*/,
                              const double /*dxdot_dx*/, const double /*dx_dx*/,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& /*local_rhs_data*/,
                              std::vector<double>& /*local_Jac_data*/) override
    {
        OGS_FATAL("not implemented");
    }

    void assembleWithJacobianForStaggeredScheme(
        double const /*t*/, std::vector<double> const& /*local_xdot*/,
        const double /*dxdot_dx*/, const double /*dx_dx*/,
        std::vector<double>& /*local_M_data*/,
        std::vector<double>& /*local_K_data*/,
        std::vector<double>& /*local_b_data*/,
        std::vector<double>& /*local_Jac_data*/,
        LocalCoupledSolutions const& /*local_coupled_solutions*/) override
    {
        OGS_FATAL("not implemented");
    }

    void preTimestepConcrete(std::vector<double> const& /*local_x*/,
                             double const /*t*/,
                             double const /*delta_t*/) override
    {
    }

    void preOutputConcrete(std::vector<double> const& local_x,
                           const double /*t*/) override;

    void postNonLinearSolverConcrete(std::vector<double> const& local_x,
                                     double const t,
                                     bool const use_monolithic_scheme) override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N_u = _ip_data[integration_point].N_v;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N_u.data(), N_u.size());
    }

private:
    TCHSStokesProcessData<VelocityDim>& _process_data;

    using BMatricesType = BMatrixPolicyType<ShapeFunctionVelocity, VelocityDim>;
    using IpData =
        IntegrationPointData<BMatricesType, ShapeMatricesTypeVelocity,
                             ShapeMatricesTypePressure, VelocityDim,
                             ShapeFunctionVelocity::NPOINTS>;
    std::vector<IpData, Eigen::aligned_allocator<IpData>> _ip_data;

    IntegrationMethod _integration_method;
    MeshLib::Element const& _element;
    bool const _is_axially_symmetric;

    static const int kelvin_vector_size =
        MathLib::KelvinVector::KelvinVectorDimensions<VelocityDim>::value;

    /*
    static const int Indices[] = {0, ShapeFunctionPressure::NPOINTS,
                                  2 * ShapeFunctionPressure::NPOINTS,
                                  3 * ShapeFunctionPressure::NPOINTS};
                             */

    struct Block
    {
    private:
        enum class BlockName : int
        {
            P_,
            T_,
            X_,
            V_
        };

    public:
        static const std::integral_constant<BlockName, BlockName::P_> P;
        static const std::integral_constant<BlockName, BlockName::T_> T;
        static const std::integral_constant<BlockName, BlockName::X_> X;
        static const std::integral_constant<BlockName, BlockName::V_> V;

        template <BlockName Row, BlockName Col, typename Matrix>
        static decltype(auto) block(Matrix&& mat,
                                    std::integral_constant<BlockName, Row>
                                        row,
                                    std::integral_constant<BlockName, Col>
                                        col)
        {
            return mat.template block<size(row), size(col)>(index(row),
                                                            index(col));
        }

        template <BlockName Row, typename Vector>
        static decltype(auto) segment(
            Vector&& vec, std::integral_constant<BlockName, Row> row)
        {
            return vec.template segment<size(row)>(index(row));
        }

        static constexpr int index(
            std::integral_constant<BlockName, BlockName::P_>)
        {
            return 0;
        }

        static constexpr int index(
            std::integral_constant<BlockName, BlockName::T_>)
        {
            return ShapeFunctionPressure::NPOINTS;
        }

        static constexpr int index(
            std::integral_constant<BlockName, BlockName::X_>)
        {
            return 2 * ShapeFunctionPressure::NPOINTS;
        }

        static constexpr int index(
            std::integral_constant<BlockName, BlockName::V_>)
        {
            return 3 * ShapeFunctionPressure::NPOINTS;
        }

        static constexpr int size(
            std::integral_constant<BlockName, BlockName::P_>)
        {
            return ShapeFunctionPressure::NPOINTS;
        }

        static constexpr int size(
            std::integral_constant<BlockName, BlockName::T_>)
        {
            return ShapeFunctionPressure::NPOINTS;
        }

        static constexpr int size(
            std::integral_constant<BlockName, BlockName::X_>)
        {
            return ShapeFunctionPressure::NPOINTS;
        }

        static constexpr int size(
            std::integral_constant<BlockName, BlockName::V_>)
        {
            return ShapeFunctionVelocity::NPOINTS * VelocityDim;
        }
    };
};

}  // namespace TCHSStokes
}  // namespace ProcessLib

#include "TCHSStokesFEM-impl.h"
