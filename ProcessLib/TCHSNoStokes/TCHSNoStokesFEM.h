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
#include "TCHSNoStokesProcessData.h"

namespace ProcessLib
{
namespace TCHSNoStokes
{
template <typename BMatricesType, typename ShapeMatrixTypeVelocity,
          typename ShapeMatricesTypePressure, int VelocityDim, int NPoints>
struct IntegrationPointData final
{
    typename ShapeMatricesTypePressure::NodalRowVectorType N_1;
    typename ShapeMatricesTypePressure::GlobalDimNodalMatrixType dNdx_1;
    double integration_weight;

    std::unique_ptr<MaterialLib::ReactiveSolidState> reactive_solid_state;
    std::unique_ptr<MaterialLib::ReactionRateData> reaction_rate_data;
    std::unique_ptr<MaterialLib::ReactiveSolidRate>
        reaction_rate;  // \hat{\rho}_{SR}

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

template <typename ShapeFunctionVelocity, typename ShapeFunctionPressure,
          typename IntegrationMethod, int VelocityDim>
class TCHSNoStokesLocalAssembler : public LocalAssemblerInterface
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

    TCHSNoStokesLocalAssembler(TCHSNoStokesLocalAssembler const&) = delete;
    TCHSNoStokesLocalAssembler(TCHSNoStokesLocalAssembler&&) = delete;

    TCHSNoStokesLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        TCHSNoStokesProcessData<VelocityDim>& process_data);

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
        auto const mat_id = _process_data.material_ids[_element.getID()];
        for (auto& ip_data : _ip_data)
        {
            if (ip_data.reactive_solid_state)
                _process_data.materials.at(mat_id).reactive_solid->preTimestep(
                    *ip_data.reactive_solid_state, *ip_data.reaction_rate);
            if (ip_data.reaction_rate_data)
                ip_data.reaction_rate_data->preTimestep(
                    *ip_data.reactive_solid_state, *ip_data.reaction_rate);
        }
    }

    void preOutputConcrete(std::vector<double> const& local_x,
                           const double /*t*/) override;

    void postNonLinearSolverConcrete(std::vector<double> const& local_x,
                                     double const t,
                                     bool const use_monolithic_scheme) override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N_1 = _ip_data[integration_point].N_1;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N_1.data(), N_1.size());
    }

    std::vector<double> const& getIntPtSolidDensity(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override;

    std::vector<double> const& getIntPtReactionRate(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override;

    std::vector<double> const& getIntPtPorosity(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::vector<double>& /*cache*/) const override;

private:
    TCHSNoStokesProcessData<VelocityDim> const& _process_data;

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

    struct Block
    {
    private:
        enum class BlockName : int
        {
            P_,
            T_,
            X_
        };

    public:
        static constexpr std::integral_constant<BlockName, BlockName::P_> P =
            std::integral_constant<BlockName, BlockName::P_>{};
        static constexpr std::integral_constant<BlockName, BlockName::T_> T =
            std::integral_constant<BlockName, BlockName::T_>{};
        static constexpr std::integral_constant<BlockName, BlockName::X_> X =
            std::integral_constant<BlockName, BlockName::X_>{};

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

        template <BlockName B>
        static decltype(auto) mapVectorSegment(
            std::vector<double> const& v,
            std::integral_constant<BlockName, B>
                b)
        {
            return Eigen::Map<Eigen::Matrix<double, size(b), 1> const>(
                v.data() + index(b), size(b));
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
    };
};

}  // namespace TCHSNoStokes
}  // namespace ProcessLib

#include "TCHSNoStokesFEM-impl.h"
