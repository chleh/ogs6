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

#include "MaterialLib/SolidModels/KelvinVector.h"
#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "IncompressibleStokesBrinkmanProcessData.h"
#include "LocalAssemblerInterface.h"

namespace ProcessLib
{
namespace IncompressibleStokesBrinkman
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
class IncompressibleStokesBrinkmanLocalAssembler
    : public LocalAssemblerInterface
{
public:
    using ShapeMatricesTypeVelocity =
        ShapeMatrixPolicyType<ShapeFunctionVelocity, VelocityDim>;

    // Types for pressure.
    using ShapeMatricesTypePressure =
        ShapeMatrixPolicyType<ShapeFunctionPressure, VelocityDim>;

    static int const KelvinVectorSize =
        ProcessLib::KelvinVectorDimensions<VelocityDim>::value;
    using Invariants = MaterialLib::SolidModels::Invariants<KelvinVectorSize>;

    IncompressibleStokesBrinkmanLocalAssembler(
        IncompressibleStokesBrinkmanLocalAssembler const&) = delete;
    IncompressibleStokesBrinkmanLocalAssembler(
        IncompressibleStokesBrinkmanLocalAssembler&&) = delete;

    IncompressibleStokesBrinkmanLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        IncompressibleStokesBrinkmanProcessData<VelocityDim>& process_data);

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

    void postTimestepConcrete(std::vector<double> const& local_x) override;

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
    IncompressibleStokesBrinkmanProcessData<VelocityDim>& _process_data;

    using BMatricesType = BMatrixPolicyType<ShapeFunctionVelocity, VelocityDim>;
    using IpData =
        IntegrationPointData<BMatricesType, ShapeMatricesTypeVelocity,
                             ShapeMatricesTypePressure, VelocityDim,
                             ShapeFunctionVelocity::NPOINTS>;
    std::vector<IpData, Eigen::aligned_allocator<IpData>> _ip_data;

    IntegrationMethod _integration_method;
    MeshLib::Element const& _element;
    bool const _is_axially_symmetric;

    static const int pressure_index = 0;
    static const int pressure_size = ShapeFunctionPressure::NPOINTS;
    static const int velocity_index = ShapeFunctionPressure::NPOINTS;
    static const int velocity_size =
        ShapeFunctionVelocity::NPOINTS * VelocityDim;
    static const int kelvin_vector_size =
        KelvinVectorDimensions<VelocityDim>::value;
};

}  // namespace IncompressibleStokesBrinkman
}  // namespace ProcessLib

#include "IncompressibleStokesBrinkmanFEM-impl.h"
