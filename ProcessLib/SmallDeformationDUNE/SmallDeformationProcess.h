/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/Extrapolation/ZienkiewiczZhuGradientRecoveryEstimator.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/VectorMatrixAssemblerDUNE.h"

#include "LocalAssemblerInterface.h"
#include "SmallDeformationProcessData.h"

namespace ProcessLib
{
namespace SmallDeformationDUNE
{
template <int DisplacementDim>
class SmallDeformationProcess final : public Process
{
    using Base = Process;

public:
    SmallDeformationProcess(
        MeshLib::FEMMesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
            process_variables,
        SmallDeformationProcessData<DisplacementDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller);

    //! \name ODESystem interface
    //! @{
    bool isLinear() const override;
    //! @}

    void setInitialConditions(const int process_id, const double t,
                              GlobalVector& x) override;

    GlobalVector const* estimateError(GlobalVector const& x,
                                      double& global_relative_error) override;
    bool refine(std::vector<char> const& elements_for_refinement) override;

private:
    using LocalAssemblerInterface =
        SmallDeformationLocalAssemblerInterface<DisplacementDim>;

    void initializeConcreteProcess(NumLib::AbstractDOFTable const& dof_table,
                                   MeshLib::FEMMesh const& mesh,
                                   unsigned const integration_order) override;

    void assembleConcreteProcess(const double t, GlobalVector const& x,
                                 GlobalMatrix& M, GlobalMatrix& K,
                                 GlobalVector& b) override;

    void assembleWithJacobianConcreteProcess(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac) override;

    void preTimestepConcreteProcess(GlobalVector const& x, double const t,
                                    double const dt,
                                    const int /*process_id*/) override;

    void computeSparsityPattern(MeshLib::DUNEMesh<DisplacementDim> const& mesh,
                                bool, const MeshLib::DUNEIdToIdxMappings&);

    void createLocalAssemblers(MeshLib::DUNEMesh<DisplacementDim> const& mesh,
                               bool global_refine,
                               const MeshLib::DUNEIdToIdxMappings&);

    void initializeExtrapolator();
    void updateDOFtables(MeshLib::DUNEMesh<DisplacementDim> const&, bool,
                         const MeshLib::DUNEIdToIdxMappings&);

    SmallDeformationProcessData<DisplacementDim> _process_data;

    VectorMatrixAssemblerDUNE _global_assembler_dune;
    std::vector<std::unique_ptr<LocalAssemblerInterface>> _local_assemblers;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap>
        _local_to_global_index_map_single_component;
    MeshLib::PropertyVector<double>* _nodal_forces = nullptr;

    struct ErrorEstimatorData
    {
        // TODO [DUNE] make const
        std::unique_ptr<NumLib::ZienkiewiczZhuGradientRecoveryEstimator>
            estimator_algorithm;
        NumLib::AbstractDOFTable* dof_table_single_component = nullptr;
    };
    ErrorEstimatorData _error_estimator_data;
};

extern template class SmallDeformationProcess<2>;
extern template class SmallDeformationProcess<3>;

}  // namespace SmallDeformationDUNE
}  // namespace ProcessLib
