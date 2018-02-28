/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cassert>

#include "MeshLib/Elements/Utils.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/TCHSStokes/CreateLocalAssemblers.h"

#include "TCHSStokesFEM.h"
#include "TCHSStokesProcess.h"
#include "TCHSStokesProcessData.h"

namespace ProcessLib
{
namespace TCHSStokes
{
template <int VelocityDim>
TCHSStokesProcess<VelocityDim>::TCHSStokesProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    TCHSStokesProcessData<VelocityDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller,
    bool const use_monolithic_scheme)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller),
              use_monolithic_scheme),
      _process_data(std::move(process_data))
{
}

template <int VelocityDim>
bool TCHSStokesProcess<VelocityDim>::isLinear() const
{
    return false;
}

template <int VelocityDim>
MathLib::MatrixSpecifications
TCHSStokesProcess<VelocityDim>::getMatrixSpecifications(
    const int process_id) const
{
    // For the monolithic scheme or the M process (deformation) in the staggered
    // scheme.
    if (_use_monolithic_scheme || process_id == 1)
    {
        auto const& l = *_local_to_global_index_map;
        return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
                &l.getGhostIndices(), &this->_sparsity_pattern};
    }

    // For staggered scheme and H process (pressure).
    auto const& l = *_local_to_global_index_map_with_base_nodes;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &_sparsity_pattern_with_linear_element};
}

template <int VelocityDim>
void TCHSStokesProcess<VelocityDim>::constructDofTable()
{
    // Create single component dof in every of the mesh's nodes.
    _mesh_subset_all_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, &_mesh.getNodes());
    // Create single component dof in the mesh's base nodes.
    _base_nodes = MeshLib::getBaseNodes(_mesh.getElements());
    _mesh_subset_base_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, &_base_nodes);

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables of stress or strain
    std::vector<MeshLib::MeshSubsets> all_mesh_subsets_single_component;
    all_mesh_subsets_single_component.emplace_back(
        _mesh_subset_all_nodes.get());
    _local_to_global_index_map_single_component =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);

    if (_use_monolithic_scheme)
    {
        // three linearly interpolated unknowns: p, T, x_mV
        std::vector<MeshLib::MeshSubsets> all_mesh_subsets;
        all_mesh_subsets.emplace_back(_mesh_subset_base_nodes.get());
        all_mesh_subsets.emplace_back(_mesh_subset_base_nodes.get());
        all_mesh_subsets.emplace_back(_mesh_subset_base_nodes.get());

        // one quadratically interpolated unknown: v_D
        const int monolithic_process_id = 0;
        std::generate_n(
            std::back_inserter(all_mesh_subsets),
            getProcessVariables(monolithic_process_id)
                .back()
                .get()
                .getNumberOfComponents(),
            [&]() {
                return MeshLib::MeshSubsets{_mesh_subset_all_nodes.get()};
            });

        std::vector<int> const vec_n_components{1, 1, 1, VelocityDim};
        _local_to_global_index_map =
            std::make_unique<NumLib::LocalToGlobalIndexMap>(
                std::move(all_mesh_subsets), vec_n_components,
                NumLib::ComponentOrder::BY_LOCATION);
        assert(_local_to_global_index_map);
    }
    else
    {
        // TODO bugs in here.
        // For velocity equation.
        const int process_id = 1;
        std::vector<MeshLib::MeshSubsets> all_mesh_subsets;
        std::generate_n(
            std::back_inserter(all_mesh_subsets),
            getProcessVariables(process_id)[0].get().getNumberOfComponents(),
            [&]() {
                return MeshLib::MeshSubsets{_mesh_subset_all_nodes.get()};
            });

        std::vector<int> const vec_n_components{VelocityDim};
        _local_to_global_index_map =
            std::make_unique<NumLib::LocalToGlobalIndexMap>(
                std::move(all_mesh_subsets), vec_n_components,
                NumLib::ComponentOrder::BY_LOCATION);

        // For pressure equation.
        // Collect the mesh subsets with base nodes in a vector.
        std::vector<MeshLib::MeshSubsets> all_mesh_subsets_base_nodes;
        all_mesh_subsets_base_nodes.emplace_back(_mesh_subset_base_nodes.get());
        _local_to_global_index_map_with_base_nodes =
            std::make_unique<NumLib::LocalToGlobalIndexMap>(
                std::move(all_mesh_subsets_base_nodes),
                // by location order is needed for output
                NumLib::ComponentOrder::BY_LOCATION);

        _sparsity_pattern_with_linear_element = NumLib::computeSparsityPattern(
            *_local_to_global_index_map_with_base_nodes, _mesh);

        assert(_local_to_global_index_map);
        assert(_local_to_global_index_map_with_base_nodes);
    }
}

template <int VelocityDim>
void TCHSStokesProcess<VelocityDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    const int mechanical_process_id = _use_monolithic_scheme ? 0 : 1;
    // const int deformation_variable_id = _use_monolithic_scheme ? 1 : 0;
    ProcessLib::TCHSStokes::createLocalAssemblers<VelocityDim,
                                                  TCHSStokesLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        // use displacement process variable to set shape function order
        getProcessVariables(mechanical_process_id)
            .back()  // [deformation_variable_id]
            .get()
            .getShapeFunctionOrder(),
        _local_assemblers, mesh.isAxiallySymmetric(), integration_order,
        _process_data);

#if 0
    Base::_secondary_variables.addSecondaryVariable(
        "sigma_xx",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigmaXX));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_yy",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigmaYY));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_zz",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigmaZZ));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_xy",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigmaXY));

    if (VelocityDim == 3)
    {
        Base::_secondary_variables.addSecondaryVariable(
            "sigma_xz",
            makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                             &LocalAssemblerInterface::getIntPtSigmaXZ));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_yz",
            makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                             &LocalAssemblerInterface::getIntPtSigmaYZ));
    }

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_xx",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsilonXX));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_yy",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsilonYY));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_zz",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsilonZZ));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_xy",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsilonXY));

    Base::_secondary_variables.addSecondaryVariable(
        "velocity",
        makeExtrapolator(mesh.getDimension(), getExtrapolator(),
                         _local_assemblers,
                         &LocalAssemblerInterface::getIntPtDarcyVelocity));
#endif

    auto mesh_prop_nodal_p = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "pressure_interpolated",
        MeshLib::MeshItemType::Node, 1);
    mesh_prop_nodal_p->resize(mesh.getNumberOfNodes());
    _process_data.mesh_prop_nodal_p = mesh_prop_nodal_p;

    auto mesh_prop_nodal_T = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "temperature_interpolated",
        MeshLib::MeshItemType::Node, 1);
    mesh_prop_nodal_T->resize(mesh.getNumberOfNodes());
    _process_data.mesh_prop_nodal_T = mesh_prop_nodal_T;

    auto mesh_prop_nodal_xmV = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(mesh), "vapour_mass_fraction_interpolated",
        MeshLib::MeshItemType::Node, 1);
    mesh_prop_nodal_xmV->resize(mesh.getNumberOfNodes());
    _process_data.mesh_prop_nodal_xmV = mesh_prop_nodal_xmV;
}

template <int VelocityDim>
void TCHSStokesProcess<VelocityDim>::initializeBoundaryConditions()
{
    if (_use_monolithic_scheme)
    {
        const int process_id_of_hydromechancs = 0;
        initializeProcessBoundaryConditionsAndSourceTerms(
            *_local_to_global_index_map, process_id_of_hydromechancs);
        return;
    }

    // Staggered scheme:
    // for the equations of pressure
    const int hydraulic_process_id = 0;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map_with_base_nodes, hydraulic_process_id);

    // for the equations of deformation.
    const int mechanical_process_id = 1;
    initializeProcessBoundaryConditionsAndSourceTerms(
        *_local_to_global_index_map, mechanical_process_id);
}

template <int VelocityDim>
void TCHSStokesProcess<VelocityDim>::assembleConcreteProcess(
    const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b)
{
    DBUG("Assemble the equations for TCHSStokes");

    // Note: This assembly function is for the Picard nonlinear solver. Since
    // only the Newton-Raphson method is employed to simulate coupled HM
    // processes in this class, this function is actually not used so far.

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        dof_table, t, x, M, K, b, _coupled_solutions);
}

template <int VelocityDim>
void TCHSStokesProcess<VelocityDim>::assembleWithJacobianConcreteProcess(
    const double t, GlobalVector const& x, GlobalVector const& xdot,
    const double dxdot_dx, const double dx_dx, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b, GlobalMatrix& Jac)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    // For the monolithic scheme
    if (_use_monolithic_scheme)
    {
        DBUG(
            "Assemble the Jacobian of TCHSStokes for "
            "the "
            "monolithic"
            " scheme.");
        dof_tables.emplace_back(*_local_to_global_index_map);
    }
    else
    {
        // For the staggered scheme
        if (_coupled_solutions->process_id == 0)
        {
            DBUG(
                "Assemble the Jacobian equations of liquid fluid process in "
                "TCHSStokes for the staggered "
                "scheme.");
        }
        else
        {
            DBUG(
                "Assemble the Jacobian equations of mechanical process in "
                "TCHSStokes for the staggered "
                "scheme.");
        }
        dof_tables.emplace_back(*_local_to_global_index_map_with_base_nodes);
        dof_tables.emplace_back(*_local_to_global_index_map);
    }

    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, dof_tables, t, x, xdot, dxdot_dx, dx_dx, M, K, b,
        Jac, _coupled_solutions);
}

template <int VelocityDim>
void TCHSStokesProcess<VelocityDim>::preTimestepConcreteProcess(
    GlobalVector const& x,
    double const t,
    double const dt,
    const int process_id)
{
    DBUG("PreTimestep TCHSStokesProcess.");

    _process_data.delta_t = dt;

    MeshLib::Location const l(_mesh.getID(), MeshLib::MeshItemType::Node,
                              _process_data.velocity_probe_node_id);

    auto const p_index =
        _local_to_global_index_map_with_base_nodes->getLocalIndex(
            l, 0 /* p */, x.getRangeBegin(), x.getRangeEnd());
    _process_data.probed_temperature = x[p_index];

    auto const T_index =
        _local_to_global_index_map_with_base_nodes->getLocalIndex(
            l, 1 /* T */, x.getRangeBegin(), x.getRangeEnd());
    _process_data.probed_temperature = x[T_index];

    Eigen::Matrix<double, VelocityDim, 1> v;
    int const global_component_offset = 3;  // p, T, x, --> v <--
    for (int component_id = 0; component_id < VelocityDim; ++component_id)
    {
        auto const global_component_id = global_component_offset + component_id;

        auto const index =
            _local_to_global_index_map_with_base_nodes->getLocalIndex(
                l, global_component_id, x.getRangeBegin(), x.getRangeEnd());

        // TODO for PETSc the global vector must be copied. Cf.
        // ProcessOutput.cpp
        v[component_id] = x[index];
    }
    _process_data.probed_velocity = v.norm();

    if (hasMechanicalProcess(process_id))
        GlobalExecutor::executeMemberOnDereferenced(
            &LocalAssemblerInterface::preTimestep, _local_assemblers,
            *_local_to_global_index_map, x, t, dt);
}

template <int VelocityDim>
void TCHSStokesProcess<VelocityDim>::preOutputConcreteProcess(
    GlobalVector const& x, const double t, const int process_id) const
{
    DBUG("PostTimestep TCHSStokesProcess.");
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::preOutput, _local_assemblers,
        getDOFTable(process_id), x, t);
}

template <int VelocityDim>
void TCHSStokesProcess<VelocityDim>::postNonLinearSolverConcreteProcess(
    GlobalVector const& x, const double t, const int process_id)
{
    if (!hasMechanicalProcess(process_id))
    {
        return;
    }

    DBUG("PostNonLinearSolver TCHSStokesProcess.");
    // Calculate strain, stress or other internal variables of mechanics.
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::postNonLinearSolver, _local_assemblers,
        getDOFTable(process_id), x, t, _use_monolithic_scheme);
}

template <int VelocityDim>
std::tuple<NumLib::LocalToGlobalIndexMap*, bool>
TCHSStokesProcess<VelocityDim>::getDOFTableForExtrapolatorData() const
{
    const bool manage_storage = false;
    return std::make_tuple(_local_to_global_index_map_single_component.get(),
                           manage_storage);
}

template <int VelocityDim>
NumLib::LocalToGlobalIndexMap const&
TCHSStokesProcess<VelocityDim>::getDOFTable(const int process_id) const
{
    if (hasMechanicalProcess(process_id))
    {
        return *_local_to_global_index_map;
    }

    // For the equation of pressure
    return *_local_to_global_index_map_with_base_nodes;
}

}  // namespace TCHSStokes
}  // namespace ProcessLib
