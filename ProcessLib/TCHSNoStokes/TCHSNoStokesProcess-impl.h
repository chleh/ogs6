/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <omp.h>
#include <cassert>

#include "MeshLib/Elements/Utils.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/TCHSNoStokes/CreateLocalAssemblers.h"

#include "TCHSNoStokesFEM.h"
#include "TCHSNoStokesProcess.h"
#include "TCHSNoStokesProcessData.h"

namespace ProcessLib
{
namespace TCHSNoStokes
{
template <int VelocityDim>
TCHSNoStokesProcess<VelocityDim>::TCHSNoStokesProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    TCHSNoStokesProcessData<VelocityDim>&& process_data,
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
bool TCHSNoStokesProcess<VelocityDim>::isLinear() const
{
    return false;
}

template <int VelocityDim>
MathLib::MatrixSpecifications
TCHSNoStokesProcess<VelocityDim>::getMatrixSpecifications(
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
    else
    {
        OGS_FATAL("wrong scheme");
    }
}

template <int VelocityDim>
void TCHSNoStokesProcess<VelocityDim>::constructDofTable()
{
    // Create single component dof in every of the mesh's nodes.
    _mesh_subset_all_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, &_mesh.getNodes());

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
        all_mesh_subsets.emplace_back(_mesh_subset_all_nodes.get());
        all_mesh_subsets.emplace_back(_mesh_subset_all_nodes.get());
        all_mesh_subsets.emplace_back(_mesh_subset_all_nodes.get());

        std::vector<int> const vec_n_components{1, 1, 1};
        _local_to_global_index_map =
            std::make_unique<NumLib::LocalToGlobalIndexMap>(
                std::move(all_mesh_subsets), vec_n_components,
                NumLib::ComponentOrder::BY_LOCATION);
        assert(_local_to_global_index_map);
    }
    else
    {
        OGS_FATAL("wrong scheme");
    }
}

template <int VelocityDim>
void TCHSNoStokesProcess<VelocityDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    const int mechanical_process_id = _use_monolithic_scheme ? 0 : 1;
    // const int deformation_variable_id = _use_monolithic_scheme ? 1 : 0;
    ProcessLib::TCHSNoStokes::createLocalAssemblers<VelocityDim,
                                                    TCHSNoStokesLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        // use displacement process variable to set shape function order
        getProcessVariables(mechanical_process_id)
            .back()  // [deformation_variable_id]
            .get()
            .getShapeFunctionOrder(),
        _local_assemblers, mesh.isAxiallySymmetric(), integration_order,
        _process_data);

    Base::_secondary_variables.addSecondaryVariable(
        "solid_density",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSolidDensity));
    Base::_secondary_variables.addSecondaryVariable(
        "reaction_rate",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtReactionRate));
    Base::_secondary_variables.addSecondaryVariable(
        "porosity",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtPorosity));

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

    auto create_mesh_prop = [&](std::string const& name,
                                MeshLib::MeshItemType type, int comp) {
        auto* mesh_prop = MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), name, type, comp);
        switch (type)
        {
            case MeshLib::MeshItemType::Node:
                mesh_prop->resize(mesh.getNumberOfNodes() * comp);
                break;
            case MeshLib::MeshItemType::Cell:
                mesh_prop->resize(mesh.getNumberOfElements() * comp);
                break;
            default:
                OGS_FATAL("Unsupported case.");
        }
        return mesh_prop;
    };

    _process_data.mesh_prop_cell_hat_rho_SR =
        create_mesh_prop("reaction_rate_cell", MeshLib::MeshItemType::Cell, 1);
    _process_data.mesh_prop_cell_reaction_enthalpy =
        create_mesh_prop("reaction_enthalpy", MeshLib::MeshItemType::Cell, 1);
    _process_data.mesh_prop_cell_rho_SR =
        create_mesh_prop("rho_SR", MeshLib::MeshItemType::Cell, 1);
    _process_data.mesh_prop_cell_rho_GR =
        create_mesh_prop("rho_GR", MeshLib::MeshItemType::Cell, 1);
    _process_data.mesh_prop_cell_lambda =
        create_mesh_prop("lambda", MeshLib::MeshItemType::Cell, VelocityDim);
    _process_data.mesh_prop_cell_cpS =
        create_mesh_prop("cpS", MeshLib::MeshItemType::Cell, 1);
    _process_data.mesh_prop_cell_cpG =
        create_mesh_prop("cpG", MeshLib::MeshItemType::Cell, 1);
    _process_data.mesh_prop_cell_v_Darcy = create_mesh_prop(
        "darcy_velocity", MeshLib::MeshItemType::Cell, VelocityDim);
    _process_data.mesh_prop_cell_mass_flux =
        create_mesh_prop("mass_flux", MeshLib::MeshItemType::Cell, VelocityDim);
    _process_data.mesh_prop_cell_vapour_mass_flux = create_mesh_prop(
        "vapour_mass_flux", MeshLib::MeshItemType::Cell, VelocityDim);
}

template <int VelocityDim>
void TCHSNoStokesProcess<VelocityDim>::initializeBoundaryConditions()
{
    if (_use_monolithic_scheme)
    {
        const int process_id_of_hydromechancs = 0;
        initializeProcessBoundaryConditionsAndSourceTerms(
            *_local_to_global_index_map, process_id_of_hydromechancs);
        return;
    }
    else
    {
        OGS_FATAL("wrong scheme");
    }
}

template <int VelocityDim>
void TCHSNoStokesProcess<VelocityDim>::assembleConcreteProcess(
    const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b)
{
    DBUG("Assemble the equations for TCHSNoStokes");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    auto const n = _local_assemblers.size();

#pragma omp parallel for
    for (std::size_t i = 0; i < n; i++)
    {
        _global_assembler.assemble(i, *_local_assemblers[i], dof_table, t, x, M,
                                   K, b, _coupled_solutions);
    }
}

template <int VelocityDim>
void TCHSNoStokesProcess<VelocityDim>::assembleWithJacobianConcreteProcess(
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
            "Assemble the Jacobian of TCHSNoStokes for "
            "the "
            "monolithic"
            " scheme.");
        dof_tables.emplace_back(*_local_to_global_index_map);
    }
    else
    {
        OGS_FATAL("wrong scheme");
    }

    auto const n = _local_assemblers.size();

#pragma omp parallel for
    for (std::size_t i = 0; i < n; i++)
    {
        _global_assembler.assembleWithJacobian(
            i, *_local_assemblers[i], dof_tables, t, x, xdot, dxdot_dx, dx_dx,
            M, K, b, Jac, _coupled_solutions);
    }
}

template <int VelocityDim>
void TCHSNoStokesProcess<VelocityDim>::preTimestepConcreteProcess(
    GlobalVector const& x,
    double const t,
    double const dt,
    const int process_id)
{
    DBUG("PreTimestep TCHSNoStokesProcess.");

    _process_data.delta_t = dt;

    MeshLib::Location const l(_mesh.getID(), MeshLib::MeshItemType::Node,
                              _process_data.velocity_probe_node_id);

    auto const p_index = _local_to_global_index_map->getLocalIndex(
        l, 0 /* p */, x.getRangeBegin(), x.getRangeEnd());
    _process_data.probed_pressure = x[p_index];

    auto const T_index = _local_to_global_index_map->getLocalIndex(
        l, 1 /* T */, x.getRangeBegin(), x.getRangeEnd());
    _process_data.probed_temperature = x[T_index];

    if (hasMechanicalProcess(process_id))
        GlobalExecutor::executeMemberOnDereferenced(
            &LocalAssemblerInterface::preTimestep, _local_assemblers,
            *_local_to_global_index_map, x, t, dt);
}

template <int VelocityDim>
void TCHSNoStokesProcess<VelocityDim>::preOutputConcreteProcess(
    GlobalVector const& x, const double t, const int process_id) const
{
    DBUG("PostTimestep TCHSNoStokesProcess.");
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::preOutput, _local_assemblers,
        getDOFTable(process_id), x, t);
}

template <int VelocityDim>
void TCHSNoStokesProcess<VelocityDim>::postNonLinearSolverConcreteProcess(
    GlobalVector const& x, const double t, const int process_id)
{
    if (!hasMechanicalProcess(process_id))
    {
        return;
    }

    DBUG("PostNonLinearSolver TCHSNoStokesProcess.");
    // Calculate strain, stress or other internal variables of mechanics.
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::postNonLinearSolver, _local_assemblers,
        getDOFTable(process_id), x, t, _use_monolithic_scheme);
}

template <int VelocityDim>
void TCHSNoStokesProcess<VelocityDim>::preIterationConcreteProcess(
    const unsigned iter, GlobalVector const& /*x*/)
{
    for (auto& mat : _process_data.materials)
    {
        mat.second.reactive_solid->preIteration(iter);
        mat.second.reaction_rate->preIteration(iter);
    }
}

template <int VelocityDim>
std::tuple<NumLib::LocalToGlobalIndexMap*, bool>
TCHSNoStokesProcess<VelocityDim>::getDOFTableForExtrapolatorData() const
{
    const bool manage_storage = false;
    return std::make_tuple(_local_to_global_index_map_single_component.get(),
                           manage_storage);
}

template <int VelocityDim>
NumLib::LocalToGlobalIndexMap const&
TCHSNoStokesProcess<VelocityDim>::getDOFTable(const int process_id) const
{
    if (hasMechanicalProcess(process_id))
    {
        return *_local_to_global_index_map;
    }
    else
    {
        OGS_FATAL("wrong scheme");
    }
}

}  // namespace TCHSNoStokes
}  // namespace ProcessLib
