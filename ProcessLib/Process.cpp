/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Process.h"

#include "BaseLib/Functional.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"
#include "ProcessVariable.h"

namespace ProcessLib
{
Process::Process(
    MeshLib::Mesh& mesh,
    NonlinearSolver& nonlinear_solver,
    std::unique_ptr<TimeDiscretization>&& time_discretization,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
    SecondaryVariableCollection&& secondary_variables,
    ProcessOutput&& process_output)
    : _mesh(mesh),
      _secondary_variables(std::move(secondary_variables)),
      _process_output(std::move(process_output)),
      _global_assembler(std::move(jacobian_assembler)),
      _nonlinear_solver(nonlinear_solver),
      _time_discretization(std::move(time_discretization)),
      _process_variables(std::move(process_variables))
{
}

void Process::output(std::string const& file_name,
                     const unsigned /*timestep*/,
                     GlobalVector const& x) const
{
    doProcessOutput(file_name, x, _mesh, *_local_to_global_index_map,
                    _process_variables, _secondary_variables, _process_output);
}

void Process::initialize()
{
    DBUG("Initialize process.");

    DBUG("Construct dof mappings.");
    constructDofTable();

    DBUG("Compute sparsity pattern");
    computeSparsityPattern();

    DBUG("Initialize the extrapolator");
    initializeExtrapolator();

    initializeConcreteProcess(*_local_to_global_index_map, _mesh,
                              _integration_order);

    DBUG("Initialize boundary conditions.");
    _boundary_conditions.addBCsForProcessVariables(
        _process_variables, *_local_to_global_index_map, _integration_order);
}

void Process::setInitialConditions(GlobalVector& x)
{
    DBUG("Set initial conditions.");
    for (int variable_id = 0;
         variable_id < static_cast<int>(_process_variables.size());
         ++variable_id)
    {
        ProcessVariable& pv = _process_variables[variable_id];
        for (int component_id = 0; component_id < pv.getNumberOfComponents();
             ++component_id)
        {
            setInitialConditions(pv, variable_id, component_id, x);
        }
    }
}

MathLib::MatrixSpecifications Process::getMatrixSpecifications() const
{
    auto const& l = *_local_to_global_index_map;
    return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
            &l.getGhostIndices(), &_sparsity_pattern};
}

void Process::assemble(const double t, GlobalVector const& x, GlobalMatrix& M,
                       GlobalMatrix& K, GlobalVector& b)
{
    assembleConcreteProcess(t, x, M, K, b);
    _boundary_conditions.apply(t, x, K, b);
}

void Process::assembleWithJacobian(const double t, GlobalVector const& x,
                                   GlobalVector const& xdot,
                                   const double dxdot_dx, const double dx_dx,
                                   GlobalMatrix& M, GlobalMatrix& K,
                                   GlobalVector& b, GlobalMatrix& Jac)
{
#if 0
    _jacobian_assembler->assembleWithJacobian(
        BaseLib::easyBind(&Process::assemble, this),
        BaseLib::easyBind(&Process::assembleWithJacobianAnalytical, this), t, x,
        xdot, dxdot_dx, dx_dx, M, K, b, Jac);
}

void Process::assembleWithJacobianAnalytical(
    const double t, GlobalVector const& x, GlobalVector const& xdot,
    const double dxdot_dx, const double dx_dx, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b, GlobalMatrix& Jac)
{
#endif
    assembleWithJacobianConcreteProcess(t, x, xdot, dxdot_dx, dx_dx, M, K, b, Jac);

    // Call global assembler for each Neumann boundary local assembler.
    for (auto const& bc : _neumann_bcs)
        bc->integrate(t, b);

    // TODO apply BCs to Jacobian.
}

void Process::constructDofTable()
{
    // Create single component dof in every of the mesh's nodes.
    _mesh_subset_all_nodes.reset(
        new MeshLib::MeshSubset(_mesh, &_mesh.getNodes()));

    // Collect the mesh subsets in a vector.
    std::vector<std::unique_ptr<MeshLib::MeshSubsets>> all_mesh_subsets;
    for (ProcessVariable const& pv : _process_variables)
    {
        std::generate_n(
            std::back_inserter(all_mesh_subsets),
            pv.getNumberOfComponents(),
            [&]() {
                return std::unique_ptr<MeshLib::MeshSubsets>{
                    new MeshLib::MeshSubsets{_mesh_subset_all_nodes.get()}};
            });
    }

    _local_to_global_index_map.reset(new NumLib::LocalToGlobalIndexMap(
        std::move(all_mesh_subsets), NumLib::ComponentOrder::BY_LOCATION));
}

void Process::initializeExtrapolator()
{
    NumLib::LocalToGlobalIndexMap const* dof_table_single_component;
    bool manage_storage;

    if (_local_to_global_index_map->getNumberOfComponents() == 1)
    {
        // For single-variable-single-component processes reuse the existing DOF
        // table.
        dof_table_single_component = _local_to_global_index_map.get();
        manage_storage = false;
    }
    else
    {
        // Otherwise construct a new DOF table.
        std::vector<std::unique_ptr<MeshLib::MeshSubsets>>
            all_mesh_subsets_single_component;
        all_mesh_subsets_single_component.emplace_back(
            new MeshLib::MeshSubsets(_mesh_subset_all_nodes.get()));

        dof_table_single_component = new NumLib::LocalToGlobalIndexMap(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);
        manage_storage = true;
    }

    std::unique_ptr<NumLib::Extrapolator> extrapolator(
        new NumLib::LocalLinearLeastSquaresExtrapolator(
            *dof_table_single_component));

    // TODO Later on the DOF table can change during the simulation!
    _extrapolator_data = ExtrapolatorData(
        std::move(extrapolator), dof_table_single_component, manage_storage);
}

void Process::setInitialConditions(ProcessVariable const& variable,
                                   int const variable_id,
                                   int const component_id,
                                   GlobalVector& x)
{
    std::size_t const n_nodes = _mesh.getNumberOfNodes();
    for (std::size_t node_id = 0; node_id < n_nodes; ++node_id)
    {
        MeshLib::Location const l(_mesh.getID(), MeshLib::MeshItemType::Node,
                                  node_id);
        auto global_index = std::abs(_local_to_global_index_map->getGlobalIndex(
            l, variable_id, component_id));
#ifdef USE_PETSC
        // The global indices of the ghost entries of the global matrix or
        // the global vectors need to be set as negative values for equation
        // assembly, however the global indices start from zero.  Therefore,
        // any ghost entry with zero index is assigned an negative value of
        // the vector size or the matrix dimension.  To assign the initial
        // value for the ghost entries, the negative indices of the ghost
        // entries are restored to zero.
        if (global_index == x.size())
            global_index = 0;
#endif
        x.set(global_index,
              variable.getInitialConditionValue(node_id, component_id));
    }
}

void Process::computeSparsityPattern()
{
    _sparsity_pattern =
        NumLib::computeSparsityPattern(*_local_to_global_index_map, _mesh);
}

ProcessVariable& findProcessVariable(
    std::vector<ProcessVariable> const& variables,
    BaseLib::ConfigTree const& pv_config, std::string const& tag)
{
    // Find process variable name in process config.
    //! \ogs_file_special
    std::string const name = pv_config.getConfigParameter<std::string>(tag);

    // Find corresponding variable by name.
    auto variable = std::find_if(
        variables.cbegin(), variables.cend(),
        [&name](ProcessVariable const& v) { return v.getName() == name; });

    if (variable == variables.end())
    {
        OGS_FATAL(
            "Could not find process variable '%s' in the provided variables "
            "list for config tag <%s>.",
            name.c_str(), tag.c_str());
    }
    DBUG("Found process variable \'%s\' for config tag <%s>.",
         variable->getName().c_str(), tag.c_str());

    // Const cast is needed because of variables argument constness.
    return const_cast<ProcessVariable&>(*variable);
}

std::vector<std::reference_wrapper<ProcessVariable>> findProcessVariables(
    std::vector<ProcessVariable> const& variables,
    BaseLib::ConfigTree const& process_config,
    std::initializer_list<std::string>
        tag_names)
{
    std::vector<std::reference_wrapper<ProcessVariable>> vars;
    vars.reserve(tag_names.size());

    //! \ogs_file_param{process__process_variables}
    auto const pv_conf = process_config.getConfigSubtree("process_variables");

    for (auto const& tag : tag_names)
    {
        vars.emplace_back(findProcessVariable(variables, pv_conf, tag));
    }

    return vars;
}

}  // namespace ProcessLib
