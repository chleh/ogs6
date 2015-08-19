/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_TESPROCESS_IMPL_H_
#define PROCESS_LIB_TESPROCESS_IMPL_H_

#include <cassert>

#include "logog/include/logog.hpp"

#include "TESProcess.h"


template <class ConfigTree>
static ProcessLib::ProcessVariable const*
find_variable(ConfigTree const& config,
              std::string const& variable_role,
              std::vector<ProcessLib::ProcessVariable> const& variables)
{
    std::string const name = config.template get<std::string>(variable_role);

    auto variable = std::find_if(variables.cbegin(), variables.cend(),
            [&name](ProcessLib::ProcessVariable const& v) {
                return v.getName() == name;
            });

    if (variable == variables.end())
    {
        ERR("Expected process variable \'%s\' not found in provided variables list.",
            name.c_str());
        std::abort();
    }


    DBUG("Found process variable %s for role %s.",
         name.c_str(), variable_role.c_str());

    return &*variable;
}

namespace ProcessLib
{

namespace TES
{

template<typename GlobalSetup>
TESProcess<GlobalSetup>::
TESProcess(MeshLib::Mesh& mesh,
        std::vector<ProcessVariable> const& variables,
        ConfigTree const& config)
    : Process(mesh)
{
    DBUG("Create TESProcess.");

    // primary variables
    {
        const std::string vars[NODAL_DOF] = { "fluid_pressure",
                                              "temperature",
                                              "vapour_mass_fraction" };

        ConfigTree proc_vars = config.get_child("process_variables");

        for (unsigned i=0; i<NODAL_DOF; ++i)
        {
            auto variable = find_variable(proc_vars, vars[i], variables);

            _process_vars[i] = const_cast<ProcessVariable*>(variable);
        }
    }

    // secondary variables
    {
        const std::string vars[NODAL_DOF_2ND] = { "solid_density" };

        ConfigTree proc_vars = config.get_child("secondary_variables");

        for (unsigned i=0; i<NODAL_DOF_2ND && i<1; ++i)
        {
            auto variable = find_variable(proc_vars, vars[i], variables);

            _secondary_process_vars[i] = const_cast<ProcessVariable*>(variable);
        }
    }

    // reactive system
    {
        std::string rsys = config.get<std::string>("reactive_system");
        DBUG("reactive_system: %s", rsys.c_str());

        _materials._adsorption = Ads::Adsorption::newInstance(rsys);
    }
}

template<typename GlobalSetup>
template <unsigned GlobalDim>
void
TESProcess<GlobalSetup>::
createLocalAssemblers()
{
    DBUG("Create local assemblers.");
    // Populate the vector of local assemblers.
    _local_assemblers.resize(_mesh.getNElements());
    // Shape matrices initializer
    using LocalDataInitializer = AssemblerLib::LocalDataInitializer<
        TES::LocalAssemblerDataInterface,
        TES::LocalAssemblerData,
        typename GlobalSetup::MatrixType,
        typename GlobalSetup::VectorType,
        GlobalDim>;

    LocalDataInitializer initializer;

    using LocalAssemblerBuilder =
        AssemblerLib::LocalAssemblerBuilder<
            MeshLib::Element,
            LocalDataInitializer>;

    LocalAssemblerBuilder local_asm_builder(
        initializer, *_local_to_global_index_map);

    DBUG("Calling local assembler builder for all mesh elements.");
    _global_setup.execute(
                local_asm_builder,
                _mesh.getElements(),
                _local_assemblers,
                _integration_order,
                this);

    DBUG("Create global assembler.");
    _global_assembler.reset(
        new GlobalAssembler(*_A, *_rhs, *_local_to_global_index_map));
    _global_assembler->setX(_x.get());
    _global_assembler->setSecondaryVariables(
                _secondary_variables.get(),
                _local_to_global_index_map_secondary.get());

    for (unsigned i=0; i<NODAL_DOF; ++i)
    {
        // TODO [CL] can't that be (partially) moved to ProcessVariable?

        MeshGeoToolsLib::MeshNodeSearcher& process_var_mesh_node_searcher =
            MeshGeoToolsLib::MeshNodeSearcher::getMeshNodeSearcher(
                _process_vars[i]->getMesh());


        DBUG("Initialize boundary conditions.");
        _process_vars[i]->initializeDirichletBCs(
                    process_var_mesh_node_searcher,
                    _dirichlet_bc.global_ids, _dirichlet_bc.values,
                    NODAL_DOF, i);


        //
        // Neumann boundary conditions.
        //
        {
            // Find mesh nodes.
            MeshGeoToolsLib::BoundaryElementsSearcher process_var_mesh_element_searcher(
                _process_vars[i]->getMesh(), process_var_mesh_node_searcher);

            // Create a neumann BC for the hydraulic head storing them in the
            // _neumann_bcs vector.
            _process_vars[i]->createNeumannBcs(
                    std::back_inserter(_neumann_bcs),
                    process_var_mesh_element_searcher,
                    _global_setup,
                    _integration_order,
                    *_local_to_global_index_map,
                    *_mesh_subset_all_nodes);
        }
    }



    for (auto bc : _neumann_bcs)
        bc->initialize(_global_setup, *_A, *_rhs, _mesh.getDimension());

}

template<typename GlobalSetup>
void
TESProcess<GlobalSetup>::
initialize()
{
    DBUG("Initialize TESProcess.");

    DBUG("Construct dof mappings.");
    // Create single component dof in every of the mesh's nodes.
    _mesh_subset_all_nodes = new MeshLib::MeshSubset(_mesh, &_mesh.getNodes());

    // Collect the mesh subsets in a vector.
    for (unsigned i=0; i<NODAL_DOF; ++i)
    {
        _all_mesh_subsets.push_back(new MeshLib::MeshSubsets(_mesh_subset_all_nodes));
    }

    _local_to_global_index_map.reset(
        new AssemblerLib::LocalToGlobalIndexMap(_all_mesh_subsets));

    DBUG("Compute sparsity pattern");
    _node_adjacency_table.createTable(_mesh.getNodes());

    DBUG("Allocate global matrix, vectors, and linear solver.");
    _A.reset(_global_setup.createMatrix(_local_to_global_index_map->dofSize()));
    _x.reset(_global_setup.createVector(_local_to_global_index_map->dofSize()));
    _rhs.reset(_global_setup.createVector(_local_to_global_index_map->dofSize()));


    for (unsigned i=0; i<NODAL_DOF_2ND; ++i)
    {
        _all_mesh_subsets_secondary.push_back(new MeshLib::MeshSubsets(_mesh_subset_all_nodes));
    }

    _local_to_global_index_map_secondary.reset(
        new AssemblerLib::LocalToGlobalIndexMap(_all_mesh_subsets_secondary));

    _secondary_variables.reset(
                _global_setup.createVector(
                    _local_to_global_index_map_secondary->dofSize()));

    if (_mesh.getDimension()==1)
        createLocalAssemblers<1>();
    else if (_mesh.getDimension()==2)
        createLocalAssemblers<2>();
    else if (_mesh.getDimension()==3)
        createLocalAssemblers<3>();
    else
        assert(false);
}

template<typename GlobalSetup>
void
TESProcess<GlobalSetup>::
solve()
{
    DBUG("Solve TESProcess.");

    _A->setZero();
    MathLib::setMatrixSparsity(*_A, _node_adjacency_table);
    *_rhs = 0;   // This resets the whole vector.

    // set initial values
    for (unsigned i=0; i<NODAL_DOF; ++i)
    {
        _process_vars[i]->setIC(*_x, NODAL_DOF, i);
    }

    // set initial values for secondary variables
    for (unsigned i=0; i<NODAL_DOF_2ND; ++i)
    {
        if (_secondary_process_vars[i])
            _secondary_process_vars[i]->setIC(*_secondary_variables, NODAL_DOF_2ND, i);
    }


    // std::cerr << __FUNCTION__ << ":" << __LINE__ << ":_x\n" << _x->getRawVector() << std::endl;

    // Call global assembler for each local assembly item.
    _global_setup.execute(*_global_assembler, _local_assemblers);

    printGlobalMatrix(_A->getRawMatrix());

    // Call global assembler for each Neumann boundary local assembler.
    for (auto bc : _neumann_bcs)
        bc->integrate(_global_setup);

    // Apply known values from the Dirichlet boundary conditions.
    MathLib::applyKnownSolution(*_A, *_rhs, _dirichlet_bc.global_ids, _dirichlet_bc.values);

    typename GlobalSetup::LinearSolver linearSolver(*_A);
    linearSolver.solve(*_rhs, *_x);
}

template<typename GlobalSetup>
void
TESProcess<GlobalSetup>::
post(std::string const& file_name)
{
    DBUG("Postprocessing TESProcess.");
    std::string const property_name = "Result";

    // Get or create a property vector for results.
    boost::optional<MeshLib::PropertyVector<double>&> result;
    if (_mesh.getProperties().hasPropertyVector(property_name))
    {
        result = _mesh.getProperties().template
            getPropertyVector<double>(property_name);
    }
    else
    {
        result = _mesh.getProperties().template
            createNewPropertyVector<double>(property_name,
                MeshLib::MeshItemType::Node);
        result->resize(_x->size());
    }
    assert(result && result->size() == _x->size());

    // Copy result
    for (std::size_t i = 0; i < _x->size(); ++i)
        (*result)[i] = (*_x)[i];

    // Write output file
    FileIO::VtuInterface vtu_interface(&_mesh, vtkXMLWriter::Binary, true);
    vtu_interface.writeToFile(file_name);
}

template<typename GlobalSetup>
TESProcess<GlobalSetup>::
~TESProcess()
{
    for (auto p : _neumann_bcs)
        delete p;

    for (auto p : _local_assemblers)
        delete p;

    for (auto p : _all_mesh_subsets)
        delete p;

    delete _mesh_subset_all_nodes;
}

} // namespace TES

} // namespace ProcessLib

#endif  // PROCESS_LIB_TESPROCESS_IMPL_H_
