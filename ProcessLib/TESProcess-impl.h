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
#include <cstdio>

#include "AssemblerLib/LocalToGlobalIndexMap.h"

#include "logog/include/logog.hpp"

#include "NumLib/TimeStepping/Algorithms/FixedTimeStepping.h"

#include "LocalLinearLeastSquaresExtrapolator.h"
#include "GlobalLinearLeastSquaresExtrapolator.h"

#include "TESProcess.h"


namespace
{

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

} // anonymous namespace


#if 0
template<typename Conf, typename InT, typename OutT>
void parseParameter(Conf const& config,
                    std::string const& param_name, OutT* (*builder)(InT const&), OutT*& target)
{
    InT in_value = config.template get<InT>(param_name);
    // DBUG("reactive_system: %s", in_value.c_str());

    target = builder(in_value);
    // _assembly_params._adsorption = Ads::Adsorption::newInstance(rsys);
}
#endif


namespace ProcessLib
{

namespace TES
{

template<typename GlobalSetup>
TESProcess<GlobalSetup>::
TESProcess(MeshLib::Mesh& mesh,
           std::vector<ProcessVariable> const& variables,
           std::vector<std::unique_ptr<ParameterBase>> const& /*parameters*/,
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
        auto const& proc_vars = config.get_child_optional("secondary_variables");
        if (proc_vars)
        {
            auto add_secondary_variable =
                    [this, &proc_vars](
                    std::string const& var, SecondaryVariables type, unsigned num_components)
            {
                auto variable = proc_vars->get_optional<std::string>(var);
                if (variable)
                {
                    _secondary_process_vars.emplace_back(type, *variable, num_components);
                }
            };

            add_secondary_variable("solid_density", SecondaryVariables::SOLID_DENSITY, 1);
            add_secondary_variable("reaction_rate", SecondaryVariables::REACTION_RATE, 1);
            add_secondary_variable("velocity_x",    SecondaryVariables::VELOCITY_X,    1);
            if (_mesh.getDimension() >= 2) add_secondary_variable("velocity_y",    SecondaryVariables::VELOCITY_Y,    1);
            if (_mesh.getDimension() >= 3) add_secondary_variable("velocity_z",    SecondaryVariables::VELOCITY_Z,    1);
        }
    }

    // variables for output
    {
        auto const& out_vars = config.get_child_optional("output.variables");
        if (out_vars)
        {
            auto const& out_vars_range = out_vars->equal_range("variable");
            for (auto it = out_vars_range.first; it!=out_vars_range.second; ++it)
            {
                // auto const& out_var = it->first; //->second.get<std::string>("variable");
                auto const& out_var = it->second.get_value<std::string>();

                if (_output_variables.find(out_var) != _output_variables.cend())
                {
                    ERR("output variable `%s' specified twice.", out_var.c_str());
                    std::abort();
                }

                auto pred = [&out_var](ProcessVariable const* pv) {
                    return pv->getName() == out_var;
                };

                // check if process variable
                auto const& pcs_var = std::find_if(
                    _process_vars.cbegin(), _process_vars.cend(),
                    pred);

                if (pcs_var == _process_vars.cend())
                {
                    auto pred2 = [&out_var](std::tuple<SecondaryVariables, std::string, unsigned> const& p) {
                        return std::get<1>(p) == out_var;
                    };

                    // check if secondary variable
                    auto const& pcs_var2 = std::find_if(
                        _secondary_process_vars.cbegin(), _secondary_process_vars.cend(),
                        pred2);

                    if (pcs_var2 == _secondary_process_vars.cend())
                    {
                        ERR("Output variable `%s' is neither a process variable nor a"
                            " secondary variable", out_var.c_str());
                        std::abort();
                    }
                }

                DBUG("adding output variable `%s'", out_var.c_str());
                _output_variables.insert(out_var);
            }
        }
    }

#if 1
    // reactive system
    {
        auto rsys = config.get<std::string>("reactive_system");
        DBUG("reactive_system: %s", rsys.c_str());

        _assembly_params._adsorption = Ads::Adsorption::newInstance(rsys);
    }
#else
    parseParameter(config, "reactive_system", Ads::Adsorption::newInstance, _assembly_params._adsorption);
#endif

    // matrix order
    {
        auto order = config.get<std::string>("global_matrix_order");
        DBUG("global_matrix_order: %s", order.c_str());

        if (order == "BY_COMPONENT")
            _global_matrix_order = AssemblerLib::ComponentOrder::BY_COMPONENT;
        else if (order == "BY_LOCATION")
            _global_matrix_order = AssemblerLib::ComponentOrder::BY_LOCATION;
        else {
            ERR("unknown global matrix order `%s'", order.c_str());
            std::abort();
        }
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

    for (unsigned i=0; i<NODAL_DOF; ++i)
    {
        // TODO [CL] can't that be (partially) moved to ProcessVariable?

        MeshGeoToolsLib::MeshNodeSearcher& process_var_mesh_node_searcher =
            MeshGeoToolsLib::MeshNodeSearcher::getMeshNodeSearcher(
                _process_vars[i]->getMesh());


        DBUG("Initialize boundary conditions.");
        _process_vars[i]->initializeDirichletBCs(
                    process_var_mesh_node_searcher,
                    _local_to_global_index_map->getMeshComponentMap(), i,
                    _dirichlet_bc.global_ids, _dirichlet_bc.values);


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
        new AssemblerLib::LocalToGlobalIndexMap(_all_mesh_subsets, _global_matrix_order));

    DBUG("Compute sparsity pattern");
    _node_adjacency_table.createTable(_mesh.getNodes());

    DBUG("Allocate global matrix, vectors, and linear solver.");
    _A.reset(_global_setup.createMatrix(_local_to_global_index_map->dofSize()));
    _x.reset(_global_setup.createVector(_local_to_global_index_map->dofSize()));
    _rhs.reset(_global_setup.createVector(_local_to_global_index_map->dofSize()));

    _x_prev_ts.reset(_global_setup.createVector(_local_to_global_index_map->dofSize()));


    // for extrapolation of secondary variables
    _all_mesh_subsets_single_component.push_back(new MeshLib::MeshSubsets(_mesh_subset_all_nodes));
    _local_to_global_index_map_single_component.reset(
                new AssemblerLib::LocalToGlobalIndexMap(_all_mesh_subsets_single_component));

    if (_mesh.getDimension()==1)
        createLocalAssemblers<1>();
    else if (_mesh.getDimension()==2)
        createLocalAssemblers<2>();
    else if (_mesh.getDimension()==3)
        createLocalAssemblers<3>();
    else
        assert(false);


    // set initial values
    for (unsigned i=0; i<NODAL_DOF; ++i)
    {
        _process_vars[i]->setIC(*_x, _local_to_global_index_map->getMeshComponentMap(), i);
    }

    std::puts("------ initial values ----------");
    printGlobalVector(_x->getRawVector());

    _picard.reset(new MathLib::Nonlinear::Picard);
    _picard->setAbsTolerance(1e-1);
    _picard->setRelTolerance(1e-6);
    _picard->setMaxIterations(100);
}

template<typename GlobalSetup>
bool TESProcess<GlobalSetup>::solve(const double delta_t)
{
    DBUG("Solve TESProcess.");

    _assembly_params._time_step = delta_t;
    _assembly_params._is_new_timestep = true;

    *_x_prev_ts = *_x;

    auto cb = [this](typename GlobalSetup::VectorType& x_prev_iter,
                             typename GlobalSetup::VectorType& x_curr)
    {
        singlePicardIteration(x_prev_iter, x_curr);
    };

    return _picard->solve(cb, *_x_prev_ts, *_x);
}


template<typename GlobalSetup>
void
TESProcess<GlobalSetup>::
postTimestep(const std::string& file_name, const unsigned timestep)
// TODO [CL] remove second parameter
{
    INFO("postprocessing timestep %i", timestep);

    std::puts("---- solution ----");
    printGlobalVector(_x->getRawVector());

    auto add_primary_var = [this](const unsigned vi)
    {
        std::string const& property_name = _process_vars[vi]->getName();
        if (_output_variables.find(property_name) == _output_variables.cend())
            return;

        DBUG("  process var %s", property_name.c_str());

        // Get or create a property vector for results.
        boost::optional<MeshLib::PropertyVector<double>&> result;

        assert(_x->size() == NODAL_DOF * _mesh.getNNodes());
        auto const N = _mesh.getNNodes();

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
            result->resize(N);
        }
        assert(result && result->size() == N);

        auto const& mcmap = _local_to_global_index_map->getMeshComponentMap();
        // Copy result
        for (std::size_t i = 0; i < N; ++i)
        {
            MeshLib::Location loc(_mesh.getID(), MeshLib::MeshItemType::Node, i);
            auto const idx = mcmap.getGlobalIndex(loc, vi);
            (*result)[i] = (*_x)[idx];
        }
    };

    for (unsigned vi=0; vi!=NODAL_DOF; ++vi)
    {
        add_primary_var(vi);
    }


#if 0
    GlobalLinearLeastSquaresExtrapolator<
            typename GlobalSetup::MatrixType,
            typename GlobalSetup::VectorType, SecondaryVariables,
            LocalAssembler>
            extrapolator(*_local_to_global_index_map_single_component);
#else
    LocalLinearLeastSquaresExtrapolator<typename GlobalSetup::VectorType, SecondaryVariables,
            LocalAssembler>
            extrapolator(*_local_to_global_index_map_single_component);
#endif

    auto add_secondary_var = [this, &extrapolator]
                             (SecondaryVariables const property,
                             std::string const& property_name,
                             const unsigned num_components
                             )
    {
        {
            if (_output_variables.find(property_name) == _output_variables.cend())
                return;

            DBUG("  process var %s", property_name.c_str());

            // Get or create a property vector for results.
            boost::optional<MeshLib::PropertyVector<double>&> result;

            auto const N = _mesh.getNNodes();

            if (_mesh.getProperties().hasPropertyVector(property_name))
            {
                result = _mesh.getProperties().template
                    getPropertyVector<double>(property_name);
            }
            else
            {
                result = _mesh.getProperties().template
                    createNewPropertyVector<double>(property_name,
                        MeshLib::MeshItemType::Node,
                        num_components);
                result->resize(N);
            }
            assert(result && result->size() == N);


            extrapolator.extrapolate(_local_assemblers, property);
            auto const& nodal_values = extrapolator.getNodalValues().getRawVector();

            // Copy result
            for (std::size_t i = 0; i < N; ++i)
            {
                (*result)[i] = nodal_values[i];
            }
        }

        if (_output_residuals) {
            DBUG("  process var %s residual", property_name.c_str());
            auto const& property_name_res = property_name + "_residual";

            // Get or create a property vector for results.
            boost::optional<MeshLib::PropertyVector<double>&> result;

            auto const num_cells = _mesh.getNElements();

            if (_mesh.getProperties().hasPropertyVector(property_name_res))
            {
                result = _mesh.getProperties().template
                    getPropertyVector<double>(property_name_res);
            }
            else
            {
                result = _mesh.getProperties().template
                    createNewPropertyVector<double>(property_name_res,
                        MeshLib::MeshItemType::Cell);
                result->resize(num_cells);
            }
            assert(result && result->size() == num_cells);

            extrapolator.calculateResiduals(_local_assemblers, property);
            auto const& residuals = extrapolator.getElementResiduals();

            // Copy result
            for (std::size_t i = 0; i < num_cells; ++i)
            {
                (*result)[i] = residuals[i];
            }
        }
    };

    for (auto const& p : _secondary_process_vars)
    {
        add_secondary_var(std::get<0>(p), std::get<1>(p), std::get<2>(p));
    }


    // Write output file
    FileIO::VtuInterface vtu_interface(&_mesh, vtkXMLWriter::Binary, true);
    vtu_interface.writeToFile(file_name);
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



template<typename GlobalSetup>
void
TESProcess<GlobalSetup>::
singlePicardIteration(typename GlobalSetup::VectorType& /*x_prev_iter*/,
                      typename GlobalSetup::VectorType& x_curr)
{
    _global_assembler->setX(&x_curr, _x_prev_ts.get());

    _A->setZero();
    MathLib::setMatrixSparsity(*_A, _node_adjacency_table); // TODO [CL] always required?
    *_rhs = 0;   // This resets the whole vector.

    // Call global assembler for each local assembly item.
    _global_setup.execute(*_global_assembler, _local_assemblers);

    // Call global assembler for each Neumann boundary local assembler.
    for (auto bc : _neumann_bcs)
        bc->integrate(_global_setup);

    // Apply known values from the Dirichlet boundary conditions.
    MathLib::applyKnownSolution(*_A, *_rhs, _dirichlet_bc.global_ids, _dirichlet_bc.values);

    typename GlobalSetup::LinearSolver linearSolver(*_A);
    linearSolver.solve(*_rhs, x_curr);

    _assembly_params._is_new_timestep = false;
}


} // namespace TES

} // namespace ProcessLib

#endif  // PROCESS_LIB_TESPROCESS_IMPL_H_
