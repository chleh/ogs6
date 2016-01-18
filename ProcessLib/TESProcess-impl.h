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

#include "MathLib/LinAlg/VectorNorms.h"

#include "TESProcess.h"

#include "BaseLib/Timing.h"


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
                    DBUG("proc var `%s'", v.getName().c_str());
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

            add_secondary_variable("reaction_kinetic_indicator", SecondaryVariables::REACTION_KINETIC_INDICATOR, 1);

            add_secondary_variable("vapour_partial_pressure", SecondaryVariables::VAPOUR_PARTIAL_PRESSURE, 1);
            add_secondary_variable("relative_humidity",       SecondaryVariables::RELATIVE_HUMIDITY,       1);
            add_secondary_variable("loading",                 SecondaryVariables::LOADING,                 1);
            add_secondary_variable("equilibrium_loading",     SecondaryVariables::EQUILIBRIUM_LOADING,     1);
            add_secondary_variable("reaction_damping_factor", SecondaryVariables::REACTION_DAMPING_FACTOR, 1);
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

            auto const& out_resid = config.get_optional<bool>("output.output_extrapolation_residuals");
            if (out_resid)
            {
                _output_residuals = *out_resid;
            }
        }
    }

    std::vector<std::pair<const std::string, double*> > params{
        { "fluid_specific_heat_source",            &_assembly_params._fluid_specific_heat_source },
        { "fluid_specific_isobaric_heat_capacity", &_assembly_params._cpG },
        // {  "solid_hydraulic_permeability",          &_assembly_params._solid_perm_tensor },
        { "solid_specific_heat_source",            &_assembly_params._solid_specific_heat_source },
        { "solid_heat_conductivity",               &_assembly_params._solid_heat_cond },
        { "solid_specific_isobaric_heat_capacity", &_assembly_params._cpS },
        { "tortuosity",                            &_assembly_params._tortuosity },
        { "diffusion_coefficient",                 &_assembly_params._diffusion_coefficient_component },
        { "porosity",                              &_assembly_params._poro },
        { "solid_density_dry",                     &_assembly_params._rho_SR_dry },
        { "solid_density_initial",                 &_assembly_params._initial_solid_density }
    };

    for (auto const& p : params)
    {
        auto const par = config.get_optional<double>(p.first);
        if (par) {
            DBUG("setting parameter `%s' to value `%g'", p.first.c_str(), *par);
            *p.second = *par;
        }
    }

    // permeability
    {
        auto const par = config.get_optional<double>("solid_hydraulic_permeability");
        if (par)
        {
            DBUG("setting parameter `solid_hydraulic_permeability' to isotropic value `%g'", *par);
            const auto dim = _mesh.getDimension();
            _assembly_params._solid_perm_tensor
                    = Eigen::MatrixXd::Identity(dim, dim) * (*par);
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

    // linear solver
    {
        auto const par = config.get_child_optional("linear_solver");

        if (par)
        {
            _linear_solver_options.reset(new BaseLib::ConfigTree(*par));
        }
        else
        {
            ERR("no linear solver configuration present.");
            std::abort();
        }
    }

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

    // debug output
    {
        auto param = config.get_optional<bool>("output_element_matrices");
        if (param)
        {
            DBUG("output_element_matrices: %s", (*param) ? "true" : "false");

            _assembly_params._output_element_matrices = *param;
        }

        param = config.get_optional<bool>("output_iteration_results");
        if (param)
        {
            DBUG("output_iteration_results: %s", (*param) ? "true" : "false");

            _output_iteration_results = *param;
        }

        param = config.get_optional<bool>("output_global_matrix");
        if (param)
        {
            DBUG("output_global_matrix: %s", (*param) ? "true" : "false");

            _output_global_matrix = *param;
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
                    *_local_to_global_index_map, i,
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
                    i,
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

    // TODO [CL]: Warning message
    /*
    if (LOGARITHMIC_VAPOUR_MASS_FRACTION)
        INFO("Vapour mass fraction is taken logarithmically!"
             " Consider that for initial and boundary conditions as well as for output.");
             */

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
    _sparsity_pattern = std::move(AssemblerLib::computeSparsityPattern(
                *_local_to_global_index_map, _mesh));

    DBUG("Allocate global matrix, vectors, and linear solver.");
    _A.reset(_global_setup.createMatrix(_local_to_global_index_map->dofSize()));
    _x.reset(_global_setup.createVector(_local_to_global_index_map->dofSize()));
    _rhs.reset(_global_setup.createVector(_local_to_global_index_map->dofSize()));

    _x_prev_ts.reset(_global_setup.createVector(_local_to_global_index_map->dofSize()));

    // for extrapolation of secondary variables
    _all_mesh_subsets_single_component.push_back(new MeshLib::MeshSubsets(_mesh_subset_all_nodes));
    _local_to_global_index_map_single_component.reset(
                new AssemblerLib::LocalToGlobalIndexMap(_all_mesh_subsets_single_component, _global_matrix_order)
                );

    _linear_solver.reset(new typename GlobalSetup::LinearSolver(
        *_A, "", _linear_solver_options.get()));
    _extrapolator.reset(new ExtrapolatorImpl(*_local_to_global_index_map_single_component));

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
        setInitialConditions(*_process_vars[i], i);
    }

    /*
    std::puts("------ initial values ----------");
    printGlobalVector(_x->getRawVector());
    */

    _picard.reset(new MathLib::Nonlinear::Picard);
    _picard->setAbsTolerance(1e-1);
    _picard->setRelTolerance(1e-6);
    _picard->setMaxIterations(100);
    _picard->printErrors(true);
}

template<typename GlobalSetup>
void TESProcess<GlobalSetup>::
setInitialConditions(ProcessVariable const& variable, std::size_t component_id)
{
    std::size_t const n = _mesh.getNNodes();
    for (std::size_t i = 0; i < n; ++i)
    {
        MeshLib::Location const l(_mesh.getID(),
                                  MeshLib::MeshItemType::Node, i);
        std::size_t const global_index =
            _local_to_global_index_map->getGlobalIndex(
                l, component_id);
        _x->set(global_index,
               variable.getInitialConditionValue(*_mesh.getNode(i)));
    }
}

template<typename GlobalSetup>
bool TESProcess<GlobalSetup>::solve(const double delta_t)
{
    DBUG("Solve TESProcess.");

#if 0
    auto tmp = *_x;
    if (false && _timestep != 0)
    {
        // this probably cannot be applied with the current reaction scheme!
        // the reaction rate is extrapolated from the solution of the last timestep
        // so this solution should not be extrapolated itself here.
        // the proper treatment would be to use previos timestep values in the local assembler.

        // this is at the beginning of a new timestep t+dt.
        // _x         contains the solution of timstep t
        // _x_prev_ts contains the solution of timstep t-dt
        // get first estimate of _x at timestep t+dt by extrapolation using the last two timesteps.
        *_x *= 1.0 + delta_t / _assembly_params._delta_t; // _assembly_params._delta_t is the old timestep
        *_x_prev_ts *= delta_t / _assembly_params._delta_t;
        *_x -= *_x_prev_ts;
    }
    *_x_prev_ts = tmp;
#else
    *_x_prev_ts = *_x;
#endif

    _assembly_params._delta_t = delta_t;
    _assembly_params._iteration_in_current_timestep = 0;
    ++ _timestep;

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
postTimestep(const std::string& file_name, const unsigned /*timestep*/)
// TODO [CL] remove second parameter
{
    INFO("postprocessing timestep");

    /*
    std::puts("---- solution ----");
    printGlobalVector(_x->getRawVector());
    */

    auto count = [](MeshLib::Mesh const& mesh, MeshLib::MeshItemType type) -> std::size_t
    {
        switch (type) {
        case MeshLib::MeshItemType::Cell: return mesh.getNElements();
        case MeshLib::MeshItemType::Node: return mesh.getNNodes();
        default: break;
        }
        return 0;
    };

    auto get_or_create_mesh_property = [this, &count](std::string const& property_name, MeshLib::MeshItemType type)
    {
        // Get or create a property vector for results.
        boost::optional<MeshLib::PropertyVector<double>&> result;

        auto const N = count(_mesh, type);

        if (_mesh.getProperties().hasPropertyVector(property_name))
        {
            result = _mesh.getProperties().template
                getPropertyVector<double>(property_name);
        }
        else
        {
            result = _mesh.getProperties().template
                createNewPropertyVector<double>(property_name, type);
            result->resize(N);
        }
        assert(result && result->size() == N);

        return result;
    };

    auto add_primary_var = [this, &get_or_create_mesh_property](const unsigned vi)
    {
        std::string const& property_name = _process_vars[vi]->getName();
        if (_output_variables.find(property_name) == _output_variables.cend())
            return;

        DBUG("  process var %s", property_name.c_str());

        auto result = get_or_create_mesh_property(property_name, MeshLib::MeshItemType::Node);
        assert(result->size() == _mesh.getNNodes());

        // Copy result
        for (std::size_t i = 0; i < _mesh.getNNodes(); ++i)
        {
            MeshLib::Location loc(_mesh.getID(), MeshLib::MeshItemType::Node, i);
            auto const idx = _local_to_global_index_map->getGlobalIndex(loc, vi);
            (*result)[i] = (*_x)[idx];
        }
    };

    assert(_x->size() == NODAL_DOF * _mesh.getNNodes());
    for (unsigned vi=0; vi!=NODAL_DOF; ++vi)
    {
        add_primary_var(vi);
    }


    auto add_secondary_var = [this, &get_or_create_mesh_property]
                             (SecondaryVariables const property,
                             std::string const& property_name,
                         #ifndef NDEBUG
                             const unsigned num_components
                         #else
                             const unsigned /*num_components*/
                         #endif
                             )
    {
        assert(num_components == 1); // TODO [CL] implement other cases

        {
            if (_output_variables.find(property_name) == _output_variables.cend())
                return;

            DBUG("  process var %s", property_name.c_str());

            auto result = get_or_create_mesh_property(property_name, MeshLib::MeshItemType::Node);
            assert(result->size() == _mesh.getNNodes());

            _extrapolator->extrapolate(*_x, *_local_to_global_index_map, _local_assemblers, property);
            auto const& nodal_values = _extrapolator->getNodalValues();

            // Copy result
            for (std::size_t i = 0; i < _mesh.getNNodes(); ++i)
            {
                (*result)[i] = nodal_values[i];
            }
        }

        if (_output_residuals) {
            DBUG("  process var %s residual", property_name.c_str());
            auto const& property_name_res = property_name + "_residual";

            auto result = get_or_create_mesh_property(property_name_res, MeshLib::MeshItemType::Cell);
            assert(result->size() == _mesh.getNElements());

            _extrapolator->calculateResiduals(*_x, *_local_to_global_index_map, _local_assemblers, property);
            auto const& residuals = _extrapolator->getElementResiduals();

            // Copy result
            for (std::size_t i = 0; i < _mesh.getNElements(); ++i)
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
singlePicardIteration(GlobalVector& x_prev_iter,
                      GlobalVector& x_curr)
{
    bool iteration_accepted = false;
    unsigned num_try = 0;

    do
    {
        INFO("-> TES process try number %u in current picard iteration", num_try);
        _assembly_params._number_of_try_of_iteration = num_try;

        _global_assembler->setX(&x_curr, _x_prev_ts.get());

        _A->setZero();
        MathLib::setMatrixSparsity(*_A, _sparsity_pattern);
        *_rhs = 0;   // This resets the whole vector.

        // Call global assembler for each local assembly item.
        _global_setup.execute(*_global_assembler, _local_assemblers);

        {
        BaseLib::TimingOneShot timing{"assembly"};
        // Call global assembler for each Neumann boundary local assembler.
        for (auto bc : _neumann_bcs)
            bc->integrate(_global_setup, &x_curr, _x_prev_ts.get());
        timing.stop();
        }

        // Apply known values from the Dirichlet boundary conditions.

        INFO("size of known values: %li", _dirichlet_bc.global_ids.size());

        {
        BaseLib::TimingOneShot timing{"apply known solutions"};
        MathLib::applyKnownSolution(*_A, *_rhs, _dirichlet_bc.global_ids, _dirichlet_bc.values);
        timing.stop();
        }

#if !defined(USE_LIS)
        // double residual = MathLib::norm((*_A) * x_curr - (*_rhs), MathLib::VecNormType::INFINITY_N);
        GlobalVector res_vec;
        _A->multiply(x_curr, res_vec);
        res_vec -= *_rhs;

        double residual = MathLib::norm(res_vec, MathLib::VecNormType::INFINITY_N);
        DBUG("residual of old solution with new matrix: %g", residual);
#endif

        MathLib::scaleDiagonal(*_A, *_rhs);

#ifndef USE_LIS
        // _A->getRawMatrix().rowwise() /= diag; //  = invDiag * _A->getRawMatrix();
        // .cwiseQuotient(diag); //  = invDiag * _rhs->getRawVector();

        _A->multiply(x_curr, res_vec);
        res_vec -= *_rhs;
        residual = MathLib::norm(res_vec, MathLib::VecNormType::INFINITY_N);
        DBUG("residual of new solution with new matrix: %g", residual);
#endif

#ifndef NDEBUG
        if (_total_iteration == 0 && num_try == 0 && _output_global_matrix)
        {
#ifdef USE_LIS
        MathLib::finalizeMatrixAssembly(*_A);
#endif
            // TODO [CL] Those files will be written to the working directory.
            //           Relative path needed.
            _A->write("global_matrix.txt");
            _rhs->write("global_rhs.txt");
        }
#endif

        {
        BaseLib::TimingOneShot timing{"linear solver"};
        _linear_solver->solve(*_rhs, x_curr);
        timing.stop();
        }

#ifndef NDEBUG
        if (_total_iteration == 0 && num_try == 0 && _output_global_matrix)
        {
            // TODO [CL] Those files will be written to the working directory.
            //           Relative path needed.
            _A->write("global_matrix_post.txt");
            _rhs->write("global_rhs_post.txt");
        }
#endif

#if 0
        // assert that no component is smaller than some lower bound

        const std::array<double, 3> bounds = { 1.0, 100.0, 1e-6 }; // lower bounds for p, T, x

        auto alpha = [](double xnew, double xold, double xmin) // computes the damping coefficient
        {
            // assert (xold >= xmin);
            return std::max(0.0, (xold - xmin) / (xold - xnew));
        };

        std::array<double, 3> damping_coeffs = { 1.0, 1.0, 1.0 };
        bool do_damping = false;
        const double pre_damp = 0.75;
        const double min_damp = 0.1;

        switch(_global_matrix_order)
        {
        case AssemblerLib::ComponentOrder::BY_COMPONENT:
        {
            for (std::size_t i=0; i<x_curr.size(); i+=3)
            {
                for (std::size_t d=0; d<NODAL_DOF; ++d)
                {
                    if (x_curr[i+d] < bounds[d]) {
                        damping_coeffs[d] = std::min(damping_coeffs[d], alpha(x_curr[i+d], x_prev_iter[i+d], bounds[d]));
                        do_damping = true;
                    }
                }
            }

            if (do_damping)
            {
                auto const damping_coeff =
                        std::max(min_damp, pre_damp * *std::min_element(damping_coeffs.cbegin(), damping_coeffs.cend()));
                DBUG("doing damping with coeff %g", damping_coeff);

                for (unsigned i=0; i<x_curr.size(); ++i)
                {
                    x_curr[i] = (1.0-damping_coeff) * x_prev_iter[i] + damping_coeff * x_curr[i];
                }
            }
            break;
        }
        case AssemblerLib::ComponentOrder::BY_LOCATION:
        {
            for (std::size_t d=0; d<NODAL_DOF; ++d)
            {
                for (std::size_t i=0; i<x_curr.size(); i+=3)
                {
                    if (x_curr[i+d] < bounds[d]) {
                        damping_coeffs[d] = std::min(damping_coeffs[d], alpha(x_curr[i+d], x_prev_iter[i+d], bounds[d]));
                        do_damping = true;
                    }
                }
            }

            if (do_damping)
            {
                auto const damping_coeff =
                        std::max(min_damp, pre_damp * *std::min_element(damping_coeffs.cbegin(), damping_coeffs.cend()));
                DBUG("doing damping with coeff %g", damping_coeff);

                for (unsigned i=0; i<x_curr.size(); ++i)
                {
                    x_curr[i] = (1.0-damping_coeff) * x_prev_iter[i] + damping_coeff * x_curr[i];
                }
            }
            break;
        }
        }
#endif

        if (_output_iteration_results)
        {
            DBUG("output results of iteration %li", _total_iteration);
            std::string fn = "tes_iter_" + std::to_string(_total_iteration) +
                             + "_ts_" + std::to_string(_timestep)
                             + "_" +    std::to_string(_assembly_params._iteration_in_current_timestep)
                             + "_" +    std::to_string(num_try)
                             + ".vtu";
            postTimestep(fn, 0);
        }

        bool check_passed = true;

        if (!Trafo::constrained)
        {
            // bounds checking only has to happen if the vapour mass fraction is non-logarithmic.

            auto& ga = *_global_assembler;

            auto check_variable_bounds
            = [&ga, &check_passed](
              std::size_t id, LocalAssembler* const loc_asm)
            {
                // DBUG("%lu", id);

                std::vector<double> const* localX;
                std::vector<double> const* localX_pts;
                ga.getLocalNodalValues(id, localX, localX_pts);

                if (!loc_asm->checkBounds(*localX, *localX_pts)) check_passed = false;
            };

            _global_setup.execute(check_variable_bounds, _local_assemblers);

            if (!check_passed)
            {
                x_curr = x_prev_iter;
            }
        }

        iteration_accepted = check_passed;

        ++num_try;
    }
    while(! iteration_accepted);

    DBUG("ts %lu iteration %lu (%lu) try %u accepted", _timestep, _total_iteration,
         _assembly_params._iteration_in_current_timestep, num_try-1);

    ++ _assembly_params._iteration_in_current_timestep;
    ++_total_iteration;
}


} // namespace TES

} // namespace ProcessLib

#endif  // PROCESS_LIB_TESPROCESS_IMPL_H_
