/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TESProcess.h"

#include "BaseLib/Functional.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "TESOGS5MaterialModels.h"

namespace ProcessLib
{
namespace TES
{
TESProcess::TESProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller,
    const BaseLib::ConfigTree& config)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller))
{
    DBUG("Create TESProcess.");

    // physical parameters for local assembly
    {
        std::vector<std::pair<std::string, double*>> params{
            //! \ogs_file_param_special{prj__processes__process__TES__fluid_specific_heat_source}
            {"fluid_specific_heat_source",
             &_assembly_params.fluid_specific_heat_source},
            //! \ogs_file_param_special{prj__processes__process__TES__fluid_specific_isobaric_heat_capacity}
            {"fluid_specific_isobaric_heat_capacity", &_assembly_params.cpG},
            //! \ogs_file_param_special{prj__processes__process__TES__solid_specific_heat_source}
            {"solid_specific_heat_source",
             &_assembly_params.solid_specific_heat_source},
            //! \ogs_file_param_special{prj__processes__process__TES__solid_heat_conductivity}
            {"solid_heat_conductivity", &_assembly_params.solid_heat_cond},
            //! \ogs_file_param_special{prj__processes__process__TES__solid_specific_isobaric_heat_capacity}
            {"solid_specific_isobaric_heat_capacity", &_assembly_params.cpS},
            //! \ogs_file_param_special{prj__processes__process__TES__tortuosity}
            {"tortuosity", &_assembly_params.tortuosity},
            //! \ogs_file_param_special{prj__processes__process__TES__solid_density_dry}
            {"solid_density_dry", &_assembly_params.rho_SR_dry},
            //! \ogs_file_param_special{prj__processes__process__TES__solid_density_lower_limit}
            {"solid_density_lower_limit", &_assembly_params.rho_SR_lower},
            //! \ogs_file_param_special{prj__processes__process__TES__solid_density_upper_limit}
            {"solid_density_upper_limit", &_assembly_params.rho_SR_upper},
            //! \ogs_file_param_special{prj__processes__process__TES__solid_density_initial}
            {"solid_density_initial", &_assembly_params.initial_solid_density}};

        for (auto const& p : params)
        {
            auto const par =
                //! \ogs_file_special
                config.getConfigParameter<double>(p.first);
            DBUG("setting parameter `%s' to value `%g'", p.first.c_str(), par);
            *p.second = par;
        }
    }

    //! \ogs_file_param_special{process__TES__porosity}
    _assembly_params.poro =
        &findParameter<double>(config, "porosity", parameters, 1);

    //! \ogs_file_param{process__TES__permeability}
    _assembly_params.permeability =
        &findParameter<double>(config, "permeability", parameters, 1);

    _assembly_params.dielectric_heating_term_enabled =
        config.getConfigParameter<bool>("dielectric_heating_term_enabled");

    if (auto const heat_loss =
            config.getConfigSubtreeOptional("volumetric_heat_loss"))
    {
        _assembly_params.volumetric_heat_loss =
            createVolumetricHeatLoss(*heat_loss);
    }

    _assembly_params.diffusion_coefficient_component =
        createDiffusionCoefficient(
            //! \ogs_file_param_special{process__TES__diffusion_coefficient}
            config.getConfigSubtree("diffusion_coefficient"));

    if (_assembly_params.dielectric_heating_term_enabled)
    {
        auto const hps = config.getConfigSubtree("heating_power_scaling");
        _assembly_params.heating_power_scaling =
            MathLib::PiecewiseLinearInterpolation(
                hps.getConfigParameter<std::vector<double>>("times"),
                hps.getConfigParameter<std::vector<double>>("scalings"), false);
    }
    else
    {
        config.ignoreConfigParameter("heating_power_scaling");
    }

    if (auto prop = config.getConfigParameterOptional<std::string>(
            "initial_solid_density_mesh_property"))
    {
        assert(!prop->empty());
        _assembly_params.initial_solid_density_mesh_property = *prop;
    }

    _assembly_params.reaction_rate = MaterialLib::createReactionRate(
        //! \ogs_file_param{process__TES__reaction_rate}
        config.getConfigSubtree("reaction_rate"));

    _assembly_params.reactive_solid = MaterialLib::createReactiveSolidModel(
        //! \ogs_file_param{process__TES__reactive_solid}
        config.getConfigSubtree("reactive_solid"), parameters);

    // debug output
    if (auto const param =
            //! \ogs_file_param{prj__processes__process__TES__output_element_matrices}
            config.getConfigParameterOptional<bool>("output_element_matrices"))
    {
        DBUG("output_element_matrices: %s", (*param) ? "true" : "false");

        _assembly_params.output_element_matrices = *param;
    }

    // TODO somewhere else
    /*
    if (auto const param =
    //! \ogs_file_param{prj__processes__process__TES__output_global_matrix}
    config.getConfigParameterOptional<bool>("output_global_matrix"))
    {
        DBUG("output_global_matrix: %s", (*param) ? "true" : "false");

        this->_process_output.output_global_matrix = *param;
    }
    */
}

void TESProcess::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh, unsigned const integration_order)
{
    _assembly_params.dof_table = _local_to_global_index_map.get();

    ProcessLib::ProcessVariable const& pv = getProcessVariables(0)[0];
    ProcessLib::createLocalAssemblers<TESLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        pv.getShapeFunctionOrder(), _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _assembly_params);

    initializeSecondaryVariables();

    // set initial solid density from mesh property
    if (!_assembly_params.initial_solid_density_mesh_property.empty())
    {
        auto prop = mesh.getProperties().template getPropertyVector<double>(
            _assembly_params.initial_solid_density_mesh_property);
        if (!prop)
            OGS_FATAL(
                "property `%s' not found.",
                _assembly_params.initial_solid_density_mesh_property.c_str());
        if (prop->getNumberOfComponents() != 1)
            OGS_FATAL(
                "property `%s' not found has %d components instead of the "
                "expected value of one.",
                _assembly_params.initial_solid_density_mesh_property.c_str(),
                prop->getNumberOfComponents());

        switch (prop->getMeshItemType())
        {
            case MeshLib::MeshItemType::Cell:
            {
                auto init_solid_density =
                    [&prop](std::size_t id,
                            TESLocalAssemblerInterface& loc_asm) {
                        // TODO loc_asm_id is assumed to be the mesh element id.
                        loc_asm.initializeSolidDensity(
                            MeshLib::MeshItemType::Cell, {{(*prop)[id]}});
                    };

                NumLib::SerialExecutor::executeDereferenced(init_solid_density,
                                                            _local_assemblers);
                break;
            }
            case MeshLib::MeshItemType::Node:
            {
                std::vector<GlobalIndexType> indices;
                std::vector<double> values;

                auto init_solid_density =
                    [&](std::size_t id, TESLocalAssemblerInterface& loc_asm) {
                        NumLib::getRowColumnIndices(
                            id, getSingleComponentDOFTable(), indices);
                        values.clear();
                        for (auto i : indices)
                            values.push_back((*prop)[i]);

                        loc_asm.initializeSolidDensity(
                            MeshLib::MeshItemType::Node, values);
                    };

                NumLib::SerialExecutor::executeDereferenced(init_solid_density,
                                                            _local_assemblers);
                break;
            }
            default:
                ERR("Unhandled mesh item type for initialization of secondary "
                    "variable.");
                std::abort();
        }
    }
}

void TESProcess::initializeSecondaryVariables()
{
    // adds a secondary variables to the collection of all secondary variables.
    auto add2nd = [&](std::string const& var_name,
                      SecondaryVariableFunctions&& fcts) {
        _secondary_variables.addSecondaryVariable(var_name, std::move(fcts));
    };

    // named functions: vapour partial pressure ////////////////////////////////
    auto p_V_fct = [=](const double p, const double x_mV) {
        const double x_nV = Adsorption::AdsorptionReaction::getMolarFraction(
            x_mV, _assembly_params.M_react, _assembly_params.M_inert);
        return p * x_nV;
    };
    _named_function_caller.addNamedFunction(
        {"vapour_partial_pressure",
         {"pressure", "vapour_mass_fraction"},
         BaseLib::easyBind(std::move(p_V_fct))});
    _named_function_caller.plug("vapour_partial_pressure", "pressure",
                                "TES_pressure");
    _named_function_caller.plug("vapour_partial_pressure",
                                "vapour_mass_fraction",
                                "TES_vapour_mass_fraction");
    // /////////////////////////////////////////////////////////////////////////

    // named functions: solid density //////////////////////////////////////////
    auto solid_density = std::make_unique<CachedSecondaryVariable>(
        "TES_solid_density", getExtrapolator(), _local_assemblers,
        &TESLocalAssemblerInterface::getIntPtSolidDensity,
        _secondary_variable_context);

    for (auto&& fct : solid_density->getNamedFunctions())
        _named_function_caller.addNamedFunction(std::move(fct));

    add2nd("solid_density", solid_density->getExtrapolator());

    _cached_secondary_variables.emplace_back(std::move(solid_density));
    // /////////////////////////////////////////////////////////////////////////

    // named functions: fluid density //////////////////////////////////////////
    auto rho_GR_fct = [](const double p, const double T, const double x_mV) {
        return fluid_density(p, T, x_mV);
    };
    _named_function_caller.addNamedFunction(
        {"fluid_density",
         {"pressure", "temperature", "vapour_mass_fraction"},
         BaseLib::easyBind(std::move(rho_GR_fct))});
    _named_function_caller.plug("fluid_density", "pressure", "TES_pressure");
    _named_function_caller.plug("fluid_density", "temperature",
                                "TES_temperature");
    _named_function_caller.plug("fluid_density", "vapour_mass_fraction",
                                "TES_vapour_mass_fraction");
    // //////////////////////////////////////////////////////////////////////////

    // named functions: diffusion coefficient
    // //////////////////////////////////////////
    auto D_fct = [this](const double p, const double T, const double p_V) {
        return _assembly_params.diffusion_coefficient_component
            ->getDiffusionCoefficient(p, T, p_V);
    };
    _named_function_caller.addNamedFunction(
        {"diffusion_coefficient",
         {"pressure", "temperature", "vapour_partial_pressure"},
         BaseLib::easyBind(std::move(D_fct))});
    _named_function_caller.plug("diffusion_coefficient", "pressure",
                                "TES_pressure");
    _named_function_caller.plug("diffusion_coefficient", "temperature",
                                "TES_temperature");
    _named_function_caller.plug("diffusion_coefficient",
                                "vapour_partial_pressure",
                                "vapour_partial_pressure");
    // //////////////////////////////////////////////////////////////////////////

    // named functions: volumetric heating power
    // ////////////////////////////////
    if (_assembly_params.dielectric_heating_term_enabled)
    {
        auto heat_power_fct = [this](const double T,
                                     const double rho_SR) -> double {
            auto const loading = Adsorption::AdsorptionReaction::getLoading(
                rho_SR, _assembly_params.rho_SR_dry);
            auto const t = _assembly_params.current_time;
            // TODO check if heating_power_scaling is provided.
            return _assembly_params.heating_power_scaling.getValue(t) *
                   getVolumetricJouleHeatingPower(T, loading);
        };
        _named_function_caller.addNamedFunction(
            {"volumetric_joule_heating_power",
             {"temperature", "solid_density"},
             BaseLib::easyBind(std::move(heat_power_fct))});
        _named_function_caller.plug("volumetric_joule_heating_power",
                                    "temperature", "TES_temperature");
        _named_function_caller.plug("volumetric_joule_heating_power",
                                    "solid_density", "TES_solid_density");
    }
    // /////////////////////////////////////////////////////////////////////////

    // named functions: from kinetics
    for (auto&& fct : _assembly_params.reactive_solid->getNamedFunctions())
        _named_function_caller.addNamedFunction(std::move(fct));
    for (auto&& fct : _assembly_params.reaction_rate->getNamedFunctions())
        _named_function_caller.addNamedFunction(std::move(fct));

    // creates an extrapolator
    auto makeEx =
        [&](unsigned const n_components,
            std::vector<double> const& (TESLocalAssemblerInterface::*method)(
                const double /*t*/,
                GlobalVector const& /*current_solution*/,
                NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
                std::vector<double>& /*cache*/)
                const) -> SecondaryVariableFunctions {
        return ProcessLib::makeExtrapolator(n_components, getExtrapolator(),
                                            _local_assemblers, method);
    };

    add2nd("reaction_rate",
           makeEx(1, &TESLocalAssemblerInterface::getIntPtReactionRate));

    add2nd("darcy_velocity",
           makeEx(_mesh.getDimension(),
                  &TESLocalAssemblerInterface::getIntPtDarcyVelocity));

    add2nd("mass_flux",
           makeEx(_mesh.getDimension(),
                  &TESLocalAssemblerInterface::getIntPtMassFlux));

    add2nd("conductive_heat_flux",
           makeEx(_mesh.getDimension(),
                  &TESLocalAssemblerInterface::getIntPtConductiveHeatFlux));

    add2nd("relative_humidity",
           {1, BaseLib::easyBind(&TESProcess::computeRelativeHumidity, this),
            nullptr});

    add2nd("fluid_viscosity",
           {1, BaseLib::easyBind(&TESProcess::computeFluidViscosity, this),
            nullptr});

    // /// reactive solid state ////////////////////////////////////////
    for (auto& internal_state :
         _assembly_params.reactive_solid->getInternalStateVariables())
    {
        auto eval = std::move(internal_state.eval);
        auto const dim = internal_state.dim;

        auto f = [eval, dim](
                     TESLocalAssemblerInterface const& loc_asm,
                     const double /*t*/,
                     GlobalVector const& /*current_solution*/,
                     NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
                     std::vector<double>& cache) -> std::vector<double> const& {

            auto const& solid_states =
                loc_asm.getLocalAssemblerData().reactive_solid_state;
            auto const n = solid_states.size();
            cache.resize(dim * n);
            std::vector<double> cache2(dim);

            for (std::size_t i = 0; i < n; ++i)
            {
                auto& eval_result = eval(*solid_states[i], cache2);
                for (std::size_t d = 0; d < dim; ++d)
                {
                    cache[i + n * d] = eval_result[d];
                }
            }

            return cache;
        };

        auto ex = makeExtrapolator(internal_state.dim, getExtrapolator(),
                                   _local_assemblers,
                                   BaseLib::easyBind(std::move(f)));
        add2nd(internal_state.name, std::move(ex));
    }

    // /// reaction rate ///////////////////////////////////////////////
    for (auto& reaction_rate_variable :
         _assembly_params.reactive_solid->getRateVariables())
    {
        auto eval = std::move(reaction_rate_variable.eval);
        auto const dim = reaction_rate_variable.dim;

        auto f = [eval, dim](
                     TESLocalAssemblerInterface const& loc_asm,
                     const double /*t*/,
                     GlobalVector const& /*current_solution*/,
                     NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
                     std::vector<double>& cache) -> std::vector<double> const& {

            auto const& reaction_rates =
                loc_asm.getLocalAssemblerData().reaction_rate;
            auto const n = reaction_rates.size();
            cache.resize(dim * n);
            std::vector<double> cache2(dim);

            for (std::size_t i = 0; i < n; ++i)
            {
                auto& eval_result = eval(*reaction_rates[i], cache2);
                for (std::size_t d = 0; d < dim; ++d)
                {
                    cache[i + n * d] = eval_result[d];
                }
            }

            return cache;
        };

        auto ex = makeExtrapolator(reaction_rate_variable.dim,
                                   getExtrapolator(), _local_assemblers,
                                   BaseLib::easyBind(std::move(f)));
        add2nd(reaction_rate_variable.name, std::move(ex));
    }
}

void TESProcess::assembleConcreteProcess(const double t,
                                         GlobalVector const& x,
                                         GlobalMatrix& M,
                                         GlobalMatrix& K,
                                         GlobalVector& b)
{
    DBUG("Assemble TESProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        dof_table, t, x, M, K, b, _coupled_solutions);
}

void TESProcess::assembleWithJacobianConcreteProcess(
    const double t, GlobalVector const& x, GlobalVector const& xdot,
    const double dxdot_dx, const double dx_dx, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b, GlobalMatrix& Jac)
{
    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_table = {std::ref(*_local_to_global_index_map)};
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, dof_table, t, x, xdot, dxdot_dx, dx_dx, M, K, b, Jac,
        _coupled_solutions);
}

void TESProcess::preTimestepConcreteProcess(GlobalVector const& x,
                                            const double t,
                                            const double delta_t,
                                            const int /*process_id*/)
{
    DBUG("new timestep");

    _assembly_params.delta_t = delta_t;
    _assembly_params.current_time = t;

    GlobalExecutor::executeMemberOnDereferenced(
        &TESLocalAssemblerInterface::preTimestep, _local_assemblers,
        *_local_to_global_index_map, x, t, delta_t);
}

void TESProcess::preIterationConcreteProcess(const unsigned iter,
                                             GlobalVector const& /*x*/)
{
    _assembly_params.reactive_solid->preIteration(iter);
    _assembly_params.reaction_rate->preIteration(iter);
}

NumLib::IterationResult TESProcess::postIterationConcreteProcess(
    GlobalVector const& /*x*/)
{
    _assembly_params.reactive_solid->postIteration();
    return NumLib::IterationResult::SUCCESS;
}

GlobalVector const& TESProcess::computeRelativeHumidity(
    double const /*t*/,
    GlobalVector const& x,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::unique_ptr<GlobalVector>& result_cache)
{
    assert(&dof_table == _local_to_global_index_map.get());

    auto const& dof_table_single = getSingleComponentDOFTable();
    result_cache = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
        {dof_table_single.dofSizeWithoutGhosts(),
         dof_table_single.dofSizeWithoutGhosts(),
         &dof_table_single.getGhostIndices(), nullptr});

    GlobalIndexType const nnodes = _mesh.getNumberOfNodes();

    for (GlobalIndexType node_id = 0; node_id < nnodes; ++node_id)
    {
        auto const p = NumLib::getNodalValue(x, _mesh, dof_table, node_id,
                                             COMPONENT_ID_PRESSURE);
        auto const T = NumLib::getNodalValue(x, _mesh, dof_table, node_id,
                                             COMPONENT_ID_TEMPERATURE);
        auto const x_mV = NumLib::getNodalValue(x, _mesh, dof_table, node_id,
                                                COMPONENT_ID_MASS_FRACTION);

        auto const x_nV = Adsorption::AdsorptionReaction::getMolarFraction(
            x_mV, _assembly_params.M_react, _assembly_params.M_inert);

        auto const p_S =
            Adsorption::AdsorptionReaction::getEquilibriumVapourPressure(T);

        // TODO Problems with PETSc? (local vs. global index)
        result_cache->set(node_id, p * x_nV / p_S);
    }

    return *result_cache;
}

GlobalVector const& TESProcess::computeFluidViscosity(
    const double /*t*/,
    GlobalVector const& x,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::unique_ptr<GlobalVector>& result_cache)
{
    assert(&dof_table == _local_to_global_index_map.get());

    auto const& dof_table_single = getSingleComponentDOFTable();
    result_cache = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
        {dof_table_single.dofSizeWithoutGhosts(),
         dof_table_single.dofSizeWithoutGhosts(),
         &dof_table_single.getGhostIndices(), nullptr});

    GlobalIndexType const nnodes = _mesh.getNumberOfNodes();

    for (GlobalIndexType node_id = 0; node_id < nnodes; ++node_id)
    {
        auto const p = NumLib::getNodalValue(x, _mesh, dof_table, node_id,
                                             COMPONENT_ID_PRESSURE);
        auto const T = NumLib::getNodalValue(x, _mesh, dof_table, node_id,
                                             COMPONENT_ID_TEMPERATURE);
        auto const x_mV = NumLib::getNodalValue(x, _mesh, dof_table, node_id,
                                                COMPONENT_ID_MASS_FRACTION);

        auto const eta_GR = fluid_viscosity(p, T, x_mV);

        // TODO Problems with PETSc? (local vs. global index)
        result_cache->set(node_id, eta_GR);
    }

    return *result_cache;
}

}  // namespace TES

}  // namespace ProcessLib
