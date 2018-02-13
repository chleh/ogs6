/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateIncompressibleStokesBrinkmanProcess.h"

#include <cassert>

#include "MaterialLib/SolidModels/CreateLinearElasticIsotropic.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"

#include "IncompressibleStokesBrinkmanProcess.h"
#include "IncompressibleStokesBrinkmanProcessData.h"

namespace ProcessLib
{
namespace IncompressibleStokesBrinkman
{
template <int DisplacementDim>
std::unique_ptr<Process> createIncompressibleStokesBrinkmanProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "INCOMPRESSIBLE_STOKES_BRINKMAN");
    DBUG("Create IncompressibleStokesBrinkmanProcess.");

    auto const staggered_scheme =
        //! \ogs_file_param{prj__processes__process__INCOMPRESSIBLE_STOKES_BRINKMAN__coupling_scheme}
        config.getConfigParameterOptional<std::string>("coupling_scheme");
    const bool use_monolithic_scheme =
        !(staggered_scheme && (*staggered_scheme == "staggered"));

    // Process variable.

    //! \ogs_file_param{prj__processes__process__INCOMPRESSIBLE_STOKES_BRINKMAN__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    ProcessVariable* variable_p;
    ProcessVariable* variable_u;
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;
    if (use_monolithic_scheme)  // monolithic scheme.
    {
        auto per_process_variables = findProcessVariables(
            variables, pv_config,
            {//! \ogs_file_param_special{prj__processes__process__INCOMPRESSIBLE_STOKES_BRINKMAN__process_variables__pressure}
             "pressure",
             //! \ogs_file_param_special{prj__processes__process__INCOMPRESSIBLE_STOKES_BRINKMAN__process_variables__velocity}
             "velocity"});
        variable_p = &per_process_variables[0].get();
        variable_u = &per_process_variables[1].get();
        process_variables.push_back(std::move(per_process_variables));
    }
    else  // staggered scheme.
    {
        using namespace std::string_literals;
        for (auto const& variable_name : {"pressure"s, "velocity"s})
        {
            auto per_process_variables =
                findProcessVariables(variables, pv_config, {variable_name});
            process_variables.push_back(std::move(per_process_variables));
        }
        variable_p = &process_variables[0][0].get();
        variable_u = &process_variables[1][0].get();
    }

    DBUG("Associate velocity with process variable \'%s\'.",
         variable_u->getName().c_str());

    if (variable_u->getNumberOfComponents() != DisplacementDim)
    {
        OGS_FATAL(
            "Number of components of the process variable '%s' is different "
            "from the velocity dimension: got %d, expected %d",
            variable_u->getName().c_str(),
            variable_u->getNumberOfComponents(),
            DisplacementDim);
    }

    DBUG("Associate pressure with process variable \'%s\'.",
         variable_p->getName().c_str());
    if (variable_p->getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Pressure process variable '%s' is not a scalar variable but has "
            "%d components.",
            variable_p->getName().c_str(),
            variable_p->getNumberOfComponents());
    }

    auto const* material_ids =
        mesh.getProperties().getPropertyVector<int>("MaterialIDs");

    if (!(material_ids &&
          material_ids->getMeshItemType() == MeshLib::MeshItemType::Cell) &&
        material_ids->getNumberOfComponents() == 1)
    {
        OGS_FATAL("Field `MaterialIDs' is not set up properly.");
    }

    auto const pellet_diameter =
        config.getConfigParameter<double>("pellet_diameter");
    auto const bed_radius = config.getConfigParameter<double>("bed_radius");

    auto const average_darcy_velocity =
        config.getConfigParameter<double>("average_darcy_velocity");

    auto const homogeneous_porosity =
        config.getConfigParameter<double>("homogeneous_porosity");

    auto const& fluid_density =
        findParameter<double>(config, "fluid_density", parameters, 1);

    auto const& fluid_viscosity =
        findParameter<double>(config, "fluid_viscosity", parameters, 1);

    auto effective_fluid_viscosity = createEffectiveFluidViscosity(
        config.getConfigSubtree("effective_fluid_viscosity"));

    IncompressibleStokesBrinkmanProcessData<DisplacementDim> process_data{
        *material_ids,        pellet_diameter,
        bed_radius,           average_darcy_velocity,
        homogeneous_porosity, fluid_density,
        fluid_viscosity,      std::move(effective_fluid_viscosity)};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"IncompressibleStokesBrinkman_pressure",
         "IncompressibleStokesBrinkman_velocity"});

    ProcessLib::createSecondaryVariables(config, secondary_variables,
                                         named_function_caller);

    return std::make_unique<
        IncompressibleStokesBrinkmanProcess<DisplacementDim>>(
        mesh, std::move(jacobian_assembler), parameters, integration_order,
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(named_function_caller),
        use_monolithic_scheme);
}

template std::unique_ptr<Process> createIncompressibleStokesBrinkmanProcess<2>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

template std::unique_ptr<Process> createIncompressibleStokesBrinkmanProcess<3>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

}  // namespace IncompressibleStokesBrinkman
}  // namespace ProcessLib
