/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateTESProcess.h"
#include "CreateJacobianAssembler.h"
#include "TESProcess.h"

namespace ProcessLib
{
namespace TES
{

std::unique_ptr<Process> createTESProcess(
    MeshLib::Mesh& mesh,
    Process::NonlinearSolver& nonlinear_solver,
    std::unique_ptr<Process::TimeDiscretization>&& time_discretization,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& /*parameters*/,
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "TES");

    DBUG("Create TESProcess.");

    auto process_variables = findProcessVariables(
        variables, config,
        {"fluid_pressure", "temperature", "vapour_mass_fraction"});

    SecondaryVariableCollection secondary_variables{
        config.getConfigSubtreeOptional("secondary_variables"),
        {"solid_density", "reaction_rate", "velocity_x", "velocity_y",
         "velocity_z", "loading", "reaction_damping_factor",
         "vapour_partial_pressure", "relative_humidity",
         "equilibrium_loading"}};

    ProcessOutput process_output{config.getConfigSubtree("output"),
                                 process_variables, secondary_variables};

    std::unique_ptr<ProcessLib::AbstractJacobianAssembler> jacobian_assembler;
    if (auto jac_asm_config =
            config.getConfigSubtreeOptional("jacobian_assembler")) {
        jacobian_assembler =
            ProcessLib::createJacobianAssembler(*jac_asm_config);
    }

    return std::unique_ptr<Process>{new TESProcess{
        mesh, nonlinear_solver, std::move(time_discretization),
        std::move(process_variables), std::move(secondary_variables),
        std::move(process_output), std::move(jacobian_assembler), config}};
}

}  // namespace TES
}  // namespace ProcessLib
