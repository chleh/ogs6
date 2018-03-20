/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateTCHSNoStokesProcess.h"

#include <cassert>

#include "Material/CreateTCHSNoStokesMaterial.h"
#include "MeshGeoToolsLib/CreateSearchLength.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "ProcessLib/Output/CreateSecondaryVariables.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "TCHSNoStokesProcess.h"
#include "TCHSNoStokesProcessData.h"

namespace ProcessLib
{
namespace TCHSNoStokes
{
template <int VelocityDim>
std::unique_ptr<Process> createTCHSNoStokesProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "TCHS_NOSTOKES");
    DBUG("Create TCHSNoStokesProcess.");

    auto const staggered_scheme =
        //! \ogs_file_param{prj__processes__process__INCOMPRESSIBLE_STOKES_BRINKMAN__coupling_scheme}
        config.getConfigParameterOptional<std::string>("coupling_scheme");
    const bool use_monolithic_scheme =
        !(staggered_scheme && (*staggered_scheme == "staggered"));

    // Process variable.

    //! \ogs_file_param{prj__processes__process__INCOMPRESSIBLE_STOKES_BRINKMAN__process_variables}
    auto const pv_config = config.getConfigSubtree("process_variables");

    ProcessVariable* variable_p;
    ProcessVariable* variable_T;
    ProcessVariable* variable_xmV;

    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>
        process_variables;

    if (use_monolithic_scheme)  // monolithic scheme.
    {
        auto per_process_variables = findProcessVariables(
            variables, pv_config,
            {//! \ogs_file_param_special{prj__processes__process__TCHS_STOKES__process_variables__pressure}
             "pressure",
             //! \ogs_file_param_special{prj__processes__process__TCHS_STOKES__process_variables__temperature}
             "temperature",
             //! \ogs_file_param_special{prj__processes__process__TCHS_STOKES__process_variables__mass_fraction}
             "mass_fraction"});

        variable_p = &per_process_variables[0].get();
        variable_T = &per_process_variables[1].get();
        variable_xmV = &per_process_variables[2].get();
        process_variables.push_back(std::move(per_process_variables));
    }
    else  // staggered scheme.
    {
        using namespace std::string_literals;
        for (auto const& variable_name :
             {"pressure"s, "temperature"s, "mass_fraction"s})
        {
            auto per_process_variables =
                findProcessVariables(variables, pv_config, {variable_name});
            process_variables.push_back(std::move(per_process_variables));
        }

        variable_p = &process_variables[0][0].get();
        variable_T = &process_variables[1][0].get();
        variable_xmV = &process_variables[2][0].get();
    }

    DBUG("Associate pressure with process variable `%s'.",
         variable_p->getName().c_str());
    if (variable_p->getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Pressure process variable '%s' is not a scalar variable but has "
            "%d components.",
            variable_p->getName().c_str(),
            variable_p->getNumberOfComponents());
    }

    DBUG("Associate temperature with process variable `%s'.",
         variable_T->getName().c_str());
    if (variable_T->getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Temperature process variable '%s' is not a scalar variable but "
            "has %d components.",
            variable_T->getName().c_str(),
            variable_T->getNumberOfComponents());
    }

    DBUG("Associate mass fraction with process variable `%s'.",
         variable_xmV->getName().c_str());
    if (variable_xmV->getNumberOfComponents() != 1)
    {
        OGS_FATAL(
            "Mass fraction process variable '%s' is not a scalar variable but "
            "has %d components.",
            variable_xmV->getName().c_str(),
            variable_xmV->getNumberOfComponents());
    }

    auto const* material_ids =
        mesh.getProperties().getPropertyVector<int>("MaterialIDs");

    if (!(material_ids &&
          material_ids->getMeshItemType() == MeshLib::MeshItemType::Cell &&
          material_ids->getNumberOfComponents() == 1))
    {
        OGS_FATAL("Field `MaterialIDs' is not set up properly.");
    }

    auto materials = Material::createTCHSNoStokesMaterials(
        config.getConfigSubtree("materials"), parameters);

    auto const probe_config = config.getConfigSubtree("darcy_velocity_probe");
    auto const probe_coords =
        probe_config.getConfigParameter<std::vector<double>>("coords");
    if (probe_coords.size() != 3)
        OGS_FATAL("Wrong number of coordinates.");

    auto search_length_algorithm =
        MeshGeoToolsLib::createSearchLengthAlgorithm(probe_config, mesh);

    auto const& mesh_node_searcher =
        MeshGeoToolsLib::MeshNodeSearcher::getMeshNodeSearcher(
            mesh, std::move(search_length_algorithm));
    GeoLib::Point pnt(probe_coords[0], probe_coords[1], probe_coords[2]);
    auto const probe_node_ids = mesh_node_searcher.getMeshNodeIDsForPoint(pnt);
    if (probe_node_ids.size() != 1)
    {
        OGS_FATAL(
            "The velocity probe could not be associated to exactly one mesh "
            "node. Instead %d nodes were found.",
            probe_node_ids.size());
    }
    auto const probe_node_id = probe_node_ids.front();
    auto const probe_coords_found = *mesh.getNode(probe_node_id);
    INFO(
        "Found vertex (%g, %g, %g) id %d for velocity probe. Requested "
        "coordinates: (%g, %g, %g).",
        probe_coords_found[0], probe_coords_found[1], probe_coords_found[2],
        probe_node_id, probe_coords[0], probe_coords[1], probe_coords[2]);

    auto const darcy_velocity_center =
        config.getConfigParameter<double>("darcy_velocity_center");

    TCHSNoStokesProcessData<VelocityDim> process_data{
        *material_ids, std::move(materials), probe_node_id,
        darcy_velocity_center};

    SecondaryVariableCollection secondary_variables;

    NumLib::NamedFunctionCaller named_function_caller(
        {"TCHSNoStokes_pressure", "TCHSNoStokes_temperature",
         "TCHSNoStokes_mass_fraction"});

    ProcessLib::createSecondaryVariables(config, secondary_variables,
                                         named_function_caller);

    return std::make_unique<TCHSNoStokesProcess<VelocityDim>>(
        mesh, std::move(jacobian_assembler), parameters, integration_order,
        std::move(process_variables), std::move(process_data),
        std::move(secondary_variables), std::move(named_function_caller),
        use_monolithic_scheme);
}

template std::unique_ptr<Process> createTCHSNoStokesProcess<2>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);

// TODO saving some more compile time.
#if 0
template std::unique_ptr<Process> createTCHSNoStokesProcess<3>(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    BaseLib::ConfigTree const& config);
#endif

}  // namespace TCHSNoStokes
}  // namespace ProcessLib
