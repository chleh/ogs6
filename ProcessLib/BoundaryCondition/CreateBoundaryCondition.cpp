/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateBoundaryCondition.h"

#include "BoundaryCondition.h"
#include "BoundaryConditionConfig.h"
#include "ConstraintDirichletBoundaryCondition.h"
#include "DirichletBoundaryCondition.h"
#include "NeumannBoundaryCondition.h"
#include "NonuniformDirichletBoundaryCondition.h"
#include "NonuniformNeumannBoundaryCondition.h"
#include "NormalTractionBoundaryCondition.h"
#include "PhaseFieldIrreversibleDamageOracleBoundaryCondition.h"
#include "RobinBoundaryCondition.h"
#ifdef OGS_USE_PYTHON
#include "PythonBoundaryCondition.h"
#endif

namespace ProcessLib
{
std::unique_ptr<BoundaryCondition> createBoundaryCondition(
    const BoundaryConditionConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table, const MeshLib::Mesh& mesh,
    const int variable_id, const unsigned integration_order,
    const unsigned shapefunction_order,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters,
    const Process& process)
{
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    auto const type = config.config.peekConfigParameter<std::string>("type");

    if (type == "Dirichlet")
    {
        return ProcessLib::createDirichletBoundaryCondition(
            config.config, config.mesh, dof_table, mesh.getID(), variable_id,
            *config.component_id, parameters);
    }
    if (type == "Neumann")
    {
        return ProcessLib::createNeumannBoundaryCondition(
            config.config, config.mesh, dof_table, variable_id,
            *config.component_id, mesh.isAxiallySymmetric(), integration_order,
            shapefunction_order, mesh.getDimension(), parameters);
    }
    if (type == "Robin")
    {
        return ProcessLib::createRobinBoundaryCondition(
            config.config, config.mesh, dof_table, variable_id,
            *config.component_id, mesh.isAxiallySymmetric(), integration_order,
            shapefunction_order, mesh.getDimension(), parameters);
    }
    if (type == "NonuniformDirichlet")
    {
        return ProcessLib::createNonuniformDirichletBoundaryCondition(
            config.config, dof_table, variable_id, *config.component_id, mesh);
    }
    if (type == "NonuniformNeumann")
    {
        return ProcessLib::createNonuniformNeumannBoundaryCondition(
            config.config, dof_table, variable_id, *config.component_id,
            integration_order, shapefunction_order, mesh);
    }
#ifdef OGS_USE_PYTHON
    if (type == "Python")
    {
        return nullptr;
#if 0
        std::unique_ptr<MeshGeoToolsLib::SearchLength> search_length_algorithm =
            MeshGeoToolsLib::createSearchLengthAlgorithm(config.config, mesh);

        MeshGeoToolsLib::MeshNodeSearcher const& mesh_node_searcher =
            MeshGeoToolsLib::MeshNodeSearcher::getMeshNodeSearcher(
                mesh, std::move(search_length_algorithm));

        MeshGeoToolsLib::BoundaryElementsSearcher boundary_element_searcher(
            mesh, mesh_node_searcher);

        return ProcessLib::createPythonBoundaryCondition(
            config.config, mesh_node_searcher.getMeshNodeIDs(config.geometry),
            getClonedElements(boundary_element_searcher, config.geometry),
            dof_table, variable_id, *config.component_id, mesh,
            integration_order, shapefunction_order, parameters);
#endif
    }
#endif

    //
    // Special boundary conditions
    //
    if (type == "ConstraintDirichlet")
    {
        return createConstraintDirichletBoundaryCondition(
            config.config, config.mesh, dof_table, variable_id,
            integration_order, *config.component_id, parameters, process);
    }
    if (type == "NormalTraction")
    {
        return ProcessLib::NormalTractionBoundaryCondition::
            createNormalTractionBoundaryCondition(
                config.config, config.mesh, dof_table, variable_id,
                mesh.isAxiallySymmetric(), integration_order,
                shapefunction_order, mesh.getDimension(), parameters);
    }
    if (type == "PhaseFieldIrreversibleDamageOracleBoundaryCondition")
    {
        return ProcessLib::
            createPhaseFieldIrreversibleDamageOracleBoundaryCondition(
                config.config, dof_table, mesh, variable_id,
                *config.component_id);
    }
    OGS_FATAL("Unknown boundary condition type: `%s'.", type.c_str());
}

}  // namespace ProcessLib
