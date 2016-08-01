/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ProcessVariable.h"

#include <utility>

#include <logog/include/logog.hpp>

#include "GeoLib/GEOObjects.h"
#include "MeshLib/Mesh.h"

#include "BaseLib/ConfigTree.h"

namespace ProcessLib
{
ProcessVariable::ProcessVariable(BaseLib::ConfigTree const& config,
                                 MeshLib::Mesh& mesh,
                                 GeoLib::GEOObjects const& geometries)
    :
    //! \ogs_file_param{prj__process_variables__process_variable__name}
    _name(config.getConfigParameter<std::string>("name")),
    _mesh(mesh),
    //! \ogs_file_param{prj__process_variables__process_variable__components}
    _n_components(config.getConfigParameter<int>("components"))
{
    DBUG("Constructing process variable %s", this->_name.c_str());

    // Initial condition
    //! \ogs_file_param{prj__process_variables__process_variable__initial_condition}
    if (auto ic_config = config.getConfigSubtreeOptional("initial_condition"))
    {
        //! \ogs_file_param{initial_condition__type}
        auto const type = ic_config->peekConfigParameter<std::string>("type");
        if (type == "Uniform")
        {
            _initial_condition =
                createUniformInitialCondition(*ic_config, _n_components);
        }
        else if (type == "MeshProperty")
        {
            _initial_condition =
                createMeshPropertyInitialCondition(*ic_config, _mesh, _n_components);
        }
        else
        {
            ERR("Unknown type of the initial condition.");
        }
    }
    else
    {
        INFO("No initial condition found.");
    }

    // Boundary conditions
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions}
    if (auto bcs_config = config.getConfigSubtreeOptional("boundary_conditions"))
    {
        for (auto bc_config :
             //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition}
             bcs_config->getConfigSubtreeList("boundary_condition"))
        {
            auto const geometrical_set_name =
                    //! \ogs_file_param{boundary_condition__geometrical_set}
                    bc_config.getConfigParameter<std::string>("geometrical_set");
            auto const geometry_name =
                    //! \ogs_file_param{boundary_condition__geometry}
                    bc_config.getConfigParameter<std::string>("geometry");

            GeoLib::GeoObject const* const geometry_ptr =
                geometries.getGeoObject(geometrical_set_name, geometry_name);
            assert(geometry_ptr != nullptr);
            GeoLib::GeoObject const& geometry = *geometry_ptr;

            DBUG(
                "Found geometry type \"%s\"",
                GeoLib::convertGeoTypeToString(geometry.getGeoType()).c_str());

            // Construct type dependent boundary condition
            //! \ogs_file_param{boundary_condition__type}
            auto const type = bc_config.peekConfigParameter<std::string>("type");

            if (type == "UniformDirichlet")
            {
                _dirichlet_bc_configs.emplace_back(std::make_pair(
                    std::unique_ptr<UniformDirichletBoundaryCondition>(
                        new UniformDirichletBoundaryCondition(geometry,
                                                              bc_config)),
                    0));  // TODO, the 0 stands for component_id. Need parser.
            }
            else if (type == "UniformNeumann")
            {
                _neumann_bc_configs.emplace_back(std::make_pair(
                    std::unique_ptr<NeumannBcConfig>{
                        new NeumannBcConfig(geometry, bc_config)},
                    0));  // TODO, the 0 stands for component_id. Need parser.
            }
            else
            {
                ERR("Unknown type \'%s\' of the boundary condition.",
                    type.c_str());
            }
        }
    } else {
        INFO("No boundary conditions found.");
    }

    // Source Terms
    config.ignoreConfigParameter("source_terms");
}

ProcessVariable::ProcessVariable(ProcessVariable&& other)
    : _name(std::move(other._name)),
      _mesh(other._mesh),
      _n_components(other._n_components),
      _initial_condition(std::move(other._initial_condition)),
      _dirichlet_bc_configs(std::move(other._dirichlet_bc_configs)),
      _neumann_bc_configs(std::move(other._neumann_bc_configs))
{
}

std::string const& ProcessVariable::getName() const
{
    return _name;
}

MeshLib::Mesh const& ProcessVariable::getMesh() const
{
    return _mesh;
}

MeshLib::PropertyVector<double>& ProcessVariable::getOrCreateMeshProperty()
{
    if (_mesh.getProperties().hasPropertyVector(_name))
    {
        auto result =
            _mesh.getProperties().template getPropertyVector<double>(_name);
        assert(result);
        assert(result->size() == _mesh.getNumberOfNodes() * _n_components);
        return *result;
    }
    else
    {
        auto result = _mesh.getProperties().template createNewPropertyVector<double>(
            _name, MeshLib::MeshItemType::Node, _n_components);
        assert(result);
        result->resize(_mesh.getNumberOfNodes() * _n_components);
        return *result;
    }
}

std::vector<std::unique_ptr<BoundaryCondition>>
ProcessVariable::getBoundaryConditions(
    const NumLib::LocalToGlobalIndexMap& dof_table,
    const int variable_id,
    unsigned const integration_order)
{
    std::vector<std::unique_ptr<BoundaryCondition>> bcs;

    if (_dirichlet_bc_configs.empty() && _neumann_bc_configs.empty())
        return bcs;

    MeshGeoToolsLib::MeshNodeSearcher& mesh_node_searcher =
        MeshGeoToolsLib::MeshNodeSearcher::getMeshNodeSearcher(_mesh);

    for (auto& bc_config : _dirichlet_bc_configs)
    {
        auto const component_id = bc_config.second;
        bcs.emplace_back(bc_config.first->getDirichletBoundaryCondition(
            mesh_node_searcher, dof_table, variable_id, component_id));
    }

    if (_neumann_bc_configs.empty())
        return bcs;

    MeshGeoToolsLib::BoundaryElementsSearcher boundary_element_searcher(
        _mesh, mesh_node_searcher);

    for (auto& bc_config : _neumann_bc_configs)
    {
        auto const component_id = bc_config.second;

        bc_config.first->initialize(boundary_element_searcher);
        bcs.emplace_back(new NeumannBc(*bc_config.first, integration_order,
                                       dof_table, variable_id, component_id,
                                       _mesh.getDimension()));
    }

    return bcs;
}

}  // namespace ProcessLib
