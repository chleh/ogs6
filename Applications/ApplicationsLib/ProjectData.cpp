/**
 * \author Karsten Rink
 * \date   2010-08-25
 * \brief  Implementation of the project data class.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ProjectData.h"

#include <algorithm>

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "BaseLib/FileTools.h"

#include "MeshLib/Mesh.h"

#include "NumLib/TimeStepping/Algorithms/FixedTimeStepping.h"

// FileIO
#include "FileIO/XmlIO/Boost/BoostXmlGmlInterface.h"
#include "FileIO/readMeshFromFile.h"

#include "BaseLib/ConfigTreeNew.h"

namespace detail
{
static
void readGeometry(std::string const& fname, GeoLib::GEOObjects & geo_objects)
{
	DBUG("Reading geometry file \'%s\'.", fname.c_str());
	FileIO::BoostXmlGmlInterface gml_reader(geo_objects);
	gml_reader.readFile(fname);
}

}

ProjectData::ProjectData(BaseLib::ConfigTreeNew const& project_config,
	std::string const& path)
{
	// geometry
	std::string const geometry_file = BaseLib::copyPathToFileName(
			project_config.getConfParam<std::string>("geometry"), path
		);
	detail::readGeometry(geometry_file, *_geoObjects);

	// mesh
	std::string const mesh_file = BaseLib::copyPathToFileName(
			project_config.getConfParam<std::string>("mesh"), path
		);

	MeshLib::Mesh* const mesh = FileIO::readMeshFromFile(mesh_file);
	if (!mesh)
		ERR("Could not read mesh from \'%s\' file. No mesh added.",
			mesh_file.c_str());
	_mesh_vec.push_back(mesh);

	// process variables

	parseProcessVariables(project_config.getConfSubtree("process_variables"));

	// parameters
	parseParameters(project_config.getConfSubtree("parameters"));

	// processes
	parseProcesses(project_config.getConfSubtree("processes"));

	// output
	parseOutput(project_config.getConfSubtree("output"), path);

	// timestepping
	parseTimeStepping(project_config.getConfSubtree("time_stepping"));
}

ProjectData::~ProjectData()
{
	delete _geoObjects;

	for (MeshLib::Mesh* m : _mesh_vec)
		delete m;
}

void ProjectData::addMesh(MeshLib::Mesh* mesh)
{
	std::string name = mesh->getName();
	isMeshNameUniqueAndProvideUniqueName(name);
	mesh->setName(name);
	_mesh_vec.push_back(mesh);
}

std::vector<MeshLib::Mesh*>::const_iterator ProjectData::findMeshByName(
		std::string const& name) const
{
	return const_cast<ProjectData&>(*this).findMeshByName(name);
}

std::vector<MeshLib::Mesh*>::iterator ProjectData::findMeshByName(
		std::string const& name)
{
	return std::find_if(_mesh_vec.begin(), _mesh_vec.end(),
			[&name](MeshLib::Mesh* mesh)
			{
				return mesh && (name == mesh->getName());
			});
}

const MeshLib::Mesh* ProjectData::getMesh(const std::string &name) const
{
	std::vector<MeshLib::Mesh*>::const_iterator it = findMeshByName(name);
	return (it == _mesh_vec.end() ? nullptr : *it);
}

bool ProjectData::removeMesh(const std::string &name)
{
	bool mesh_found = false;
	std::vector<MeshLib::Mesh*>::iterator it = findMeshByName(name);
	while (it != _mesh_vec.end())
	{
		delete *it;
		*it = nullptr;
		it = findMeshByName(name);
		mesh_found = true;
	}

	_mesh_vec.erase(std::remove(_mesh_vec.begin(), _mesh_vec.end(), nullptr),
			_mesh_vec.end());
	return mesh_found;
}

bool ProjectData::meshExists(const std::string &name) const
{
	return findMeshByName(name) != _mesh_vec.end();
}

bool ProjectData::isMeshNameUniqueAndProvideUniqueName(std::string &name) const
{
	int count(0);
	bool isUnique(false);
	std::string cpName;

	while (!isUnique)
	{
		isUnique = true;
		cpName = name;

		count++;
		// If the original name already exists we start to add numbers to name for
		// as long as it takes to make the name unique.
		if (count > 1)
			cpName = cpName + "-" + std::to_string(count);

		for (std::vector<MeshLib::Mesh*>::const_iterator it = _mesh_vec.begin();
				it != _mesh_vec.end(); ++it)
			if ( cpName.compare((*it)->getName()) == 0 )
				isUnique = false;
	}

	// At this point cpName is a unique name and isUnique is true.
	// If cpName is not the original name, "name" is changed and isUnique is set to false,
	// indicating that a vector with the original name already exists.
	if (count > 1)
	{
		isUnique = false;
		name = cpName;
	}
	return isUnique;
}

void ProjectData::parseProcessVariables(
	BaseLib::ConfigTreeNew const& process_variables_config)
{
	DBUG("Parse process variables:")
	if (_geoObjects == nullptr) {
		ERR("Geometric objects are required to define process variables.");
		ERR("No geometric objects present.");
		return;
	}

	// TODO at the moment we have only one mesh, later there
	// can be several meshes. Then we have to check for correct mesh here and
	// assign the referenced mesh below.
	if (_mesh_vec.empty() || _mesh_vec[0] == nullptr) {
		ERR("A mesh is required to define process variables.");
		return;
	}

	// _process_variables.reserve(process_variables_config.size());

	for (auto var_config
		 : process_variables_config.getConfSubtreeList("process_variable")) {
		// TODO Extend to referenced meshes.
		_process_variables.emplace_back(var_config, *_mesh_vec[0], *_geoObjects);
	}
}

void ProjectData::parseParameters(BaseLib::ConfigTreeNew const& parameters_config)
{
	using namespace ProcessLib;

	DBUG("Reading parameters:");
	for (auto parameter_config : parameters_config.getConfSubtreeList("parameter"))
	{
		auto name = parameter_config.getConfParam<std::string>("name");
		auto type = parameter_config.peekConfParam<std::string>("type");

		// Create parameter based on the provided type.
		if (type == "Constant")
		{
			INFO("ConstantParameter: %s.", name.c_str());
			_parameters.push_back(createConstParameter(parameter_config));
			_parameters.back()->name = name;
		}
		else if (type == "MeshProperty")
		{
			INFO("MeshPropertyParameter: %s", name.c_str());
			_parameters.push_back(
			    createMeshPropertyParameter(parameter_config, *_mesh_vec[0]));
			_parameters.back()->name = name;
		}
		else
		{
			ERR("Cannot construct property of given type \'%s\'.",
			    type.c_str());
			std::abort();
		}
	}
}

void ProjectData::parseProcesses(BaseLib::ConfigTreeNew const& processes_config)
{
	DBUG("Reading processes:");
	for (auto process_config : processes_config.getConfSubtreeList("process")) {
		// process type must be specified.
		process_config.peekConfParam<std::string>("type");
		process_config.ignoreConfParam("name");
		_process_configs.push_back(std::move(process_config));
	}
}

void ProjectData::parseOutput(BaseLib::ConfigTreeNew const& output_config,
	std::string const& path)
{
	output_config.checkConfParam("type", "VTK");
	DBUG("Parse output configuration:");

	_output.reset(ProcessLib::Output::newInstance(output_config, path));
}

void ProjectData::parseTimeStepping(BaseLib::ConfigTreeNew const& timestepping_config)
{
	using namespace ProcessLib;

	DBUG("Reading timestepping configuration.");

	auto const type = timestepping_config.peekConfParam<std::string>("type");

	if (type == "FixedTimeStepping")
	{
		_time_stepper = NumLib::FixedTimeStepping::newInstance(timestepping_config);
	}
	else if (type == "SingleStep")
	{
		timestepping_config.ignoreConfParam("type");
		_time_stepper.reset(new NumLib::FixedTimeStepping(0.0, 1.0, 1.0));
	}
	else
	{
		ERR("Unknown timestepper type: `%s'.", type.c_str());
		std::abort();
	}

	if (!_time_stepper)
	{
		ERR("Initialization of timestepper failed.");
		std::abort();
	}
}



