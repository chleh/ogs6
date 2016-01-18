/**
 * \author Karsten Rink
 * \date   2010-08-25
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROJECTDATA_H_
#define PROJECTDATA_H_

#include <memory>

#include "BaseLib/ConfigTreeNew.h"

#include "GeoLib/GEOObjects.h"

#include "NumLib/TimeStepping/Algorithms/ITimeStepAlgorithm.h"

#include "ProcessLib/ProcessVariable.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/Parameter.h"
#include "ProcessLib/Output.h"
#include "ProcessLib/GroundwaterFlowProcess-fwd.h"

namespace MeshLib {
	class Mesh;
}

/**
 * The ProjectData Object contains all the data needed for a certain project, i.e. all
 * geometric data (stored in a GEOObjects-object), all the meshes, processes,
 * and process variables.
 */
class ProjectData
{
public:
	/// The empty constructor used in the gui, for example, when the project's
	/// configuration is not loaded yet.
	ProjectData() = default;

	/// Constructs project data by parsing provided configuration.
	/// The additional  path is used to find files referenced in the
	/// configuration.
	ProjectData(BaseLib::ConfigTreeNew const& config_tree,
	            std::string const& path);

	ProjectData(ProjectData&) = delete;
	virtual ~ProjectData();

	/// Returns the GEOObjects containing all points, polylines and surfaces.
	GeoLib::GEOObjects* getGEOObjects()
	{
		return _geoObjects;
	}

	/// Adds a new mesh under a (possibly new) unique name.
	/// \attention This might change the given mesh's name.
	void addMesh(MeshLib::Mesh* mesh);

	/// Returns the mesh with the given name or a \c nullptr if the mesh was not
	/// found.
	const MeshLib::Mesh* getMesh(const std::string &name) const;

	/// Returns all the meshes with their respective names
	/// \attention This method should be used only by methods garanteeing
	/// read-only access to the meshes.
	/// \todo This method breaks encapsulation.
	const std::vector<MeshLib::Mesh*>& getMeshObjects() const
	{
		return _mesh_vec;
	}

	/// Deletes all meshes with the given name and removes them from the list of
	/// saved meshes. If any mesh was found for removal, true is returned and
	/// false otherwise.
	bool removeMesh(const std::string &name);

	//
	// Process interface
	//

	/// Builder for processes. The supplied template defines the types of global
	/// vectors and matrices, and the global executor. These types are passed to
	/// every of the constructed processes.
	template <typename GlobalSetupType>
	void buildProcesses()
	{
		for (auto const& pc : _process_configs)
		{
			auto const type = pc.peekConfParam<std::string>("type");
			if (type == "GROUNDWATER_FLOW") {
				// The existence check of the in the configuration referenced
				// process variables is checked in the physical process.
				// TODO at the moment we have only one mesh, later there can be
				// several meshes. Then we have to assign the referenced mesh
				// here.
				_processes.emplace_back(
					new ProcessLib::GroundwaterFlowProcess<GlobalSetupType>(
						*_mesh_vec[0], _process_variables, _parameters, pc));
			}
			else
			{
				ERR("Unknown process type: %s\n", type.c_str());
			}
		}

		// process configs are not needed anymore, so clear the storage
		// in order to trigger config tree checks
		_process_configs.clear();
	}

	/// Iterator access for processes.
	/// Provides read access to the process container.
	std::vector<
	    std::unique_ptr<ProcessLib::Process<GlobalSetupType>>>::const_iterator
	processesBegin() const
	{
		return _processes.begin();
	}
	std::vector<std::unique_ptr<ProcessLib::Process<GlobalSetupType>>>::iterator
	processesBegin()
	{
		return _processes.begin();
	}

	/// Iterator access for processes as in processesBegin().
	std::vector<
	    std::unique_ptr<ProcessLib::Process<GlobalSetupType>>>::const_iterator
	processesEnd() const
	{
		return _processes.end();
	}
	std::vector<std::unique_ptr<ProcessLib::Process<GlobalSetupType>>>::iterator
	processesEnd()
	{
		return _processes.end();
	}

	std::string const&
	getOutputFilePrefix() const
	{
		return _output_file_prefix;
	}

	NumLib::ITimeStepAlgorithm const& getTimeStepper() const
	{
		return *_time_stepper;
	}
	NumLib::ITimeStepAlgorithm& getTimeStepper()
	{
		return *_time_stepper;
	}

private:
	/// Checks if a mesh with the same name exists and provides a unique name in
	/// case of already existing mesh. Returns true if the mesh name is unique.
	/// Returns false and changes the provided name to a unique name otherwise.
	bool isMeshNameUniqueAndProvideUniqueName(std::string &name) const;

	/// Returns true if a mesh with the same name exists and false otherwise.
	bool meshExists(const std::string &name) const;


	/// Returns an iterator to the first found mesh with the given name.
	std::vector<MeshLib::Mesh*>::const_iterator findMeshByName(
		std::string const& name) const;
	std::vector<MeshLib::Mesh*>::iterator findMeshByName(
		std::string const& name);

	/// Parses the process variables configuration and creates new variables for
	/// each variable entry passing the corresponding subtree to the process
	/// variable constructor.
	void parseProcessVariables(BaseLib::ConfigTreeNew const& process_variables_config);

	/// Parses the parameters configuration and saves them in a list.
	/// Checks if a parameter has name tag.
	void parseParameters(BaseLib::ConfigTreeNew const& parameters_config);

	/// Parses the processes configuration and creates new processes for each
	/// process entry passing the corresponding subtree to the process
	/// constructor.
	void parseProcesses(BaseLib::ConfigTreeNew const& process_config);

	/// Parses the output configuration.
	/// Parses the file tag and sets output file prefix.
	void parseOutput(BaseLib::ConfigTreeNew const& output_config, std::string const& path);

	void parseTimeStepping(BaseLib::ConfigTreeNew const& timestepping_config);

private:
	GeoLib::GEOObjects *_geoObjects = new GeoLib::GEOObjects();
	std::vector<MeshLib::Mesh*> _mesh_vec;
	std::vector<std::unique_ptr<ProcessLib::Process<GlobalSetupType>>>
	    _processes;
	std::vector<ProcessLib::ProcessVariable> _process_variables;

	/// Buffer for each process' config used in the process building function.
	std::vector<BaseLib::ConfigTreeNew> _process_configs;

	/// Buffer for each parameter config passed to the process.
	std::vector<std::unique_ptr<ProcessLib::ParameterBase>> _parameters;

	/// Timestepper
	std::unique_ptr<NumLib::ITimeStepAlgorithm> _time_stepper;

	std::unique_ptr<ProcessLib::Output> _output;
};

#endif //PROJECTDATA_H_
