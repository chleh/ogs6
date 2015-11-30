/**
 * \date   2014-08-04
 * \brief  Implementation of OpenGeoSys simulation application
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <memory>

#ifndef NDEBUG
#include <fenv.h>
#endif

// ThirdParty/tclap
#include "tclap/CmdLine.h"

// BaseLib
#include "BaseLib/BuildInfo.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/ConfigTreeNew.h"
#include "BaseLib/FileTools.h"

#include "FileIO/VtkIO/PVDFile.h"

#include "Applications/ApplicationsLib/LinearSolverLibrarySetup.h"
#include "Applications/ApplicationsLib/LogogSetup.h"
#include "Applications/ApplicationsLib/ProjectData.h"

#include "ProcessLib/NumericsConfig.h"

void solveProcesses(ProjectData &project)
{
	INFO("Solve processes.");

	auto &time_stepper = project.getTimeStepper();
	auto &output_control = project.getOutputControl();

	std::string const out_pref = output_control.getFilePrefix();

	auto do_output = [out_pref](unsigned pcs, ProcessLib::Process* p,
					 unsigned ts, double current_time, FileIO::PVDFile& pvdf)
	{
		std::string const output_file_name =
				out_pref + "_pcs_" + std::to_string(pcs)
				+ "_ts_" + std::to_string(ts)
				+ "_t_"  + std::to_string(current_time)
				+ ".vtu";
		DBUG("output to %s", output_file_name.c_str());
		// auto const& p = _processes[pcs];
		p->postTimestep(output_file_name, ts);
		pvdf.addVTUFile(output_file_name, current_time);
	};

	std::vector<std::unique_ptr<FileIO::PVDFile> > pvd_files;

	{
		unsigned i = 0;  // process counter, used to distinguish output files
		for (auto p = project.processesBegin(); p != project.processesEnd();
			 ++p)
		{
			std::string const fn = out_pref + "_pcs_" + std::to_string(i)
								   + ".pvd";
			pvd_files.emplace_back(new FileIO::PVDFile(fn));

			do_output(i, *p, 0, 0.0, *pvd_files.back());
		}
	}

	bool output_timestep = false;

	while (time_stepper.next()) // skips zeroth timestep, but OK since end of first timestep is after first delta t
	{
		BaseLib::TimingOneShot timing("timestep");

		const auto dt = time_stepper.getTimeStep().dt();
		const auto current_time = time_stepper.getTimeStep().current();
		const auto timestep = time_stepper.getTimeStep().steps();

		output_timestep = output_control.doOutput(timestep);

		INFO("=================== timestep %i === %g s ===================",
			 timestep, current_time);

		bool accepted = true;

		unsigned i = 0;  // process counter, used to distinguish output files
		for (auto p = project.processesBegin(); p != project.processesEnd();
		     ++p)
		{
			accepted = accepted && (*p)->solve(current_time, dt);

			if (!accepted) {
				ERR("Timestep has not been accepted. Aborting.");
				break;
			}

			if (output_timestep)
			{
				do_output(i, *p, timestep, current_time, *pvd_files[i]);
			}

			++i;
		}

		timing.stop();

		if (!accepted) break;
	}

	// output result of last timestep if not done automatically
	if (!output_timestep)
	{
		const auto timestep = time_stepper.getTimeStep().steps();
		const auto current_time = time_stepper.getTimeStep().current();

		unsigned i = 0;  // process counter, used to distinguish output files
		for (auto p = project.processesBegin(); p != project.processesEnd();
			 ++p)
		{
			do_output(i, *p, timestep, current_time, *pvd_files[i]);
		}
	}
}

int main(int argc, char *argv[])
{
#ifndef NDEBUG
	feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

	// Parse CLI arguments.
	TCLAP::CmdLine cmd("OpenGeoSys-6 software.\n"
			"Copyright (c) 2012-2016, OpenGeoSys Community "
			"(http://www.opengeosys.org) "
			"Distributed under a Modified BSD License. "
			"See accompanying file LICENSE.txt or "
			"http://www.opengeosys.org/project/license",
		' ',
		BaseLib::BuildInfo::git_describe);

	TCLAP::UnlabeledValueArg<std::string> project_arg(
		"project-file",
		"Path to the ogs6 project file.",
		true,
		"",
		"PROJECT FILE");
	cmd.add(project_arg);

	TCLAP::SwitchArg nonfatal_arg("",
		"config-warnings-nonfatal",
		"warnings from parsing the configuration file will not trigger program abortion");
	cmd.add(nonfatal_arg);

	cmd.parse(argc, argv);

	ApplicationsLib::LogogSetup logog_setup;
	ApplicationsLib::LinearSolverLibrarySetup linear_solver_library_setup(
	    argc, argv);

	// Project's configuration
	BaseLib::ConfigTree project_config =
	    BaseLib::read_xml_config(project_arg.getValue());

	std::unique_ptr<ProjectData> project;
	{
		// Nested scope in order to trigger config tree checks early.
		// Caution: The top level config tree must not be saved inside
		//          ProjectData and the boost::property_tree must not be
		//          created inside this same scope!
		using Conf = BaseLib::ConfigTreeNew;
		Conf conf(project_config.get_child("OpenGeoSysProject"),
		          Conf::onerror,
		          nonfatal_arg.getValue() ? Conf::onwarning : Conf::onerror);

		project.reset(new ProjectData(
		              conf, BaseLib::extractPath(project_arg.getValue())));
	}

	// Create processes.
	project->buildProcesses<GlobalSetupType>();

	INFO("Initialize processes.");
	for (auto p_it = project->processesBegin(); p_it != project->processesEnd(); ++p_it)
	{
		(*p_it)->initialize();
	}

	solveProcesses(*project);

	return 0;
}
