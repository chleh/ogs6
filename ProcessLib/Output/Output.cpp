/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Output.h"

#include <cassert>
#include <fstream>
#include <vector>

#include <logog/include/logog.hpp>

#include "Applications/InSituLib/Adaptor.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/RunTime.h"
#include "ProcessLib/Process.h"

namespace
{
//! Determines if there should be output at the given \c timestep.
template <typename CountsSteps>
bool shallDoOutput(unsigned timestep, CountsSteps const& repeats_each_steps)
{
    unsigned each_steps = 1;

    for (auto const& pair : repeats_each_steps)
    {
        each_steps = pair.each_steps;

        if (timestep > pair.repeat * each_steps)
        {
            timestep -= pair.repeat * each_steps;
        }
        else
        {
            break;
        }
    }

    return timestep % each_steps == 0;
}

//! Converts a vtkXMLWriter's data mode string to an int. See
/// Output::_output_file_data_mode.
int convertVtkDataMode(std::string const& data_mode)
{
    if (data_mode == "Ascii")
    {
        return 0;
    }
    if (data_mode == "Binary")
    {
        return 1;
    }
    if (data_mode == "Appended")
    {
        return 2;
    }
    OGS_FATAL(
        "Unsupported vtk output file data mode \"%s\". Expected Ascii, "
        "Binary, or Appended.",
        data_mode.c_str());
}

double getTNonlinearIteration(double t, double t_prev, unsigned iteration)
{
    auto const dt = t - t_prev;
    // four digits are reserved to represent the iteration number
    return t_prev + iteration * std::pow(10, std::floor(std::log10(dt)) - 4);
}
}  // namespace

namespace ProcessLib
{
Output::Output(std::string output_directory, std::string prefix,
               bool const compress_output, std::string const& data_mode,
               bool const output_nonlinear_iteration_results,
               std::vector<PairRepeatEachSteps> repeats_each_steps)
    : _output_directory(std::move(output_directory)),
      _output_file_prefix(std::move(prefix)),
      _output_file_compression(compress_output),
      _output_file_data_mode(convertVtkDataMode(data_mode)),
      _output_nonlinear_iteration_results(output_nonlinear_iteration_results),
      _repeats_each_steps(std::move(repeats_each_steps))
{
}

void Output::addProcess(ProcessLib::Process const& process,
                        const int process_id)
{
    auto const filename = BaseLib::joinPaths(
        _output_directory,
        _output_file_prefix + "_pcs_" + std::to_string(process_id) + ".pvd");
    auto const filename_nonlinear_iterations =
        BaseLib::joinPaths(_output_directory,
                           _output_file_prefix + "_pcs_" +
                               std::to_string(process_id) + "_nliter.pvd");
    _process_to_process_data.emplace(
        std::piecewise_construct,
        std::forward_as_tuple(&process),
        std::forward_as_tuple(filename, _output_nonlinear_iteration_results
                                            ? &filename_nonlinear_iterations
                                            : nullptr));
}

// TODO return a reference.
Output::ProcessData* Output::findProcessData(Process const& process,
                                             const int process_id)
{
    auto range = _process_to_process_data.equal_range(&process);
    int counter = 0;
    ProcessData* process_data = nullptr;
    for (auto spd_it = range.first; spd_it != range.second; ++spd_it)
    {
        if (counter == process_id)
        {
            process_data = &spd_it->second;
            break;
        }
        counter++;
    }
    if (process_data == nullptr)
    {
        OGS_FATAL(
            "The given process is not contained in the output"
            " configuration. Aborting.");
    }

    return process_data;
}

void Output::doOutputAlways(Process const& process,
                            const int process_id,
                            ProcessOutput const& process_output,
                            unsigned timestep,
                            const double t,
                            GlobalVector const& x)
{
    updateTPrev(t);
    BaseLib::RunTime time_output;
    time_output.start();

    process.preOutput(x, t, process_id);

    // Need to add variables of process to vtu even no output takes place.
    processOutputData(t, x, process.getMesh(), process.getDOFTable(process_id),
                      process.getProcessVariables(process_id),
                      process.getSecondaryVariables(), process_output);

    // For the staggered scheme for the coupling, only the last process, which
    // gives the latest solution within a coupling loop, is allowed to make
    // output.
    if (!(process_id == static_cast<int>(_process_to_process_data.size()) - 1 ||
          process.isMonolithicSchemeUsed()))
        return;

    std::string const output_file_name =
        _output_file_prefix + "_pcs_" + std::to_string(process_id) + "_ts_" +
        std::to_string(timestep) + "_t_" + std::to_string(t) + ".vtu";
    std::string const output_file_path =
        BaseLib::joinPaths(_output_directory, output_file_name);

    DBUG("output to %s", output_file_path.c_str());

    makeOutput(output_file_path, process.getMesh(), _output_file_compression,
               _output_file_data_mode);

    ProcessData* process_data = findProcessData(process, process_id);
    process_data->pvd_file.addVTUFile(output_file_name, t);
    if (_output_nonlinear_iteration_results)
        process_data->pvd_file_nonlinear_iterations->addVTUFile(
            output_file_name, t);

    INFO("[time] Output of timestep %d took %g s.", timestep,
         time_output.elapsed());
}

void Output::doOutput(Process const& process,
                      const int process_id,
                      ProcessOutput const& process_output,
                      unsigned timestep,
                      const double t,
                      GlobalVector const& x)
{
    updateTPrev(t);
    if (shallDoOutput(timestep, _repeats_each_steps))
    {
        doOutputAlways(process, process_id, process_output, timestep, t, x);
    }
#ifdef USE_INSITU
    // Note: last time step may be output twice: here and in
    // doOutputLastTimestep() which throws a warning.
    InSituLib::CoProcess(process.getMesh(), t, timestep, false);
#endif
}

void Output::doOutputLastTimestep(Process const& process,
                                  const int process_id,
                                  ProcessOutput const& process_output,
                                  unsigned timestep,
                                  const double t,
                                  GlobalVector const& x)
{
    updateTPrev(t);
    if (!shallDoOutput(timestep, _repeats_each_steps))
    {
        doOutputAlways(process, process_id, process_output, timestep, t, x);
    }
#ifdef USE_INSITU
    InSituLib::CoProcess(process.getMesh(), t, timestep, true);
#endif
}

void Output::doOutputNonlinearIteration(Process const& process,
                                        const int process_id,
                                        ProcessOutput const& process_output,
                                        const unsigned timestep, const double t,
                                        GlobalVector const& x,
                                        const unsigned iteration)
{
    updateTPrev(t);
    if (!_output_nonlinear_iteration_results)
    {
        return;
    }

    BaseLib::RunTime time_output;
    time_output.start();

    process.preOutput(x, t, process_id);

    processOutputData(t, x, process.getMesh(), process.getDOFTable(process_id),
                      process.getProcessVariables(process_id),
                      process.getSecondaryVariables(), process_output);

    // For the staggered scheme for the coupling, only the last process, which
    // gives the latest solution within a coupling loop, is allowed to make
    // output.
    if (!(process_id == static_cast<int>(_process_to_process_data.size()) - 1 ||
          process.isMonolithicSchemeUsed()))
        return;

    // Only check whether a process data is available for output.
    auto* process_data = findProcessData(process, process_id);

    std::string const output_file_name =
        _output_file_prefix + "_pcs_" + std::to_string(process_id) + "_ts_" +
        std::to_string(timestep) + "_t_" + std::to_string(t) + "_nliter_" +
        std::to_string(iteration) + ".vtu";
    std::string const output_file_path =
        BaseLib::joinPaths(_output_directory, output_file_name);

    process_data->pvd_file_nonlinear_iterations->addVTUFile(
        output_file_name, getTNonlinearIteration(t, _t_prev, iteration));

    DBUG("output iteration results to %s", output_file_path.c_str());

    makeOutput(output_file_path, process.getMesh(), _output_file_compression,
               _output_file_data_mode);

    INFO("[time] Output took %g s.", time_output.elapsed());
}

void Output::updateTPrev(double t)
{
    if (t != _t_last_call)
    {
        _t_prev = _t_last_call;
    }

    _t_last_call = t;
}

}  // namespace ProcessLib
