/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ProcessOutputDUNE.h"

#if OGS_USE_DUNE

#include <vtk/vtkXMLWriter.h>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include "MeshLib/DUNEMesh.h"
#include "NumLib/DOF/AbstractDOFTable.h"

namespace
{
template <int GlobalDim>
static void doProcessOutputDUNE(
    std::string const& file_name,
    bool const /*compress_output*/,
    int const data_mode,
    const double t,
    GlobalVector const& x,
    MeshLib::DUNEMesh<GlobalDim>& mesh,
    NumLib::AbstractDOFTable const& dof_table,
    std::vector<std::reference_wrapper<ProcessLib::ProcessVariable>> const&
        process_variables,
    ProcessLib::SecondaryVariableCollection secondary_variables,
    ProcessLib::ProcessOutput const& process_output)
{
    Dune::VTK::OutputType output_type;
    switch (data_mode)
    {
        case vtkXMLWriter::Ascii:
            output_type = Dune::VTK::OutputType::ascii;
            break;
        case vtkXMLWriter::Binary:
            output_type = Dune::VTK::OutputType::base64;
            break;
        case vtkXMLWriter::Appended:
            output_type = Dune::VTK::OutputType::appendedbase64;
            break;
    }

    auto const& output_variables = process_output.output_variables;
    std::set<std::string> already_output;

    auto gridView = mesh.getMesh().leafGridView();
    auto const num_nodes = gridView.size(GlobalDim);
    auto const num_elements = gridView.size(0);

    Dune::VTKWriter<decltype(gridView)> vtkWriter(gridView);

    // TODO [DUNE] implement all
    OGS_ALWAYS_ASSERT(process_variables.size() == 1);

    auto const& pv = process_variables.front().get();
    OGS_ALWAYS_ASSERT(output_variables.find(pv.getName()) !=
                      output_variables.cend());

    already_output.insert(pv.getName());

    OGS_ALWAYS_ASSERT(x.size() == pv.getNumberOfComponents() * num_nodes);

    // TODO [DUNE] assumes by-location ordering of x!
    vtkWriter.addVertexData(x, pv.getName(), pv.getNumberOfComponents());

#ifndef USE_PETSC
    std::vector<std::unique_ptr<std::vector<double>>> value_cache;
    auto add_secondary_var = [&](ProcessLib::SecondaryVariable const& var,
                                 std::string const& output_name) {
        {
            DBUG("  secondary variable %s", output_name.c_str());

            auto nodal_values_mesh = std::make_unique<std::vector<double>>(
                num_nodes * var.fcts.num_components);

            std::unique_ptr<GlobalVector> result_cache;
            auto const& nodal_values =
                var.fcts.eval_field(t, x, dof_table, result_cache);
            if (nodal_values_mesh->size() !=
                static_cast<std::size_t>(nodal_values.size()))
            {
                OGS_FATAL(
                    "Secondary variable `%s' did not evaluate to the right "
                    "number of components. Expected: %d, actual: %d.",
                    var.name.c_str(), nodal_values_mesh->size(),
                    nodal_values.size());
            }

            // Copy result
            for (GlobalIndexType i = 0; i < nodal_values.size(); ++i)
            {
                assert(!std::isnan(nodal_values[i]));
                (*nodal_values_mesh)[i] = nodal_values[i];
            }

            value_cache.push_back(std::move(nodal_values_mesh));

            vtkWriter.addVertexData(*value_cache.back(), output_name,
                                    var.fcts.num_components);
        }

        // output extrapolation residuals
        if (process_output.output_residuals && var.fcts.eval_residuals)
        {
            DBUG("  secondary variable %s residual", output_name.c_str());
            auto const& property_name_res = output_name + "_residual";

            auto residuals_mesh = std::make_unique<std::vector<double>>(
                num_elements * var.fcts.num_components);

            std::unique_ptr<GlobalVector> result_cache;
            auto const& residuals =
                var.fcts.eval_residuals(t, x, dof_table, result_cache);
            if (residuals_mesh->size() !=
                static_cast<std::size_t>(residuals.size()))
            {
                OGS_FATAL(
                    "Thee residual of secondary variable `%s' did not evaluate "
                    "to the right number of components. Expected: %d, actual: "
                    "%d.",
                    var.name.c_str(), residuals_mesh->size(), residuals.size());
            }

            // Copy result
            for (GlobalIndexType i = 0; i < residuals.size(); ++i)
            {
                assert(!std::isnan(residuals[i]));
                (*residuals_mesh)[i] = residuals[i];
            }

            value_cache.push_back(std::move(residuals_mesh));

            vtkWriter.addCellData(*value_cache.back(), property_name_res,
                                  var.fcts.num_components);
        }
    };

    for (auto const& external_variable_name : output_variables)
    {
        if (!already_output.insert(external_variable_name).second)
        {
            // no insertion took place, output already done
            continue;
        }

        add_secondary_var(secondary_variables.get(external_variable_name),
                          external_variable_name);
    }
#else
    (void)secondary_variables;
    (void)t;
#endif  // USE_PETSC

    vtkWriter.write(file_name, output_type);
}
}  // anonymous namespace

namespace ProcessLib
{
void doProcessOutputDUNE(
    std::string const& file_name,
    bool const compress_output,
    int const data_mode,
    const double t,
    GlobalVector const& x,
    MeshLib::FEMMesh& fem_mesh,
    NumLib::AbstractDOFTable const& dof_table,
    std::vector<std::reference_wrapper<ProcessVariable>> const&
        process_variables,
    SecondaryVariableCollection secondary_variables,
    ProcessOutput const& process_output)
{
    DBUG("Process output DUNE.");

    if (auto* mesh = dynamic_cast<MeshLib::DUNEMesh<2>*>(&fem_mesh))
    {
        ::doProcessOutputDUNE(file_name, compress_output, data_mode, t, x,
                              *mesh, dof_table, process_variables,
                              secondary_variables, process_output);
    }
    else if (auto* mesh = dynamic_cast<MeshLib::DUNEMesh<3>*>(&fem_mesh))
    {
        ::doProcessOutputDUNE(file_name, compress_output, data_mode, t, x,
                              *mesh, dof_table, process_variables,
                              secondary_variables, process_output);
    }
    else
    {
        OGS_FATAL("unsupported mesh type");
    }
}

}  // namespace ProcessLib

#endif
