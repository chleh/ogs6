/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateSourceTerm.h"

#include "CreateNodalSourceTerm.h"
// TODO [DUNE] fixme
#if 0
#include "CreateVolumetricSourceTerm.h"
#include "SourceTerm.h"
#endif
#include "SourceTermConfig.h"

namespace ProcessLib
{
std::unique_ptr<SourceTerm> createSourceTerm(
    const SourceTermConfig& config, const NumLib::AbstractDOFTable& dof_table,
    const MeshLib::FEMMesh& mesh, const int variable_id,
    const unsigned integration_order, const unsigned shapefunction_order,
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters)
{
    //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__type}
    auto const type = config.config.peekConfigParameter<std::string>("type");

    // TODO [DUNE] re-enable
#if 0
    if (type == "Nodal")
    {
        return ProcessLib::createNodalSourceTerm(
            config.config, config.mesh, dof_table, mesh.getID(), variable_id,
            *config.component_id, parameters);
    }

    if (type == "Volumetric")
    {
        return ProcessLib::createVolumetricSourceTerm(
            config.config, config.mesh, dof_table, parameters,
            integration_order, shapefunction_order, variable_id,
            *config.component_id);
    }
#endif

    OGS_FATAL("Unknown source term type: `%s'.", type.c_str());
}
}  // namespace ProcessLib
