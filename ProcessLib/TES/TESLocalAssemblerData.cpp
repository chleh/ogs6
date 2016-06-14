/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "TESLocalAssemblerData.h"
#include "TESReactionAdaptor.h"


// TODO move somwhere else

std::vector<double>
initializeIntegrationPointValues(
        std::size_t number_of_integration_points,
        ProcessLib::Parameter<double, MeshLib::Element const&> const& parameter,
        MeshLib::Element const& e)
{
    std::vector<double> integration_point_values;

    // set initial solid density from mesh property
    if (!_assembly_params.initial_solid_density_mesh_property.empty())
    {
        auto prop = mesh.getProperties().template getPropertyVector<double>(
            _assembly_params.initial_solid_density_mesh_property);
        assert(prop->getNumberOfComponents() == 1);

        switch (prop->getMeshItemType())
        {
            case MeshLib::MeshItemType::Cell:
            {
                auto init_solid_density = [&prop](std::size_t id,
                                                  LocalAssembler& loc_asm) {
                    // TODO loc_asm_id is assumed to be the mesh element id.
                    loc_asm.initializeSolidDensity(MeshLib::MeshItemType::Cell,
                                                   {{(*prop)[id]}});
                };

                GlobalSetup::executeDereferenced(init_solid_density,
                                                 _local_assemblers);
                break;
            }
            case MeshLib::MeshItemType::Node:
            {
                std::vector<GlobalIndexType> indices;
                std::vector<double> values;

                auto init_solid_density = [&](std::size_t id,
                                              LocalAssembler& loc_asm) {
                    getRowColumnIndices_(
                        id, *_local_to_global_index_map_single_component,
                        indices);
                    values.clear();
                    for (auto i : indices)
                        values.push_back((*prop)[i]);

                    loc_asm.initializeSolidDensity(MeshLib::MeshItemType::Node,
                                                   values);
                };

                GlobalSetup::executeDereferenced(init_solid_density,
                                                 _local_assemblers);
                break;
            }
            default:
                ERR("Unhandled mesh item type for initialization of secondary "
                    "variable.");
                std::abort();
        }
    }
}


namespace ProcessLib
{
namespace TES
{
TESLocalAssemblerData::TESLocalAssemblerData(MeshLib::Element const& e,
                                             AssemblyParams const& ap_,
                                             const unsigned num_int_pts,
                                             const unsigned dimension)
    : ap(ap_),
      solid_density(initializeIntegrationPointValues(
          num_int_pts, ap_.initial_solid_density, e)),
      reaction_rate(num_int_pts),
      velocity(dimension, std::vector<double>(num_int_pts)),
      reaction_adaptor(TESFEMReactionAdaptor::newInstance(*this)),
      solid_density_prev_ts(solid_density),
      reaction_rate_prev_ts(num_int_pts)
{
}

TESLocalAssemblerData::~TESLocalAssemblerData() = default;
}
}  // namespaces
