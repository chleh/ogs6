/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "TESLocalAssemblerData.h"

#include "BaseLib/Error.h"

namespace ProcessLib
{
namespace TES
{
TESLocalAssemblerData::TESLocalAssemblerData(AssemblyParams const& ap_,
                                             const std::size_t element_id_,
                                             const unsigned num_int_pts,
                                             const unsigned /*dimension*/)
    : ap(ap_),
      element_id(element_id_),
      reactive_solid_state(num_int_pts),
      reaction_rate_data(num_int_pts),
      reaction_rate(num_int_pts)
{
    SpatialPosition pos;
    pos.setElementID(element_id);

    for (unsigned ip = 0; ip < reactive_solid_state.size(); ++ip)
    {
        pos.setIntegrationPoint(ip);
        auto& d = reactive_solid_state[ip];
        // TODO warning: 0.0 is the time!
        d = ap.reactive_solid->createReactiveSolidState(0.0, pos);
        if (!ap.reaction_rate->isStateCompatible(*d))
            OGS_FATAL(
                "reaction rate and reactive solid state are incompatible.");
    }

    for (auto& d : reaction_rate)
    {
        d = ap.reactive_solid->createReactiveSolidRate();
    }

    for (auto& d : reaction_rate_data)
    {
        d = ap.reaction_rate->createReactionRateData();
    }
}

}  // namespace TES
}  // namespace ProcessLib
