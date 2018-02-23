/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MaterialLib/ReactionKinetics/ReactionRateData.h"
#include "MaterialLib/ReactionKinetics/ReactiveSolidState.h"
#include "TESAssemblyParams.h"

namespace ProcessLib
{
namespace TES
{
struct TESLocalAssemblerData
{
    TESLocalAssemblerData(AssemblyParams const& ap_,
                          const std::size_t element_id_,
                          const unsigned num_int_pts, const unsigned dimension);

    ~TESLocalAssemblerData() = default;

    AssemblyParams const& ap;

    std::size_t const element_id;

    // integration point quantities
    std::vector<std::unique_ptr<MaterialLib::ReactiveSolidState>>
        reactive_solid_state;
    std::vector<std::unique_ptr<MaterialLib::ReactionRateData>>
        reaction_rate_data;
    std::vector<std::unique_ptr<MaterialLib::ReactiveSolidRate>>
        reaction_rate;  // \hat{\rho}_{SR}

    // integration point values of unknowns -- temporary storage
    double p = std::numeric_limits<double>::quiet_NaN();  // gas pressure
    double T = std::numeric_limits<double>::quiet_NaN();  // temperature
    double vapour_mass_fraction =
        std::numeric_limits<double>::quiet_NaN();  // fluid mass fraction of the
                                                   // second component

    // temporary storage for some properties
    // values do not change during the assembly of one integration point
    double rho_GR = std::numeric_limits<double>::quiet_NaN();
    double p_V =
        std::numeric_limits<double>::quiet_NaN();  // vapour partial pressure

    double poro = std::numeric_limits<double>::quiet_NaN();
};
}
}  // namespaces
