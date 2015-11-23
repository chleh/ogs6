/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "reaction.h"

namespace Ads
{

class ReactionSinusoidal final : public Reaction
{
public:
    double get_enthalpy(const double /*p_Ads*/, const double /*T_Ads*/,
                        const double /*M_Ads*/) const {
        return 0.0;
    }

    double get_reaction_rate(const double /*p_Ads*/, const double /*T_Ads*/, const double /*M_Ads*/,
                             const double /*loading*/) const {
        return 0.0;
    }
};

}
