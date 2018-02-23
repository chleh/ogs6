/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/NamedFunctionProvider.h"

namespace MaterialLib
{
namespace Enerchem
{
struct LoadingLoading
{
    double loading_zeolite;
    double loading_salt_solution;
};

struct HeatOfReactionTwoParts
{
    double heat_of_reaction_zeolite;
    double heat_of_reaction_salt_solution;
};

class ReactionEquilibriumTwoPartsLoadingLoadingModel
    : public NumLib::NamedFunctionProvider
{
public:
    virtual LoadingLoading getEquilibriumLoading(const double p,
                                                 const double T) = 0;
    virtual HeatOfReactionTwoParts getHeatsOfReaction(
        const double p, const double T, const LoadingLoading& loadings) = 0;
};

}  // namespace Enerchem
}  // namespace MaterialLib
