/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALLIB_REACTIONEQUILIBRIUM_H
#define MATERIALLIB_REACTIONEQUILIBRIUM_H

#include <memory>

#include "BaseLib/Error.h"
#include "NumLib/NamedFunctionProvider.h"

#include "ReactionHeatModel.h"
#include "ReactionKinetics.h"

namespace MaterialLib
{
class ReactionEquilibrium : public NumLib::NamedFunctionProvider,
                            public ReactionHeatModel
{
public:
    virtual double getEquilibriumDensity(const double p, const double T) = 0;

    virtual double getDEquilibriumDensityDp(const double p, const double T) = 0;
    virtual double getDEquilibriumDensityDT(const double /*p*/,
                                            const double /*T*/)
    {
        OGS_FATAL("not implemented");
    }
};

std::unique_ptr<ReactionEquilibrium> createReactionEquilibrium(
    BaseLib::ConfigTree const& config);

}  // namespace MaterialLib

#endif  // MATERIALLIB_REACTIONEQUILIBRIUM_H
