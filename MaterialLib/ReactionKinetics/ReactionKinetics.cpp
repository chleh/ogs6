/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ReactionKinetics.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "LinearDrivingForceConstantCoefficient.h"
#include "LinearDrivingForceGorbachCoefficient.h"
#include "LinearDrivingForceMetteCoefficient.h"

#include "Enerchem/LDFTwoComponentTest2.h"
#include "LinearDrivingForceEnerchemTest1.h"

#include "Monod.h"
#include "ReactionCaOH2.h"
#include "ReactionInert.h"
#include "ReactionMnO.h"

namespace MaterialLib
{
std::unique_ptr<ReactionKinetics> createReactionKinetics(
    BaseLib::ConfigTree const& config)
{
    auto const type = config.peekConfigParameter<std::string>("type");

    if (type == "LinearDrivingForceConstantCoefficient")
        return createLinearDrivingForceConstantCoefficient(config);
    // TODO, this is a quick hack! fix it properly!
    else if (type == "LinearDrivingForceCoefficientKast")
        return createLinearDrivingForceMetteCoefficient(config);
    else if (type == "LinearDrivingForceGorbachCoefficient")
        return createLinearDrivingForceGorbachCoefficient(config);
    else if (type == "LinearDrivingForceEnerchemTest1")
        return createLinearDrivingForceEnerchemTest1(config);
    else if (type == "LDFTwoComponentTest2")
        return createLDFTwoComponentTest2(config);
    else if (type == "Monod")
        return createMonodKinetics(config);
    else if (type == "Inert")
    {
        config.checkConfigParameter("type", "Inert");
        return std::unique_ptr<ReactionKinetics>(new ReactionInert);
    }
    else if (type == "CaOH2_Schaube")
        return createReactionCaOH2(config);
    else if (type == "MnO")
        return createReactionMnO(config);

    OGS_FATAL("There is no reaction kinetics of type `%s'.", type.c_str());
}

}  // namespace MaterialLib
