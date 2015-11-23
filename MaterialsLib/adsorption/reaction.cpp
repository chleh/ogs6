/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "logog/include/logog.hpp"

#include "reaction.h"

#include "density_legacy.h"
#include "density_100MPa.h"
#include "density_const.h"
#include "density_cook.h"
#include "density_dubinin.h"
#include "density_hauer.h"
#include "density_mette.h"
#include "density_nunez.h"

#include "reaction_inert.h"
#include "reaction_sinusoidal.h"


namespace Ads
{

std::unique_ptr<Reaction>
Reaction::
newInstance(BaseLib::ConfigTree const& conf)
{
	auto const type = conf.get<std::string>("type", "");

	if (type == "Z13XBF")
		return std::unique_ptr<Reaction>(new DensityLegacy);
	else if (type == "Z13XBF_100MPa")
		return std::unique_ptr<Reaction>(new Density100MPa);
	else if (type == "Z13XBF_Const")
		return std::unique_ptr<Reaction>(new DensityConst);
	else if (type == "Z13XBF_Cook")
		return std::unique_ptr<Reaction>(new DensityCook);
	else if (type == "Z13XBF_Dubinin")
		return std::unique_ptr<Reaction>(new DensityDubinin);
	else if (type == "Z13XBF_Hauer")
		return std::unique_ptr<Reaction>(new DensityHauer);
	else if (type == "Z13XBF_Mette")
		return std::unique_ptr<Reaction>(new DensityMette);
	else if (type == "Z13XBF_Nunez")
		return std::unique_ptr<Reaction>(new DensityNunez);
	else if (type == "Inert")
		return std::unique_ptr<Reaction>(new ReactionInert);
	else if (type == "Sinusoidal")
		return std::unique_ptr<Reaction>(new ReactionSinusoidal(conf));

	if (type.empty()) {
		ERR("No reactive system specified.");
	} else {
		ERR("Unknown reactive system: %s.", type.c_str());
	}
	return nullptr;
}

} // namespace Ads
