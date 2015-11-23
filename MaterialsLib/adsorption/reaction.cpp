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

bool
stringToReactiveSystem(std::string const& name, SolidReactiveSystem& rsys_out)
{
	if (name == "Z13XBF")
		rsys_out = SolidReactiveSystem::Z13XBF;
	else if (name == "Z13XBF_100MPa")
		rsys_out = SolidReactiveSystem::Z13XBF_100MPa;
	else if (name == "Z13XBF_Const")
		rsys_out = SolidReactiveSystem::Z13XBF_Const;
	else if (name == "Z13XBF_Cook")
		rsys_out = SolidReactiveSystem::Z13XBF_Cook;
	else if (name == "Z13XBF_Dubinin")
		rsys_out = SolidReactiveSystem::Z13XBF_Dubinin;
	else if (name == "Z13XBF_Hauer")
		rsys_out = SolidReactiveSystem::Z13XBF_Hauer;
	else if (name == "Z13XBF_Mette")
		rsys_out = SolidReactiveSystem::Z13XBF_Mette;
	else if (name == "Z13XBF_Nunez")
		rsys_out = SolidReactiveSystem::Z13XBF_Nunez;
	else if (name == "Inert")
		rsys_out = SolidReactiveSystem::Inert;
	else if (name == "Sinusoidal")
		rsys_out = SolidReactiveSystem::Sinusoidal;
	else return false;

	return true;
}

std::unique_ptr<Reaction>
Reaction::
newInstance(std::string const& rsys)
{
	SolidReactiveSystem r;
	if (stringToReactiveSystem(rsys, r)) {
		return newInstance(r);
	} else {
		ERR("unknown reactive system: %s", rsys.c_str());
		return nullptr;
	}
}

std::unique_ptr<Reaction>
Reaction::
newInstance(const SolidReactiveSystem rsys)
{
	switch (rsys)
	{
	case SolidReactiveSystem::Z13XBF:
		return std::unique_ptr<Reaction>(new DensityLegacy);
	case SolidReactiveSystem::Z13XBF_100MPa:
		return std::unique_ptr<Reaction>(new Density100MPa);
	case SolidReactiveSystem::Z13XBF_Const:
		return std::unique_ptr<Reaction>(new DensityConst);
	case SolidReactiveSystem::Z13XBF_Cook:
		return std::unique_ptr<Reaction>(new DensityCook);
	case SolidReactiveSystem::Z13XBF_Dubinin:
		return std::unique_ptr<Reaction>(new DensityDubinin);
	case SolidReactiveSystem::Z13XBF_Hauer:
		return std::unique_ptr<Reaction>(new DensityHauer);
	case SolidReactiveSystem::Z13XBF_Mette:
		return std::unique_ptr<Reaction>(new DensityMette);
	case SolidReactiveSystem::Z13XBF_Nunez:
		return std::unique_ptr<Reaction>(new DensityNunez);
	case SolidReactiveSystem::Inert:
		return std::unique_ptr<Reaction>(new ReactionInert);
	case SolidReactiveSystem::Sinusoidal:
		return std::unique_ptr<Reaction>(new ReactionSinusoidal);
	}

	return nullptr;
}

} // namespace Ads
