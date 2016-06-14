/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <logog/include/logog.hpp>

#include "Reaction.h"

#include "Density100MPa.h"
#include "DensityConst.h"
#include "DensityCook.h"
#include "DensityDubinin.h"
#include "DensityHauer.h"
#include "DensityLegacy.h"
#include "DensityMette.h"
#include "DensityNunez.h"

#include "DensityHauerNaYStach.h"

#include "ReactionCaOH2.h"
#include "ReactionInert.h"
#include "ReactionSinusoidal.h"


namespace Adsorption
{

std::unique_ptr<Reaction>
Reaction::
newInstance(BaseLib::ConfigTree const& conf)
{
    //! \ogs_file_param{materials__adsorption__reaction__type}
    auto const type = conf.getConfigParameter<std::string>("type");

    if (type.find("Z13XBF") == 0)
    {
        auto const k_rate = conf.getConfigParameter<double>("k_rate");

        if (type == "Z13XBF")
            return std::unique_ptr<Reaction>(new DensityLegacy(k_rate));
        else if (type == "Z13XBF_100MPa")
            return std::unique_ptr<Reaction>(new Density100MPa(k_rate));
        else if (type == "Z13XBF_Const")
            return std::unique_ptr<Reaction>(new DensityConst(k_rate));
        else if (type == "Z13XBF_Cook")
            return std::unique_ptr<Reaction>(new DensityCook(k_rate));
        else if (type == "Z13XBF_Dubinin")
            return std::unique_ptr<Reaction>(new DensityDubinin(k_rate));
        else if (type == "Z13XBF_Hauer")
            return std::unique_ptr<Reaction>(new DensityHauer(k_rate));
        else if (type == "Z13XBF_Mette")
            return std::unique_ptr<Reaction>(new DensityMette(k_rate));
        else if (type == "Z13XBF_Nunez")
            return std::unique_ptr<Reaction>(new DensityNunez(k_rate));
    }
    else if (type == "Inert")
        return std::unique_ptr<Reaction>(new ReactionInert);
    else if (type == "Sinusoidal")
        return std::unique_ptr<Reaction>(new ReactionSinusoidal(conf));
    else if (type == "CaOH2")
        return std::unique_ptr<Reaction>(new ReactionCaOH2(conf));
    else if (type == "NaYStach")
        return std::unique_ptr<Reaction>(new DensityHauerNaYStach(
            conf.getConfigParameter<double>("k_rate")));

    ERR("Unknown reactive system: %s.", type.c_str());
    std::abort();

    return nullptr;
}

} // namespace Adsorption
