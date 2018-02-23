/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DubininPolanyiOrTPD.h"
#include <cmath>
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "BaseLib/Functional.h"
#include "MaterialLib/PhysicalConstant.h"
#include "ReactionEquilibrium.h"

namespace MaterialLib
{
double DubininPolanyiOrTPD::getEquilibriumLoading(const double p_Ads,
                                                  const double T_Ads) const
{
    const double C_eq_Dub = DubininPolanyi::getEquilibriumLoading(p_Ads, T_Ads);

    // Cf. Kraus' Diss. Fig. 38 p. 71
    const double A1 = 0.079, A2 = 27.28, A3 = 388.5, A4 = 33.72;
    const double C_eq_TPD =
        A1 + (A2 - A1) / (1.0 + std::exp((T_Ads - A3) / A4));

    return std::min(C_eq_Dub, C_eq_TPD);
}

std::unique_ptr<DubininPolanyiOrTPD> createDubininPolanyiOrTPD(
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "DubininPolanyiOrTPD");

    auto const rho_SR_dry =
        config.getConfigParameter<double>("adsorbent_density_dry");
    auto data =
        createDubininPolanyiData(config.getConfigSubtree("equilibrium_data"));
    return std::unique_ptr<DubininPolanyiOrTPD>(
        new DubininPolanyiOrTPD(rho_SR_dry, std::move(data)));
}

}  // namespace MaterialLib
