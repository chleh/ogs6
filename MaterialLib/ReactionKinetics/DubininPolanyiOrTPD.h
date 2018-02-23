/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALLIB_DUBININPOLANYIORTPD_H
#define MATERIALLIB_DUBININPOLANYIORTPD_H

#include "DubininPolanyi.h"

namespace MaterialLib
{
class DubininPolanyiOrTPD final : public DubininPolanyi
{
public:
    using DubininPolanyi::DubininPolanyi;

    double getEquilibriumLoading(const double p_Ads,
                                 const double T_Ads) const override;
};

std::unique_ptr<DubininPolanyiOrTPD> createDubininPolanyiOrTPD(
    BaseLib::ConfigTree const& config);

}  // namespace MaterialLib

#endif  // MATERIALLIB_DUBININPOLANYIORTPD_H
