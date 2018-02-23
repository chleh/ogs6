/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ReactiveSolidModel.h"

namespace MaterialLib
{
class ReactiveSolidModelLoadingBased : public ReactiveSolidModel
{
public:
    explicit ReactiveSolidModelLoadingBased(
        const double sorbent_density_dry,
        unsigned num_components,
        ProcessLib::Parameter<double> const& initial_conversion)
        : _sorbent_density_dry{sorbent_density_dry},
          _num_components(num_components),
          _initial_conversion(initial_conversion)
    {
    }

    std::vector<ReactiveSolidStateInternalVariable> getInternalStateVariables()
        const override;

    std::vector<ReactiveSolidRateVariable> getRateVariables() const override;

    std::unique_ptr<ReactiveSolidState> createReactiveSolidState(
        double const t, ProcessLib::SpatialPosition const& pos) const override;

    std::unique_ptr<ReactiveSolidRate> createReactiveSolidRate() const override;

    double getOverallRate(ReactiveSolidRate const& rate) const override;

    double getSolidDensity(ReactiveSolidState const& state) const override;

    double getHeatingRate(HeatOfReactionData const& heat_of_reaction,
                          ReactiveSolidRate const& rate) override;

private:
    const double _sorbent_density_dry;
    unsigned _num_components;
    ProcessLib::Parameter<double> const& _initial_conversion;
};

std::unique_ptr<ReactiveSolidModelLoadingBased>
createReactiveSolidModelLoadingBased(
    BaseLib::ConfigTree const& config,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters);

}  // namespace MaterialLib
