/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <functional>
#include <memory>

#include "NumLib/NamedFunctionProvider.h"
#include "ProcessLib/Parameter/Parameter.h"

#include "ReactionHeatModel.h"
#include "ReactionRate.h"
#include "ReactiveSolidRate.h"

namespace BaseLib
{
class ConfigTree;
}  // namespace BaseLib

namespace MaterialLib
{
struct ReactiveSolidStateInternalVariable
{
    using Fct = std::function<std::vector<double> const&(
        ReactiveSolidState const& /*solid_state*/, std::vector<double>&
        /*cache*/)>;

    ReactiveSolidStateInternalVariable(std::string const& name_,
                                       unsigned const dim_, Fct&& eval_)
        : name{name_}, dim{dim_}, eval{std::move(eval_)}
    {
    }

    std::string name;
    unsigned dim;
    Fct eval;
};

struct ReactiveSolidRateVariable
{
    using Fct = std::function<std::vector<double> const&(
        ReactiveSolidRate const& /*solid_rate*/, std::vector<double>&
        /*cache*/)>;

    ReactiveSolidRateVariable(std::string const& name_, unsigned const dim_,
                              Fct&& eval_)
        : name{name_}, dim{dim_}, eval{std::move(eval_)}
    {
    }

    std::string name;
    unsigned dim;
    Fct eval;
};

class ReactiveSolidModel : public NumLib::NamedFunctionProvider
{
public:
    virtual void preIteration(const unsigned iter) { (void)iter; }
    virtual void postIteration() {}
    virtual void preTimestep(ReactiveSolidState& solid_state,
                             ReactiveSolidRate const& reaction_rate)
    {
        (void)solid_state;
        (void)reaction_rate;
    }

    virtual std::vector<ReactiveSolidStateInternalVariable>
    getInternalStateVariables() const = 0;

    virtual std::vector<ReactiveSolidRateVariable> getRateVariables() const = 0;

    virtual std::unique_ptr<ReactiveSolidState> createReactiveSolidState(
        double const t, ProcessLib::SpatialPosition const& pos) const = 0;

    virtual std::unique_ptr<ReactiveSolidRate> createReactiveSolidRate()
        const = 0;

    virtual double getOverallRate(ReactiveSolidRate const& rate) const = 0;

    virtual double getSolidDensity(ReactiveSolidState const& state) const = 0;

    virtual double getHeatingRate(HeatOfReactionData const& heat_of_reaction,
                                  ReactiveSolidRate const& rate) = 0;
};

std::unique_ptr<ReactiveSolidModel> createReactiveSolidModel(
    BaseLib::ConfigTree const& config,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters);

}  // namespace MaterialLib
