/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ReactiveSolidModelLoadingBased.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "BaseLib/Functional.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "ProcessLib/Parameter/SpatialPosition.h"
#include "ProcessLib/Utils/ProcessUtils.h"

#include "ReactiveSolidStateOneComponent.h"
#include "ReactiveSolidStateTwoComponents.h"

namespace MaterialLib
{
std::unique_ptr<ReactiveSolidState>
ReactiveSolidModelLoadingBased::createReactiveSolidState(
    double const t, ProcessLib::SpatialPosition const& pos) const
{
    if (_num_components == 1)
    {
        return std::unique_ptr<ReactiveSolidStateOneComponent>{
            new ReactiveSolidStateOneComponent{_initial_conversion(t, pos)}};
    }
    else if (_num_components == 2)
    {
        return std::unique_ptr<ReactiveSolidStateTwoComponents>{
            new ReactiveSolidStateTwoComponents{_initial_conversion(t, pos)}};
    }
    else
    {
        OGS_FATAL("Unsupported number of components: %u.", _num_components);
    }
}

std::vector<ReactiveSolidStateInternalVariable>
ReactiveSolidModelLoadingBased::getInternalStateVariables() const
{
    std::vector<ReactiveSolidStateInternalVariable> vars;

    auto f = [](ReactiveSolidState const& solid_state,
                std::vector<double>& cache) -> std::vector<double> const& {
        auto const& conv = solid_state.conversion();
        cache.resize(conv.size());
        MathLib::toVector(cache) = conv;
        return cache;
    };

    vars.emplace_back(
        "conversion", _num_components, BaseLib::easyBind(std::move(f)));

    return vars;
}

std::vector<ReactiveSolidRateVariable>
ReactiveSolidModelLoadingBased::getRateVariables() const
{
    std::vector<ReactiveSolidRateVariable> vars;

    auto f = [](ReactiveSolidRate const& rate,
                std::vector<double>& cache) -> std::vector<double> const& {
        cache.resize(rate.size());
        MathLib::toVector(cache) = rate;
        return cache;
    };

    vars.emplace_back(
        "conversion_rate", _num_components, BaseLib::easyBind(std::move(f)));

    return vars;
}

std::unique_ptr<ReactiveSolidRate>
ReactiveSolidModelLoadingBased::createReactiveSolidRate() const
{
    auto qR = std::unique_ptr<ReactiveSolidRate>{
        new Eigen::VectorXd{_num_components}};
    qR->setZero();
    return qR;
}

double ReactiveSolidModelLoadingBased::getOverallRate(
    ReactiveSolidRate const& rate) const
{
    // TODO warning: All existing kinetics have to be fixed!!!
    return rate.sum() * _sorbent_density_dry;
}

double ReactiveSolidModelLoadingBased::getSolidDensity(
    const ReactiveSolidState& state) const
{
    // TODO warning: All existing kinetics have to be fixed!!!
    // TODO adapt for different kinds of state
    return (1.0 + state.conversion().sum()) * _sorbent_density_dry;
}

double ReactiveSolidModelLoadingBased::getHeatingRate(
    const HeatOfReactionData& heat_of_reaction, const ReactiveSolidRate& rate)
{
    // TODO warning: All existing kinetics have to be fixed!!!
    return _sorbent_density_dry * heat_of_reaction.dot(rate);
}

std::unique_ptr<ReactiveSolidModelLoadingBased>
createReactiveSolidModelLoadingBased(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters)
{
    config.checkConfigParameter("type", "LoadingBased");

    auto const num_components =
        config.getConfigParameter<unsigned>("solid_state_num_components");
    if (num_components > 2)
        OGS_FATAL("Too many components: %u.", num_components);

    auto const& init_conv_name =
        config.getConfigParameter<std::string>("initial_conversion");
    auto const& init_conv_param = ProcessLib::findParameter<double>(
        init_conv_name, parameters, num_components);

    auto const sorbent_density_dry =
        config.getConfigParameter<double>("sorbent_density_dry");

    return std::unique_ptr<ReactiveSolidModelLoadingBased>(
        new ReactiveSolidModelLoadingBased{sorbent_density_dry, num_components,
                                           init_conv_param});
}

}  // namespace MaterialLib
