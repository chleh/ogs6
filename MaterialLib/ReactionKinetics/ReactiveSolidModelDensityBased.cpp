/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ReactiveSolidModelDensityBased.h"

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
ReactiveSolidModelDensityBased::createReactiveSolidState(
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
ReactiveSolidModelDensityBased::getInternalStateVariables() const
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
ReactiveSolidModelDensityBased::getRateVariables() const
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
ReactiveSolidModelDensityBased::createReactiveSolidRate() const
{
    auto qR = std::unique_ptr<ReactiveSolidRate>{
        new Eigen::VectorXd{_num_components}};
    qR->setZero();
    return qR;
}

double ReactiveSolidModelDensityBased::getOverallRate(
    ReactiveSolidRate const& rate) const
{
    return rate.sum();
}

double ReactiveSolidModelDensityBased::getSolidDensity(
    const ReactiveSolidState& state) const
{
    return state.conversion().sum();
}

double ReactiveSolidModelDensityBased::getHeatingRate(
    const HeatOfReactionData& heat_of_reaction, const ReactiveSolidRate& rate)
{
    // TODO warning: All existing kinetics have to be fixed!!!
    return heat_of_reaction.dot(rate);
}

std::unique_ptr<ReactiveSolidModelDensityBased>
createReactiveSolidModelDensityBased(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters)
{
    config.checkConfigParameter("type", "DensityBased");

    auto const num_components =
        config.getConfigParameter<unsigned>("solid_state_num_components");
    if (num_components > 2)
        OGS_FATAL("Too many components: %u.", num_components);

    auto const& init_conv_name =
        config.getConfigParameter<std::string>("initial_conversion");
    auto const& init_conv_param = ProcessLib::findParameter<double>(
        init_conv_name, parameters, num_components);

    return std::unique_ptr<ReactiveSolidModelDensityBased>(
        new ReactiveSolidModelDensityBased{num_components, init_conv_param});
}

}  // namespace MaterialLib
