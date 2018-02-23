#pragma once

#include "MathLib/ODE/ODESolver.h"
#include "ReactionKinetics.h"
#include "ReactionRate.h"

namespace MaterialLib
{
class ReactionRateRawWithODESolver final : public ReactionRate
{
public:
    ReactionRateRawWithODESolver(
        std::unique_ptr<ReactionKinetics>&& react_kin,
        std::unique_ptr<MathLib::ODE::ODESolver<1>>&& ode_solver,
        bool const only_once_per_ts);

    void preIteration(const unsigned iter) override;

    bool computeReactionRate(double const delta_t, double const p,
                             double const T, double& p_V, ReactiveSolidRate& qR,
                             ReactiveSolidState* solid_state,
                             ReactionRateData* reaction_rate_data) override;

    HeatOfReactionData getHeatOfReaction(
        const double p_V, const double T,
        const ReactiveSolidState* const state) const
    {
        return _react_kin->getHeatOfReaction(p_V, T, state);
    }

    std::vector<NumLib::NamedFunction> getNamedFunctions() const override
    {
        return _react_kin->getNamedFunctions();
    }

    bool isStateCompatible(ReactiveSolidState& state) const override;

    std::unique_ptr<ReactionRateData> createReactionRateData() const override;

private:
    std::unique_ptr<ReactionKinetics> _react_kin;
    unsigned _iter;
    std::unique_ptr<MathLib::ODE::ODESolver<1>> _ode_solver;
    const bool _only_once_per_ts;
};

std::unique_ptr<ReactionRateRawWithODESolver>
createReactionRateRawWithODESolver(BaseLib::ConfigTree const& config);

}  // namespace MaterialLib
