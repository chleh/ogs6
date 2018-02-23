#include "ReactionRateRawWithODESolver.h"

#include "BaseLib/ConfigTree.h"
#include "MathLib/ODE/ODESolverBuilder.h"

#include "ReactionRateDataRawWithODESolver.h"
#include "ReactiveSolidStateOneComponent.h"

namespace MaterialLib
{
ReactionRateRawWithODESolver::ReactionRateRawWithODESolver(
    std::unique_ptr<ReactionKinetics>&& react_kin,
    std::unique_ptr<MathLib::ODE::ODESolver<1>>&& ode_solver,
    bool const only_once_per_ts)
    : _react_kin(std::move(react_kin)),
      _ode_solver(std::move(ode_solver)),
      _only_once_per_ts(only_once_per_ts)
{
    // TODO: make configurable
    const double abs_tol = 1e-10;
    const double rel_tol = 1e-10;

    _ode_solver->setTolerance(abs_tol, rel_tol);
}

bool ReactionRateRawWithODESolver::computeReactionRate(
    const double delta_t, double const p, const double T, double& p_V,
    ReactiveSolidRate& qR, ReactiveSolidState* solid_state,
    ReactionRateData* reaction_rate_data)
{
    auto& rhoSR = solid_state->conversion();
    assert(rhoSR.size() == 1);

    if (_only_once_per_ts && _iter != 1)
        return false;

    auto f_lambda = [p_V, p, T, this](
                        double /*t*/,
                        MathLib::ODE::MappedConstVector<1> const& y,
                        MathLib::ODE::MappedVector<1>& ydot) {
        ReactiveSolidStateOneComponent s{{y[0]}};
        auto const qR = _react_kin->getReactionRate(p_V, p, T, s);
        assert(qR.size() == 1);
        ydot[0] = qR[0];
        return true;
    };

    auto* d =
        dynamic_cast<ReactionRateDataRawWithODESolver*>(reaction_rate_data);
    if (!d)
        OGS_FATAL("data is not of type ReactionRateDataRawWithODESolver*.");
    const double rhoSR_prev_ts = d->solid_density_prev_ts;

    _ode_solver->setFunction(f_lambda, nullptr);
    _ode_solver->setIC(0.0, {rhoSR_prev_ts});
    _ode_solver->preSolve();
    if (!_ode_solver->solve(delta_t))
        OGS_FATAL("The ODE solver failed.");
    rhoSR[0] = _ode_solver->getSolution()[0];
    qR = _react_kin->getReactionRate(p_V, p, T, *solid_state);

    return false;
}

void ReactionRateRawWithODESolver::preIteration(const unsigned iter)
{
    _iter = iter;
}

std::unique_ptr<ReactionRateData>
ReactionRateRawWithODESolver::createReactionRateData() const
{
    return std::unique_ptr<ReactionRateData>(
        new ReactionRateDataRawWithODESolver);
}

bool ReactionRateRawWithODESolver::isStateCompatible(
    ReactiveSolidState& state) const
{
    return dynamic_cast<ReactiveSolidStateOneComponent*>(&state) != nullptr;
}

std::unique_ptr<ReactionRateRawWithODESolver>
createReactionRateRawWithODESolver(BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "ReactionRateRawWithODESolver");

    auto only_once_per_ts =
        config.getConfigParameter<bool>("compute_only_once_per_timestep");

    auto react_kin =
        createReactionKinetics(config.getConfigSubtree("reaction_kinetics"));

    auto ode_solver =
        MathLib::ODE::createODESolver<1>(config.getConfigSubtree("ode_solver"));

    return std::unique_ptr<ReactionRateRawWithODESolver>(
        new ReactionRateRawWithODESolver(
            std::move(react_kin), std::move(ode_solver), only_once_per_ts));
}

}  // namespace MaterialLib
