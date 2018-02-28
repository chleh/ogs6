/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include <logog/include/logog.hpp>

#include "BaseLib/Functional.h"

#include "MaterialLib/PhysicalConstant.h"
#include "ProcessLib/TES/Material/DiffusionCoefficient.h"

#include "ReactionEquilibrium.h"
#include "ReactionKinetics.h"

namespace MaterialLib
{
class LinearDrivingForceMetteCoefficient final : public ReactionKinetics
{
public:
    LinearDrivingForceMetteCoefficient(
        const double pellet_radius, const double pellet_porosity,
        const double pellet_tortuosity,
        const double molar_mass_reactive_component,
        std::unique_ptr<ProcessLib::TES::DiffusionCoefficient>&&
            diffusion_coefficient,
        std::unique_ptr<ReactionEquilibrium>&& equil)
        : _prefactor(15.0 / pellet_radius / pellet_radius * pellet_porosity /
                     pellet_tortuosity * molar_mass_reactive_component /
                     PhysicalConstant::IdealGasConstant),
          _diffusion_coefficient(std::move(diffusion_coefficient)),
          _equil(std::move(equil))
    {
    }

    ReactiveSolidRate getReactionRate(
        const double p_V, const double p, const double T,
        ReactiveSolidState const& solid_state) override
    {
        assert(solid_state.conversion().size() == 1);
        auto const rho_SR = solid_state.conversion()[0];

        auto const D = _diffusion_coefficient->getDiffusionCoefficient(p, T);
        auto drhodp = _equil->getDEquilibriumDensityDp(p_V, T);

        // Limit isotherm gradient^{-1} to maximum observed for Ca(90)-X.
        // For Ca(90)-X (dC/dp_V)^{-1} < 10^5.
        if (drhodp < 0.01)
            drhodp = 0.01;  // ==> 1/drhodp < 100 always

        // INFO("drhodp = %g, p_V = %g, T = %g", drhodp, p_V, T);
        auto const k_LDF = _prefactor * D / T / drhodp;
        return Eigen::VectorXd::Constant(
            1, k_LDF * (_equil->getEquilibriumDensity(p_V, T) - rho_SR));
    }

    HeatOfReactionData getHeatOfReaction(
        const double p_V, const double T,
        const ReactiveSolidState* const state) const
    {
        return _equil->getHeatOfReaction(p_V, T, state);
    }

    std::vector<NumLib::NamedFunction> getNamedFunctions() const override
    {
        auto fcts = _equil->getNamedFunctions();
        auto cb = [&](const double p_V, const double p,
                      const double T) -> double {
            auto const D =
                _diffusion_coefficient->getDiffusionCoefficient(p, T);
            auto const drhodp = _equil->getDEquilibriumDensityDp(p_V, T);
            INFO("drhodp = %g, p_V = %g, T = %g", drhodp, p_V, T);
            auto const k_LDF = _prefactor * D / T / drhodp;
            return k_LDF;
        };

        fcts.emplace_back(NumLib::NamedFunction{
            "LDF_prefactor",
            {"vapour_partial_pressure", "pressure", "temperature"},
            BaseLib::easyBind(std::move(cb))});

        return fcts;
    }

    bool isStateCompatible(ReactiveSolidState& state) const override
    {
        return state.conversion().size() == 1;
    }

private:
    const double _prefactor;
    std::unique_ptr<ProcessLib::TES::DiffusionCoefficient> const
        _diffusion_coefficient;
    std::unique_ptr<ReactionEquilibrium> const _equil;
};

std::unique_ptr<ReactionKinetics> createLinearDrivingForceMetteCoefficient(
    BaseLib::ConfigTree const& config);

}  // namespace MaterialLib
