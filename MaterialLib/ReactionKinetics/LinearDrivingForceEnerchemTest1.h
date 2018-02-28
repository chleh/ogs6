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
#include "ReactionKineticsByEquilibrium.h"

namespace MaterialLib
{
class LinearDrivingForceEnerchemTest1 final
    : public ReactionKineticsByEquilibrium
{
public:
    LinearDrivingForceEnerchemTest1(
        const double pellet_radius, const double pellet_porosity,
        const double pellet_tortuosity, const double pellet_density_dry,
        const double loading_dependent_damping_factor,
        const double molar_mass_reactive_component,
        std::unique_ptr<ProcessLib::TES::DiffusionCoefficient>&&
            diffusion_coefficient,
        std::unique_ptr<ReactionEquilibrium>&& equil)
        : _prefactor(15.0 / pellet_radius / pellet_radius * pellet_porosity /
                     pellet_tortuosity * molar_mass_reactive_component /
                     PhysicalConstant::IdealGasConstant),
          _pellet_density_dry(pellet_density_dry),
          _loading_dependent_damping_factor(loading_dependent_damping_factor),
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
        auto const rho_SR_eq = _equil->getEquilibriumDensity(p_V, T);

        // For Ca(90)XBF roughly equivalent to 0.1 * (dC_eq/dp_V)^{-1}
        auto const isotherm_factor =
            std::min(100.0, p_V / (rho_SR_eq - _pellet_density_dry));

        auto const loading = rho_SR / _pellet_density_dry - 1.0;
        auto const loading_dependent_damping =
            std::max(1e-3, 1.0 - _loading_dependent_damping_factor * loading);

        auto const k_LDF =
            _prefactor * D * loading_dependent_damping / T * isotherm_factor;
        return Eigen::VectorXd::Constant(1, k_LDF * (rho_SR_eq - rho_SR));
    }

    HeatOfReactionData getHeatOfReaction(
        const double p_V, const double T,
        const ReactiveSolidState* const state) const
    {
        return _equil->getHeatOfReaction(p_V, T, state);
    }

    std::vector<NumLib::NamedFunction> getNamedFunctions() const override
    {
        return _equil->getNamedFunctions();
    }

    bool isStateCompatible(ReactiveSolidState& state) const override
    {
        return state.conversion().size() == 1;
    }

private:
    const double _prefactor;
    const double _pellet_density_dry;
    const double _loading_dependent_damping_factor;
    std::unique_ptr<ProcessLib::TES::DiffusionCoefficient> const
        _diffusion_coefficient;
    std::unique_ptr<ReactionEquilibrium> const _equil;
};

std::unique_ptr<ReactionKinetics> createLinearDrivingForceEnerchemTest1(
    BaseLib::ConfigTree const& config);

}  // namespace MaterialLib
