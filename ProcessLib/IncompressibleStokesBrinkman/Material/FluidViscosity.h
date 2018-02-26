#pragma once

#include <cmath>
#include <memory>
#include "BaseLib/ConfigTree.h"
#include "MathLib/Curve/CreatePiecewiseLinearCurve.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

namespace ProcessLib
{
class EffectiveFluidViscosity
{
public:
    virtual double operator()(double const t, const double fluid_viscosity,
                              double const fluid_density) = 0;
    virtual ~EffectiveFluidViscosity() = default;
};

class EffectiveFluidViscosityIdentity : public EffectiveFluidViscosity
{
public:
    double operator()(double const /*t*/,
                      double const fluid_viscosity,
                      double const /*fluid_density*/) override
    {
        return fluid_viscosity;
    }
};

class EffectiveFluidViscosityGiese : public EffectiveFluidViscosity
{
public:
    EffectiveFluidViscosityGiese(
        std::unique_ptr<MathLib::PiecewiseLinearInterpolation>&&
            average_darcy_velocity,
        double const pellet_diameter)
        : _average_darcy_velocity(std::move(average_darcy_velocity)),
          _pellet_diameter(pellet_diameter)
    {
    }
    double operator()(double const t,
                      double const fluid_viscosity,
                      double const fluid_density) override
    {
        auto const Re0 = _average_darcy_velocity->getValue(t) *
                         _pellet_diameter * fluid_density / fluid_viscosity;
        return 2.0 * fluid_viscosity * std::exp(2e-3 * Re0);
    }

private:
    std::unique_ptr<MathLib::PiecewiseLinearInterpolation> const
        _average_darcy_velocity;
    const double _pellet_diameter;
};

inline std::unique_ptr<EffectiveFluidViscosity> createEffectiveFluidViscosity(
    BaseLib::ConfigTree const& config)
{
    auto const type = config.getConfigParameter<std::string>("type");
    if (type == "Identity")
        return std::make_unique<EffectiveFluidViscosityIdentity>();
    if (type == "Giese")
    {
        auto v = MathLib::createPiecewiseLinearCurve<
            MathLib::PiecewiseLinearInterpolation>(
            config.getConfigSubtree("average_darcy_velocity"));
        auto const d = config.getConfigParameter<double>("pellet_diameter");
        return std::make_unique<EffectiveFluidViscosityGiese>(std::move(v), d);
    }

    OGS_FATAL("Unknown effective viscosity model: %s.", type.c_str());
}

}  // namespace ProcessLib
