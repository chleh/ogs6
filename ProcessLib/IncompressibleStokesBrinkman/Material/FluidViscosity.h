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
    virtual double getViscosity(const double fluid_viscosity,
                                double const Re_0) const = 0;
    virtual ~EffectiveFluidViscosity() = default;
};

class EffectiveFluidViscosityIdentity : public EffectiveFluidViscosity
{
public:
    double getViscosity(double const fluid_viscosity,
                        double const /*Re_0*/) const override
    {
        return fluid_viscosity;
    }
};

class EffectiveFluidViscosityGiese : public EffectiveFluidViscosity
{
public:
    double getViscosity(double const fluid_viscosity,
                        double const Re_0) const override
    {
        return 2.0 * fluid_viscosity * std::exp(2e-3 * Re_0);
    }
};

class EffectiveFluidViscosityGieseCustomCoefficient
    : public EffectiveFluidViscosity
{
public:
    explicit EffectiveFluidViscosityGieseCustomCoefficient(double coeff)
        : _coeff(coeff)
    {
    }

    double getViscosity(double const fluid_viscosity,
                        double const Re_0) const override
    {
        return _coeff * fluid_viscosity * std::exp(2e-3 * Re_0);
    }

private:
    double _coeff;
};

class EffectiveFluidViscosityConstantFactor : public EffectiveFluidViscosity
{
public:
    explicit EffectiveFluidViscosityConstantFactor(double factor)
        : _factor(factor)
    {
    }

    double getViscosity(double const fluid_viscosity,
                        double const /*Re_0*/) const override
    {
        return _factor * fluid_viscosity;
    }

private:
    double _factor;
};

class ReynoldsNumber
{
public:
    virtual double getRe(double t,
                         double density,
                         double velocity,
                         double viscosity) const = 0;
    virtual ~ReynoldsNumber() = default;
};

class ReynoldsNumberLocal final : public ReynoldsNumber
{
public:
    ReynoldsNumberLocal(double characteristic_length)
        : _characteristic_length(characteristic_length)
    {
    }

    double getRe(double /*t*/,
                 double density,
                 double velocity,
                 double viscosity) const override
    {
        return density * velocity * _characteristic_length / viscosity;
    }

private:
    double const _characteristic_length;
};

class ReynoldsNumberNonLocal final : public ReynoldsNumber
{
public:
    ReynoldsNumberNonLocal(
        double characteristic_length,
        std::unique_ptr<MathLib::PiecewiseLinearInterpolation>&&
            average_velocity)
        : _characteristic_length(characteristic_length),
          _average_velocity(std::move(average_velocity))
    {
    }

    double getRe(double t,
                 double density,
                 double /*velocity*/,
                 double viscosity) const override
    {
        return density * _average_velocity->getValue(t) *
               _characteristic_length / viscosity;
    }

private:
    double const _characteristic_length;
    std::unique_ptr<MathLib::PiecewiseLinearInterpolation> const
        _average_velocity;
};

inline std::unique_ptr<EffectiveFluidViscosity> createEffectiveFluidViscosity(
    BaseLib::ConfigTree const& config)
{
    auto const type = config.getConfigParameter<std::string>("type");
    if (type == "Identity")
        return std::make_unique<EffectiveFluidViscosityIdentity>();
    if (type == "ConstantFactor")
    {
        auto const factor = config.getConfigParameter<double>("factor");
        return std::make_unique<EffectiveFluidViscosityConstantFactor>(factor);
    }
    if (type == "Giese")
    {
        return std::make_unique<EffectiveFluidViscosityGiese>();
    }
    if (type == "GieseCustomCoefficient")
    {
        auto const coeff = config.getConfigParameter<double>("coefficient");
        return std::make_unique<EffectiveFluidViscosityGieseCustomCoefficient>(
            coeff);
    }

    OGS_FATAL("Unknown effective viscosity model: %s.", type.c_str());
}

inline std::unique_ptr<ReynoldsNumber> createReynoldsNumber(
    BaseLib::ConfigTree const& config)
{
    auto const type = config.getConfigParameter<std::string>("type");
    if (type == "Local")
    {
        auto const l =
            config.getConfigParameter<double>("characteristic_length");
        return std::make_unique<ReynoldsNumberLocal>(l);
    }
    if (type == "NonLocal")
    {
        auto const l =
            config.getConfigParameter<double>("characteristic_length");
        auto v = MathLib::createPiecewiseLinearCurve<
            MathLib::PiecewiseLinearInterpolation>(
            config.getConfigSubtree("average_velocity"));
        return std::make_unique<ReynoldsNumberNonLocal>(l, std::move(v));
    }

    OGS_FATAL("Unknown effective viscosity model: %s.", type.c_str());
}

}  // namespace ProcessLib
