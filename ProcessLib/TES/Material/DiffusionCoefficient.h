#pragma once

#include <cmath>
#include <memory>

#include <boost/math/constants/constants.hpp>
#include "MaterialLib/PhysicalConstant.h"

namespace BaseLib
{
class ConfigTree;
}  // namespace BaseLib

namespace ProcessLib
{
namespace TES
{
class DiffusionCoefficient
{
public:
    // p -- pressure of the gas mixture
    // T -- temperature of the gas mixture
    // p_V -- vapour mass fraction
    virtual double getDiffusionCoefficient(const double p, const double T,
                                           const double p_V) const = 0;

    virtual ~DiffusionCoefficient() = default;
};

class DiffusionCoefficientConstant final : public DiffusionCoefficient
{
public:
    explicit DiffusionCoefficientConstant(const double D) : _d(D) {}
    virtual double getDiffusionCoefficient(const double /*p*/,
                                           const double /*T*/,
                                           const double /*p_V*/) const override
    {
        return _d;
    }

private:
    const double _d;
};

class DiffusionCoefficientWaterNitrogenMarrero final
    : public DiffusionCoefficient
{
public:
    virtual double getDiffusionCoefficient(const double p, const double T,
                                           const double /*p_V*/) const override
    {
        // Marrero, T.R., Mason, E.A., 1972. Gaseous Diffusion Coefficients.
        // Journal of Physical and Chemical Reference Data 1, 3–118.
        // doi:10.1063/1.3253094

        // Here we don't account for changes of D with gas composition.
        return 1.87e-6 * std::pow(T, 2.072) * 1e-4 * 101325 / p;
    }
};

class DiffusionCoefficientWaterAirFuller final : public DiffusionCoefficient
{
public:
    virtual double getDiffusionCoefficient(const double p, const double T,
                                           const double /*p_V*/) const override
    {
        // Poling, B.E., Prausnitz, J.M., O’Connell, J.P., 2001. The properties
        // of gases and liquids, 5th ed. ed. McGraw-Hill, New York.
        // Eq. (11-4.4)

        // Original:
        // Fuller, E.N., Ensley, K., Giddings, J.C., 1969. Diffusion of
        // halogenated hydrocarbons in helium. The effect of structure on
        // collision cross sections. The Journal of Physical Chemistry 73,
        // 3679–3685. doi:10.1021/j100845a020

        double const diff_vol_water = 13.1;
        double const diff_vol_air = 19.7;
        double const diff_vol_factor =
            std::cbrt(diff_vol_water) + std::cbrt(diff_vol_air);

        double const M_air =
            MaterialLib::PhysicalConstant::MolarMass::Air * 1000.0;  // in g/mol
        double const M_water = MaterialLib::PhysicalConstant::MolarMass::Water *
                               1000.0;  // in g/mol
        double const M_AB = 2.0 / (1.0 / M_air + 1.0 / M_water);

        double const p_bar = p * 1e-5;  // pressure in bar

        double const D_AB = 0.00143 * std::pow(T, 1.75) /
                            (p_bar * std::sqrt(M_AB) * diff_vol_factor *
                             diff_vol_factor);  // in cm^2/s
        return D_AB / 1e4;                      // in m^2/s
    }
};

class DiffusionCoefficientKnudsen final : public DiffusionCoefficient
{
public:
    DiffusionCoefficientKnudsen(const double pore_diam, const double molar_mass)
        : _pore_diam(pore_diam), _molar_mass(molar_mass)
    {
    }

    virtual double getDiffusionCoefficient(const double /*p*/, const double T,
                                           const double /*p_V*/) const override
    {
        auto const R = MaterialLib::PhysicalConstant::IdealGasConstant;
        return 4.0 / 3.0 * _pore_diam *
               std::sqrt(R * T / 2.0 / boost::math::constants::pi<double>() /
                         _molar_mass);
    }

private:
    const double _pore_diam;   //!< in m
    const double _molar_mass;  //!< in kg/mol
};

class DiffusionCoefficientGasAndKnudsen final : public DiffusionCoefficient
{
public:
    DiffusionCoefficientGasAndKnudsen(
        std::unique_ptr<DiffusionCoefficient>&& D_gas,
        std::unique_ptr<DiffusionCoefficient>&& D_Knudsen)
        : _diff_coeff_gas(std::move(D_gas)),
          _diff_coeff_Knudsen(std::move(D_Knudsen))
    {
    }

    virtual double getDiffusionCoefficient(const double p, const double T,
                                           const double p_V) const override
    {
        auto const D_gas = _diff_coeff_gas->getDiffusionCoefficient(p, T, p_V);
        auto const D_Knudsen =
            _diff_coeff_Knudsen->getDiffusionCoefficient(p, T, p_V);

        return 1.0 / (1.0 / D_gas + 1.0 / D_Knudsen);
    }

private:
    std::unique_ptr<DiffusionCoefficient> const _diff_coeff_gas;
    std::unique_ptr<DiffusionCoefficient> const _diff_coeff_Knudsen;
};

std::unique_ptr<DiffusionCoefficient> createDiffusionCoefficient(
    BaseLib::ConfigTree const& config);

}  // namespace TES

}  // namespace ProcessLib
