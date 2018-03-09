#pragma once

#include <cassert>
#include <memory>

#include "MaterialLib/ReactionKinetics/ReactionRate.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "ProcessLib/IncompressibleStokesBrinkman/Material/FluidViscosity.h"
#include "ProcessLib/TES/Material/DiffusionCoefficient.h"
#include "ProcessLib/TES/TESOGS5MaterialModels.h"

namespace detail
{
inline double polynomial(double x, double const* coeffs, int n_coeffs)
{
    double x_pow = 1.0;
    double res = 0.0;
    for (int i = 0; i < n_coeffs; ++i)
    {
        res += x_pow * coeffs[i];
        x_pow *= x;
    }
    return res;
}
}  // namespace detail

namespace ProcessLib
{
namespace TCHSStokes
{
namespace Material
{
class FluidDensity
{
public:
    virtual double getDensity(const double p,
                              const double T,
                              const double x_mV) const = 0;
    virtual double getThermalExpansionCoefficient(const double p,
                                                  const double T,
                                                  const double x_mV) const = 0;
    virtual double getCompressibility(const double p,
                                      const double T,
                                      const double x_mV) const = 0;
    virtual double getDensityChangeWithComposition(const double p,
                                                   const double T,
                                                   const double x_mV) const = 0;

    virtual ~FluidDensity() = default;
};

class FluidDensityMixtureWaterNitrogen final : public FluidDensity
{
public:
    double getDensity(const double p,
                      const double T,
                      const double x_mV) const override
    {
        return ProcessLib::TES::fluid_density(p, T, x_mV);
    }

    double getThermalExpansionCoefficient(const double /*p*/,
                                          const double T,
                                          const double /*x_mV*/) const override
    {
        return 1.0 / T;
    }

    double getCompressibility(const double p,
                              const double /*T*/,
                              const double /*x_mV*/) const override
    {
        return 1.0 / p;
    }

    double getDensityChangeWithComposition(const double /*p*/,
                                           const double /*T*/,
                                           const double x_mV) const override
    {
        double const M_R = MaterialLib::PhysicalConstant::MolarMass::Water;
        double const M_I = MaterialLib::PhysicalConstant::MolarMass::N2;

        // 1 / M_G * dM_G/dx_mV
        double const dMG_dxmV_over_MG =
            (M_R - M_I) / (M_I * x_mV + M_R * (1.0 - x_mV));

        return dMG_dxmV_over_MG;
    }
};

class FluidDensityConstant final : public FluidDensity
{
public:
    explicit FluidDensityConstant(double value) : _value(value) {}

    double getDensity(const double /*p*/,
                      const double /*T*/,
                      const double /*x_mV*/) const override
    {
        return _value;
    }

    double getThermalExpansionCoefficient(const double /*p*/,
                                          const double /*T*/,
                                          const double /*x_mV*/) const override
    {
        return 0.0;
    }

    double getCompressibility(const double /*p*/,
                              const double /*T*/,
                              const double /*x_mV*/) const override
    {
        return 0.0;
    }

    double getDensityChangeWithComposition(const double /*p*/,
                                           const double /*T*/,
                                           const double /*x_mV*/) const override
    {
        return 0.0;
    }

private:
    double const _value;
};

class FluidViscosity
{
public:
    virtual double getViscosity(const double p,
                                const double T,
                                const double x) const = 0;
    virtual ~FluidViscosity() = default;
};

class FluidViscosityMixtureWaterNitrogen final : public FluidViscosity
{
public:
    double getViscosity(const double p,
                        const double T,
                        const double x) const override
    {
        return ProcessLib::TES::fluid_viscosity(p, T, x);
    }
};

class FluidViscosityConstant final : public FluidViscosity
{
public:
    explicit FluidViscosityConstant(double value) : _value(value) {}

    double getViscosity(const double /*p*/,
                        const double /*T*/,
                        const double /*x*/) const override
    {
        return _value;
    }

private:
    double const _value;
};

class FluidHeatCapacity
{
public:
    virtual double getHeatCapacity(double T, double x_mV) const = 0;
    virtual ~FluidHeatCapacity() = default;
};

class FluidHeatCapacityMixtureWaterNitrogen final : public FluidHeatCapacity
{
public:
    double getHeatCapacity(double T, double x_mV) const override
    {
        // Poling, B.E., Prausnitz, J.M., O’Connell, J.P., 2001. The properties
        // of gases and liquids, 5th ed. ed. McGraw-Hill, New York. Appendix
        // Section C
        const double a_N2[5] = {3.539, -0.261e-3, 0.007e-5, 0.157e-8,
                                -0.099e-11};
        const double a_H2O[5] = {4.395, -4.186e-3, 1.405e-5, -1.564e-8,
                                 0.632e-11};

        const double cp_N2 = MaterialLib::PhysicalConstant::IdealGasConstant *
                             ::detail::polynomial(T, a_N2, 5);
        const double cp_H2O = MaterialLib::PhysicalConstant::IdealGasConstant *
                              ::detail::polynomial(T, a_H2O, 5);

        return x_mV * cp_H2O + (1.0 - x_mV) * cp_N2;
    }
};

class FluidHeatCapacityConstant final : public FluidHeatCapacity
{
public:
    explicit FluidHeatCapacityConstant(double value) : _value(value) {}

    double getHeatCapacity(double /*T*/, double /*x_mV*/) const override
    {
        return _value;
    }

private:
    double const _value;
};

class PecletNumberHeat
{
public:
    virtual double getPe(double t,
                         double density,
                         double velocity,
                         double heat_capacity,
                         double heat_conductivity) const = 0;
    virtual ~PecletNumberHeat() = default;
};

class PecletNumberHeatLocal final : public PecletNumberHeat
{
public:
    PecletNumberHeatLocal(double characteristic_length)
        : _characteristic_length(characteristic_length)
    {
    }

    double getPe(double /*t*/,
                 double density,
                 double velocity,
                 double heat_capacity,
                 double heat_conductivity) const override
    {
        return density * velocity * heat_capacity * _characteristic_length /
               heat_conductivity;
    }

private:
    double const _characteristic_length;
};

class PecletNumberHeatNonLocal final : public PecletNumberHeat
{
public:
    PecletNumberHeatNonLocal(
        double characteristic_length,
        std::unique_ptr<MathLib::PiecewiseLinearInterpolation>&&
            average_velocity)
        : _characteristic_length(characteristic_length),
          _average_velocity(std::move(average_velocity))
    {
    }

    double getPe(double t,
                 double density,
                 double /*velocity*/,
                 double heat_capacity,
                 double heat_conductivity) const override
    {
        return density * _average_velocity->getValue(t) * heat_capacity *
               _characteristic_length / heat_conductivity;
    }

private:
    double const _characteristic_length;
    std::unique_ptr<MathLib::PiecewiseLinearInterpolation> const
        _average_velocity;
};

class HeatConductivity
{
public:
    virtual Eigen::DiagonalMatrix<double, 2> getHeatConductivity(
        double const t, double const p, double const T, double const x_mV,
        double const r, double const porosity, double const rho_GR,
        double const c_pG, double const Re_0, double const v_Darcy,
        double const v_Darcy_center) const = 0;

    virtual ~HeatConductivity() = default;
};

class HeatConductivityLambdaRNoRadiation final : public HeatConductivity
{
public:
    HeatConductivityLambdaRNoRadiation(
        double bed_radius,
        double pellet_diameter,
        std::unique_ptr<HeatConductivity>&& lambda_fluid,
        std::unique_ptr<HeatConductivity>&& lambda_pellet,
        std::unique_ptr<PecletNumberHeat>&& peclet_number_heat)
        : _bed_radius(bed_radius),
          _pellet_diameter(pellet_diameter),
          _lambda_fluid(std::move(lambda_fluid)),
          _lambda_pellet(std::move(lambda_pellet)),
          _peclet_number_heat(std::move(peclet_number_heat))
    {
    }

    Eigen::DiagonalMatrix<double, 2> getHeatConductivity(
        double const t, double const p, double const T, double const x_mV,
        double const r, double const porosity, double const rho_GR,
        double const c_pG, double const Re_0, double const v_Darcy,
        double const v_Darcy_center) const override
    {
        auto const lambda_f_mat = _lambda_fluid->getHeatConductivity(
            t, p, T, x_mV, r, porosity, rho_GR, c_pG, Re_0, v_Darcy,
            v_Darcy_center);
        assert(lambda_f_mat.diagonal()[0] == lambda_f_mat.diagonal()[1]);
        double const lambda_f = lambda_f_mat.diagonal()[0];

        auto const lambda_p_mat = _lambda_pellet->getHeatConductivity(
            t, p, T, x_mV, r, porosity, rho_GR, c_pG, Re_0, v_Darcy,
            v_Darcy_center);
        assert(lambda_p_mat.diagonal()[0] == lambda_p_mat.diagonal()[1]);
        double const lambda_p = lambda_p_mat.diagonal()[0];

        double const Pe_0 =
            _peclet_number_heat->getPe(t, rho_GR, v_Darcy, c_pG, lambda_f);

        double const B = 1.25 * std::pow(1.0 / porosity - 1.0, 10.0 / 9.0);
        double const k_p = lambda_p / lambda_f;
        double const N = 1.0 - B / k_p;
        double const k_c = 2.0 / N *
                           (B / N / N * (1.0 - 1.0 / k_p) * std::log(k_p / B) -
                            0.5 * (B + 1.0) - (B - 1.0) / N);
        double const k_bed =
            1 - std::sqrt(1.0 - porosity) + std::sqrt(1.0 - porosity) * k_c;
        double const lambda_bed = k_bed * lambda_f;

        double const K_ax = 2.0;
        double const lambda_ax = lambda_bed + Pe_0 / K_ax * lambda_f;

        auto f = [this, Re_0](double x) {
            double const K_2 = 0.44 + 4.0 * std::exp(-Re_0 / 70.0);
            if (x < K_2 * _pellet_diameter)
                return boost::math::pow<2 /* n */>(x / K_2 / _pellet_diameter);
            return 1.0;
        };
        double const K_1 = 0.125;
        double const lambda_r = lambda_bed + K_1 * Pe_0 * v_Darcy_center /
                                                 v_Darcy * f(_bed_radius - r) *
                                                 lambda_f;

        return Eigen::DiagonalMatrix<double, 2>(lambda_r, lambda_ax);
    }

private:
    double const _bed_radius;
    double const _pellet_diameter;
    std::unique_ptr<HeatConductivity> _lambda_fluid;
    std::unique_ptr<HeatConductivity> _lambda_pellet;
    std::unique_ptr<PecletNumberHeat> _peclet_number_heat;
};

class HeatConductivityMixtureWaterNitrogen final : public HeatConductivity
{
public:
    Eigen::DiagonalMatrix<double, 2> getHeatConductivity(
        double const /*t*/, double const p, double const T, double const x_mV,
        double const /*r*/, double const /*porosity*/, double const /*rho_GR*/,
        double const /*c_pG*/, double const /*Re_0*/, double const /*v_Darcy*/,
        double const /*v_Darcy_center*/) const override
    {
        double const lambda =
            ProcessLib::TES::fluid_heat_conductivity(p, T, x_mV);
        return Eigen::DiagonalMatrix<double, 2>(lambda, lambda);
    }
};

class HeatConductivityConstant final : public HeatConductivity
{
public:
    HeatConductivityConstant(double value) : _value(value) {}

    Eigen::DiagonalMatrix<double, 2> getHeatConductivity(
        double const /*t*/, double const /*p*/, double const /*T*/,
        double const /*x_mV*/, double const /*r*/, double const /*porosity*/,
        double const /*rho_GR*/, double const /*c_pG*/, double const /*Re_0*/,
        double const /*v_Darcy*/,
        double const /*v_Darcy_center*/) const override
    {
        return Eigen::DiagonalMatrix<double, 2>(_value, _value);
    }

private:
    double const _value;
};

class PecletNumberMass
{
public:
    virtual double getPe(double t,
                         double velocity,
                         double diffusion_coefficient) const = 0;
    virtual ~PecletNumberMass() = default;
};

class PecletNumberMassLocal final : public PecletNumberMass
{
public:
    explicit PecletNumberMassLocal(double characteristic_length)
        : _characteristic_length(characteristic_length)
    {
    }

    double getPe(double /*t*/,
                 double velocity,
                 double diffusion_coefficient) const override
    {
        return velocity * _characteristic_length / diffusion_coefficient;
    }

private:
    double const _characteristic_length;
};

class PecletNumberMassNonLocal final : public PecletNumberMass
{
public:
    PecletNumberMassNonLocal(
        double characteristic_length,
        std::unique_ptr<MathLib::PiecewiseLinearInterpolation>&&
            average_velocity)
        : _characteristic_length(characteristic_length),
          _average_velocity(std::move(average_velocity))
    {
    }

    double getPe(double t,
                 double /*velocity*/,
                 double diffusion_coefficient) const override
    {
        return _average_velocity->getValue(t) * _characteristic_length /
               diffusion_coefficient;
    }

private:
    double const _characteristic_length;
    std::unique_ptr<MathLib::PiecewiseLinearInterpolation> const
        _average_velocity;
};

class MassDispersion
{
public:
    virtual Eigen::DiagonalMatrix<double, 2> getMassDispersion(
        double const t, double const p, double const T, double const v_Darcy,
        double const r, double const porosity, double const p_center,
        double const T_center, double const v_Darcy_center) const = 0;

    virtual ~MassDispersion() = default;
};

class MassDispersionLambdaRModel final : public MassDispersion
{
public:
    MassDispersionLambdaRModel(
        double bed_radius,
        double pellet_diameter,
        std::unique_ptr<ProcessLib::TES::DiffusionCoefficient>&&
            diffusion_coefficient,
        std::unique_ptr<PecletNumberMass>&& peclet_number_mass)
        : _bed_radius(bed_radius),
          _pellet_diameter(pellet_diameter),
          _diffusion_coefficient(std::move(diffusion_coefficient)),
          _peclet_number_mass(std::move(peclet_number_mass)),
          _peclet_number_mass_center(_pellet_diameter)
    {
    }

    Eigen::DiagonalMatrix<double, 2> getMassDispersion(
        double const t, double const p, double const T, double const v_Darcy,
        double const r, double const porosity, double const p_center,
        double const T_center, double const v_Darcy_center) const override
    {
        double const diff =
            _diffusion_coefficient->getDiffusionCoefficient(p, T);
        double const diff_bed = diff * (1.0 - std::sqrt(1.0 - porosity));
        double const Pe_0 = _peclet_number_mass->getPe(t, v_Darcy, diff);
        double const K_ax = 2.0;
        double const D_ax = diff_bed + Pe_0 / K_ax * diff;

        auto f = [this](double x) {
            double const K_2 = 0.44;
            if (x < K_2 * _pellet_diameter)
                return boost::math::pow<2 /* n */>(x / K_2 / _pellet_diameter);
            return 1.0;
        };
        double const diff_center =
            _diffusion_coefficient->getDiffusionCoefficient(p_center, T_center);
        double const sqrt_Pe_0_center = std::sqrt(
            _peclet_number_mass_center.getPe(t, v_Darcy_center, diff_center));
        double const K_1 = sqrt_Pe_0_center * 0.125 / (3.0 + sqrt_Pe_0_center);
        double D_r = diff_bed;
        if (std::abs(v_Darcy) > std::numeric_limits<double>::epsilon())
            D_r += K_1 * Pe_0 * v_Darcy_center / v_Darcy * f(_bed_radius - r) *
                   diff;

        return {D_r, D_ax};
    }

private:
    double const _bed_radius;
    double const _pellet_diameter;
    std::unique_ptr<ProcessLib::TES::DiffusionCoefficient> const
        _diffusion_coefficient;
    std::unique_ptr<PecletNumberMass> const _peclet_number_mass;
    PecletNumberMassLocal const _peclet_number_mass_center;
};

class MassDispersionIdentity final : public MassDispersion
{
public:
    MassDispersionIdentity(
        std::unique_ptr<ProcessLib::TES::DiffusionCoefficient>&&
            diffusion_coefficient)
        : _diffusion_coefficient(std::move(diffusion_coefficient))
    {
    }

    Eigen::DiagonalMatrix<double, 2> getMassDispersion(
        double const /*t*/, double const p, double const T,
        double const /*v_Darcy*/, double const /*r*/, double const /*porosity*/,
        double const /*p_center*/, double const /*T_center*/,
        double const /*v_Darcy_center*/) const override
    {
        double const diff =
            _diffusion_coefficient->getDiffusionCoefficient(p, T);

        return {diff, diff};
    }

private:
    std::unique_ptr<ProcessLib::TES::DiffusionCoefficient> const
        _diffusion_coefficient;
};

class FluidMomentumProductionCoefficient
{
public:
    virtual double getCoeffOfV(double porosity, double viscosity) const = 0;
    virtual double getCoeffOfVSquared(double porosity,
                                      double fluid_density) const = 0;

    virtual ~FluidMomentumProductionCoefficient() = default;
};

class FluidMomentumProductionCoefficientErgun final
    : public FluidMomentumProductionCoefficient
{
public:
    FluidMomentumProductionCoefficientErgun(const double pellet_diameter)
        : _pellet_diameter(pellet_diameter)
    {
    }

    double getCoeffOfV(double porosity, double viscosity) const override
    {
        auto const poro3 = boost::math::pow<3>(porosity);
        return 150.0 * boost::math::pow<2>(1.0 - porosity) / poro3 * viscosity /
               _pellet_diameter / _pellet_diameter;
    }
    double getCoeffOfVSquared(double porosity,
                              double fluid_density) const override
    {
        auto const poro3 = boost::math::pow<3>(porosity);
        return 1.75 * (1.0 - porosity) / poro3 * fluid_density /
               _pellet_diameter;
    }

private:
    double const _pellet_diameter;
};

class FluidMomentumProductionCoefficientZero final
    : public FluidMomentumProductionCoefficient
{
    double getCoeffOfV(double /*porosity*/, double /*viscosity*/) const override
    {
        return 0;
    }
    double getCoeffOfVSquared(double /*porosity*/,
                              double /*fluid_density*/) const override
    {
        return 0;
    }
};

class SolidHeatCapacity
{
public:
    virtual double getSpecificHeatCapacity(double rho_SR, double T) const = 0;
    virtual ~SolidHeatCapacity() = default;
};

class SolidHeatCapacityZeoliteNaAWaterVucelic final : public SolidHeatCapacity
{
public:
    SolidHeatCapacityZeoliteNaAWaterVucelic(double cp_zeo_dry,
                                            double rho_SR_dry)
        : _cp_zeo_dry(cp_zeo_dry), _rho_SR_dry(rho_SR_dry)
    {
    }

    double getSpecificHeatCapacity(double rho_SR, double T) const override
    {
        // Vučelić, V., Vučelić, D., 1985. Heat Capacities of Water on Zeolites,
        // in: Studies in Surface Science and Catalysis. Elsevier, pp. 475–480.

        double const cp_min_water = 836.0;  // J/kg/K
        double const A = 3.76e1;
        double const B = 3.976e-4;
        double const sqrtB = std::sqrt(B);
        double const T_0 = 335;    // K
        double const T_min = 220;  // K
        double const sqrtpi = boost::math::constants::root_pi<double>();

        double const loading = rho_SR / _rho_SR_dry - 1.0;

        double const cp_water =
            cp_min_water +
            0.5 * sqrtpi * A / sqrtB *
                (std::erf(sqrtB * (T - T_0)) - std::erf(sqrtB * (T_min - T_0)));

        return 1.0 / (1.0 + loading) * _cp_zeo_dry +
               loading / (1.0 + loading) * cp_water;
    }

private:
    double const _cp_zeo_dry;
    double const _rho_SR_dry;
};

class SolidHeatCapacityZeoliteCaAWaterVucelic final : public SolidHeatCapacity
{
public:
    SolidHeatCapacityZeoliteCaAWaterVucelic(double cp_zeo_dry,
                                            double rho_SR_dry)
        : _cp_zeo_dry(cp_zeo_dry), _rho_SR_dry(rho_SR_dry)
    {
    }

    double getSpecificHeatCapacity(double rho_SR, double T) const override
    {
        // Vučelić, V., Vučelić, D., 1985. Heat Capacities of Water on Zeolites,
        // in: Studies in Surface Science and Catalysis. Elsevier, pp. 475–480.

        double const cp_min_water = 1331;  // J/kg/K
        double const A = 7.74e1;
        double const B = 27.027e-4;
        double const sqrtB = std::sqrt(B);
        double const T_0 = 368;    // K
        double const T_min = 220;  // K
        double const sqrtpi = boost::math::constants::root_pi<double>();

        double const loading = rho_SR / _rho_SR_dry - 1.0;

        double const cp_water =
            cp_min_water +
            0.5 * sqrtpi * A / sqrtB *
                (std::erf(sqrtB * (T - T_0)) - std::erf(sqrtB * (T_min - T_0)));

        return 1.0 / (1.0 + loading) * _cp_zeo_dry +
               loading / (1.0 + loading) * cp_water;
    }

private:
    double const _cp_zeo_dry;
    double const _rho_SR_dry;
};

class SolidHeatCapacityConstant final : public SolidHeatCapacity
{
public:
    SolidHeatCapacityConstant(double value) : _value(value) {}

    double getSpecificHeatCapacity(double /*rho_SR*/,
                                   double /*T*/) const override
    {
        return _value;
    }

private:
    double const _value;
};

class Porosity
{
public:
    virtual double getPorosity(double r) const = 0;
    virtual double getDPorosityDr(double r) const = 0;
    ~Porosity() = default;
};

class PorosityConstant : public Porosity
{
public:
    PorosityConstant(double value) : _value(value) {}

    double getPorosity(double /*r*/) const override { return _value; }
    double getDPorosityDr(double /*r*/) const override { return 0.0; }

private:
    double const _value;
};

class PorosityRadialProfileGiese : public Porosity
{
public:
    PorosityRadialProfileGiese(const double bed_radius,
                               const double pellet_diameter,
                               const double homogeneous_porosity)
        : _bed_radius(bed_radius),
          _pellet_diameter(pellet_diameter),
          _homogeneous_porosity(homogeneous_porosity)
    {
    }

    double getPorosity(double r) const override
    {
        auto const exp_term =
            std::exp(-5.0 * (_bed_radius - r) / _pellet_diameter);
        return _homogeneous_porosity + _homogeneous_porosity * 1.36 * exp_term;
    }

    double getDPorosityDr(double r) const override
    {
        auto const exp_term =
            std::exp(-5.0 * (_bed_radius - r) / _pellet_diameter);
        return _homogeneous_porosity * 1.36 * 5.0 / _pellet_diameter * exp_term;
    }

private:
    const double _bed_radius;
    const double _pellet_diameter;
    const double _homogeneous_porosity;
};

struct TCHSStokesMaterial
{
    std::unique_ptr<FluidDensity> fluid_density;
    std::unique_ptr<FluidViscosity> fluid_viscosity;
    std::unique_ptr<EffectiveFluidViscosity> effective_fluid_viscosity;
    std::unique_ptr<FluidHeatCapacity> fluid_heat_capacity;
    std::unique_ptr<HeatConductivity> heat_conductivity;
    std::unique_ptr<MassDispersion> mass_dispersion;
    std::unique_ptr<FluidMomentumProductionCoefficient>
        fluid_momentum_production_coefficient;
    std::unique_ptr<SolidHeatCapacity> solid_heat_capacity;
    std::unique_ptr<MaterialLib::ReactiveSolidModel> reactive_solid;
    std::unique_ptr<MaterialLib::ReactionRate> reaction_rate;
    std::unique_ptr<Porosity> porosity;

    std::unique_ptr<ProcessLib::ReynoldsNumber> reynolds_number;
    std::unique_ptr<PecletNumberHeat> peclet_number_heat;

    double molar_mass_reactive = std::numeric_limits<double>::quiet_NaN();
    double molar_mass_inert = std::numeric_limits<double>::quiet_NaN();
};

}  // namespace Material

}  // namespace TCHSStokes

}  // namespace ProcessLib
