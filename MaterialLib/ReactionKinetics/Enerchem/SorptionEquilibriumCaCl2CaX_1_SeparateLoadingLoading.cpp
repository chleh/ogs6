/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SorptionEquilibriumCaCl2CaX_1_SeparateLoadingLoading.h"

#include <cmath>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "BaseLib/Functional.h"
#include "MaterialLib/Adsorption/Adsorption.h"
#include "MaterialLib/ReactionKinetics/DubininPolanyi.h"
#include "MaterialLib/ReactionKinetics/ReactiveSolidStateTwoComponents.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

namespace
{
//! Converts p_V/Pa to X/(g/kg).
static double p_V_to_X(const double p_V)
{
    auto const M_H2O = 18.02;
    auto const M_Air = 28.96;
    auto const p0 = 1.01325e5;

    auto const xnV = p_V / p0;
    auto const xmV =
        Adsorption::AdsorptionReaction::getMassFraction(xnV, M_H2O, M_Air);
    auto const X = 1000.0 * xmV / (1.0 - xmV);
    return X;
}

//! CaCl2 deliquescence line
static const MathLib::PiecewiseLinearInterpolation CaCl2_deliquescence_line(
    {// deliquescence line, T
     3.174869999999999663e+02,  3.194769999999999754e+02,
     3.218899999999999864e+02,  3.241819999999999595e+02,
     3.265949999999999704e+02,  3.288870000000000005e+02,
     3.309379999999999882e+02,  3.334719999999999800e+02,
     3.356429999999999723e+02,  3.382969999999999686e+02,
     3.405899999999999750e+02,  3.422789999999999964e+02,
     3.452949999999999591e+02,  3.474660000000000082e+02,
     3.498789999999999623e+02,  3.522919999999999732e+02,
     3.538599999999999568e+02,  3.565149999999999864e+02,
     3.589269999999999641e+02,  3.610989999999999895e+02,
     3.636329999999999814e+02,  3.667690000000000055e+02,
     3.697849999999999682e+02,  3.720779999999999745e+02,
     3.750939999999999941e+02,  3.783509999999999991e+02,
     3.823319999999999936e+02,  3.857099999999999795e+02,
     3.887259999999999991e+02,  3.918629999999999427e+02,
     3.952409999999999854e+02,  3.986189999999999714e+02,
     4.025999999999999659e+02,  4.062199999999999704e+02,
     4.088739999999999668e+02,  4.118899999999999864e+02,
     4.149059999999999491e+02,  4.180430000000000064e+02,
     4.221449999999999818e+02,  4.258849999999999909e+02,
     4.291419999999999391e+02,  4.325199999999999818e+02,
     4.348120000000000118e+02,  4.379489999999999554e+02,
     4.414479999999999791e+02,  4.436189999999999714e+02,
     4.459109999999999445e+02,  4.479619999999999891e+02,
     3.0337189999999999573e+02, 3.047589999999999577e+02,
     3.064479999999999791e+02,  3.078949999999999818e+02,
     3.093429999999999609e+02,  3.110319999999999823e+02,
     3.132039999999999509e+02,  3.146509999999999536e+02,
     3.163399999999999750e+02,  3.174868999999999663e+02,
     2.733909999999999627e+02,  2.753220000000000027e+02,
     2.774929999999999950e+02,  2.796649999999999636e+02,
     2.820779999999999745e+02,  2.844899999999999523e+02,
     2.865409999999999968e+02,  2.888340000000000032e+02,
     2.907639999999999532e+02,  2.928149999999999977e+02,
     2.945039999999999623e+02,  2.964339999999999691e+02,
     2.982439999999999714e+02,  2.996920000000000073e+02,
     3.021039999999999850e+02,  3.033719999999999573e+02},
    // deliquescence line, p_V
    {1.552268046000000140e+03, 1.702415282400000024e+03,
     1.925436324000000013e+03, 2.161016298000000006e+03,
     2.388330308000000059e+03, 2.743100149999999758e+03,
     2.939750100000000202e+03, 3.248923818000000210e+03,
     3.674620964000000185e+03, 4.253105121999999938e+03,
     4.664403491999999460e+03, 5.076368472000000111e+03,
     5.785774833999999828e+03, 6.443718904000000293e+03,
     7.176589938000000075e+03, 7.992787221999999929e+03,
     8.565805177999998705e+03, 9.613582776000001104e+03,
     1.062469682399999874e+04, 1.156288373800000045e+04,
     1.287783862399999998e+04, 1.456449524600000041e+04,
     1.622048780799999986e+04, 1.765316602000000057e+04,
     1.996496949999999924e+04, 2.223410994000000210e+04,
     2.553649587999999858e+04, 2.888021164000000135e+04,
     3.094936908000000039e+04, 3.554631163999999990e+04,
     3.868204507999999623e+04, 4.308033785999999964e+04,
     4.872119168000000354e+04, 5.342879149999999936e+04,
     5.770042838000000484e+04, 6.231070314000000508e+04,
     6.729027984000000288e+04, 7.266715609999999288e+04,
     7.968922583999999915e+04, 8.539407421999999497e+04,
     9.221749418000000878e+04, 9.882093284000000858e+04,
     1.050910665000000008e+05, 1.092013837600000115e+05,
     1.091853851200000063e+05, 1.050510699000000022e+05,
     9.876493759999999020e+04, 9.285610656000000017e+04,
     9.272011812000000646e+02, 1.024712892000000011e+03,
     1.132530393400000094e+03, 1.232575222200000098e+03,
     1.310808571799999982e+03, 1.437611126000000013e+03,
     1.528830038400000149e+03, 1.552454696799999965e+03,
     1.552348039199999903e+03, 1.552268046000000140e+03,
     2.605378524000000198e+02, 3.015610317999999666e+02,
     3.490503282000000240e+02, 4.134448542000000089e+02,
     4.859720221999999694e+02, 5.538862489999999070e+02,
     6.510513226000000486e+02, 7.306978854000000183e+02,
     7.891329180000000179e+02, 8.456881104000000278e+02,
     8.720858663999999862e+02, 8.923908069999999952e+02,
     9.131757068000000572e+02, 9.201751118000000815e+02,
     9.272011812000000646e+02, 9.272011812000000646e+02});

static MaterialLib::DubininPolanyi CaX_equilibrium(
    1.0 /* Using fake density here */,
    MaterialLib::createDubininPolanyiData("CaX80NoOutliers_Hauer"));

//! Computes deliquescence humidity of CaCl2.
static double CaCl2_p_V_deliquescence(const double T)
{
    return CaCl2_deliquescence_line.getValue(T);
}

static double CaCl2_water_fraction_hydrate_solution(const double p_Ads,
                                                    const double T)
{
    // TODO Daten fuer Hydrate noch nicht eingepflegt
    assert(p_Ads >= CaCl2_p_V_deliquescence(T));

    auto X_sorp = p_V_to_X(p_Ads);

    // Only data with X < 50 have been fitted.
    // Higher X was not observed in experiments
    if (X_sorp > 50.0)
        X_sorp = 50.0;

    const double a0 = 959.889873;
    const double a1 = -6.14849284;
    const double a2 = 0.00985078;

    const double b0 = 42.3974494;
    const double b1 = -0.24087759;
    const double b2 = 0.00034151;

    const double c0 = 9479.06005;
    const double c1 = -62.9458832;
    const double c2 = 0.10470995;

    const double d0 = 1228.32018;
    const double d1 = -8.08494628;
    const double d2 = 0.01332128;

    // a0:   959.889873 +/- 466.1101 (48.56%) (init= 96.14133)
    // a1:  -6.14849284 +/- 3.021169 (49.14%) (init=-0.5897181)
    // a2:   0.00985078 +/- 0.004894 (49.68%) (init= 0.00091039)
    // b0:   42.3974494 +/- 15.09632 (35.61%) (init= 52.81675)
    // b1:  -0.24087759 +/- 0.093532 (38.83%) (init=-0.3084013)
    // b2:   0.00034151 +/- 0.000145 (42.45%) (init= 0.00045097)
    // c0:   9479.06005 +/- 1.92e+05 (2028.67%) (init= 9603.371)
    // c1:  -62.9458832 +/- 1.25e+03 (1984.61%) (init=-63.07857)
    // c2:   0.10470995 +/- 2.028163 (1936.94%) (init= 0.1038204)
    // d0:   1228.32018 +/- 1.25e+04 (1020.93%) (init= 1271.209)
    // d1:  -8.08494628 +/- 81.46205 (1007.58%) (init=-8.265628)
    // d2:   0.01332128 +/- 0.132251 (992.78%) (init= 0.01345437)

    const double a = a0 + a1 * T + a2 * T * T;
    const double b = b0 + b1 * T + b2 * T * T;
    const double c = c0 + c1 * T + c2 * T * T;
    const double d = d0 + d1 * T + d2 * T * T;

    const double nu = a + b * X_sorp + std::exp((X_sorp - c) / d);
    return nu;
}

//! Water loading of the not included salt. ThNo. Eq. (4.6)
static double CaCl2_loading_not_included_salt(const double s_ni,
                                              const double p_Ads,
                                              const double T)
{
    auto const nu = CaCl2_water_fraction_hydrate_solution(p_Ads, T);
    return s_ni * nu;
}

//! Donnan-Equilibrium. ThNo Eq. (4.7), (4.23), Tab. 4.17
static double CaCl2_included_salt_Donnan(const double nu)
{
    const double e0 = 278.3;
    const double e1 = -1.264;

    return e0 * std::pow(nu, e1);
}

static double CaCl2_CaX_get_molar_mass_composite(const double salt_loading)
{
    // fitted Data for CaCl2/CaX
    return 113.529 * salt_loading + 13256.397;  // g/mol
}

//! ThNo Eqs. (4.4), (4.11)
static double CaCl2_CaX_loading_reduction_included_salt(const double s_i)
{
    auto const beta_exp = 2.8;  // CaCl_2, ThNo Tab 4.12
    return beta_exp * s_i;
}

double CaCl2_CaX_get_C_eq_zeolite(const double p_V,
                                  const double T,
                                  const double Delta_C_i)
{
    auto const rho_SR_eq = CaX_equilibrium.getEquilibriumDensity(p_V, T);
    auto const C_eq = CaX_equilibrium.getLoading(rho_SR_eq);
    auto const C0 = CaX_equilibrium.getMaximumLoading(T);

    // Salz blockiert Teil der Pore. Nur effektiv oberhalb bestimmter
    // Wasserbeladung
    if (C_eq + Delta_C_i <= C0)
    {
        return C_eq;
    }
    else
    {
        return C0 - Delta_C_i;
    }
}

static Eigen::Array2d CaCl2_CaX_get_equilibrium(const double p_Ads,
                                                const double T_Ads,
                                                const double salt_loading)
{
    if (p_Ads < 0.0)
    {
        return {0.0, 0.0};
    }

    auto const p_VD = CaCl2_p_V_deliquescence(T_Ads);

    const double s_i_max = 1e8;  // CaCl_2 : >23, Tab. 4.11

    double s_i, Delta_w_ni;

    if (p_Ads <= p_VD)
    {  // no deliquescence
        if (salt_loading <= s_i_max)
        {
            s_i = salt_loading;
            Delta_w_ni = 0.0;
        }
        else
        {
            OGS_FATAL("Fuer CaCl_2 wird s_i_max im Experiment nie erreicht");
            s_i = s_i_max;
            auto const s_ni = salt_loading - s_i_max;
            Delta_w_ni = CaCl2_loading_not_included_salt(s_ni, p_Ads, T_Ads);
        }
    }
    else
    {  // deliquescence
        auto const nu = CaCl2_water_fraction_hydrate_solution(p_Ads, T_Ads);
        s_i = CaCl2_included_salt_Donnan(nu);

        if (s_i >= salt_loading)
        {
            // no solution in mesopores
            s_i = salt_loading;
            Delta_w_ni = 0.0;
        }
        else
        {
            auto const s_ni = salt_loading - s_i;
            Delta_w_ni = CaCl2_loading_not_included_salt(s_ni, p_Ads, T_Ads);
        }
    }

    const double M_Comp = CaCl2_CaX_get_molar_mass_composite(salt_loading);
    const double M_H2O = 18.02;  // g/mol
    const double w_to_C = M_H2O / M_Comp;

    // limit loading due to not included salt
    const double Delta_C_ni = std::min(2.0, Delta_w_ni * w_to_C);

    auto const Delta_C_i =
        CaCl2_CaX_loading_reduction_included_salt(s_i) * w_to_C;
    auto const C_eq_Z = CaCl2_CaX_get_C_eq_zeolite(p_Ads, T_Ads, Delta_C_i);
    // auto const C_eq = C_eq_Z + Delta_C_ni;  // ~ ThNo Eq. (4.3)

    //    if (C_eq > 2.0)
    //        INFO("p_V = %g, T = %g, C_Z = %g, C = %g",
    //             p_Ads, T_Ads, C_eq_Z, C_eq);
    return {C_eq_Z, Delta_C_ni};
}

}  // anonymous namespace

namespace MaterialLib
{
void SorptionEquilibriumCaCl2CaX_1_SeparateLoadingLoading::getEquilibrium(
    const double p_Ads,
    const double T_Ads,
    ReactiveSolidState const& /*solid_state*/,
    ReactiveSolidState& equilibrium) const
{
    equilibrium.conversion() =
        CaCl2_CaX_get_equilibrium(p_Ads, T_Ads, _salt_loading);
}

HeatOfReactionData
SorptionEquilibriumCaCl2CaX_1_SeparateLoadingLoading::getHeatOfReaction(
    const double p_Ads,
    const double T_Ads,
    const ReactiveSolidState* const state) const
{
    // TODO use other enthalpies?
    return (Eigen::VectorXd(2)
                << CaX_equilibrium.getHeatOfReaction(p_Ads, T_Ads, state),
            waterEnthalpyOfEvaporation(T_Ads))
        .finished();
}

std::vector<NumLib::NamedFunction>
SorptionEquilibriumCaCl2CaX_1_SeparateLoadingLoading::getNamedFunctions() const
{
    auto f1 = [&](double p_V, double T) -> double {
        return CaCl2_CaX_get_equilibrium(p_V, T, _salt_loading)[0];
    };
    auto f2 = [&](double p_V, double T) -> double {
        return CaCl2_CaX_get_equilibrium(p_V, T, _salt_loading)[1];
    };
    return {{"equilibrium_loading1",
             {"adsorptive_partial_pressure", "adsorbent_temperature"},
             BaseLib::easyBind(std::move(f1))},
            {"equilibrium_loading2",
             {"adsorptive_partial_pressure", "adsorbent_temperature"},
             BaseLib::easyBind(std::move(f2))}};
}

std::unique_ptr<SorptionEquilibriumCaCl2CaX_1_SeparateLoadingLoading>
createSorptionEquilibriumCaCl2CaX_1_SeparateLoadingLoading(
    BaseLib::ConfigTree const& config)
{
    config.checkConfigParameter("type", "CaCl2CaX_1_SeparateLoadingLoading");

    auto const salt_loading = config.getConfigParameter<double>("salt_loading");

    return std::unique_ptr<
        SorptionEquilibriumCaCl2CaX_1_SeparateLoadingLoading>{
        new SorptionEquilibriumCaCl2CaX_1_SeparateLoadingLoading{salt_loading}};
}

}  // namespace MaterialLib
