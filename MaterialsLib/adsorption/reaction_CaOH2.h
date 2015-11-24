/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include<Eigen/Core>

#include "reaction.h"
#include "adsorption.h"


namespace ProcessLib
{
template<typename>
class TESFEMReactionAdaptorCaOH2;
}

namespace Ads
{

class ReactionCaOH2 final : public Reaction
{
public:
    explicit ReactionCaOH2(BaseLib::ConfigTree const& conf)
        : ode_solver_config{conf.get_child("ode_solver_config", BaseLib::ConfigTree{})}
    {
        /*auto const param = conf.get_optional<double>("reaction_enthalpy");
        if (param) {
            _enthalpy = *param;
        } else {
            ERR("<reaction_enthalpy> not specified.");
            std::abort();
        }*/
    }

    double get_enthalpy(const double /*p_Ads*/, const double /*T_Ads*/,
                        const double /*M_Ads*/) const override;

    double get_reaction_rate(const double /*p_Ads*/, const double /*T_Ads*/, const double /*M_Ads*/,
                             const double /*loading*/) const override;


private:
    void eval(double /*t*/, Eigen::VectorXd const& y, Eigen::VectorXd &dydx);
    void calculate_qR();
    void set_chemical_equilibrium();
    double Ca_hydration();


    static constexpr double R = Ads::GAS_CONST;  // [J/mol/K]
    double rho_s;    // solid phase density
    double rho_s_0;  // initial solid phase density
    double phi_solid; //solid volume fraction
    double p_gas;    // gas phase pressure;
    double p_r_g;    // pressure of H2O on gas phase;
    double p_eq;     // equilibrium pressure; // unit in bar
    double T_eq;     // equilibrium temperature;
    double T_s;      // solid phase temperature;
    double T;        // gas phase temperature;
    double qR;       // rate of solid density change;
    double x_react;    // mass fraction of water in gas phase;
    double X_D;      // mass fraction of dehydration (CaO) in the solid phase;
    double X_H;      // mass fraction of hydration in the solid phase;
    double dt;       // time step size;
    double rho_low; //lower density limit
    double rho_up; //upper density limit
    double reaction_enthalpy;
    double reaction_entropy;
    double M_carrier; //inert component molar mass
    double M_react; //reactive component molar mass
    double W0; //maximum specific adsorbed volume
    double C_eq; //equilibrium loading of sorbens
    double p_min; //minimum pressure for which we find a valid adsorption potential

    double tol_l;
    double tol_u;
    double tol_rho;

    double lower_solid_density_limit;
    double upper_solid_density_limit;

    const BaseLib::ConfigTree ode_solver_config;

    template<typename>
    friend class ProcessLib::TESFEMReactionAdaptorCaOH2;
};

}
