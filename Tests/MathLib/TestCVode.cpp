/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <gtest/gtest.h>

#include <array>

#include "MaterialsLib/adsorption/adsorption.h"
#include "ProcessLib/OdeSolver.h"

#include <cstdio>

using namespace ProcessLib::Ode;
using namespace Ads;

const unsigned NEQ = 2; // number of equations



const double T = 303.15; // K
const double R = 8.314;  // J/mol/K
const double M = 0.018;  // kg/mol
const double phi = 0.4;
const double rho_SR0 = 1150.0; // kg/m^3
// const double k = 6e-3;   // s^-1

Adsorption* ads;


void f(const double, double const*const y, double *const ydot)
{
    const double pV = y[0];
    const double C  = y[1];
    ydot[1] = ads->get_reaction_rate(pV, T, M, C);
    ydot[0] = -R*T/M/phi * (1.0-phi) * rho_SR0 * ydot[1];
}

TEST(MathLibCVodeTest, ZeoliteAdsorption)
{
    // initial values
    const double C  = 0.0;
    const double pV = 100.0; // Pa

    ads = Adsorption::newInstance(SolidReactiveSystem::Z13XBF_Hauer);

    ConcreteOdeSolver<2, CVodeSolverInternal> ode_solver;

    ode_solver.init();
    ode_solver.setTolerance(1e-8, 1e-6);

    // double const ic[] = { pV, C };
    ode_solver.setIC(0.0, { pV, C });

    ode_solver.solve(f, 0.00001);

    double const* y = ode_solver.getSolution();

    std::printf("pV: %g, C: %g\n", y[0], y[1]);
}
