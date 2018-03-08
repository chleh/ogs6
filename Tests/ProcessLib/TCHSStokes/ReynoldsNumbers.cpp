
#include <gtest/gtest.h>

#include <iostream>
#include "ProcessLib/TCHSStokes/Material/TCHSStokesMaterial.h"

TEST(ProcessLibTCHSStokes, ReynoldsNumbers)
{
    using namespace ProcessLib::TCHSStokes::Material;

    double const p = 1e5;
    double const xmV = 0.02;
    double const d_pel = 0.002;

    FluidViscosityMixtureWaterNitrogen mu_law;
    FluidDensityMixtureWaterNitrogen rho_law;

    auto Re = [=](double rho, double mu, double v) {
        return rho * v * d_pel / mu;
    };

    for (auto mass_flux : {0.3, 0.445})
    {
        for (double T = 20 + 273.15; T < 180 + 275.15; T += 10.0)
        {
            auto const rho = rho_law.getDensity(p, T, xmV);
            auto const mu = mu_law.getViscosity(p, T, xmV);
            auto const v = mass_flux / rho;

            std::cout << mass_flux << "\t" << T << "\t" << Re(rho, mu, v)
                      << '\n';
        }
    }
}
