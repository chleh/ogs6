
#include <gtest/gtest.h>

#include <iostream>
#include "ProcessLib/TCHSStokes/Material/TCHSStokesMaterial.h"

TEST(ProcessLibTCHSStokes, SolidHeatCapacityZeoliteCaAWaterVucelic)
{
    using namespace ProcessLib::TCHSStokes::Material;

    double const rho_SR_dry = 1.0;
    double const cp_zeo_dry = 0.0;
    SolidHeatCapacityZeoliteCaAWaterVucelic CaA(cp_zeo_dry, rho_SR_dry);

    double const rho_SR = 2.0 * rho_SR_dry;  // loading == 1

    std::cout << std::erf(3e-2 * (300.0 - 220.0)) << " "
              << std::erf(-3e-2 * (300.0 - 220.0)) << "\n\n";

    for (double T = 220; T < 500; T += 10.0)
    {
        std::cout << T << "\t" << 2.0 * CaA.getSpecificHeatCapacity(rho_SR, T)
                  << '\n';
    }
}
