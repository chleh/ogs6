/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <limits>
#include <gtest/gtest.h>
#include <logog/include/logog.hpp>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/AnalyticalJacobianAssembler.h"
#include "ProcessLib/CentralDifferencesJacobianAssembler.h"

class LocalAssemblerM : public ProcessLib::LocalAssemblerInterface
{
public:
    void assemble(
        double const /*t*/, std::vector<double> const& local_x,
        std::vector<double>& local_M_data, std::vector<double>& /*local_K_data*/,
        std::vector<double>& /*local_b_data*/) override
    {
        auto local_M = MathLib::toZeroedMatrix(local_M_data, local_x.size(),
                                               local_x.size());
        local_M.diagonal().noalias() =
            MathLib::toVector(local_x, local_x.size());
    }

    void assembleWithJacobian(double const t,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_xdot,
                              const double dxdot_dx, const double /*dx_dx*/,
                              std::vector<double>& local_M_data,
                              std::vector<double>& local_K_data,
                              std::vector<double>& local_b_data,
                              std::vector<double>& local_Jac_data) override
    {
        assemble(t, local_x, local_M_data, local_K_data, local_b_data);

        auto local_Jac = MathLib::toZeroedMatrix(local_Jac_data, local_x.size(),
                                                 local_x.size());
        local_Jac.diagonal().noalias() =
            MathLib::toVector(local_xdot, local_xdot.size());

        auto local_M =
            MathLib::toMatrix(local_M_data, local_x.size(), local_x.size());

        local_Jac.noalias() += dxdot_dx * local_M;
    }
};

TEST(ProcessLib, CentralDifferencesJacobianAssembler)
{
    ProcessLib::AnalyticalJacobianAssembler jac_asm_ana;
    ProcessLib::CentralDifferencesJacobianAssembler jac_asm_cd;
    LocalAssemblerM loc_asm;

    double const eps = std::numeric_limits<double>::epsilon();
    double const eps_cd = 1e-8;

    std::vector<double> x, xdot, M_data_cd, K_data_cd, b_data_cd, Jac_data_cd,
        M_data_ana, K_data_ana, b_data_ana, Jac_data_ana;
    double const dxdot_dx = 0.25, dx_dx = 0.75,
            t = 0.0;

    x = { 1.0, 2.0, 4.0 };
    xdot = { 8.0, 16.0, 32.0 };

    jac_asm_cd.assembleWithJacobian(loc_asm, t, x, xdot, dxdot_dx, dx_dx, M_data_cd,
                                 K_data_cd, b_data_cd, Jac_data_cd);

    ASSERT_EQ(x.size()*x.size(), M_data_cd.size());
    ASSERT_EQ(x.size()*x.size(), Jac_data_cd.size());

    jac_asm_ana.assembleWithJacobian(loc_asm, t, x, xdot, dxdot_dx, dx_dx,
                                     M_data_ana, K_data_ana, b_data_ana,
                                     Jac_data_ana);

    ASSERT_EQ(x.size()*x.size(), M_data_ana.size());
    ASSERT_EQ(x.size()*x.size(), Jac_data_ana.size());

    for (std::size_t i=0; i<x.size()*x.size(); ++i) {
        EXPECT_NEAR(M_data_ana[i], M_data_cd[i], eps);
        EXPECT_NEAR(Jac_data_ana[i], Jac_data_cd[i], eps_cd*xdot[i/x.size()]);
    }
}
