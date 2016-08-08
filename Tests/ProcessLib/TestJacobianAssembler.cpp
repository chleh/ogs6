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

struct MatDiagX
{
    // M = diag(x1, x2, x3, ...)
    static void getMat(std::vector<double> const& x_data,
                       std::vector<double>& mat_data)
    {
        auto mat =
            MathLib::toZeroedMatrix(mat_data, x_data.size(), x_data.size());
        mat.diagonal().noalias() = MathLib::toVector(x_data, x_data.size());
    }

    // dM/dx * y = diag(y1, y2, y3, ...)
    static void getDMatDxTimesY(std::vector<double> const& x_data,
                                std::vector<double> const& y_data,
                                std::vector<double>& dMdxTy_data)
    {
        auto dMdxTy =
            MathLib::toZeroedMatrix(dMdxTy_data, x_data.size(), x_data.size());
        dMdxTy.diagonal().noalias() = MathLib::toVector(y_data, y_data.size());
    }
};

struct VecX
{
    // v = (x1, x2, x3, ...)
    static void getVec(std::vector<double> const& x_data,
                       std::vector<double>& vec_data)
    {
        vec_data = x_data;
    }

    // v = Identity
    static void getDVecDx(std::vector<double> const& x_data,
                          std::vector<double>& mat_data)
    {
        auto mat =
            MathLib::toZeroedMatrix(mat_data, x_data.size(), x_data.size());
        mat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                            Eigen::RowMajor>::Identity(x_data.size(),
                                                       x_data.size());
    }
};

template <typename Mat>
class LocalAssemblerM final : public ProcessLib::LocalAssemblerInterface
{
public:
    void assemble(double const /*t*/, std::vector<double> const& local_x,
                  std::vector<double>& local_M_data,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_b_data*/) override
    {
        Mat::getMat(local_x, local_M_data);
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

        Mat::getDMatDxTimesY(local_x, local_xdot, local_Jac_data);

        auto local_Jac =
            MathLib::toMatrix(local_Jac_data, local_x.size(), local_x.size());
        auto local_M =
            MathLib::toMatrix(local_M_data, local_x.size(), local_x.size());
        local_Jac.noalias() += dxdot_dx * local_M;
    }
};

template <typename Mat>
class LocalAssemblerK final : public ProcessLib::LocalAssemblerInterface
{
public:
    void assemble(double const /*t*/, std::vector<double> const& local_x,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& local_K_data,
                  std::vector<double>& /*local_b_data*/) override
    {
        Mat::getMat(local_x, local_K_data);
    }

    void assembleWithJacobian(double const t,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_xdot,
                              const double /*dxdot_dx*/, const double dx_dx,
                              std::vector<double>& local_M_data,
                              std::vector<double>& local_K_data,
                              std::vector<double>& local_b_data,
                              std::vector<double>& local_Jac_data) override
    {
        assemble(t, local_x, local_M_data, local_K_data, local_b_data);

        Mat::getDMatDxTimesY(local_x, local_xdot, local_Jac_data);

        auto local_Jac =
            MathLib::toMatrix(local_Jac_data, local_x.size(), local_x.size());
        auto local_K =
            MathLib::toMatrix(local_K_data, local_x.size(), local_x.size());
        local_Jac.noalias() += dx_dx * local_K;
    }
};

template <typename Vec>
class LocalAssemblerB final : public ProcessLib::LocalAssemblerInterface
{
public:
    void assemble(double const /*t*/, std::vector<double> const& local_x,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& local_b_data) override
    {
        Vec::getVec(local_x, local_b_data);
    }

    void assembleWithJacobian(double const t,
                              std::vector<double> const& local_x,
                              std::vector<double> const& /*local_xdot*/,
                              const double /*dxdot_dx*/, const double /*dx_dx*/,
                              std::vector<double>& local_M_data,
                              std::vector<double>& local_K_data,
                              std::vector<double>& local_b_data,
                              std::vector<double>& local_Jac_data) override
    {
        assemble(t, local_x, local_M_data, local_K_data, local_b_data);

        Vec::getDVecDx(local_x, local_Jac_data);
    }
};

TEST(ProcessLib, CentralDifferencesJacobianAssembler)
{
    ProcessLib::AnalyticalJacobianAssembler jac_asm_ana;
    ProcessLib::CentralDifferencesJacobianAssembler jac_asm_cd;
    LocalAssemblerM<MatDiagX> loc_asm;

    double const eps = std::numeric_limits<double>::epsilon();
    double const eps_cd = 1e-8;

    std::vector<double> x, xdot, M_data_cd, K_data_cd, b_data_cd, Jac_data_cd,
        M_data_ana, K_data_ana, b_data_ana, Jac_data_ana;
    double const dxdot_dx = 0.25, dx_dx = 0.75, t = 0.0;

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
