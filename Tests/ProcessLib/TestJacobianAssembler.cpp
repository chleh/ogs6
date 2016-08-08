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

//! Fills a vector with values whose absolute value is between \c abs_min and
//! \c abs_max.
void fillRandomlyConstrainedAbsoluteValues(std::vector<double>& xs,
                                           double const abs_min,
                                           double const abs_max)
{
    std::random_device rd;
    std::mt19937 random_number_generator(rd());
    double const abs_range = abs_max - abs_min;
    std::uniform_real_distribution<double> rnd(abs_min,
                                               abs_min + 2.0 * abs_range);

    for (auto& x : xs) {
        // v in [ abs_min, abs_min + 2 abs_range ]
        auto v = rnd(random_number_generator);
        if (v > abs_max) {
            // produce negative values
            // (v - abs_range) in [ abs_min, abs_max ]
            v = -(v - abs_range);
        }
        x = v;
    }
}

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

struct MatVecDiagX
{
    using Mat = MatDiagX;
    using Vec = VecX;
};

template <typename MatVec>
class LocalAssemblerM final : public ProcessLib::LocalAssemblerInterface
{
public:
    void assemble(double const /*t*/, std::vector<double> const& local_x,
                  std::vector<double>& local_M_data,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_b_data*/) override
    {
        MatVec::Mat::getMat(local_x, local_M_data);
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

        // dM/dx * xdot
        MatVec::Mat::getDMatDxTimesY(local_x, local_xdot, local_Jac_data);

        auto local_Jac =
            MathLib::toMatrix(local_Jac_data, local_x.size(), local_x.size());
        auto local_M =
            MathLib::toMatrix(local_M_data, local_x.size(), local_x.size());
        local_Jac.noalias() += dxdot_dx * local_M;
    }

    static const bool asmM = true;
    static const bool asmK = false;
    static const bool asmb = false;
};

template <typename MatVec>
class LocalAssemblerK final : public ProcessLib::LocalAssemblerInterface
{
public:
    void assemble(double const /*t*/, std::vector<double> const& local_x,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& local_K_data,
                  std::vector<double>& /*local_b_data*/) override
    {
        MatVec::Mat::getMat(local_x, local_K_data);
    }

    void assembleWithJacobian(double const t,
                              std::vector<double> const& local_x,
                              std::vector<double> const& /*local_xdot*/,
                              const double /*dxdot_dx*/, const double dx_dx,
                              std::vector<double>& local_M_data,
                              std::vector<double>& local_K_data,
                              std::vector<double>& local_b_data,
                              std::vector<double>& local_Jac_data) override
    {
        assemble(t, local_x, local_M_data, local_K_data, local_b_data);

        // dK/dx * x
        MatVec::Mat::getDMatDxTimesY(local_x, local_x, local_Jac_data);

        auto local_Jac =
            MathLib::toMatrix(local_Jac_data, local_x.size(), local_x.size());
        auto local_K =
            MathLib::toMatrix(local_K_data, local_x.size(), local_x.size());
        local_Jac.noalias() += dx_dx * local_K;
    }

    static const bool asmM = false;
    static const bool asmK = true;
    static const bool asmb = false;
};

template <typename MatVec>
class LocalAssemblerB final : public ProcessLib::LocalAssemblerInterface
{
public:
    void assemble(double const /*t*/, std::vector<double> const& local_x,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& local_b_data) override
    {
        MatVec::Vec::getVec(local_x, local_b_data);
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

        // db/dx
        MatVec::Vec::getDVecDx(local_x, local_Jac_data);

        auto local_Jac =
            MathLib::toMatrix(local_Jac_data, local_x.size(), local_x.size());
        // J = -db/dx !!
        local_Jac = -local_Jac;
    }

    static const bool asmM = false;
    static const bool asmK = false;
    static const bool asmb = true;
};

template<class LocAsm>
struct ProcessLibCentralDifferencesJacobianAssembler : public ::testing::Test
{
    static void test()
    {
        ProcessLib::AnalyticalJacobianAssembler jac_asm_ana;
        ProcessLib::CentralDifferencesJacobianAssembler jac_asm_cd;
        LocAsm loc_asm;

        double const eps = std::numeric_limits<double>::epsilon();
        double const eps_cd = 1e-8;

        std::vector<double> x, xdot, M_data_cd, K_data_cd, b_data_cd, Jac_data_cd,
            M_data_ana, K_data_ana, b_data_ana, Jac_data_ana;
        double const dxdot_dx = 0.25, dx_dx = 0.75, t = 0.0;


        {
            std::random_device rd;
            std::mt19937 random_number_generator(rd());
            std::uniform_int_distribution<std::size_t> rnd(3, 64);

            auto const size = rnd(random_number_generator);
            x.resize(size);
            xdot.resize(size);

            // all components will be of order of magnitude one
            fillRandomlyConstrainedAbsoluteValues(x, 0.5, 1.5);
            fillRandomlyConstrainedAbsoluteValues(xdot, 0.5, 1.5);
        }

        jac_asm_cd.assembleWithJacobian(loc_asm, t, x, xdot, dxdot_dx, dx_dx, M_data_cd,
                                     K_data_cd, b_data_cd, Jac_data_cd);


        jac_asm_ana.assembleWithJacobian(loc_asm, t, x, xdot, dxdot_dx, dx_dx,
                                         M_data_ana, K_data_ana, b_data_ana,
                                         Jac_data_ana);

        if (LocAsm::asmM) {
            ASSERT_EQ(x.size()*x.size(), M_data_cd.size());
            ASSERT_EQ(x.size()*x.size(), M_data_ana.size());
            for (std::size_t i=0; i<x.size()*x.size(); ++i)
                EXPECT_NEAR(M_data_ana[i], M_data_cd[i], eps);
        }

        if (LocAsm::asmK) {
            ASSERT_EQ(x.size()*x.size(), K_data_cd.size());
            ASSERT_EQ(x.size()*x.size(), K_data_ana.size());
            for (std::size_t i=0; i<x.size()*x.size(); ++i)
                EXPECT_NEAR(K_data_ana[i], K_data_cd[i], eps);
        }

        if (LocAsm::asmb) {
            ASSERT_EQ(x.size(), b_data_cd.size());
            ASSERT_EQ(x.size(), b_data_ana.size());
            for (std::size_t i=0; i<x.size(); ++i)
                EXPECT_NEAR(b_data_ana[i], b_data_cd[i], eps);
        }

        ASSERT_EQ(x.size()*x.size(), Jac_data_cd.size());
        ASSERT_EQ(x.size()*x.size(), Jac_data_ana.size());
        for (std::size_t i=0; i<x.size()*x.size(); ++i)
            EXPECT_NEAR(Jac_data_ana[i], Jac_data_cd[i], eps_cd);
    }
};

typedef ::testing::Types<LocalAssemblerM<MatVecDiagX>,
                         LocalAssemblerK<MatVecDiagX>,
                         LocalAssemblerB<MatVecDiagX>>
    TestCases;

TYPED_TEST_CASE(ProcessLibCentralDifferencesJacobianAssembler, TestCases);

TYPED_TEST(ProcessLibCentralDifferencesJacobianAssembler, Test)
{
    TestFixture::test();
}
