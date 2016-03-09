#include "MathLib/LinAlg/Sparse/LOLMatrix.h"
#include "MathLib/LinAlg/Sparse/LOLCSRTools.h"

#include "MathLib/LinAlg/Dense/DenseVector.h"

#include "MathLib/LinAlg/ApplyKnownSolution.h"
#include "MathLib/LinAlg/Scaling.h"

#include <gtest/gtest.h>

#include <random>
#include <numeric>

#include <iostream>



TEST(MathLibSparseLOLMatrix, BasicTest)
{
    const unsigned nrows = 20;
    const unsigned ncols = 10;

    double values[nrows][ncols];

    std::random_device random_dev;
    std::mt19937 shuffle_gen(random_dev());

    std::uniform_real_distribution<double> values_dist(1.0, 11.0);

    for (unsigned r=0; r<nrows; ++r) {
        for (unsigned c=0; c<ncols; ++c) {
            values[r][c] = values_dist(shuffle_gen);
        }
    }


    std::vector<unsigned> row_idcs(nrows);
    std::vector<unsigned> col_idcs(ncols);
    std::iota(row_idcs.begin(), row_idcs.end(), 0);
    std::iota(col_idcs.begin(), col_idcs.end(), 0);


    MathLib::LOLMatrix lolmat(nrows, ncols);

    ASSERT_EQ(nrows, lolmat.getNRows());
    ASSERT_EQ(ncols, lolmat.getNCols());


    // fill coo matrix in random order
    std::shuffle(row_idcs.begin(), row_idcs.end(), shuffle_gen);

    for (unsigned r : row_idcs) {
        std::shuffle(col_idcs.begin(), col_idcs.end(), shuffle_gen);

        for (unsigned c : col_idcs) {
            lolmat.add(r, c, values[r][c]);
        }
    }

    // get values from coo matrix in random order
    std::shuffle(row_idcs.begin(), row_idcs.end(), shuffle_gen);

    for (unsigned r : row_idcs) {
        std::shuffle(col_idcs.begin(), col_idcs.end(), shuffle_gen);

        for (unsigned c : col_idcs) {
            EXPECT_EQ(values[r][c], lolmat.get(r, c));
        }
    }



    MathLib::CSRMatrix csrmat = toCSR(lolmat);

    ASSERT_EQ(lolmat.getNCols(), csrmat.getNCols());
    ASSERT_EQ(lolmat.getNRows(), csrmat.getNRows());

    ASSERT_EQ(lolmat, toLOL(csrmat));

}

#ifdef OGS_USE_MKL
TEST(MathLibSparseLOLMatrix, KnownSolution)
{
    const unsigned nrows = 3;
    const unsigned ncols = 3;

    double values[3][3] = {
        { 1.0, 2.0, 3.0 },
        { 4.0, 5.0, 7.0 },
        { 6.0, 8.0, 9.0 }
    };

    double rvalues[3] = { 1.0, 1.0, 1.0 };

    MathLib::LOLMatrix lolmat(3, 3);
    MathLib::DenseVector<double> rhs(3);

    for (unsigned r=0; r<nrows; ++r) {
        for (unsigned c=0; c<ncols; ++c) {
            lolmat.setValue(r, c, values[r][c]);
        }
        rhs[r] = rvalues[r];
    }

    std::vector<std::size_t> known_ids{{1}};
    std::vector<double> known_vals{{1.0}};

    MathLib::applyKnownSolution(lolmat, rhs, known_ids, known_vals);

    double values_result[3][3] = {
        { 1.0, 0.0, 3.0 },
        { 0.0, 1.0, 0.0 },
        { 6.0, 0.0, 9.0 }
    };

    double rvalues_result[3] = { -1.0, 1.0, -7.0 };

    for (unsigned r=0; r<nrows; ++r) {
        for (unsigned c=0; c<ncols; ++c) {
            EXPECT_EQ(values_result[r][c], lolmat.get(r, c));
        }
        EXPECT_EQ(rvalues_result[r], rhs[r]);
    }
}

TEST(MathLibSparseLOLMatrix, Scaling)
{
    const unsigned nrows = 3;
    const unsigned ncols = 3;

    double values[3][3] = {
        { 2.0, 2.0, 3.0 },
        { 4.0, 4.0, 7.0 },
        { 6.0, 8.0, 8.0 }
    };

    double rvalues[3] = { 1.0, 1.0, 1.0 };

    MathLib::LOLMatrix lolmat(3, 3);
    MathLib::DenseVector<double> rhs(3);

    for (unsigned r=0; r<nrows; ++r) {
        for (unsigned c=0; c<ncols; ++c) {
            lolmat.setValue(r, c, values[r][c]);
        }
        rhs[r] = rvalues[r];
    }

    MathLib::scaleDiagonal(lolmat, rhs);

    double values_result[3][3] = {
        { 1.0,  1.0, 1.5  },
        { 1.0,  1.0, 1.75 },
        { 0.75, 1.0, 1.0  }
    };

    double rvalues_result[3] = { 0.5, 0.25, 0.125 };

    for (unsigned r=0; r<nrows; ++r) {
        for (unsigned c=0; c<ncols; ++c) {
            EXPECT_EQ(values_result[r][c], lolmat.get(r, c));
        }
        EXPECT_EQ(rvalues_result[r], rhs[r]);
    }
}
#endif
