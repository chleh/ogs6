#include "MathLib/LinAlg/Sparse/LOLMatrix.h"
#include "MathLib/LinAlg/Sparse/LOLCSRTools.h"

#include <gtest/gtest.h>

#include <random>
#include <numeric>



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

