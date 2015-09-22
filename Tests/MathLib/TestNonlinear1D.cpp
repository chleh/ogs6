#include <gtest/gtest.h>

#include <cstdio>

#include "MathLib/Nonlinear/Root1D.h"

namespace
{

double f(double x)
{
    return x;
}

double f2(double x)
{
    return x*x-1;
}

} // anonymous namespace


TEST(MathLibNonlinear1DTest, Bisect)
{

    MathLib::Nonlinear::Bisect bi(f2, -0.1, 1.1);

    std::printf("x ~ %14.7g, range = %14.7g\n", bi.get_result(), bi.get_range());

    for (unsigned n=0; n<10; ++n)
    {
        bi.step(1);
        std::printf("x ~ %14.7g, range = %14.7g\n", bi.get_result(), bi.get_range());
    }
}


template <typename T>
class RegulaFalsiTest : public ::testing::Test
{
};


using namespace MathLib::Nonlinear;
typedef ::testing::Types<Unmodified, Illinois, Pegasus, AndersonBjorck> RegulaFalsiTypes;
TYPED_TEST_CASE(RegulaFalsiTest, RegulaFalsiTypes);



TYPED_TEST(RegulaFalsiTest, QuadraticFunction)
{
    RegulaFalsi<TypeParam> rf(f2, -0.1, 1.1);

    std::printf(" 0 -- x ~ %14.7g, range = %14.7g\n", rf.get_result(), rf.get_range());

    for (unsigned n=0; n<10; ++n)
    {
        rf.step(1);
        std::printf("%2i --  x ~ %14.7g, range = %14.7g\n", n+1, rf.get_result(), rf.get_range());
    }
}


