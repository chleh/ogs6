/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "logog/include/logog.hpp"

#include "PardisoLinearSolver.h"


namespace MathLib
{

PardisoLinearSolver::PardisoLinearSolver(const BaseLib::ConfigTree &/*option*/)
{
    /*
    boost::optional<std::string> solver_type
            = option.get_optional<std::string>("solver_type");
    if (solver_type)
    {

    }
    else
    {
        ERR("option <solver_type> not given");
    }
    */
}

void PardisoLinearSolver
::solve(LOLMatrix &A, DenseVector<double> &b, DenseVector<double> &x)
{
    //
}

}
