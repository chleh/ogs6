/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "BaseLib/ConfigTree.h"
#include "MathLib/LinAlg/Sparse/LOLMatrix.h"
#include "MathLib/LinAlg/Dense/DenseVector.h"

namespace MathLib
{

class PardisoLinearSolverImpl;

class PardisoLinearSolver final
{
public:
    /**
     * Constructor
     * @param A           Coefficient matrix object
     * @param solver_name A name used as a prefix for command line options
     *                    if there are such options available.
     * @param option      A pointer to a linear solver option. In case you omit
     *                    this argument, default settings follow those of
     *                    LisOption struct.
     */
    explicit PardisoLinearSolver(BaseLib::ConfigTree const& option);


    /**
     * solve a given linear equations
     *
     * @param b     RHS vector
     * @param x     Solution vector
     */
    void solve(LOLMatrix& A, DenseVector<double> &b, DenseVector<double> &x);

    ~PardisoLinearSolver();

private:
    std::unique_ptr<PardisoLinearSolverImpl> _data;
};

}
