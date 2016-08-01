/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_BOUNDARYCONDITION_H
#define PROCESSLIB_BOUNDARYCONDITION_H

#include "NumLib/NumericsConfig.h"

namespace ProcessLib
{

class BoundaryCondition
{
public:
    virtual void apply(const double t, GlobalVector const& x, GlobalMatrix& K,
                       GlobalVector& b) = 0;
    virtual ~BoundaryCondition() = default;
};

}  // ProcessLib

#endif  // PROCESSLIB_BOUNDARYCONDITION_H
