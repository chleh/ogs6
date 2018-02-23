/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Dense>

namespace MaterialLib
{
/*
class HeatOfReactionData
{
public:
    virtual ~HeatOfReactionData() = default;
};
*/
using HeatOfReactionData = Eigen::VectorXd;

}  // namespace MaterialLib
