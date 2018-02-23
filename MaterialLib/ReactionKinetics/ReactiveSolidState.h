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

#include "BaseLib/Error.h"
#include "MaterialLib/ReactionKinetics/ReactiveSolidRate.h"

namespace MaterialLib
{
class ReactiveSolidState
{
public:
    virtual ReactiveSolidState& operator+=(Eigen::VectorXd const&) = 0;
    /*{
        OGS_FATAL("not implemented");
    }*/
    virtual ReactiveSolidState& operator+=(
        Eigen::Ref<Eigen::VectorXd> const&) = 0;
    /*{
        OGS_FATAL("not implemented");
    }*/

    virtual Eigen::VectorXd const& conversion() const = 0;
    virtual Eigen::VectorXd& conversion() = 0;

    virtual ~ReactiveSolidState() = default;
};

}  // namespace MaterialLib
