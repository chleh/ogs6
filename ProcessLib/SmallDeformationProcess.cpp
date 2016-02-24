/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SmallDeformationProcess-fwd.h"
#include "SmallDeformationProcess.h"

namespace ProcessLib
{

template class SmallDeformationProcess<GlobalSetupType, 1>;
template class SmallDeformationProcess<GlobalSetupType, 2>;
template class SmallDeformationProcess<GlobalSetupType, 3>;

}   // namespace ProcessLib
