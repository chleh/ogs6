/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "IncompressibleStokesBrinkmanProcess.h"
#include "IncompressibleStokesBrinkmanProcess-impl.h"

namespace ProcessLib
{
namespace IncompressibleStokesBrinkman
{
template class IncompressibleStokesBrinkmanProcess<2>;
template class IncompressibleStokesBrinkmanProcess<3>;

}  // namespace IncompressibleStokesBrinkman
}  // namespace ProcessLib
