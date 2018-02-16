/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "IncompressibleStokesBrinkmanModifiedProcess.h"
#include "IncompressibleStokesBrinkmanModifiedProcess-impl.h"

namespace ProcessLib
{
namespace IncompressibleStokesBrinkmanModified
{
template class IncompressibleStokesBrinkmanModifiedProcess<2>;
template class IncompressibleStokesBrinkmanModifiedProcess<3>;

}  // namespace IncompressibleStokesBrinkmanModified
}  // namespace ProcessLib
