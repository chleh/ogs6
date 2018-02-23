/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TCHSStokesProcess.h"
#include "TCHSStokesProcess-impl.h"

namespace ProcessLib
{
namespace TCHSStokes
{
template class TCHSStokesProcess<2>;
template class TCHSStokesProcess<3>;

}  // namespace TCHSStokes
}  // namespace ProcessLib
