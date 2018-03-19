/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TCHSNoStokesProcess.h"
#include "TCHSNoStokesProcess-impl.h"

namespace ProcessLib
{
namespace TCHSNoStokes
{
template class TCHSNoStokesProcess<2>;
// TODO saving compilation time
// template class TCHSNoStokesProcess<3>;

}  // namespace TCHSNoStokes
}  // namespace ProcessLib
