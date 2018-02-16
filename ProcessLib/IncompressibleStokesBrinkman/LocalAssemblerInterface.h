/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "ProcessLib/LocalAssemblerInterface.h"

namespace ProcessLib
{
namespace IncompressibleStokesBrinkman
{
struct LocalAssemblerInterface : public ProcessLib::LocalAssemblerInterface,
                                 public NumLib::ExtrapolatableElement
{
};

}  // namespace IncompressibleStokesBrinkman
}  // namespace ProcessLib
