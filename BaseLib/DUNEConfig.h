#pragma once

#include "Config.h"

#if OGS_USE_DUNE
#define ENABLE_UG 1
#include "dune-config.h"
#if HAVE_UG
#pragma message "HAVE_UG set"
#else
#pragma message "HAVE_UG not set"
#endif

namespace Dune
{
template <int d>
class UGGrid;
}

namespace BaseLib
{
template <int d>
using DUNEGridType = Dune::UGGrid<d>;

}  // namespace BaseLib

#endif  // OGS_USE_DUNE
