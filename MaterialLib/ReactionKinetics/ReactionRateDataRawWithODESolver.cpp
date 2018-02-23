#include "ReactionRateDataRawWithODESolver.h"

namespace MaterialLib
{
void ReactionRateDataRawWithODESolver::preTimestep(
    ReactiveSolidState const& solid_state,
    ReactiveSolidRate const& /*reaction_rate*/)
{
    solid_density_prev_ts = solid_state.conversion()[0];
}

}  // namespace MaterialLib
