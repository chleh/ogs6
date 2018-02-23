#pragma once

#include "ReactionRateData.h"

namespace MaterialLib
{
struct ReactionRateDataRawWithODESolver final : public ReactionRateData
{
public:
    void preTimestep(ReactiveSolidState const& solid_state,
                     ReactiveSolidRate const& reaction_rate) override;

    double solid_density_prev_ts;
};

}  // namespace MaterialLib
