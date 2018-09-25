#pragma once

#include <memory>

#include "BaseLib/Error.h"
#include "EquationSystem.h"

namespace BaseLib
{
class ConfigTree;
}  // namespace BaseLib

namespace NumLib
{
//! Decides which elements should be adaptively refined or coarsened.
class AdaptiveRefinementCriterion
{
public:
    AdaptiveRefinementCriterion(double reltol, double refine_fraction)
        : _reltol(reltol), _refine_fraction(refine_fraction)
    {
        OGS_ALWAYS_ASSERT(0 < _refine_fraction && _refine_fraction <= 1);
    }

    bool refine(EquationSystem& sys, GlobalVector const& x);

private:
    //! Global relative error tolerance. If the global tolerance is met no mesh
    //! adaptation will be triggered.
    double const _reltol;

    //! Only adapt this fraction of the elements that are candidates for either
    //! coarsening or refinement.
    double const _refine_fraction;
};

//! Creates a new instance of \c AdaptiveRefinementCriterion.
std::unique_ptr<AdaptiveRefinementCriterion> createAdaptiveRefinementCriterion(
    BaseLib::ConfigTree const& config);

}  // namespace NumLib
