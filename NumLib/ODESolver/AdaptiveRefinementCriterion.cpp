#include "AdaptiveRefinementCriterion.h"

#include <numeric>

#include "BaseLib/ConfigTree.h"

#include "EquationSystemWithRefinementSupport.h"

namespace NumLib
{
bool AdaptiveRefinementCriterion::refine(EquationSystem& sys,
                                         const GlobalVector& x)
{
    auto sys_ = dynamic_cast<EquationSystemWithRefinementSupport*>(&sys);
    if (sys_ == nullptr)
    {
        DBUG("Equation system not suitable for refinement.");
        return false;
    }

    double global_relative_error;
    auto const& element_errors = *sys_->estimateError(x, global_relative_error);

    if (global_relative_error <= _reltol)
    {
        // TODO [DUNE] maybe check anyway if we would only want to coarsen.
        INFO("no adaptation necessary. error (%g) < tolerance (%g).",
             global_relative_error, _reltol);
        return false;
    }
    auto const num_elements = static_cast<std::size_t>(element_errors.size());
    auto const sqrt_num_elements = std::sqrt(num_elements);

    std::vector<std::pair<double, std::size_t>> indicators_and_element_ids;
    indicators_and_element_ids.reserve(num_elements);
    for (std::size_t i = 0; i < num_elements; ++i)
    {
        // Cf. Grätsch, T., Bathe, K.-J., 2005. A posteriori error estimation
        // techniques in practical finite element analysis. Computers &
        // Structures 83, 235–265. doi:10.1016/j.compstruc.2004.08.011
        // and
        // Zienkiewicz, O.C., Zhu, J.Z., 1987. A simple error estimator and
        // adaptive procedure for practical engineering analysis. International
        // Journal for Numerical Methods in Engineering 24, 337–357.
        // doi:10.1002/nme.1620240206
        auto const indicator = sqrt_num_elements * element_errors[i] / _reltol;
        indicators_and_element_ids.emplace_back(indicator, i);
    }

    std::sort(indicators_and_element_ids.begin(),
              indicators_and_element_ids.end(),
              [](auto a, auto b) { return a.first > b.first; });

    std::vector<char> mark_for_refinement(num_elements);

    std::size_t num_refine = 0;
    std::size_t num_coarsen = 0;
    for (std::size_t i = 0; i < num_elements; ++i)
    {
        auto const& ei = indicators_and_element_ids[i];
        if (ei.first > 1)
        {
            mark_for_refinement[ei.second] = 1;
            ++num_refine;
        }
        else if (ei.first < 1)
        {
            mark_for_refinement[ei.second] = -1;
            ++num_coarsen;
        }
    }

    if (num_refine == 0 && num_coarsen == 0)
    {
        DBUG("No elements found to refine/coarsen.");
        return false;
    }

    if (num_refine != 0)
    {
        auto const num_elements_to_mark_for_refinement =
            static_cast<std::size_t>(std::ceil(_refine_fraction * num_refine));
        DBUG("marking %d (= %g * %d) out of %d elements for refinement",
             num_elements_to_mark_for_refinement, _refine_fraction, num_refine,
             num_elements);

        // unmark elements marked in excess
        for (std::size_t ii = num_elements_to_mark_for_refinement;
             ii < num_refine;
             ++ii)
        {
            auto const& ei = indicators_and_element_ids[ii];
            mark_for_refinement[ei.second] = 0;
        }
    }
    if (num_coarsen != 0)
    {
        // TODO [DUNE] maybe drop the percentage for coarsening or use a
        // separate setting?
        auto const num_elements_to_mark_for_coarsening =
            static_cast<std::size_t>(std::ceil(_refine_fraction * num_coarsen));
        DBUG("marking %d (= %g * %d) out of %d elements for coarsening",
             num_elements_to_mark_for_coarsening, _refine_fraction, num_coarsen,
             num_elements);

        // unmark elements marked in excess
        auto num_coarsen_eff = num_coarsen;
        for (std::size_t ii = 0;
             ii < num_coarsen - num_elements_to_mark_for_coarsening;
             ++ii)
        {
            auto const iii = num_elements - num_coarsen + ii;
            auto const& ei = indicators_and_element_ids[iii];
            mark_for_refinement[ei.second] = 0;
            --num_coarsen_eff;
        }
        assert(num_coarsen_eff == num_elements_to_mark_for_coarsening);
    }

    // TODO [DUNE] maybe add setting to not adapt the mesh if only coarsening
    return sys_->refine(mark_for_refinement);
}

std::unique_ptr<AdaptiveRefinementCriterion> createAdaptiveRefinementCriterion(
    BaseLib::ConfigTree const& config)
{
    auto const reltol = config.getConfigParameter<double>("reltol");
    auto const refine_fraction =
        config.getConfigParameter<double>("refine_fraction");

    return std::make_unique<AdaptiveRefinementCriterion>(reltol,
                                                         refine_fraction);
}

}  // namespace NumLib
