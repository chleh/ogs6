#pragma once

#include "EquationSystem.h"

namespace NumLib
{
//! Adds refinement support to the \c EquationSystem interface.
class EquationSystemWithRefinementSupport : public EquationSystem
{
public:
    //! Performs the error estimation.
    //! \param x The solution vector.
    //! \param global_relative_error The estimated global relative error of the
    //! solution vector.
    //! \return A pointer to the vector of local error estimates, one entry per
    //! finite element.
    virtual GlobalVector const* estimateError(GlobalVector const& x,
                                              double& global_relative_error)
    {
        (void)x;
        (void)global_relative_error;
        return nullptr;
    }

    //! Trigger mesh adaptation.
    //! \param elements_for_refinement A flag for each mesh element if it should
    //! be kept (0), refined (1) or coarsened (-1)
    //! \return true if the mesh has been adapted, false otherwise.
    virtual bool refine(std::vector<char> const& elements_for_refinement)
    {
        (void)elements_for_refinement;
        return false;
    }
};

}  // namespace NumLib
