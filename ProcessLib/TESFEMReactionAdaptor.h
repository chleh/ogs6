/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>
#include <memory>
#include "TESFEM-data.h"



namespace ProcessLib
{
namespace TES
{

template<typename Traits>
class LADataNoTpl;

template<typename Traits>
class TESFEMReactionAdaptor
{
public:
    virtual bool
    checkBounds(std::vector<double> const& localX,
                std::vector<double> const& localX_pts) = 0;

    virtual void
    initReaction(const unsigned int_pt) = 0;

    virtual void
    preZerothTryAssemble() = 0;

    // TODO: remove
    virtual double getReactionDampingFactor() const = 0;

    virtual ~TESFEMReactionAdaptor() = default;

    static std::unique_ptr<TESFEMReactionAdaptor<Traits> >
    newInstance(LADataNoTpl<Traits>& data);
};


template<typename Traits>
class TESFEMReactionAdaptorAdsorption : public TESFEMReactionAdaptor<Traits>
{
public:
    TESFEMReactionAdaptorAdsorption(LADataNoTpl<Traits>& data);

    bool checkBounds(const std::vector<double> &localX, const std::vector<double> &localX_pts)
    override;

    void initReaction(const unsigned int_pt) override {
        initReaction_slowDownUndershootStrategy(int_pt);
    }

    void preZerothTryAssemble() override;

private:
    void initReaction_slowDownUndershootStrategy(const unsigned int_pt);

    /// returns estimated equilibrium vapour pressure
    /// based on a local (i.e. no diffusion/advection) balance
    /// of adsorbate loading and vapour partial pressure
    double estimateAdsorptionEquilibrium(const double p_V0, const double C0) const;

    double getReactionDampingFactor() const override {
        return _reaction_damping_factor;
    }

    double _reaction_damping_factor = 1.0;
    std::vector<bool> _bounds_violation;

    LADataNoTpl<Traits>& _data;
};

}
}

#include "TESFEMReactionAdaptor-impl.h"
