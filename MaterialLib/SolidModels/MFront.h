/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MechanicsBase.h"

namespace aster
{
// Cf. mfront/include/MFront/Aster/Aster.hxx and
// https://www.viva64.com/en/a/0050/ and
// https://en.cppreference.com/w/cpp/language/types
using AsterInt = std::ptrdiff_t;
using AsterReal = double;
}  // namespace aster

namespace MaterialLib
{
namespace Solids
{
namespace MFront
{
template <int DisplacementDim>
class MFront : public MechanicsBase<DisplacementDim>
{
public:
    struct MaterialStateVariables
        : public MechanicsBase<DisplacementDim>::MaterialStateVariables
    {
        void pushBackState() override {}

        MaterialStateVariables& operator=(MaterialStateVariables const&) =
            default;

        typename MechanicsBase<DisplacementDim>::MaterialStateVariables&
        operator=(typename MechanicsBase<DisplacementDim>::
                      MaterialStateVariables const& state) noexcept override
        {
            return operator=(static_cast<MaterialStateVariables const&>(state));
        }
    };

    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    using Callback = void (*)(
        aster::AsterReal* const STRESS, aster::AsterReal* const STATEV,
        aster::AsterReal* const DDSOE, const aster::AsterReal* const STRAN,
        const aster::AsterReal* const DSTRAN,
        const aster::AsterReal* const DTIME, const aster::AsterReal* const TEMP,
        const aster::AsterReal* const DTEMP,
        const aster::AsterReal* const PREDEF,
        const aster::AsterReal* const DPRED, const aster::AsterInt* const NTENS,
        const aster::AsterInt* const NSTATV,
        const aster::AsterReal* const PROPS,
        const aster::AsterInt* const NPROPS, const aster::AsterReal* const DROT,
        aster::AsterReal* const PNEWDT, const aster::AsterInt* const NUMMOD);

    explicit MFront(Callback callback);

    std::unique_ptr<
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables>
    createMaterialStateVariables() const override;

    boost::optional<std::tuple<KelvinVector,
                               std::unique_ptr<typename MechanicsBase<
                                   DisplacementDim>::MaterialStateVariables>,
                               KelvinMatrix>>
    integrateStress(
        double const t,
        ProcessLib::SpatialPosition const& x,
        double const dt,
        KelvinVector const& eps_prev,
        KelvinVector const& eps,
        KelvinVector const& sigma_prev,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables,
        double const T) const override;

    double computeFreeEnergyDensity(
        double const t,
        ProcessLib::SpatialPosition const& x,
        double const dt,
        KelvinVector const& eps,
        KelvinVector const& sigma,
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
            material_state_variables) const override;

private:
    Callback _callback;
};

extern template class MFront<2>;
extern template class MFront<3>;

}  // namespace MFront
}  // namespace Solids
}  // namespace MaterialLib
