/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MFront.h"

namespace MaterialLib
{
namespace Solids
{
namespace MFront
{
template <int DisplacementDim>
MFront<DisplacementDim>::MFront(Callback callback) : _callback(callback)
{
}

template <int DisplacementDim>
std::unique_ptr<typename MechanicsBase<DisplacementDim>::MaterialStateVariables>
MFront<DisplacementDim>::createMaterialStateVariables() const
{
    return std::make_unique<MaterialStateVariables>();
}

template <int DisplacementDim>
boost::optional<std::tuple<typename MFront<DisplacementDim>::KelvinVector,
                           std::unique_ptr<typename MechanicsBase<
                               DisplacementDim>::MaterialStateVariables>,
                           typename MFront<DisplacementDim>::KelvinMatrix>>
MFront<DisplacementDim>::integrateStress(
    double const t,
    ProcessLib::SpatialPosition const& x,
    double const dt,
    KelvinVector const& eps_prev,
    KelvinVector const& eps,
    KelvinVector const& sigma_prev,
    typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
        material_state_variables,
    double const T) const
{
    DBUG("Integrate stress.");
    _callback(nullptr /* aster::AsterReal* const STRESS */,
              nullptr /* aster::AsterReal* const STATEV */,
              nullptr /* aster::AsterReal* const DDSOE */,
              nullptr /* const aster::AsterReal* const STRAN */,
              nullptr /* const aster::AsterReal* const DSTRAN */,
              nullptr /* const aster::AsterReal* const DTIME */,
              nullptr /* const aster::AsterReal* const TEMP */,
              nullptr /* const aster::AsterReal* const DTEMP */,
              nullptr /* const aster::AsterReal* const PREDEF */,
              nullptr /* const aster::AsterReal* const DPRED */,
              nullptr /* const aster::AsterInt* const NTENS */,
              nullptr /* const aster::AsterInt* const NSTATV */,
              nullptr /* const aster::AsterReal* const PROPS */,
              nullptr /* const aster::AsterInt* const NPROPS */,
              nullptr /* const aster::AsterReal* const DROT */,
              nullptr /* aster::AsterReal* const PNEWDT */,
              nullptr /* const aster::AsterInt* const NUMMOD */);
    return boost::none;
}

template <int DisplacementDim>
double MFront<DisplacementDim>::computeFreeEnergyDensity(
    double const t,
    ProcessLib::SpatialPosition const& x,
    double const dt,
    KelvinVector const& eps,
    KelvinVector const& sigma,
    typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
        material_state_variables) const
{
    return std::numeric_limits<double>::quiet_NaN();
}

template class MFront<2>;
template class MFront<3>;

}  // namespace MFront
}  // namespace Solids
}  // namespace MaterialLib
