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
    using MaterialStateVariables =
        typename MechanicsBase<DisplacementDim>::MaterialStateVariables;

    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    std::unique_ptr<MaterialStateVariables> createMaterialStateVariables()
        const override
    {
        return nullptr;
    }

    boost::optional<std::tuple<KelvinVector,
                               std::unique_ptr<MaterialStateVariables>,
                               KelvinMatrix>>
    integrateStress(double const t,
                    ProcessLib::SpatialPosition const& x,
                    double const dt,
                    KelvinVector const& eps_prev,
                    KelvinVector const& eps,
                    KelvinVector const& sigma_prev,
                    MaterialStateVariables const& material_state_variables,
                    double const T) const override
    {
        return boost::none;
    }

    double computeFreeEnergyDensity(
        double const t,
        ProcessLib::SpatialPosition const& x,
        double const dt,
        KelvinVector const& eps,
        KelvinVector const& sigma,
        MaterialStateVariables const& material_state_variables) const override
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
};
}  // namespace MFront
}  // namespace Solids
}  // namespace MaterialLib
