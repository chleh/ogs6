/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <boost/optional.hpp>
#include <functional>
#include <memory>
#include <tuple>
#include <vector>

#include "ProcessLib/Deformation/BMatrixPolicy.h"

#include "reflect-lib/reflect-macros.h"

namespace ProcessLib
{
class SpatialPosition;
}

namespace MeshLib
{
class Element;
}

namespace MaterialLib
{
namespace Solids
{
/// Interface for mechanical solid material models. Provides updates of the
/// stress for a given current state and also a tangent at that position. If the
/// implemented material model stores an internal state, the nested
/// MaterialStateVariables class should be used; it's only responsibility is to
/// provide state's push back possibility.
template <int DisplacementDim>
struct MechanicsBase
{
    /// The MaterialStateVariables may store material model specific state
    /// (other than sigma and eps), which are usually material history
    /// dependent. The objects are stored by the user (usually in assembly per
    /// integration point) and are created via \ref
    /// createMaterialStateVariables().
    struct MaterialStateVariables
    {
        virtual ~MaterialStateVariables() = default;
        virtual MaterialStateVariables& operator=(
            MaterialStateVariables const&) = default;

        virtual void pushBackState() = 0;

        REFLECT((MaterialStateVariables), FIELDS(), METHODS(pushBackState))
    };

    /// Polymorphic creator for MaterialStateVariables objects specific for a
    /// material model.
    virtual std::unique_ptr<MaterialStateVariables>
    createMaterialStateVariables() = 0;

    using KelvinVector = ProcessLib::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix = ProcessLib::KelvinMatrixType<DisplacementDim>;

    /// Dynamic size Kelvin vector and matrix wrapper for the polymorphic
    /// constitutive relation compute function.
    /// Returns nothing in case of errors in the computation if Newton
    /// iterations did not converge, for example.
    boost::optional<std::tuple<KelvinVector,
                               std::unique_ptr<MaterialStateVariables>,
                               KelvinMatrix>>
    integrateStress(double const t,
                    ProcessLib::SpatialPosition const& x,
                    double const dt,
                    Eigen::Matrix<double, Eigen::Dynamic, 1> const& eps_prev,
                    Eigen::Matrix<double, Eigen::Dynamic, 1> const& eps,
                    Eigen::Matrix<double, Eigen::Dynamic, 1> const& sigma_prev,
                    MaterialStateVariables const& material_state_variables)
    {
        // TODO Avoid copies of data:
        // Using MatrixBase<Derived> not possible because template functions
        // cannot be virtual. Maybe there is a workaround for this.  Using
        // Map<Matrix<double, ...>> makes the interface (for the material model
        // implementation) unnecessary difficult.
        KelvinVector const eps_prev_{eps_prev};
        KelvinVector const eps_{eps};
        KelvinVector const sigma_prev_{sigma_prev};

        return integrateStress(
            t, x, dt, eps_prev_, eps_, sigma_prev_, material_state_variables);
    }
    bool integrateStressPythonUniqueName(
        double const t,
        ProcessLib::SpatialPosition const& x,
        double const dt,
        Eigen::Matrix<double, Eigen::Dynamic, 1> const& eps_prev,
        Eigen::Matrix<double, Eigen::Dynamic, 1> const& eps,
        Eigen::Matrix<double, Eigen::Dynamic, 1> const& sigma_prev,
        MaterialStateVariables const& material_state_variables_in,
        Eigen::Ref<KelvinVector>& stress_out,
        MaterialStateVariables& material_state_variables_out,
        Eigen::Ref<KelvinMatrix>& tangent_modulus_out)
    {
#if 1
        auto res = integrateStress(
            t, x, dt, eps_prev, eps, sigma_prev, material_state_variables_in);
        if (!res)
            return false;

        stress_out = std::get<0>(*res);
        material_state_variables_out = *std::get<1>(*res);
        tangent_modulus_out = std::get<2>(*res);
        return true;
#else
        auto res = integrateStress(
            t, x, dt, eps_prev, eps, sigma_prev, material_state_variables);
        if (!res)
            return boost::none;
        return {{std::get<0>(*res), std::get<2>(*res)}};
#endif
    }

    /// Computation of the constitutive relation for specific material model.
    /// This should be implemented in the derived model. Fixed Kelvin vector and
    /// matrix size version; for dynamic size arguments there is an overloaded
    /// wrapper function.
    /// Returns nothing in case of errors in the computation if Newton
    /// iterations did not converge, for example.
    virtual boost::optional<std::tuple<KelvinVector,
                                       std::unique_ptr<MaterialStateVariables>,
                                       KelvinMatrix>>
    integrateStress(double const t,
                    ProcessLib::SpatialPosition const& x,
                    double const dt,
                    KelvinVector const& eps_prev,
                    KelvinVector const& eps,
                    KelvinVector const& sigma_prev,
                    MaterialStateVariables const& material_state_variables) = 0;

    /// Helper type for providing access to internal variables.
    struct InternalVariable
    {
        using Getter = std::function<std::vector<double> const&(
            MaterialStateVariables const&, std::vector<double>& /*cache*/)>;

        /// name of the internal variable
        std::string const name;

        /// number of components of the internal variable
        unsigned const num_components;

        /// function accessing the internal variable
        Getter const getter;
    };

    /// Returns internal variables defined by the specific material model, if
    /// any.
    virtual std::vector<InternalVariable> getInternalVariables() const
    {
        return {};
    }

    virtual ~MechanicsBase() = default;

    REFLECT((MechanicsBase<DisplacementDim>),
            FIELDS(),
            METHODS(createMaterialStateVariables,
                    integrateStressPythonUniqueName))
};

}  // namespace Solids
}  // namespace MaterialLib
