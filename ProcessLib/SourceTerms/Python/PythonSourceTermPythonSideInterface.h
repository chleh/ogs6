/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

namespace ProcessLib
{
//! Base class for source terms.
//! This class will get Python bindings and is intended to be to be derived in
//! Python.
class PythonSourceTermPythonSideInterface
{
public:
    /*!
     * Computes source term values for the provided arguments
     * (time, position of the node, node id, primary variables at the node).
     *
     * \return value value of the source term at that node.
     */
    virtual double getSourceTermValue(
        double /*t*/, std::array<double, 3> /*x*/, std::size_t /*node_id*/,
        std::vector<double> const& /*primary_variables*/) const
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    //! Tells if getSourceTermValue() has been overridden in the derived class
    //! in Python.
    //!
    //! \pre getSourceTermValue() must already have been called once.
    bool isOverriddenSourceTermValue() const
    {
        return _overridden_source_term_value;
    }

    /*!
     * Computes the flux for the provided arguments (time, position of the node,
     * node id, primary variables at the node).
     *
     * \return flux Flux of the source term at that node and derivative of the
     * flux w.r.t. all primary variables.
     */
    virtual std::pair<double, std::array<double, 3>> getFlux(
        double /*t*/, std::array<double, 3> /*x*/, std::size_t /*node_id*/,
        std::vector<double> const& /*primary_variables*/) const
    {
        return {std::numeric_limits<double>::quiet_NaN(), {}};
    }

    //! Tells if getFlux() has been overridden in the derived class in Python.
    //!
    //! \pre getFlux() must already have been called once.
    bool isOverriddenGetFlux() const
    {
        return _overridden_get_flux;
    }

    virtual ~PythonSourceTermPythonSideInterface() = default;

private:
    //! Tells if getSourceTermValue() has been overridden in the derived class
    //! in Python.
    mutable bool _overridden_source_term_value = true;
    //! Tells if getFlux() has been overridden in the derived class in Python.
    mutable bool _overridden_get_flux = true;
};
}  // namespace ProcessLib
