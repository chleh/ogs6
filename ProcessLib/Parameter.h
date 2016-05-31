/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_PARAMETER_H_
#define PROCESS_LIB_PARAMETER_H_

#include <memory>

#include <logog/include/logog.hpp>
#include <boost/optional.hpp>

#include "BaseLib/ConfigTree.h"
#include "MeshLib/Elements/Element.h"

namespace MeshLib
{
template <typename T>
class PropertyVector;
}

namespace ProcessLib
{
/// Base class for all parameters, not an interface class. This avoids using of
/// void* when storing parameters and convenient destruction.
/// Its property name helps addressing the right parameter.
class ParameterBase
{
public:
    explicit ParameterBase(std::string const& name)
        : _name(name)
    {}

    std::string const& getName() const { return _name; }

    virtual ~ParameterBase() = default;

private:
    std::string const _name;
};

/// A parameter is representing a value or function of any type.
/// The ReturnType can represent an underlying type of an aggregate type like
/// tuple or matrix (\see tupleSize()).
/// The total number of stored tuples is provided.
template <typename ReturnType, typename... Args>
class Parameter : public ParameterBase
{
public:
    explicit Parameter(std::string const& name)
        : ParameterBase(name)
    {}

    virtual ReturnType operator()(Args&&... args) const = 0;
};

/// Single, constant value parameter.
template <typename ReturnType>
class ConstParameter final
    : public Parameter<ReturnType, MeshLib::Element const&>
{
public:
    ConstParameter(std::string const& name, ReturnType value)
        : Parameter<ReturnType, MeshLib::Element const&>(name)
        , _value(value)
    {
    }

    ReturnType operator()(MeshLib::Element const&) const override
    {
        return _value;
    }

private:
    ReturnType _value;
};

std::unique_ptr<ParameterBase> createConstParameter(
        std::string const& name, BaseLib::ConfigTree const& config);

/// A parameter represented by a mesh property vector.
template <typename ReturnType>
class MeshPropertyParameter final
    : public Parameter<ReturnType, MeshLib::Element const&>
{
public:
    MeshPropertyParameter(std::string const& name,
                          MeshLib::PropertyVector<ReturnType> const& property)
        : Parameter<ReturnType, MeshLib::Element const&>(name)
        , _property(property)
    {
    }

    ReturnType operator()(MeshLib::Element const& e) const override
    {
        return _property[e.getID()];
    }

private:
    MeshLib::PropertyVector<ReturnType> const& _property;
};

std::unique_ptr<ParameterBase> createMeshPropertyParameter(
        std::string const& name,
        BaseLib::ConfigTree const& config, MeshLib::Mesh const& mesh);

}  // namespace ProcessLib

#endif  // PROCESS_LIB_PARAMETER_H_
