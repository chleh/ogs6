/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <functional>
#include <map>
#include <memory>
#include <type_traits>
#include <typeindex>
#include <typeinfo>

#include "BaseLib/DUNEConfig.h"

#include <dune/geometry/type.hh>

#include "MeshLib/Elements/Elements.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/Fem/Integration/GaussLegendreIntegrationPolicy.h"

#ifndef OGS_MAX_ELEMENT_DIM
static_assert(false, "The macro OGS_MAX_ELEMENT_DIM is undefined.");
#endif

#ifndef OGS_MAX_ELEMENT_ORDER
static_assert(false, "The macro OGS_MAX_ELEMENT_ORDER is undefined.");
#endif

// The following macros decide which element types will be compiled, i.e.
// which element types will be available for use in simulations.

#ifdef OGS_ENABLE_ELEMENT_SIMPLEX
#define ENABLED_ELEMENT_TYPE_SIMPLEX 1u
#else
#define ENABLED_ELEMENT_TYPE_SIMPLEX 0u
#endif

#ifdef OGS_ENABLE_ELEMENT_CUBOID
#define ENABLED_ELEMENT_TYPE_CUBOID 1u << 1
#else
#define ENABLED_ELEMENT_TYPE_CUBOID 0u
#endif

#ifdef OGS_ENABLE_ELEMENT_PRISM
#define ENABLED_ELEMENT_TYPE_PRISM 1u << 2
#else
#define ENABLED_ELEMENT_TYPE_PRISM 0u
#endif

#ifdef OGS_ENABLE_ELEMENT_PYRAMID
#define ENABLED_ELEMENT_TYPE_PYRAMID 1u << 3
#else
#define ENABLED_ELEMENT_TYPE_PYRAMID 0u
#endif

// Dependent element types.
// Faces of tets, pyramids and prisms are triangles
#define ENABLED_ELEMENT_TYPE_TRI                                       \
    ((ENABLED_ELEMENT_TYPE_SIMPLEX) | (ENABLED_ELEMENT_TYPE_PYRAMID) | \
     (ENABLED_ELEMENT_TYPE_PRISM))
// Faces of hexes, pyramids and prisms are quads
#define ENABLED_ELEMENT_TYPE_QUAD                                     \
    ((ENABLED_ELEMENT_TYPE_CUBOID) | (ENABLED_ELEMENT_TYPE_PYRAMID) | \
     (ENABLED_ELEMENT_TYPE_PRISM))

// All enabled element types
#define OGS_ENABLED_ELEMENTS                                          \
    ((ENABLED_ELEMENT_TYPE_SIMPLEX) | (ENABLED_ELEMENT_TYPE_CUBOID) | \
     (ENABLED_ELEMENT_TYPE_PYRAMID) | (ENABLED_ELEMENT_TYPE_PRISM))

// Include only what is needed (Well, the conditions are not sharp).
#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_SIMPLEX) != 0
#include "NumLib/Fem/ShapeFunction/ShapeTet10.h"
#include "NumLib/Fem/ShapeFunction/ShapeTet4.h"
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_TRI) != 0
#include "NumLib/Fem/ShapeFunction/ShapeTri3.h"
#include "NumLib/Fem/ShapeFunction/ShapeTri6.h"
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_CUBOID) != 0
#include "NumLib/Fem/ShapeFunction/ShapeHex20.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex8.h"
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_QUAD) != 0
#include "NumLib/Fem/ShapeFunction/ShapeQuad4.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad8.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad9.h"
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_PRISM) != 0
#include "NumLib/Fem/ShapeFunction/ShapePrism15.h"
#include "NumLib/Fem/ShapeFunction/ShapePrism6.h"
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_PYRAMID) != 0
#include "NumLib/Fem/ShapeFunction/ShapePyra13.h"
#include "NumLib/Fem/ShapeFunction/ShapePyra5.h"
#endif

// TODO [DUNE] quadratic element support with DUNE
#undef OGS_MAX_ELEMENT_ORDER
#define OGS_MAX_ELEMENT_ORDER 1

namespace ProcessLib
{
/// The LocalDataInitializer is a functor creating a local assembler data with
/// corresponding to the mesh element type shape functions and calling
/// initialization of the new local assembler data.
/// For example for MeshLib::Quad a local assembler data with template argument
/// NumLib::ShapeQuad4 is created.
template <typename LocalAssemblerInterface,
          template <typename, typename, int, typename>
          class SmallDeformationLocalAssembler,
          int GlobalDim, typename Basis, typename... ConstructorArgs>
class LocalDataInitializer final
{
public:
    using LADataIntfPtr = std::unique_ptr<LocalAssemblerInterface>;

    explicit LocalDataInitializer(
        /*NumLib::LocalToGlobalIndexMap const& dof_table*/)
    // : _dof_table(dof_table)
    {
        Dune::GeometryType dgt;

        // /// Quads and Hexahedra ///////////////////////////////////

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_QUAD) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 2 && OGS_MAX_ELEMENT_ORDER >= 1
        dgt.makeCube(2);
        _builder[dgt] = makeLocalAssemblerBuilder<NumLib::ShapeQuad4>();
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_CUBOID) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 1
        dgt.makeCube(3);
        _builder[dgt] = makeLocalAssemblerBuilder<NumLib::ShapeHex8>();
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_QUAD) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 2 && OGS_MAX_ELEMENT_ORDER >= 2
        _builder[std::type_index(typeid(MeshLib::Quad8))] =
            makeLocalAssemblerBuilder<NumLib::ShapeQuad8>();
        _builder[std::type_index(typeid(MeshLib::Quad9))] =
            makeLocalAssemblerBuilder<NumLib::ShapeQuad9>();
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_CUBOID) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 2
        _builder[std::type_index(typeid(MeshLib::Hex20))] =
            makeLocalAssemblerBuilder<NumLib::ShapeHex20>();
#endif

        // /// Simplices ////////////////////////////////////////////////

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_TRI) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 2 && OGS_MAX_ELEMENT_ORDER >= 1
        dgt.makeSimplex(2);
        _builder[dgt] = makeLocalAssemblerBuilder<NumLib::ShapeTri3>();
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_SIMPLEX) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 1
        dgt.makeSimplex(3);
        _builder[dgt] = makeLocalAssemblerBuilder<NumLib::ShapeTet4>();
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_TRI) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 2 && OGS_MAX_ELEMENT_ORDER >= 2
        _builder[std::type_index(typeid(MeshLib::Tri6))] =
            makeLocalAssemblerBuilder<NumLib::ShapeTri6>();
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_SIMPLEX) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 2
        _builder[std::type_index(typeid(MeshLib::Tet10))] =
            makeLocalAssemblerBuilder<NumLib::ShapeTet10>();
#endif

        // /// Prisms ////////////////////////////////////////////////////

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_PRISM) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 1
        dgt.makePrism();
        _builder[dgt] = makeLocalAssemblerBuilder<NumLib::ShapePrism6>();
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_PRISM) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 2
        _builder[std::type_index(typeid(MeshLib::Prism15))] =
            makeLocalAssemblerBuilder<NumLib::ShapePrism15>();
#endif

        // /// Pyramids //////////////////////////////////////////////////

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_PYRAMID) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 1
        dgt.makePyramid();
        _builder[dgt] = makeLocalAssemblerBuilder<NumLib::ShapePyra5>();
#endif

#if (OGS_ENABLED_ELEMENTS & ENABLED_ELEMENT_TYPE_PYRAMID) != 0 && \
    OGS_MAX_ELEMENT_DIM >= 3 && OGS_MAX_ELEMENT_ORDER >= 2
        _builder[std::type_index(typeid(MeshLib::Pyramid13))] =
            makeLocalAssemblerBuilder<NumLib::ShapePyra13>();
#endif
    }

    /// Sets the provided \c data_ptr to the newly created local assembler data.
    ///
    /// \attention
    /// The index \c id is not necessarily the mesh item's id. Especially when
    /// having multiple meshes it will differ from the latter.
    LADataIntfPtr operator()(typename BaseLib::DUNEGridType<GlobalDim>::Traits::
                                 template Codim<0>::Entity const& mesh_item,
                             Basis const& basis,
                             ConstructorArgs&&... args) const
    {
        auto const geom_type = mesh_item.type();
        auto const it = _builder.find(geom_type);

        if (it != _builder.end())
        {
            // auto const num_local_dof = _dof_table.getNumberOfElementDOF(id);
            return it->second(mesh_item, basis,
                              std::forward<ConstructorArgs>(args)...);
        }
        else
        {
            OGS_FATAL(
                "You are trying to build a local assembler for an unknown mesh "
                "element type (dim %d, is none %d, topolody id %d)."
                " Maybe you have disabled this mesh element type in your build "
                "configuration or this process requires higher order elements.",
                geom_type.dim(), geom_type.isNone(), geom_type.id());
        }
    }

private:
    using LADataBuilder = std::function<LADataIntfPtr(
        typename BaseLib::DUNEGridType<GlobalDim>::Traits::template Codim<
            0>::Entity const& e,
        Basis const& basis,
        ConstructorArgs&&...)>;

    template <typename ShapeFunction>
    using IntegrationMethod = typename NumLib::GaussLegendreIntegrationPolicy<
        typename ShapeFunction::MeshElement>::IntegrationMethod;

    template <typename ShapeFunction>
    using LAData = SmallDeformationLocalAssembler<
        ShapeFunction, IntegrationMethod<ShapeFunction>, GlobalDim, Basis>;

    /// A helper forwarding to the correct version of makeLocalAssemblerBuilder
    /// depending whether the global dimension is less than the shape function's
    /// dimension or not.
    template <typename ShapeFunction>
    static LADataBuilder makeLocalAssemblerBuilder()
    {
        return makeLocalAssemblerBuilder<ShapeFunction>(
            static_cast<std::integral_constant<
                bool, (GlobalDim >= ShapeFunction::DIM)>*>(nullptr));
    }

    /// Mapping of element types to local assembler constructors.
    std::map<Dune::GeometryType, LADataBuilder> _builder;

    // NumLib::LocalToGlobalIndexMap const& _dof_table;

    // local assembler builder implementations.
private:
    /// Generates a function that creates a new LocalAssembler of type
    /// LAData<ShapeFunction>. Only functions with shape function's dimension
    /// less or equal to the global dimension are instantiated, e.g.  following
    /// combinations of shape functions and global dimensions: (Line2, 1),
    /// (Line2, 2), (Line2, 3), (Hex20, 3) but not (Hex20, 2) or (Hex20, 1).
    template <typename ShapeFunction>
    static LADataBuilder makeLocalAssemblerBuilder(std::true_type*)
    {
        return [](typename BaseLib::DUNEGridType<
                      GlobalDim>::Traits::template Codim<0>::Entity const& e,
                  Basis const& basis,
                  ConstructorArgs&&... args) {
            return LADataIntfPtr{new LAData<ShapeFunction>{
                e, basis, std::forward<ConstructorArgs>(args)...}};
        };
    }

    /// Returns nullptr for shape functions whose dimensions are less than the
    /// global dimension.
    template <typename ShapeFunction>
    static LADataBuilder makeLocalAssemblerBuilder(std::false_type*)
    {
        return nullptr;
    }
};

}  // namespace ProcessLib

#undef ENABLED_ELEMENT_TYPE_SIMPLEX
#undef ENABLED_ELEMENT_TYPE_CUBOID
#undef ENABLED_ELEMENT_TYPE_PYRAMID
#undef ENABLED_ELEMENT_TYPE_PRISM
#undef ENABLED_ELEMENT_TYPE_TRI
#undef ENABLED_ELEMENT_TYPE_QUAD
#undef OGS_ENABLED_ELEMENTS
