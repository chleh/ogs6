/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <vector>

#include <logog/include/logog.hpp>

#include "BaseLib/DUNEConfig.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"

#include "LocalDataInitializer.h"

namespace ProcessLib
{
namespace SmallDeformationDUNE
{
namespace detail
{
template <int GlobalDim, typename Basis,
          template <typename, typename, int, typename>
          class LocalAssemblerImplementation,
          typename LocalAssemblerInterface, typename... ExtraCtorArgs>
void createLocalAssemblers(
    Basis const& basis,
    // NumLib::LocalToGlobalIndexMap const& dof_table,
    // std::vector<MeshLib::Element*> const& mesh_elements,
    std::vector<std::unique_ptr<LocalAssemblerInterface>>& local_assemblers,
    ExtraCtorArgs&&... extra_ctor_args)
{
    // Shape matrices initializer
    using LDI = LocalDataInitializer<LocalAssemblerInterface,
                                     LocalAssemblerImplementation, GlobalDim,
                                     Basis, ExtraCtorArgs...>;

    // Populate the vector of local assemblers.
    auto const& gridView = basis.gridView();
    auto const num_elements = gridView.size(0 /* codim */);
    DBUG("num elements: %i.", num_elements);
    local_assemblers.clear();
    local_assemblers.reserve(num_elements);

    LDI initializer;

    for (auto const& element : Dune::elements(gridView))
    {
        // DBUG("type: %s", Dune::className(element).c_str());
        local_assemblers.push_back(initializer(
            element, basis, std::forward<ExtraCtorArgs>(extra_ctor_args)...));
    }
}

}  // namespace detail

/*! Creates local assemblers for each element of the given \c mesh.
 *
 * \tparam LocalAssemblerImplementation the individual local assembler type
 * \tparam LocalAssemblerInterface the general local assembler interface
 * \tparam ExtraCtorArgs types of additional constructor arguments.
 *         Those arguments will be passed to the constructor of
 *         \c LocalAssemblerImplementation.
 *
 * The first two template parameters cannot be deduced from the arguments.
 * Therefore they always have to be provided manually.
 */
template <int GlobalDim, typename Basis,
          template <typename, typename, int, typename>
          class LocalAssemblerImplementation,
          typename LocalAssemblerInterface, typename... ExtraCtorArgs>
void createLocalAssemblers(
    Basis const& basis,
    // NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<std::unique_ptr<LocalAssemblerInterface>>& local_assemblers,
    ExtraCtorArgs&&... extra_ctor_args)
{
    DBUG("Create local assemblers.");

    detail::createLocalAssemblers<GlobalDim, Basis,
                                  LocalAssemblerImplementation>(
        basis, local_assemblers,
        std::forward<ExtraCtorArgs>(extra_ctor_args)...);
}
}  // namespace SmallDeformationDUNE

}  // namespace ProcessLib
