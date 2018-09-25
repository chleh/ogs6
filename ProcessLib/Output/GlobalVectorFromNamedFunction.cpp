/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GlobalVectorFromNamedFunction.h"
#include "MathLib/LinAlg/FinalizeVectorAssembly.h"
#include "MathLib/LinAlg/MatrixVectorTraits.h"
#include "NumLib/DOF/DOFTableUtil.h"

namespace ProcessLib
{
GlobalVectorFromNamedFunction::GlobalVectorFromNamedFunction(
    NumLib::SpecificFunctionCaller&& function_caller,
    MeshLib::FEMMesh const& mesh,
    NumLib::AbstractDOFTable const& dof_table_single,
    SecondaryVariableContext& context)
    : _function_caller(std::move(function_caller)),
      _mesh(mesh),
      _dof_table_single(dof_table_single),
      _context(context)
{
    assert(dof_table_single.getNumberOfComponents() == 1);
}

GlobalVector const& GlobalVectorFromNamedFunction::call(
    const double /*t*/,
    GlobalVector const& x,
    NumLib::AbstractDOFTable const& dof_table,
    std::unique_ptr<GlobalVector>& result)
{
    result = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(
        {_dof_table_single.dofSizeWithoutGhosts(),
         _dof_table_single.dofSizeWithoutGhosts(),
         /*&_dof_table_single.getGhostIndices()*/ nullptr, nullptr, nullptr,
         nullptr});

    auto* mesh = dynamic_cast<MeshLib::Mesh const*>(&_mesh);
    auto* dof_table_ =
        dynamic_cast<NumLib::LocalToGlobalIndexMap const*>(&dof_table);
    if (mesh != nullptr)
    {
        GlobalIndexType nnodes = mesh->getNumberOfNodes();

        auto const n_args = _function_caller.getNumberOfUnboundArguments();
        assert(static_cast<std::size_t>(dof_table_->getNumberOfComponents()) ==
               n_args);
        std::vector<double> args(n_args);

        for (GlobalIndexType node_id = 0; node_id < nnodes; ++node_id)
        {
            // TODO maybe fill args via callback mechanism or remove this class
            // entirely. Caution: The order of args will be the same as the
            // order of the components in the global vector!
            for (std::size_t i = 0; i < n_args; ++i)
            {
                args[i] =
                    NumLib::getNodalValue(x, *mesh, *dof_table_, node_id, i);
            }

            _context.index = node_id;
            auto const value = _function_caller.call(args);

            result->set(node_id, value);
        }
    }

    MathLib::finalizeVectorAssembly(*result);
    return *result;
}

}  // namespace ProcessLib
