#pragma once

#include <vector>

#include "AssemblerLib/LocalToGlobalIndexMap.h"

#include "LocalNodalDOF.h"

namespace NumLib
{

template <typename GlobalVector>
class LocalNodalDOFImpl : public NumLib::LocalNodalDOF
{
public:
    LocalNodalDOFImpl(
            std::size_t index,
            GlobalVector const& global_nodal_values,
            AssemblerLib::LocalToGlobalIndexMap const& index_map)
        : _index{index}, _index_map{index_map}, _global_nodal_values{global_nodal_values}
    {}

    // this function is very similar to my [CL] current implementation of VectorMatrixAssembler
    std::vector<double> getElementNodalValues() override
    {
        std::vector<double> localX;
        auto const num_comp = _index_map.getNumComponents();

        // The local matrix will always be ordered by component,
        // no matter what the order of the global matrix is.
        for (unsigned c=0; c<num_comp; ++c)
        {
            auto const idcs = _index_map(_index, c).rows;
            localX.reserve(localX.size() + idcs.size());
            for (auto ip : idcs)
            {
                localX.emplace_back(_global_nodal_values[ip]);
            }
        }

        return localX;
    }

    std::vector<double> getElementNodalValues(unsigned component) override
    {
        std::vector<double> localX;

        auto const idcs = _index_map(_index, component).rows;
        localX.reserve(idcs.size());
        for (auto ip : idcs)
        {
            localX.emplace_back(_global_nodal_values[ip]);
        }

        return localX;
    }

private:
    const std::size_t _index;
    AssemblerLib::LocalToGlobalIndexMap const& _index_map;
    GlobalVector const& _global_nodal_values;
};

}
