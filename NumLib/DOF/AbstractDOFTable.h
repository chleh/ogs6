#pragma once

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "MathLib/LinAlg/RowColumnIndices.h"
#include "NumLib/NumericsConfig.h"

namespace NumLib
{
class AbstractDOFTable
{
public:
    using RowColumnIndices = MathLib::RowColumnIndices<GlobalIndexType>;

    virtual std::size_t dofSizeWithoutGhosts() const = 0;
    virtual std::size_t size() const = 0;
    virtual int getNumberOfVariables() const = 0;
    virtual int getNumberOfVariableComponents(int variable_id) const = 0;
    virtual int getNumberOfComponents() const = 0;
    virtual RowColumnIndices operator()(std::size_t const mesh_item_id,
                                        const int component_id) const = 0;
    virtual std::vector<GlobalIndexType> const& getGhostIndices() const = 0;
    virtual ~AbstractDOFTable() = default;
};

}  // namespace NumLib
