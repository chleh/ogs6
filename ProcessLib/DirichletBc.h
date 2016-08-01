/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_DIRICHLETBC_H
#define PROCESS_LIB_DIRICHLETBC_H

#include <vector>

#include "NumLib/IndexValueVector.h"
#include "BoundaryCondition.h"

namespace ProcessLib
{

/// A dirichlet boundary condition is represented by a list of global indices
/// with corresponding values.
template <typename IndexType>
using DirichletBc = NumLib::IndexValueVector<IndexType>;

class DirichletBoundaryCondition : public BoundaryCondition
{
public:
    DirichletBoundaryCondition(NumLib::IndexValueVector<GlobalIndexType>&& bc)
        : _bc(std::move(bc))
    {
    }

    void apply(const double /*t*/,
               GlobalVector const& /*x*/,
               GlobalMatrix& /*K*/,
               GlobalVector& /*b*/) override
    {
        // Do nothing. Dirichlet BCs are handled specially.
    }

private:
    NumLib::IndexValueVector<GlobalIndexType> _bc;
};

}  // namespace ProcessLib

#endif  // PROCESS_LIB_DIRICHLETBC_H
