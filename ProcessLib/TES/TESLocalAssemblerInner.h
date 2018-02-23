/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * The code of this file is used to decouple the evaluation of matrix elements
 * from the rest of OGS6,
 * not all of OGS6 has to be recompiled every time a small change is done.
 */

#pragma once

#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/LocalAssemblerTraits.h"

#include "TESLocalAssemblerData.h"

namespace ProcessLib
{
namespace TES
{
template <typename Traits>
class TESLocalAssemblerInner
{
public:
    explicit TESLocalAssemblerInner(AssemblyParams const& ap,
                                    const std::size_t element_id,
                                    const unsigned num_int_pts,
                                    const unsigned dimension);

    void assembleIntegrationPoint(
        unsigned integration_point,
        const double t,
        std::vector<double> const& localX,
        typename Traits::ShapeMatrices const& sm,
        const double weight,
        Eigen::Map<typename Traits::LocalMatrix>& local_M,
        Eigen::Map<typename Traits::LocalMatrix>& local_K,
        Eigen::Map<typename Traits::LocalVector>& local_b);

    void preTimestep();

    // TODO better encapsulation
    AssemblyParams const& getAssemblyParameters() const { return _d.ap; }
    TESLocalAssemblerData const& getData() const { return _d; }
    TESLocalAssemblerData& getData() { return _d; }

private:
    Eigen::Matrix3d getMassCoeffMatrix(const unsigned int_pt);
    typename Traits::LaplaceMatrix getLaplaceCoeffMatrix(const double t,
                                                         const unsigned int_pt,
                                                         const unsigned dim);
    Eigen::Matrix3d getAdvectionCoeffMatrix(const unsigned int_pt);
    Eigen::Matrix3d getContentCoeffMatrix(const unsigned int_pt);
    Eigen::Vector3d getRHSCoeffVector(const unsigned int_pt);

    void preEachAssembleIntegrationPoint(
        const unsigned int_pt,
        const double t,
        std::vector<double> const& localX,
        typename Traits::ShapeMatrices const& sm);

    TESLocalAssemblerData _d;
};

}  // namespace TES

}  // namespace ProcessLib

#include "TESLocalAssemblerInner-impl.h"
