/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * The code of this file is used to decouple the evaluation of matrix elements
 * from the rest of OGS6,
 * not all of OGS6 has to be recompiled every time a small change is done.
 */

#ifndef PROCESS_LIB_TES_FEM_NOTPL_H_
#define PROCESS_LIB_TES_FEM_NOTPL_H_

#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/VariableTransformation.h"

#include "TESLocalAssemblerData.h"

namespace ProcessLib
{
namespace TES
{
enum class TESIntPtVariables : unsigned
{
    SOLID_DENSITY,
    REACTION_RATE,
    VELOCITY_X,
    VELOCITY_Y,
    VELOCITY_Z,
    LOADING,
    REACTION_DAMPING_FACTOR
};

template <typename Traits>
class TESLocalAssemblerInner
{
public:
    explicit TESLocalAssemblerInner(AssemblyParams const& ap,
                                    const unsigned num_int_pts,
                                    const unsigned dimension);

    void assembleIntegrationPoint(
        unsigned integration_point,
        std::vector<double> const& localX,
        typename Traits::ShapeMatrices::ShapeType const& smN,
        typename Traits::ShapeMatrices::DxShapeType const& smDNdx,
        typename Traits::ShapeMatrices::JacobianType const& smJ,
        const double smDetJ,
        const double weight,
        typename Traits::LocalMatrix& local_M,
        typename Traits::LocalMatrix& local_K,
        typename Traits::LocalVector& local_b);

    void preEachAssemble();

    std::vector<double> const& getIntegrationPointValues(
        TESIntPtVariables const var, std::vector<double>& cache) const;

    // TODO better encapsulation
    AssemblyParams const& getAssemblyParameters() const { return _d.ap; }
    TESFEMReactionAdaptor const& getReactionAdaptor() const
    {
        return *_d.reaction_adaptor;
    }
    TESFEMReactionAdaptor& getReactionAdaptor() {
        return *_d.reaction_adaptor;
    }
    TESLocalAssemblerData const& getData() const {
        return _d;
    }
    TESLocalAssemblerData& getData() {
        return _d;
    }

private:
    Eigen::Matrix3d getMassCoeffMatrix(const unsigned int_pt);
    typename Traits::LaplaceMatrix getLaplaceCoeffMatrix(const unsigned int_pt,
                                                         const unsigned dim);
    Eigen::Matrix3d getAdvectionCoeffMatrix(const unsigned int_pt);
    Eigen::Matrix3d getContentCoeffMatrix(const unsigned int_pt);
    Eigen::Vector3d getRHSCoeffVector(const unsigned int_pt);

    void preEachAssembleIntegrationPoint(
        const unsigned int_pt,
        std::vector<double> const& localX,
        typename Traits::ShapeMatrices::ShapeType const& smN,
        typename Traits::ShapeMatrices::DxShapeType const& smDNdx,
        typename Traits::ShapeMatrices::JacobianType const& smJ,
        const double smDetJ);

    void initReaction(const unsigned int_pt);

    TESLocalAssemblerData _d;
};

}  // namespace TES

}  // namespace ProcessLib

// tricking cmake dependency checker
#include "TESLocalAssemblerInner-impl-incl.h"

#endif  // PROCESS_LIB_TES_FEM_NOTPL_H_
