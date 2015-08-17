/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_TES_FEM_H_
#define PROCESS_LIB_TES_FEM_H_

#include <memory>
#include <vector>

#include "TESFEM-notpl.h"
#include "TESProcess-notpl.h"

namespace ProcessLib
{

namespace TES
{

template <typename GlobalMatrix, typename GlobalVector>
class LocalAssemblerDataInterface
{
public:
    virtual ~LocalAssemblerDataInterface() = default;

    virtual void init(MeshLib::Element const& e,
                      std::size_t const local_matrix_size,
                      unsigned const integration_order,
                      TESProcessInterface const* process) = 0;

    virtual void assemble(std::vector<double> const& localX,
                          std::vector<double> const& localSecondaryVariables) = 0;

    virtual void addToGlobal(GlobalMatrix& A, GlobalVector& rhs,
                             AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const&) const = 0;
};



template <typename ShapeFunction_,
          typename IntegrationMethod_,
          typename GlobalMatrix,
          typename GlobalVector,
          unsigned GlobalDim>
class LocalAssemblerData
        : public LocalAssemblerDataInterface<GlobalMatrix, GlobalVector>
{
public:
    using ShapeFunction = ShapeFunction_;
    // using NodalMatrixType   = typename ShapeMatrixPolicyType<ShapeFunction, GlobalDim>::NodalMatrixType;
    // using NodalVectorType   = typename ShapeMatrixPolicyType<ShapeFunction, GlobalDim>::NodalVectorType;

    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices     = typename ShapeMatricesType::ShapeMatrices;


    void
    init(MeshLib::Element const& e,
         std::size_t const local_matrix_size,
         unsigned const integration_order,
         TESProcessInterface const* process) override;

    void assemble(std::vector<double> const& localX,
                  std::vector<double> const& localSecondaryVariables) override;

    void addToGlobal(GlobalMatrix& A, GlobalVector& rhs,
                     AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const& indices) const override;

private:
    std::vector<ShapeMatrices> _shape_matrices;
    LADataNoTpl _data;

    static const unsigned NODAL_DOF = 3;
    static const unsigned MAT_SIZE = ShapeFunction::NPOINTS * NODAL_DOF;
    using NodalMatrixType = Eigen::Matrix<double, MAT_SIZE, MAT_SIZE>;
    using NodalVectorType = Eigen::Matrix<double, MAT_SIZE, 1>;

    // std::unique_ptr<NodalMatrixType> _localA;
    // std::unique_ptr<NodalVectorType> _localRhs;
    Eigen::MatrixXd _localA;
    Eigen::VectorXd _localRhs;

    unsigned _integration_order = 2;
};


}   // namespace TES
}   // namespace ProcessLib


#include "TESFEM-impl.h"


#endif  // PROCESS_LIB_TES_FEM_H_
