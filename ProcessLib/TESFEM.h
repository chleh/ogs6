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

#include "NumLib/Extrapolation/LocalNodalDOF.h"

namespace ProcessLib
{

namespace TES
{

class Extrapolatable
{
public:
    virtual Eigen::VectorXd const& getShapeMatrix(const unsigned integration_point) const = 0;

    virtual std::shared_ptr< const std::vector<double> >
    getIntegrationPointValues(SecondaryVariables var, NumLib::LocalNodalDOF& nodal_dof) const = 0;
};


template <typename GlobalMatrix, typename GlobalVector>
class LocalAssemblerDataInterface
        : public Extrapolatable
{
public:
    virtual ~LocalAssemblerDataInterface() = default;

    virtual void init(MeshLib::Element const& e,
                      std::size_t const local_matrix_size,
                      unsigned const integration_order,
                      TESProcessInterface* process) = 0;

    virtual void assemble(std::vector<double> const& localX,
                          std::vector<double> const& localXPrevTs) = 0;

    virtual void addToGlobal(GlobalMatrix& A, GlobalVector& rhs,
                             AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const&) const = 0;

    virtual bool checkBounds(std::vector<double> const& localX,
                             std::vector<double> const& localX_pts) = 0;
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
         TESProcessInterface* process) override;

    void assemble(std::vector<double> const& localX,
                  std::vector<double> const& localXPrevTs) override;

    void addToGlobal(GlobalMatrix& A, GlobalVector& rhs,
                     AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const& indices) const override;

    Eigen::VectorXd const& getShapeMatrix(const unsigned integration_point) const override {
        return _shape_matrices[integration_point].N;
        // const auto& shp_mats = _shape_matrices[integration_point];
        // return (_shape_matrices.data() + integration_point)->N;
        // return shp_mats.N;
    }

    bool checkBounds(std::vector<double> const& localX,
                     std::vector<double> const& localX_pts);

    std::shared_ptr<const std::vector<double> >
    getIntegrationPointValues(SecondaryVariables var, NumLib::LocalNodalDOF& nodal_dof) const override;

private:
    std::vector<ShapeMatrices> _shape_matrices;
    LADataNoTpl _data;

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
