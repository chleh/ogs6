/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "PythonBoundaryCondition.h"

#include "NumLib/DOF/DOFTableUtil.h"

#include "GenericNaturalBoundaryConditionLocalAssembler.h"
#include "PythonBoundaryConditionDetail.h"

namespace ProcessLib
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class PythonConditionLocalAssembler final
    : public GenericNaturalBoundaryConditionLocalAssembler<
          ShapeFunction, IntegrationMethod, GlobalDim>
{
    using Base = GenericNaturalBoundaryConditionLocalAssembler<
        ShapeFunction, IntegrationMethod, GlobalDim>;

public:
    /// The neumann_bc_value factor is directly integrated into the local
    /// element matrix.
    PythonConditionLocalAssembler(MeshLib::Element const& e,
                                  std::size_t const local_matrix_size,
                                  bool is_axially_symmetric,
                                  unsigned const integration_order,
                                  PythonBoundaryConditionData const& data)
        : Base(e, is_axially_symmetric, integration_order),
          _data(data),
          _local_rhs(local_matrix_size),
          _element(e)
    {
    }

    void assemble(std::size_t const id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, const GlobalVector& /*x*/,
                  GlobalMatrix& /*K*/, GlobalVector& b) override
    {
        using ShapeMatricesType =
            ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
        using FemType =
            NumLib::TemplateIsoparametric<ShapeFunction, ShapeMatricesType>;
        FemType fe(*static_cast<const typename ShapeFunction::MeshElement*>(
            &_element));

        auto* bc =
            _data.scope[_data.bc_object.c_str()].cast<::PyBoundaryCondition*>();

        _local_rhs.setZero();

        unsigned const n_integration_points =
            Base::_integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& sm = Base::_shape_matrices[ip];
            auto const coords = fe.interpolateCoordinates(sm.N);
            auto const res = bc->getFlux(t, coords);
            if (!res.first)
                return;
            if (!bc->isOverriddenNatural())
                throw PyNotOverridden{};

            auto const& wp = Base::_integration_method.getWeightedPoint(ip);
            _local_rhs.noalias() += sm.N * res.second * sm.detJ *
                                    wp.getWeight() * sm.integralMeasure;
        }

        auto const indices = NumLib::getIndices(id, dof_table_boundary);
        b.add(indices, _local_rhs);
    }

private:
    PythonBoundaryConditionData const& _data;
    typename Base::NodalVectorType _local_rhs;
    MeshLib::Element const& _element;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace ProcessLib
