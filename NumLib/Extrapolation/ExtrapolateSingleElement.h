#pragma once

#include "ExtrapolatableElement.h"

namespace NumLib
{
//! Extrapolates integration point data of a single finite element to the
//! element's nodes.
//!
//! \param element The finite element.
//! \param integration_point_values_mat Matrix of integration point values,
//! dimension (\#components x \#integration_points).
//!
//! \return Matrix of nodal values, dimension (\#components x \#nodes).
Eigen::MatrixXd extrapolateSingleElement(
    ExtrapolatableElement const& element,
    Eigen::Map<Eigen::Matrix<double,
                             Eigen::Dynamic,
                             Eigen::Dynamic,
                             Eigen::RowMajor>> const&
        integration_point_values_mat);

}  // namespace NumLib
