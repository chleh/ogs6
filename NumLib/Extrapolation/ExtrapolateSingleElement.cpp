#include "ExtrapolateSingleElement.h"

#include <map>

#include <Eigen/SVD>

#include "BaseLib/Error.h"

struct CachedData
{
    //! The matrix A.
    Eigen::MatrixXd A;

    //! Moore-Penrose pseudo-inverse of A.
    Eigen::MatrixXd A_pinv;
};

/*! Maps (\#nodes, \#int_pts) to (N_0, QR decomposition),
 * where N_0 is the shape matrix of the first integration point.
 *
 * \note It is assumed that the pair (\#nodes, \#int_pts) uniquely
 * identifies the set of all shape matrices N for a mesh element (i.e., only
 * N, not dN/dx etc.).
 *
 * \todo Add the element dimension as identifying criterion, or change to
 * typeid.
 */
static std::map<std::pair<unsigned, unsigned>, CachedData> pseudo_inverse_cache;

namespace NumLib
{
Eigen::MatrixXd /* (#components x #nodes) */ extrapolateSingleElement(
    ExtrapolatableElement const& element,
    Eigen::Map<
        Eigen::Matrix<double,
                      Eigen::Dynamic,
                      Eigen::Dynamic,
                      Eigen::RowMajor>> const& /* (#components x #int_pts) */
        integration_point_values_mat)
{
    auto const& N_0 = element.getShapeMatrix(0);
    auto const num_nodes = static_cast<unsigned>(N_0.cols());
    // auto const num_components = integration_point_values_mat.rows();
    auto const num_int_pts = integration_point_values_mat.cols();

    // number of integration points in the element

    if (static_cast<unsigned>(num_int_pts) < num_nodes)
        OGS_FATAL(
            "Least squares is not possible if there are more nodes than"
            "integration points.");

    auto const pair_it_inserted = pseudo_inverse_cache.emplace(
        std::make_pair(num_nodes, num_int_pts), CachedData{});

    auto& cached_data = pair_it_inserted.first->second;
    if (pair_it_inserted.second)
    {
        DBUG("Computing new singular value decomposition");

        // interpolation_matrix * nodal_values = integration_point_values
        // We are going to pseudo-invert this relation now using singular value
        // decomposition.
        auto& interpolation_matrix = cached_data.A;
        interpolation_matrix.resize(num_int_pts, num_nodes);

        interpolation_matrix.row(0) = N_0;
        for (unsigned int_pt = 1; int_pt < num_int_pts; ++int_pt)
        {
            auto const& shp_mat = element.getShapeMatrix(int_pt);
            assert(shp_mat.cols() == num_nodes);

            // copy shape matrix to extrapolation matrix row-wise
            interpolation_matrix.row(int_pt) = shp_mat;
        }

        // JacobiSVD is extremely reliable, but fast only for small matrices.
        // But we usually have small matrices and we don't compute very often.
        // Cf.
        // http://eigen.tuxfamily.org/dox/group__TopicLinearAlgebraDecompositions.html
        //
        // Decomposes interpolation_matrix = U S V^T.
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(
            interpolation_matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);

        auto const& S = svd.singularValues();
        auto const& U = svd.matrixU();
        auto const& V = svd.matrixV();

        // Compute and save the pseudo inverse V * S^{-1} * U^T.
        auto const rank = svd.rank();
        assert(rank == num_nodes);

        // cf. http://eigen.tuxfamily.org/dox/JacobiSVD_8h_source.html
        cached_data.A_pinv.noalias() = U.leftCols(rank) *
                                       S.head(rank).asDiagonal().inverse() *
                                       V.leftCols(rank).transpose();
    }
    else if (cached_data.A.row(0) != N_0)
    {
        OGS_FATAL("The cached and the passed shapematrices differ.");
    }

    // Apply the pre-computed pseudo-inverse.
    Eigen::MatrixXd nodal_values =
        integration_point_values_mat * cached_data.A_pinv;

    return nodal_values;
}

}  // namespace NumLib
