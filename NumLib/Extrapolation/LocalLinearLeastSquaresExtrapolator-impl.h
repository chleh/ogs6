#include "logog/include/logog.hpp"

#include "Eigen/Core"

#include "NumLib/Function/Interpolation.h"

#include "LocalLinearLeastSquaresExtrapolator.h"

#include "LocalNodalDOF-impl.h"


namespace NumLib
{

template<typename GlobalVector, typename VariableEnum, typename LocalAssembler>
void
LocalLinearLeastSquaresExtrapolator<GlobalVector, VariableEnum, LocalAssembler>::
extrapolate(
        GlobalVector const& global_nodal_values,
        AssemblerLib::LocalToGlobalIndexMap const& global_nodal_values_map,
        LocalAssemblers const& loc_asms, VariableEnum var)
{
    _nodal_values = 0.0;

    GlobalVector counts(_nodal_values.size());

    for (std::size_t i=0; i<loc_asms.size(); ++i)
    {
        extrapolateElement(i, global_nodal_values, global_nodal_values_map, loc_asms[i], var, counts);
    }

    _nodal_values.componentwiseDivide(counts);
}

template<typename GlobalVector, typename VariableEnum, typename LocalAssembler>
void
LocalLinearLeastSquaresExtrapolator<GlobalVector, VariableEnum, LocalAssembler>::
calculateResiduals(
        GlobalVector const& global_nodal_values,
        AssemblerLib::LocalToGlobalIndexMap const& global_nodal_values_map,
        LocalAssemblers const& loc_asms, VariableEnum var)
{
    _residuals.resize(loc_asms.size());

    for (std::size_t i=0; i<loc_asms.size(); ++i)
    {
        _residuals[i] = calculateResiudalElement(
                            i, global_nodal_values, global_nodal_values_map, loc_asms[i], var);
    }
}

template<typename GlobalVector, typename VariableEnum, typename LocalAssembler>
void
LocalLinearLeastSquaresExtrapolator<GlobalVector, VariableEnum, LocalAssembler>::
extrapolateElement(
        std::size_t index,
        GlobalVector const& global_nodal_values,
        AssemblerLib::LocalToGlobalIndexMap const& global_nodal_values_map,
        LocalAssembler const* loc_asm, VariableEnum var,
        GlobalVector& counts
        )
{
    NumLib::LocalNodalDOFImpl<GlobalVector> nodal_dof{index, global_nodal_values, global_nodal_values_map};

    auto gp_vals = loc_asm->getIntegrationPointValues(var, nodal_dof);

    const unsigned nn = loc_asm->getShapeMatrix(0).rows(); // number of mesh nodes
    // const unsigned nn = 5;
    const unsigned ni = gp_vals->size();        // number of gauss points


    Eigen::MatrixXd N(ni, nn);

    for (unsigned gp=0; gp<ni; ++gp)
    {
        auto const& shp_mat = loc_asm->getShapeMatrix(gp);
        assert(shp_mat.rows() == nn);

        for (unsigned n=0; n<nn; ++n)
        {
            N(gp, n) = shp_mat(n);
        }
    }

    const Eigen::Map<const Eigen::VectorXd> gpvs(gp_vals->data(), gp_vals->size());

    {
        // TODO: now always zeroth component is used
        auto const& indices = _local_to_global(index, 0).rows;

        Eigen::VectorXd elem_nodal_vals = (N.transpose() * N).ldlt().solve(N.transpose() * gpvs);
        _nodal_values.add(indices, elem_nodal_vals);
        counts.add(indices, 1.0);
    }
}

template<typename GlobalVector, typename VariableEnum, typename LocalAssembler>
double
LocalLinearLeastSquaresExtrapolator<GlobalVector, VariableEnum, LocalAssembler>::
calculateResiudalElement(
        std::size_t index,
        GlobalVector const& global_nodal_values,
        AssemblerLib::LocalToGlobalIndexMap const& global_nodal_values_map,
        LocalAssembler const* loc_asm, VariableEnum var)
{
    NumLib::LocalNodalDOFImpl<GlobalVector> nodal_dof{index, global_nodal_values, global_nodal_values_map};

    auto gp_vals = loc_asm->getIntegrationPointValues(var, nodal_dof);
    const unsigned ni = gp_vals->size();        // number of gauss points

    // TODO: now always zeroth component is used
    const auto& indices = _local_to_global(index, 0).rows;

    // filter nodal values of the current element
    std::vector<double> nodal_vals_element;
    nodal_vals_element.resize(indices.size());
    for (unsigned i=0; i<indices.size(); ++i) {
        nodal_vals_element[i] = _nodal_values[indices[i]];
    }

    double residual = 0.0;
    double gp_val_extrapol = 0.0;
    std::array<double*, 1> gp_val_extrapol2 = { &gp_val_extrapol };

    for (unsigned gp=0; gp<ni; ++gp)
    {
        NumLib::shapeFunctionInterpolate(
                    nodal_vals_element, loc_asm->getShapeMatrix(gp),
                    gp_val_extrapol2);
        auto const& ax_m_b = gp_val_extrapol - (*gp_vals)[gp];
        residual += ax_m_b * ax_m_b;
    }

    return residual / ni;
}

} // namespace ProcessLib
