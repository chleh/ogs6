#ifndef PROCESSLIB_EXTRAPOLATOR_H
#define PROCESSLIB_EXTRAPOLATOR_H

#include <vector>

#include "AssemblerLib/LocalToGlobalIndexMap.h"

namespace NumLib
{


/* // TODO: implement
// see http://eigen.tuxfamily.org/dox-devel/group__LeastSquares.html
enum class LinearLeastSquaresBy { SVD, QR, NormalEquation };
const LinearLeastSquaresBy linear_least_squares = LinearLeastSquaresBy::NormalEquation;
*/


template<typename GlobalVector, typename VariableEnum, typename LocalAssembler>
class Extrapolator
{
public:
    using LocalAssemblers = std::vector<LocalAssembler*>;

    virtual void extrapolate(
            GlobalVector const& global_nodal_values,
            AssemblerLib::LocalToGlobalIndexMap const& global_nodal_values_map,
            LocalAssemblers const& loc_asms, VariableEnum var) = 0;

    virtual void calculateResiduals(
            GlobalVector const& global_nodal_values,
            AssemblerLib::LocalToGlobalIndexMap const& global_nodal_values_map,
            LocalAssemblers const& loc_asms, VariableEnum var) = 0;

    virtual GlobalVector const& getNodalValues() const = 0;
    virtual GlobalVector const& getElementResiduals() const = 0;
};

} // namespace ProcessLib

#endif // PROCESSLIB_EXTRAPOLATOR_H
