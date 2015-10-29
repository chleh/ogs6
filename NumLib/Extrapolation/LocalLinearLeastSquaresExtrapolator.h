#ifndef PROCESSLIB_LOCAL_LLSQ_EXTRAPOLATOR_H
#define PROCESSLIB_LOCAL_LLSQ_EXTRAPOLATOR_H

#include "Extrapolator.h"

namespace NumLib
{

template<typename GlobalVector, typename VariableEnum, typename LocalAssembler>
class LocalLinearLeastSquaresExtrapolator
        : public Extrapolator<GlobalVector, VariableEnum, LocalAssembler>
{
public:
    using LocalAssemblers = typename Extrapolator<GlobalVector, VariableEnum, LocalAssembler>
                                     ::LocalAssemblers;

    explicit LocalLinearLeastSquaresExtrapolator(
            AssemblerLib::LocalToGlobalIndexMap const& local_to_global)
        : _nodal_values(local_to_global.dofSize()),
          _local_to_global(local_to_global)
    {}


    void extrapolate(
            GlobalVector const& global_nodal_values,
            AssemblerLib::LocalToGlobalIndexMap const& global_nodal_values_map,
            LocalAssemblers const& loc_asms, VariableEnum var) override;

    void calculateResiduals(
            GlobalVector const& global_nodal_values,
            AssemblerLib::LocalToGlobalIndexMap const& global_nodal_values_map,
            LocalAssemblers const& loc_asms, VariableEnum var) override;

    GlobalVector const& getNodalValues() const override { return _nodal_values; }

    GlobalVector const& getElementResiduals() const override { return _residuals; }

private:

    void extrapolateElement(
            std::size_t index,
            GlobalVector const& global_nodal_values,
            AssemblerLib::LocalToGlobalIndexMap const& global_nodal_values_map,
            LocalAssembler const& loc_asm, VariableEnum var,
            GlobalVector& counts
            );

    double calculateResiudalElement(
            std::size_t index,
            GlobalVector const& global_nodal_values,
            AssemblerLib::LocalToGlobalIndexMap const& global_nodal_values_map,
            LocalAssembler const* loc_asm, VariableEnum var);

    GlobalVector _nodal_values;
    GlobalVector _residuals;
    AssemblerLib::LocalToGlobalIndexMap const& _local_to_global;
};

}

#include "LocalLinearLeastSquaresExtrapolator-impl.h"

#endif // PROCESSLIB_LOCAL_LLSQ_EXTRAPOLATOR_H
