#ifndef PROCESS_LIB_TESPROCESS_NOTPL_H_
#define PROCESS_LIB_TESPROCESS_NOTPL_H_


#include "MaterialsLib/adsorption/adsorption.h"
#include "Eigen/Sparse"

namespace ProcessLib
{

namespace TES
{


struct Materials
{
    Ads::Adsorption* _adsorption;
    const double _time_step = 5.0;
};


class TESProcessInterface
{
public:
    Materials const& getMaterials() const {
        return _materials;
    }

    virtual ~TESProcessInterface() = default;

protected:
    Materials _materials;
};


void printGlobalMatrix(const Eigen::SparseMatrix<double>& mat);

} // namespace TES

} // namespace ProcessLib

#endif // PROCESS_LIB_TESPROCESS_NOTPL_H_
