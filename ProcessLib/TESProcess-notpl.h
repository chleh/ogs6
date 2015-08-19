#ifndef PROCESS_LIB_TESPROCESS_NOTPL_H_
#define PROCESS_LIB_TESPROCESS_NOTPL_H_


#include "MaterialsLib/adsorption/adsorption.h"
#include "Eigen/Sparse"
#include "Eigen/Eigen"

namespace ProcessLib
{

namespace TES
{


struct Materials
{
    Ads::Adsorption* _adsorption;
    double _time_step = 0.1;
    double _initial_solid_density = 1382.36248218;
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


bool calculateError(Eigen::VectorXd* current_solution,
                    const Eigen::Ref<Eigen::VectorXd>& previous_solution, Materials* materials);
// bool calculateError(const Eigen::SparseMatrix<double>& current_solution,
//                     const Eigen::SparseMatrix<double>& previous_solution);


void printGlobalMatrix(const Eigen::SparseMatrix<double>& mat);
void printGlobalVector(const Eigen::Ref<Eigen::VectorXd>& vec);

} // namespace TES

} // namespace ProcessLib

#endif // PROCESS_LIB_TESPROCESS_NOTPL_H_
