#ifndef PROCESS_LIB_TESPROCESS_NOTPL_H_
#define PROCESS_LIB_TESPROCESS_NOTPL_H_


#include "MaterialsLib/adsorption/adsorption.h"
#include "Eigen/Sparse"
#include "Eigen/Eigen"

namespace ProcessLib
{

namespace TES
{


struct ProcessParams
{
    Ads::Adsorption* _adsorption;
    const double   _time_step = 0.5;
    const double   _time_max  = 2000;
    const unsigned _output_every_nth_step = 50;
    const double   _initial_solid_density = 1382.36248218;
    bool  _is_new_timestep = true;
};


class TESProcessInterface
{
public:
    ProcessParams& getParams() {
        return _parameters;
    }
    ProcessParams const& getParams() const {
        return _parameters;
    }

    virtual ~TESProcessInterface() = default;

protected:
    ProcessParams _parameters;
};


bool calculateError(Eigen::VectorXd* current_solution,
                    const Eigen::Ref<Eigen::VectorXd>& previous_solution, ProcessParams* materials);
// bool calculateError(const Eigen::SparseMatrix<double>& current_solution,
//                     const Eigen::SparseMatrix<double>& previous_solution);


void printGlobalMatrix(const Eigen::SparseMatrix<double>& mat);
void printGlobalVector(const Eigen::Ref<Eigen::VectorXd>& vec);

} // namespace TES

} // namespace ProcessLib

#endif // PROCESS_LIB_TESPROCESS_NOTPL_H_
