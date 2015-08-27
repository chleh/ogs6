#ifndef PROCESS_LIB_TESPROCESS_NOTPL_H_
#define PROCESS_LIB_TESPROCESS_NOTPL_H_


#include "MaterialsLib/adsorption/adsorption.h"
#include "Eigen/Sparse"
#include "Eigen/Eigen"

namespace ProcessLib
{

namespace TES
{

const unsigned NODAL_DOF = 3;
const unsigned NODAL_DOF_2ND = 2; // loading or solid density, and reaction rate

const double M_N2  = 0.028013;
const double M_H2O = 0.018016;

struct AssemblyParams
{
    Ads::Adsorption* _adsorption;

    const double _fluid_specific_heat_source = 0.0;
    const double _cpG = 1012.0; // specific isobaric fluid heat capacity

    const Eigen::MatrixXd _solid_perm_tensor = Eigen::MatrixXd::Identity(3, 3) * 1.e-7; // TODO get dimensions
    const double _solid_specific_heat_source = 0.0;
    const double _solid_heat_cond = 0.4;
    const double _cpS = 880.0;    // specific isobaric solid heat capacity

    const double _tortuosity = 1.0;
    const double _diffusion_coefficient_component = 9.65e-2; // ???

    const double _poro = 0.7;

    const double _rho_SR_dry = 1150.0;

    const double _M_inert = M_N2; // N2
    const double _M_react = M_H2O;

    const double   _initial_solid_density = 1382.36248218;

    const double   _time_step = 0.5;
    bool  _is_new_timestep = true;
};


class TESProcessInterface
{
public:
    /*AssemblyParams& getAssemblyParams() {
        return _assembly_params;
    }*/
    AssemblyParams const& getAssemblyParams() const {
        return _assembly_params;
    }

    virtual ~TESProcessInterface() = default;

protected:
    AssemblyParams _assembly_params;
};


bool calculateError(Eigen::VectorXd* current_solution,
                    const Eigen::Ref<Eigen::VectorXd>& previous_solution, AssemblyParams* materials);
// bool calculateError(const Eigen::SparseMatrix<double>& current_solution,
//                     const Eigen::SparseMatrix<double>& previous_solution);


void printGlobalMatrix(const Eigen::SparseMatrix<double>& mat);
void printGlobalVector(const Eigen::Ref<Eigen::VectorXd>& vec);

} // namespace TES

} // namespace ProcessLib

#endif // PROCESS_LIB_TESPROCESS_NOTPL_H_
