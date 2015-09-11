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

    double _fluid_specific_heat_source = 0.0;
    double _cpG = 1.0; // specific isobaric fluid heat capacity

    Eigen::MatrixXd _solid_perm_tensor = Eigen::MatrixXd::Constant(3, 3, 0.0); // TODO get dimensions
    double _solid_specific_heat_source = 0.0;
    double _solid_heat_cond = 1.0;
    double _cpS = 1.0;    // specific isobaric solid heat capacity

    double _tortuosity = 1.0;
    double _diffusion_coefficient_component = 1.0; // ???

    double _poro = 0.5;

    double _rho_SR_dry = 1000.0;

    const double _M_inert = M_N2; // N2
    const double _M_react = M_H2O;

    double _initial_solid_density = 1000.0;

    double       _time_step = 1.0;
    bool         _is_new_timestep = true;

    bool _output_element_matrices = false;
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
