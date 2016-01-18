/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * The code of this file is used to decouple the evaluation of matrix elements from the rest of OGS6,
 * not all of OGS6 has to be recompiled every time a small change is done.
 */

#ifndef PROCESS_LIB_TES_FEM_NOTPL_H_
#define PROCESS_LIB_TES_FEM_NOTPL_H_

#include <memory>
#include <Eigen/Eigen>

#include "TESProcess-notpl.h"


namespace ProcessLib
{

namespace TES
{

enum class SecondaryVariables {
    SOLID_DENSITY, REACTION_RATE,
    VELOCITY_X, VELOCITY_Y, VELOCITY_Z,
    REACTION_KINETIC_INDICATOR,
    VAPOUR_PARTIAL_PRESSURE,
    RELATIVE_HUMIDITY,
    LOADING,
    EQUILIBRIUM_LOADING,
    REACTION_DAMPING_FACTOR
};

const double NONLINEAR_ERROR_TOLERANCE = 1e-6;


/**
 * y ... variable in global matrix
 * x ... "physical" process variable in local assembly
 *
 * x = exp(y), dx = dx/dy * dy
 * dx/dy = exp(y) = x
 */
struct TrafoLog
{
    static const bool constrained = true;

    /// Converts global matrix entry to "physical" variable
    /// used in local assembly.
    static double x(const double y) { return std::exp(y); }

    /// Derivative of the "physical" variable x w.r.t y.
    /// the argument is x!
    static double dxdy(const double x) { return x; }
};

struct TrafoIdentity
{
    static const bool constrained = false;

    /// Converts global matrix entry to "physical" variable
    /// used in local assembly.
    static double x(const double y) { return y; }

    /// Derivative of the "physical" variable x w.r.t y.
    /// the argument is x!
    constexpr static double dxdy(const double /*x*/) { return 1.0; }
};

struct TrafoTanh
{
    static const bool constrained = true;

    /// Converts global matrix entry to "physical" variable
    /// used in local assembly.
    static double x(const double y) { return 0.5 * std::tanh(y) + 0.5; }

    /// Derivative of the "physical" variable x w.r.t y.
    /// the argument is x!
    constexpr static double dxdy(const double x) { return 2.0*x*(1.0-x); }
};

typedef TrafoIdentity Trafo;


class LADataNoTpl
{
public:
    typedef Eigen::Ref<const Eigen::MatrixXd> MatRef;
    typedef Eigen::Ref<const Eigen::VectorXd> VecRef;
    typedef std::shared_ptr<std::vector<double> > SharedVector;

    void assembleIntegrationPoint(
            unsigned integration_point,
            Eigen::MatrixXd* localA,
            Eigen::VectorXd* localRhs,
            std::vector<double> const& localX,
            VecRef const& smN,
            MatRef const& smDNdx,
            MatRef const& smJ,
            const double smDetJ,
            const double weight
            );

    // TESProcessInterface const* _process;
    AssemblyParams const* _AP;

    void init(const unsigned num_int_pts, const unsigned dimension);

    void preEachAssemble();
    void postEachAssemble(Eigen::MatrixXd* localA, Eigen::VectorXd* localRhs,
                          const Eigen::VectorXd& oldX);

    std::shared_ptr<const std::vector<double> >
    getIntegrationPointValues(SecondaryVariables var) const;

    double reaction_damping_factor = 1.0;
    std::vector<bool> bounds_violation;

private:
    Eigen::Matrix3d getMassCoeffMatrix(const unsigned int_pt);
    Eigen::MatrixXd getLaplaceCoeffMatrix(const unsigned int_pt, const unsigned dim);
    Eigen::Matrix3d getAdvectionCoeffMatrix(const unsigned int_pt);
    Eigen::Matrix3d getContentCoeffMatrix(const unsigned int_pt);
    Eigen::Vector3d getRHSCoeffVector(const unsigned int_pt);

    void preEachAssembleIntegrationPoint(
            const unsigned int_pt,
            std::vector<double> const& localX,
            VecRef const& smN,
            MatRef const& smDNdx,
            const MatRef& smJ,
            const double smDetJ
            );

    void initReaction(
            const unsigned int_pt,
            std::vector<double> const& localX,
            const MatRef& smDNdx,
            const MatRef& smJ,
            const double smDetJ);

    void initReaction_localVapourUptakeStrategy(const unsigned int_pt);

    void initReaction_localDiffusionStrategy(
            const unsigned int_pt,
            std::vector<double> const& localX,
            const MatRef& smDNdx,
            const MatRef& smJ,
            const double smDetJ);

    void initReaction_simpleStrategy(const unsigned int_pt);

    void initReaction_readjustEquilibriumLoadingStrategy(const unsigned int_pt);

    void initReaction_slowDownUndershootStrategy(const unsigned int_pt);

    /// returns estimated equilibrium vapour pressure
    /// based on a local (i.e. no diffusion/advection) balance
    /// of adsorbate loading and vapour partial pressure
    double estimateAdsorptionEquilibrium(const double p_V0, const double C0) const;

    // nodal quantities, secondary variables
    std::vector<double> _solid_density;
    std::vector<double> _solid_density_prev_ts;

    std::vector<double> _reaction_rate; // dC/dt * _rho_SR_dry
    std::vector<double> _reaction_rate_prev_ts; // could also be calculated from _solid_density_prev_ts

    std::vector<double> _equilibrium_loading;
    std::vector<double> _equilibrium_loading_prev_ts;

    std::vector<bool>   _is_equilibrium_reaction;   ///< true if equilibrium reaction is used in this timestep

    /** the value of p_V that the equilibrium reaction estimated
     *  in the first iteration of this timestep */
    std::vector<double> _estimated_vapour_pressure;

    std::vector<std::vector<double> > _velocity;
    // Eigen::MatrixXd _velocity; // row index: gauss point, column index: dimension x/y/z

    std::vector<double> _reaction_rate_indicator; // TODO [CL] get rid of this

    bool is_var_out_of_bounds = false;

    // bool this_is_repeated = false;
    // bool last_was_repeated = false;

    // integration point values of unknowns
    double _p = -888.888; // gas pressure
    double _T = -888.888; // temperature
    double _vapour_mass_fraction = -888.888;     // fluid mass fraction of the second component

    // temporary storage for some properties
    // values do not change during the assembly of one integration point
    double _rho_GR = -888.888;
    double _p_V = -888.888; // vapour partial pressure
    double _qR = 88888.88888;  // reaction rate, use this in assembly!!!

    std::unique_ptr<Eigen::MatrixXd> _Lap;
    std::unique_ptr<Eigen::MatrixXd> _Mas;
    std::unique_ptr<Eigen::MatrixXd> _Adv;
    std::unique_ptr<Eigen::MatrixXd> _Cnt;
    std::unique_ptr<Eigen::VectorXd> _rhs;
};


void
ogs5OutVec(const LADataNoTpl::VecRef& vec);

void
ogs5OutMat(const LADataNoTpl::MatRef& vec);


} // namespace TES

} // namespace ProcessLib

#endif // PROCESS_LIB_TES_FEM_NOTPL_H_
