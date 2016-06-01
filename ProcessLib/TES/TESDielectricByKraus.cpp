/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TESDielectricByKraus.h"
#include <cmath>

const double kB = 1.38064852e-23;                 // J/K
const double electron_charge = 1.6021766208e-19;  // C
const double pi = 3.141592653589793238462643383279502884197169399;

// values taken from Kraus' PhD thesis

// properties of the NaY zeolite
const double radio_freq = 13.56e6;  // Hz

const double t_relax0_II_HF = 5e-15;   // s
const double t_relax0_III_HF = 2e-12;  // s

const double diff_eps_II_HF = 0.2;    // no unit
const double diff_eps_III_HF = 0.54;  // no unit

const double EA_II_HF = 0.74 * electron_charge;  // J

const double a_Vert = -1.18;  // (Ma.-%)^-1

// properties of the experimental setup
const double d_inner = 0.007;  // m

/*! Kraus p. 62
 *
 * \param loading zeolite loading, not a percentage!
 *
 * \return activation energy in Joule
 */
double get_EA_III_HF(const double loading)
{
    return electron_charge * (0.46 - 0.07 * loading * 100.0);
}

/*! Kraus formula (5.2)
 *
 * \param diff_eps \f$ \epsilon_{r,stat} - \epsilon_{r,\infty} \f$
 * \param omega    angular frequency of electric field.
 * \param t_relax  relaxation time
 */
double getImPartOfPermittivityDebye(const double diff_eps, const double omega,
                                    const double t_relax)
{
    return diff_eps * omega * t_relax /
           (1.0 + omega * omega * t_relax * t_relax);
}

//! Kraus formula (5.6)
double getImPartOfPermittivityHF_II(const double T)
{
    const double omega = 2 * pi * radio_freq;

    const double t_relax_II_HF = t_relax0_II_HF * std::exp(EA_II_HF / kB / T);
    const double imPart_II_HF =
        getImPartOfPermittivityDebye(diff_eps_II_HF, omega, t_relax_II_HF);

    return imPart_II_HF;
}

//! Kraus formula (5.6)
double getImPartOfPermittivityHF_III(const double T, const double loading)
{
    const double omega = 2 * pi * radio_freq;

    const double t_relax_III_HF =
        t_relax0_III_HF * std::exp(get_EA_III_HF(loading) / kB / T);
    const double imPart_III_HF =
        getImPartOfPermittivityDebye(diff_eps_III_HF, omega, t_relax_III_HF);

    return imPart_III_HF;
}

//! Kraus formula (5.6)
double getImPartOfPermittivityHF(const double T, const double loading)
{
    const double w = std::exp(a_Vert * loading * 100.0);  // weighting factor
    return w * getImPartOfPermittivityHF_II(T) +
           (1.0 - w) * getImPartOfPermittivityHF_III(T, loading);
}

//! Kraus formula (5.30)
double getImPartOfPermittivityZ(const double T)
{
    return -0.03 + 3e-4 * T;
}

//! Kraus formula (5.31)
double getImPartOfPermittivityBulk(const double loading)
{
    if (loading < 0.09)
        return 0.0;

    return 0.0075 * (100.0 * loading - 9.0);
}

//! Kraus formula (5.29)
double getImPartOfPermittivity(const double T, const double loading)
{
    return getImPartOfPermittivityHF(T, loading) + getImPartOfPermittivityZ(T) +
           getImPartOfPermittivityBulk(loading);
}

//! Kraus formula (5.28)
//! \return factor of the heating power in W/m.
double getZ(const double loading)
{
    return 433.0 + 2.1 * loading * 100.0;
}

namespace ProcessLib
{
namespace TES
{
double getVolumetricJouleHeatingPower(const double T, const double loading)
{
    // Note: The term m_ZS/rho_ZS is not needed here, sinve we return a
    // volumetric value!
    return getZ(loading) * getImPartOfPermittivity(T, loading) * 4 / pi /
           d_inner / d_inner;
}

}  // namespace TES
}  // namespace ProcessLib
