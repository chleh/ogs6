/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALLIB_DUBININPOLANYIDATA_H
#define MATERIALLIB_DUBININPOLANYIDATA_H

#include "functional"

namespace BaseLib
{
class ConfigTree;
}  // namespace BaseLib

namespace MaterialLib
{
class DubininPolanyiData final
{
public:
    using Fct = double (*)(const double);
    using EnthalpyFct = double (*)(double, double);

    DubininPolanyiData(const double M_Ads_, const Fct rho_Ads,
                       const Fct alpha_T_Ads, const Fct char_curve,
                       const Fct char_curve_deriv,
                       const EnthalpyFct custom_enthalpy_function)
        : M_Ads(M_Ads_),
          custom_enthalpy_method(custom_enthalpy_function),
          _rho_Ads(rho_Ads),
          _alpha_T_Ads(alpha_T_Ads),
          _char_curve(char_curve),
          _char_curve_deriv(char_curve_deriv)
    {
    }
    double getAdsorbateDensity(const double T_Ads) const
    {
        return _rho_Ads(T_Ads);
    }
    double getAdsorbateThermalExpansion(const double T_Ads) const
    {
        return _alpha_T_Ads(T_Ads);
    }
    double characteristicCurve(const double A) const { return _char_curve(A); }
    double characteristicCurveDeriv(const double A) const
    {
        return _char_curve_deriv(A);
    }

    const double M_Ads;
    const EnthalpyFct custom_enthalpy_method;

private:
    const Fct _rho_Ads;
    const Fct _alpha_T_Ads;
    const Fct _char_curve;
    const Fct _char_curve_deriv;
};

DubininPolanyiData createDubininPolanyiData(BaseLib::ConfigTree const& config);
DubininPolanyiData createDubininPolanyiData(std::string const& working_pair);

}  // namespace MaterialLib

#endif  // MATERIALLIB_DUBININPOLANYIDATA_H
