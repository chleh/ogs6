/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DubininPolanyiData.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

#include "MaterialLib/Adsorption/Density100MPa.h"
#include "MaterialLib/Adsorption/DensityConst.h"
#include "MaterialLib/Adsorption/DensityCook.h"
#include "MaterialLib/Adsorption/DensityDubinin.h"
#include "MaterialLib/Adsorption/DensityHauer.h"
#include "MaterialLib/Adsorption/DensityHauer15CaCl2CaX.h"
#include "MaterialLib/Adsorption/DensityHauerCaX80NoOutliers.h"
#include "MaterialLib/Adsorption/DensityHauerNaYStach.h"
#include "MaterialLib/Adsorption/DensityHauerNaYStachCubicAndExp.h"
#include "MaterialLib/Adsorption/DensityLegacy.h"
#include "MaterialLib/Adsorption/DensityMette.h"
#include "MaterialLib/Adsorption/DensityNunez.h"

// Cf. http://stackoverflow.com/a/16824239
#include <type_traits>

// Primary template with a static assertion
// for a meaningful error message
// if it ever gets instantiated.
// We could leave it undefined if we didn't care.

template <typename, typename T>
struct HasGetEnthalpy
{
    static_assert(std::integral_constant<T, false>::value,
                  "Second template parameter needs to be of function type.");
};

// specialization that does the checking

template <typename C, typename Ret, typename... Args>
struct HasGetEnthalpy<C, Ret(Args...)>
{
private:
    template <typename T>
    static constexpr auto check(T*) -> typename std::is_same<
        decltype(std::declval<T>().getEnthalpy(std::declval<Args>()...)),
        Ret       // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        >::type;  // attempt to call it and see if the return type is correct

    template <typename>
    static constexpr std::false_type check(...);

    typedef decltype(check<C>(nullptr)) type;

public:
    static constexpr bool value = type::value;
};

template <typename C, bool>
struct GetEnthalpyFunction
{
    static constexpr double (*enthalpy_function)(double, double) = nullptr;
};

template <typename C>
struct GetEnthalpyFunction<C, true>
{
    static constexpr double (*enthalpy_function)(double,
                                                 double) = &C::getEnthalpy;
};

template <typename DataClass>
MaterialLib::DubininPolanyiData createDubininPolanyiData()
{
    using Fct = MaterialLib::DubininPolanyiData::Fct;
    const Fct rho_Ads = &DataClass::getAdsorbateDensity;
    const Fct alpha_T_Ads = &DataClass::getAlphaT;
    const Fct char_curve = &DataClass::characteristicCurve;
    const Fct char_curve_deriv = &DataClass::dCharacteristicCurve;
    const double M_Ads = DataClass::M_Ads;

    using EnthalpyFct = MaterialLib::DubininPolanyiData::EnthalpyFct;
    const EnthalpyFct custom_enthalpy_function =
        GetEnthalpyFunction<DataClass,
                            HasGetEnthalpy<DataClass, double(double, double)>::
                                value>::enthalpy_function;

    return MaterialLib::DubininPolanyiData(M_Ads,
                                           rho_Ads,
                                           alpha_T_Ads,
                                           char_curve,
                                           char_curve_deriv,
                                           custom_enthalpy_function);
}

namespace MaterialLib
{
DubininPolanyiData createDubininPolanyiData(BaseLib::ConfigTree const& config)
{
    auto const working_pair = config.getValue<std::string>();
    return createDubininPolanyiData(working_pair);
}

DubininPolanyiData createDubininPolanyiData(std::string const& working_pair)
{
    if (working_pair == "Z13XBF")
        return ::createDubininPolanyiData<Adsorption::DensityLegacy>();
    else if (working_pair == "Z13XBF_100MPa")
        return ::createDubininPolanyiData<Adsorption::Density100MPa>();
    else if (working_pair == "Z13XBF_Const")
        return ::createDubininPolanyiData<Adsorption::DensityConst>();
    else if (working_pair == "Z13XBF_Cook")
        return ::createDubininPolanyiData<Adsorption::DensityCook>();
    else if (working_pair == "Z13XBF_Dubinin")
        return ::createDubininPolanyiData<Adsorption::DensityDubinin>();
    else if (working_pair == "Z13XBF_Hauer")
        return ::createDubininPolanyiData<Adsorption::DensityHauer>();
    else if (working_pair == "Z13XBF_Mette")
        return ::createDubininPolanyiData<Adsorption::DensityMette>();
    else if (working_pair == "Z13XBF_Nunez")
        return ::createDubininPolanyiData<Adsorption::DensityNunez>();
    else if (working_pair == "CaX80NoOutliers_Hauer")
        return ::createDubininPolanyiData<
            Adsorption::DensityHauerCaX80NoOutliers>();
    else if (working_pair == "15CaCl2CaX_Hauer")
        return ::createDubininPolanyiData<Adsorption::DensityHauer15CaCl2CaX>();
    else if (working_pair == "NaYStach_Hauer")
        return ::createDubininPolanyiData<Adsorption::DensityHauerNaYStach>();
    else if (working_pair == "NaYStach_Hauer_CubicAndExp")
        return ::createDubininPolanyiData<
            Adsorption::DensityHauerNaYStachCubicAndExp>();

    OGS_FATAL("There is no data for the adsorption material `%s'.",
              working_pair.c_str());
}

}  // namespace MaterialLib
