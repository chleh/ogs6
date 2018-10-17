/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateMFront.h"

#include <dlfcn.h>
#include <typeinfo>

namespace
{
void symbolInfo(std::string const& symbol, char const* value)
{
    INFO("Found symbol `%s' with value `%s'.", symbol.c_str(), value);
}

void symbolInfo(std::string const& symbol, double value)
{
    INFO("Found symbol `%s' with value %g.", symbol.c_str(), value);
}

void symbolInfo(std::string const& symbol, short value)
{
    INFO("Found symbol `%s' with value %hi.", symbol.c_str(), value);
}

void symbolInfo(std::string const& symbol, unsigned short value)
{
    INFO("Found symbol `%s' with value %hu.", symbol.c_str(), value);
}

void symbolInfo(std::string const& symbol, int value)
{
    INFO("Found symbol `%s' with value %i.", symbol.c_str(), value);
}

void symbolInfo(std::string const& symbol, unsigned value)
{
    INFO("Found symbol `%s' with value %u.", symbol.c_str(), value);
}

void symbolInfo(std::string const& symbol, long value)
{
    INFO("Found symbol `%s' with value %li.", symbol.c_str(), value);
}

void symbolInfo(std::string const& symbol, unsigned long value)
{
    INFO("Found symbol `%s' with value %lu.", symbol.c_str(), value);
}

template <typename T>
void symbolInfo(std::string const& symbol, T* value)
{
    INFO("Found symbol `%s' with value %p.", symbol.c_str(), value);
}

void* getRawSymbol(void* lib, std::string const& model,
                   std::string const& symbol)
{
    auto const symbol_name = model + "_" + symbol;
    dlerror();
    auto* sym = dlsym(lib, symbol_name.c_str());
    auto* err = dlerror();
    if (err != nullptr)
    {
        OGS_FATAL("Could not find symbol `%s'. The error is: `%s'.",
                  symbol_name.c_str(), err);
    }

    return sym;
}

template <typename T>
T getSymbol(void* lib, std::string const& model, std::string const& symbol);

template <typename T>
T getSymbol(void* lib, std::string const& model, std::string const& symbol, T*)
{
    auto* sym = getRawSymbol(lib, model, symbol);
    auto value = *static_cast<T*>(sym);
    INFO("YYY");
    symbolInfo(symbol, value);

    return value;
}

template <>
int* getSymbol<int*>(void* lib, std::string const& model,
                     std::string const& symbol, int**)
{
    auto* sym = getRawSymbol(lib, model, symbol);
    auto value = static_cast<int*>(sym);
    INFO("XXX");
    symbolInfo(symbol, value);

    return value;
}

template <typename T>
T** getSymbol(void* lib, std::string const& model, std::string const& symbol,
              T***)
{
    auto* sym = getRawSymbol(lib, model, symbol);
    auto value = static_cast<T**>(sym);
    INFO("ZZZ");
    symbolInfo(symbol, value);

    return value;
}

template <typename T>
std::vector<T> getSymbol(void* lib, std::string const& model,
                         std::string const& symbol, std::vector<T>*)
{
    auto const n = getSymbol<unsigned short>(lib, model, 'n' + symbol);
    auto const arr = getSymbol<T*>(lib, model, symbol);
    std::vector<T> vec;
    vec.reserve(n);

    for (unsigned short i = 0; i < n; ++i)
    {
        vec.push_back(arr[i]);
    }
    return vec;
}

template <typename T>
T getSymbol(void* lib, std::string const& model, std::string const& symbol)
{
    return getSymbol(lib, model, symbol, static_cast<T*>(nullptr));
}

std::pair<std::vector<const char*>, std::vector<int>> getNamesAndTypes(
    void* lib, std::string const& model, std::string const& symbol)
{
    auto const n = getSymbol<unsigned short>(lib, model, 'n' + symbol);
    auto const arr_names = getSymbol<const char**>(lib, model, symbol);
    auto const arr_types = getSymbol<int*>(lib, model, symbol + "Types");

    std::vector<const char*> vec_names;
    std::vector<int> vec_types;

    for (unsigned short i = 0; i < n; ++i)
    {
        vec_names.push_back(arr_names[i]);
        vec_types.push_back(arr_types[i]);
    }
    return {vec_names, vec_types};
}

}  // namespace

namespace MaterialLib
{
namespace Solids
{
namespace MFront
{
template <int DisplacementDim>
std::unique_ptr<MFront<DisplacementDim>> createMFront(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
    INFO("### MFRONT ########################################################");

    config.checkConfigParameter("type", "MFront");

    // TODO pass prj file dir.
    auto const lib_path = config.getConfigParameter<std::string>("library");
    auto const model_name = config.getConfigParameter<std::string>("model");

    dlerror();
    auto* lib = dlopen(lib_path.c_str(), RTLD_NOW);
    auto* err = dlerror();
    if (err != nullptr)
    {
        OGS_FATAL("Could not load shared library `%s'. The error is: `%s'.",
                  lib_path.c_str(), err);
    }

    if (getSymbol<const char*>(lib, model_name, "mfront_ept") != model_name)
    {
        OGS_FATAL("Loaded shared library is internally inconsistent.");
    }

#if 1
    getSymbol<const char*>(lib, model_name, "tfel_version");

    getSymbol<unsigned short>(lib, model_name, "mfront_mkt");

    getSymbol<const char*>(lib, model_name, "mfront_interface");

    getSymbol<const char*>(lib, model_name, "src");

    auto const hypotheses = getSymbol<std::vector<const char*>>(
        lib, model_name, "ModellingHypotheses");
    for (auto* str : hypotheses)
    {
        INFO("  --> `%s'", str);
    }

    getSymbol<unsigned short>(lib, model_name, "BehaviourType");

    getSymbol<unsigned short>(lib, model_name, "BehaviourKinematic");

    getSymbol<unsigned short>(lib, model_name, "SymmetryType");

    getSymbol<unsigned short>(lib, model_name, "ElasticSymmetryType");

    getSymbol<unsigned short>(lib, model_name, "savesTangentOperator");
    getSymbol<unsigned short>(lib, model_name,
                              "UsableInPurelyImplicitResolution");

    auto const properties = getSymbol<std::vector<const char*>>(
        lib, model_name, "MaterialProperties");
    for (auto* str : properties)
    {
        INFO("  --> `%s'", str);
    }

    /*
    auto const internal_vars = getSymbol<std::vector<const char*>>(
        lib, model_name, "InternalStateVariables");
    for (auto* str : internal_vars)
    {
        INFO("  --> `%s'", str);
    }*/
    auto const internal_vars =
        getNamesAndTypes(lib, model_name, "InternalStateVariables");
    for (std::size_t i = 0; i < internal_vars.first.size(); ++i)
    {
        INFO("  --> `%s' type: %i.", internal_vars.first[i],
             internal_vars.second[i]);
    }

    auto const external_vars = getSymbol<std::vector<const char*>>(
        lib, model_name, "ExternalStateVariables");
    for (auto* str : external_vars)
    {
        INFO("  --> `%s'", str);
    }

    auto const parameters_ =
        getSymbol<std::vector<const char*>>(lib, model_name, "Parameters");
    for (auto* str : parameters_)
    {
        INFO("  --> `%s'", str);
    }

    /*
getSymbol<int>(lib, model_name, "InternalStateVariablesTypes[] = {1, 0, 0, 0,
                                                           0, 0, 0};

getSymbol<int>(lib, model_name, "ParametersTypes[] = {0, 0, 0, 0, 0, 0, 2};
                                                                           */

    getSymbol<double>(lib, model_name, "theta_ParameterDefaultValue");

    getSymbol<double>(lib, model_name, "Tref_ParameterDefaultValue");

    getSymbol<double>(lib, model_name,
                      "minimal_time_step_scaling_factor_ParameterDefaultValue");

    getSymbol<double>(lib, model_name,
                      "maximal_time_step_scaling_factor_ParameterDefaultValue");

    getSymbol<double>(lib, model_name, "epsilon_ParameterDefaultValue");

    getSymbol<double>(lib, model_name,
                      "numerical_jacobian_epsilon_ParameterDefaultValue");

    getSymbol<unsigned short>(lib, model_name, "iterMax_ParameterDefaultValue");

    getSymbol<unsigned short>(lib, model_name, "requiresStiffnessTensor");
    getSymbol<unsigned short>(lib, model_name,
                              "requiresThermalExpansionCoefficientTensor");

    /*
getSymbol<int>(lib, model_name, "setParameter(const char* const key,
                                       const double value)
{
using tfel::material::BDTParametersInitializer;
auto& i = BDTParametersInitializer::get();
try
{
    i.set(key, value);
}
catch (std::runtime_error& e)
{
    std::cerr << e.what() << std::endl;
    return 0;
}
return 1;
}

getSymbol<int>(lib, model_name, "setUnsignedShortParameter(
const char* const key, const unsigned short value)
{
using tfel::material::BDTParametersInitializer;
auto& i = BDTParametersInitializer::get();
try
{
    i.set(key, value);
}
catch (std::runtime_error& e)
{
    std::cerr << e.what() << std::endl;
    return 0;
}
return 1;
}

getSymbol<void>(lib, model_name, "setOutOfBoundsPolicy(const int p)
{
if (p == 0)
{
   >(lib, model_name, "getOutOfBoundsPolicy() = tfel::material::None;
}
else if (p == 1)
{
   >(lib, model_name, "getOutOfBoundsPolicy() = tfel::material::Warning;
}
else if (p == 2)
{
   >(lib, model_name, "getOutOfBoundsPolicy() = tfel::material::Strict;
}
else
{
    std::cerr << "asterbdt_setOutOfBoundsPolicy : invalid argument\n";
}
}
    */

#if 0

    char*>(lib, model_name, "getIntegrationErrorMessage()
    {
#if (defined __GNUC__) && (!defined __clang__) && \
    (!defined __INTEL_COMPILER) && (!defined __PGI)
#if __GNUC__ * 10000 + __GNUC_MINOR__ * 100 < 40800
        static __thread char msg[128];
#else
        static thread_local char msg[128];
#endif
#else  /* (defined __GNUC__) ...*/
        static thread_local char msg[128];
#endif /* (defined __GNUC__) ...*/
        return msg;
    }  // end of>(lib, model_name, "getIntegrationErrorMessage

    getSymbol<void BDT(
        aster::AsterReal* const STRESS, aster::AsterReal* const STATEV,
        aster::AsterReal* const DDSOE, const aster::AsterReal* const STRAN,
        const aster::AsterReal* const DSTRAN,
        const aster::AsterReal* const DTIME, const aster::AsterReal* const TEMP,
        const aster::AsterReal* const DTEMP,
        const aster::AsterReal* const PREDEF,
        const aster::AsterReal* const DPRED, const aster::AsterInt* const NTENS,
        const aster::AsterInt* const NSTATV,
        const aster::AsterReal* const PROPS,
        const aster::AsterInt* const NPROPS, const aster::AsterReal* const DROT,
        aster::AsterReal* const PNEWDT, const aster::AsterInt* const NUMMOD)
    {
        char* msg =>(lib, model_name, "getIntegrationErrorMessage();
        if (aster::AsterInterface<tfel::material::BDT>::exe(
                msg, NTENS, DTIME, DROT, DDSOE, STRAN, DSTRAN, TEMP, DTEMP,
                PROPS, NPROPS, PREDEF, DPRED, STATEV, NSTATV, STRESS, NUMMOD,
               >(lib, model_name, "getOutOfBoundsPolicy(),
                aster::AsterStandardSmallStrainStressFreeExpansionHandler) != 0)
        {
            *PNEWDT = -1.;
            return;
        }
    }

    getSymbol<void asterbdt(
        aster::AsterReal* const STRESS, aster::AsterReal* const STATEV,
        aster::AsterReal* const DDSOE, const aster::AsterReal* const STRAN,
        const aster::AsterReal* const DSTRAN,
        const aster::AsterReal* const DTIME, const aster::AsterReal* const TEMP,
        const aster::AsterReal* const DTEMP,
        const aster::AsterReal* const PREDEF,
        const aster::AsterReal* const DPRED, const aster::AsterInt* const NTENS,
        const aster::AsterInt* const NSTATV,
        const aster::AsterReal* const PROPS,
        const aster::AsterInt* const NPROPS, const aster::AsterReal* const DROT,
        aster::AsterReal* const PNEWDT, const aster::AsterInt* const NUMMOD)
    {
        BDT(STRESS, STATEV, DDSOE, STRAN, DSTRAN, DTIME, TEMP, DTEMP, PREDEF,
            DPRED, NTENS, NSTATV, PROPS, NPROPS, DROT, PNEWDT, NUMMOD);
    }
#endif

#endif

    INFO("### MFRONT END ####################################################");

    return std::make_unique<MFront<DisplacementDim>>();
}

template std::unique_ptr<MFront<2>> createMFront<2>(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);
template std::unique_ptr<MFront<3>> createMFront<3>(
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);

}  // namespace MFront
}  // namespace Solids
}  // namespace MaterialLib
