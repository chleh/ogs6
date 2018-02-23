
#include "DiffusionCoefficient.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

namespace ProcessLib
{
namespace TES
{
std::unique_ptr<DiffusionCoefficient> createDiffusionCoefficient(
    BaseLib::ConfigTree const& config)
{
    auto const type = config.getConfigParameter<std::string>("type");

    if (type == "Constant")
        return std::unique_ptr<DiffusionCoefficient>(
            new DiffusionCoefficientConstant(
                config.getConfigParameter<double>("value")));
    else if (type == "WaterNitrogenMarrero")
        return std::unique_ptr<DiffusionCoefficient>(
            new DiffusionCoefficientWaterNitrogenMarrero());
    else if (type == "Knudsen")
    {
        auto const pore_diam =
            config.getConfigParameter<double>("pore_diameter");
        auto const molar_mass = config.getConfigParameter<double>("molar_mass");
        return std::unique_ptr<DiffusionCoefficient>(
            new DiffusionCoefficientKnudsen(pore_diam, molar_mass));
    }
    else if (type == "GasAndKnudsen")
    {
        auto D_gas = createDiffusionCoefficient(
            config.getConfigSubtree("diffusion_coefficient_gas"));
        auto D_Knudsen = createDiffusionCoefficient(
            config.getConfigSubtree("diffusion_coefficient_knudsen"));

        return std::unique_ptr<DiffusionCoefficient>(
            new DiffusionCoefficientGasAndKnudsen(std::move(D_gas),
                                                  std::move(D_Knudsen)));
    }

    OGS_FATAL("Unknown type of diffusion coefficient `%s'", type.c_str());
}
}
}
