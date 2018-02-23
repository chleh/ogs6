#include "VolumetricHeatLoss.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

namespace ProcessLib
{
namespace TES
{
std::unique_ptr<VolumetricHeatLoss> createVolumetricHeatLoss(
    BaseLib::ConfigTree const& config)
{
    auto const type = config.getConfigParameter<std::string>("type");

    if (type == "ConstantHeatTransferCoefficient")
    {
        auto const h = config.getConfigParameter<double>("h");
        auto const T_amb =
            config.getConfigParameter<double>("ambient_temperature");
        return std::unique_ptr<
            VolumetricHeatLossConstantHeatTransferCoefficient>(
            new VolumetricHeatLossConstantHeatTransferCoefficient{h, T_amb});
    }
    else if (type == "LinearHeatTransferCoefficient")
    {
        auto const h0 = config.getConfigParameter<double>("h0");
        auto const dhdT = config.getConfigParameter<double>("dhdtemp");
        auto const T_amb =
            config.getConfigParameter<double>("ambient_temperature");
        return std::unique_ptr<VolumetricHeatLossLinearHeatTransferCoefficient>(
            new VolumetricHeatLossLinearHeatTransferCoefficient{h0, dhdT,
                                                                T_amb});
    }

    OGS_FATAL("unknown heat loss type `%s'.", type.c_str());
}

}  // TES
}  // ProcessLib
