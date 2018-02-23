#pragma once

#include <memory>
namespace BaseLib
{
class ConfigTree;
}  // BaseLib

namespace ProcessLib
{
namespace TES
{
class VolumetricHeatLoss
{
public:
    virtual double getHeatLoss(const double T) const = 0;
    virtual ~VolumetricHeatLoss() = default;
};

class VolumetricHeatLossConstantHeatTransferCoefficient
    : public VolumetricHeatLoss
{
public:
    VolumetricHeatLossConstantHeatTransferCoefficient(
        const double volumetric_heat_transfer_coefficient, const double T_amb)
        : _volumetric_heat_transfer_coefficient{volumetric_heat_transfer_coefficient},
          _T_amb{T_amb}
    {
    }

    double getHeatLoss(const double T) const override
    {
        return _volumetric_heat_transfer_coefficient * (_T_amb - T);
    }

private:
    const double _volumetric_heat_transfer_coefficient;
    const double _T_amb;
};

class VolumetricHeatLossLinearHeatTransferCoefficient
    : public VolumetricHeatLoss
{
public:
    VolumetricHeatLossLinearHeatTransferCoefficient(const double h0,
                                                    const double dhdT,
                                                    const double T_amb)
        : _h0{h0}, _dhdT{dhdT}, _T_amb{T_amb}
    {
    }

    double getHeatLoss(const double T) const override
    {
        return (_h0 + T * _dhdT) * (_T_amb - T);
    }

private:
    const double _h0;
    const double _dhdT;
    const double _T_amb;
};

std::unique_ptr<VolumetricHeatLoss> createVolumetricHeatLoss(
    BaseLib::ConfigTree const& config);

}  // TES
}  // ProcessLib
