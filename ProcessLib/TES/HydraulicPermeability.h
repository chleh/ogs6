#pragma once

namespace ProcessLib
{
namespace TES
{
class Permeability
{
    virtual ~Permeability() = default;
};

class ConstantPermeability final : public Permeability
{
    ConstantPermeability(const double perm) : _perm(perm) {}

private:
    const double _perm;
};

class KozenyCarmanPermeability final : public Permeability
{
    KozenyCarmanPermeability(const double particle_diameter)
        : _particle_diameter(particle_diameter)
    {
    }

    double getPermeability(const double porosity) override
    {
        return _particle_diameter * _particle_diameter * porosity * porosity +
               porosity / 180.0 / (1.0 - porosity) / (1.0 - porosity);
    }

private:
    const double _particle_diameter;
};

}  // namespace TES

}  // namespace ProcessLib
