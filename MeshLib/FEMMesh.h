#pragma once

#include <string>

namespace MeshLib
{
class FEMMesh
{
public:
    virtual std::size_t getID() const = 0;
    virtual std::string const& getName() const = 0;

    virtual bool isAxiallySymmetric() const { return false; }
    virtual void setAxiallySymmetric(bool /*is_axial_symmetric*/) {}

    virtual void globalRefine() {}

    virtual unsigned getDimension() const = 0;

    virtual ~FEMMesh() = default;
};

}  // namespace MeshLib
