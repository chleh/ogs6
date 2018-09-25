#pragma once

#include <memory>
#include <string>

#include "BaseLib/Config.h"
#include "MeshLib/FEMMesh.h"

#if OGS_USE_DUNE
namespace MeshLib
{
namespace IO
{
std::unique_ptr<MeshLib::FEMMesh> readDUNEMeshFromFile(
    const std::string& file_name, unsigned dimension,
    unsigned global_refinement);
}  // namespace IO
}  // namespace MeshLib
#endif  // OGS_USE_DUNE
