#include "readDUNEMesh.h"

#if OGS_USE_DUNE

#include "BaseLib/DUNEConfig.h"
#include "BaseLib/Error.h"
#include "BaseLib/FileTools.h"

#include "MeshLib/DUNEMesh.h"

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/uggrid.hh>

namespace MeshLib
{
namespace IO
{
std::unique_ptr<MeshLib::FEMMesh> readDUNEMeshFromFile(
    const std::string& file_name, unsigned dimension,
    unsigned global_refinement)
{
    if (!BaseLib::hasFileExtension("msh", file_name))
        OGS_FATAL("Wrong file extension (%s).", file_name.c_str());

    auto const mesh_name = BaseLib::extractBaseNameWithoutExtension(file_name);

    switch (dimension)
    {
        /*case 1:
        {
            std::unique_ptr<BaseLib::DUNEGridType<1>> grid(
                Dune::GmshReader<BaseLib::DUNEGridType<1>>::read(file_name));
            return std::make_unique<MeshLib::DUNEMesh<1>>(std::move(grid));
        }*/
        case 2:
        {
            std::unique_ptr<BaseLib::DUNEGridType<2>> grid(
                Dune::GmshReader<BaseLib::DUNEGridType<2>>::read(file_name));
            DBUG("DUNE Mesh %p read", grid.get());
            return std::make_unique<MeshLib::DUNEMesh<2>>(
                mesh_name, std::move(grid), global_refinement);
        }
        case 3:
        {
            std::unique_ptr<BaseLib::DUNEGridType<3>> grid(
                Dune::GmshReader<BaseLib::DUNEGridType<3>>::read(file_name));
            return std::make_unique<MeshLib::DUNEMesh<3>>(
                mesh_name, std::move(grid), global_refinement);
        }
    }

    OGS_FATAL("Wrong mesh dimension");
}
}  // namespace IO
}  // namespace MeshLib
#endif  // OGS_USE_DUNE
