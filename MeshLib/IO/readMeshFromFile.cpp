/**
 * \file
 * \author Karsten Rink
 * \date   2012-09-27
 * \brief  Implementation of readMeshFromFile function.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * @file readMeshFromFile.cpp
 * @date 2012-09-27
 * @author Karsten Rink
 */

#include "readMeshFromFile.h"

#ifdef USE_PETSC
#include <petsc.h>
#endif

#include <boost/algorithm/string/erase.hpp>

#include <logog/include/logog.hpp>

#include "BaseLib/Config.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"

#include "MeshLib/Mesh.h"

#include "MeshLib/IO/DUNE_IO/readDUNEMesh.h"
#include "MeshLib/IO/Legacy/MeshIO.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"

#ifdef USE_PETSC
#include "MeshLib/IO/MPI_IO/NodePartitionedMeshReader.h"
#include "MeshLib/NodePartitionedMesh.h"
#endif

namespace MeshLib
{
namespace IO
{
std::unique_ptr<MeshLib::FEMMesh> readMeshFromFile(const std::string& file_name,
                                                   unsigned dimension,
                                                   unsigned global_refinement)
{
#ifdef USE_PETSC
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (world_size > 1)
    {
        MeshLib::IO::NodePartitionedMeshReader read_pmesh(PETSC_COMM_WORLD);
        const std::string file_name_base = BaseLib::dropFileExtension(file_name);
        return read_pmesh.read(file_name_base);
    }
    else if (world_size == 1)
    {
        MeshLib::Mesh* mesh = readMeshFromFileSerial(file_name);
        MeshLib::NodePartitionedMesh* part_mesh
                               = new MeshLib::NodePartitionedMesh(*mesh);
        delete mesh;
        return part_mesh;
    }
    return nullptr;
#else
    return readMeshFromFileSerial(file_name, dimension, global_refinement);
#endif
}

std::unique_ptr<MeshLib::FEMMesh> readMeshFromFileSerial(
    const std::string& file_name,
    unsigned dimension,
    unsigned global_refinement)
{
    if (BaseLib::hasFileExtension("msh", file_name))
    {
        return readDUNEMeshFromFile(file_name, dimension, global_refinement);
    }

    if (BaseLib::hasFileExtension("vtu", file_name))
        return std::unique_ptr<MeshLib::FEMMesh>(
            MeshLib::IO::VtuInterface::readVTUFile(file_name));

    ERR("readMeshFromFile(): Unknown mesh file format in file %s.", file_name.c_str());
    return nullptr;
}


} // end namespace IO
} // end namespace MeshLib
