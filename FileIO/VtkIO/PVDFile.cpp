/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <fstream>
#include <iomanip>
#include "PVDFile.h"

#include "logog/include/logog.hpp"

namespace FileIO
{

PVDFile::PVDFile(const std::string &fn)
    : _fh{fn.c_str()}
{
    if (!_fh) {
        ERR("could not open file `%s'", fn.c_str());
    }

    _fh << std::setprecision(16);

    _fh << "<?xml version=\"1.0\"?>\n"
           "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n"
           "  <Collection>\n";
}

void PVDFile::addVTUFile(const std::string &fn, double timestep)
{
    _fh << "    <DataSet timestep=\"" << timestep << "\" group=\"\" part=\"0\" file=\"" << fn << "\"/>\n";
}

PVDFile::~PVDFile()
{
    _fh << "  </Collection>\n</VTKFile>\n";
}

}
