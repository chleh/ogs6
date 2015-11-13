/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>
#include <iosfwd>

namespace FileIO
{

class PVDFile
{
public:
    PVDFile(std::string const& fn);

    void addVTUFile(std::string const& fn, double timestep);

    ~PVDFile();

private:
    std::ofstream _fh;
};

}
