/**
 * \author Christoph Lehmann
 * \date   2015-09-01
 * \brief  Do output every n timesteps
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_OUTPUT_H
#define PROCESSLIB_OUTPUT_H

#include "BaseLib/ConfigTreeNew.h"

namespace ProcessLib
{

class Output
{
public:
    static Output* newInstance(const BaseLib::ConfigTreeNew& config, const std::string& path);

    std::string const& getFilePrefix() const { return _output_file_prefix; }

    typedef struct _pairCountSteps {
        _pairCountSteps() = default;
        explicit _pairCountSteps(unsigned c, unsigned e)
            : count(c), each_steps(e) {}
        unsigned count;
        unsigned each_steps;
    } pairCountSteps;

    bool doOutput(std::size_t timestep) const;

private:
    Output(std::string const& prefix)
        : _output_file_prefix(prefix)
    {}

    std::string _output_file_prefix;
    std::vector<pairCountSteps> _countsSteps;
};

}

#endif // PROCESSLIB_OUTPUT_H

