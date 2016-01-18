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

#include<cassert>

#include "logog/include/logog.hpp"

#include "Output.h"

namespace ProcessLib
{

// TODO: change to unique_ptr
Output*
Output::newInstance(const BaseLib::ConfigTreeNew &config, std::string const& path)
{
	auto out = new Output(path + config.getConfParam<std::string>("file"));

	if (auto const timesteps = config.getConfSubtreeOptional("timesteps"))
	{
		for (auto pair : timesteps->getConfSubtreeList("pair"))
		{
			auto count      = pair.getConfParam<std::size_t>("count");
			auto each_steps = pair.getConfParam<std::size_t>("each-steps");

			assert(count != 0 && each_steps != 0);
			out->_countsSteps.emplace_back(count, each_steps);
		}

		// TODO: change assertion
		assert(out->_countsSteps.size() != 0);
	}
	else
	{
		out->_countsSteps.emplace_back(1, 1);
	}

	return out;
}


bool
Output::doOutput(std::size_t timestep) const
{
	std::size_t each_steps = 1;

	for (auto const& pair : _countsSteps)
	{
		each_steps = pair.each_steps;

		if (timestep > pair.count * each_steps)
		{
			timestep -= pair.count * each_steps;
		} else {
			break;
		}
	}

	return timestep % each_steps == 0;
}

}
