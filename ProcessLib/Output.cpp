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

Output*
Output::newInstance(const ConfigTree &config, std::string const& path)
{
	auto const file = config.get_optional<std::string>("file");
	if (!file) {
		ERR("The output config does not provide a file tag.");
		ERR("    Output file not set.");
		return nullptr;
	}

	auto out = new Output(path + *file);

	auto const timesteps = config.get_child_optional("timesteps");
	if (timesteps)
	{
		auto const range = timesteps->equal_range("pair");
		for (auto it=range.first; it!=range.second; ++it)
		{
			auto count      = it->second.get_optional<std::size_t>("count");
			auto each_steps = it->second.get_optional<std::size_t>("each-steps");

			if (count && each_steps) {
				assert(*count != 0 && *each_steps != 0);
				out->_countsSteps.emplace_back(*count, *each_steps);
			} else {
				ERR("<count> or <each-steps> missing in <pair>");
			}
		}

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
