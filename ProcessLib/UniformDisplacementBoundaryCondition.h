/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_UNIFORM_DISPLACEMENT_BOUNDARY_CONDITION_H_
#define PROCESS_LIB_UNIFORM_DISPLACEMENT_BOUNDARY_CONDITION_H_

#include <algorithm>
#include <vector>

#include "logog/include/logog.hpp"

#include "BaseLib/ConfigTree.h"
#include "AssemblerLib/LocalToGlobalIndexMap.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"

namespace GeoLib
{
class GeoObject;
}

namespace ProcessLib
{
/// The UniformDisplacementBoundaryCondition class describes a constant in space
/// and time Displacement boundary condition.
/// The expected parameter in the passed configuration is "value" which, when
/// not present defaults to zero.
class UniformDisplacementBoundaryCondition
{
public:
	UniformDisplacementBoundaryCondition(
	    GeoLib::GeoObject const* const geometry,
	    std::map<std::string,
	             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
	        curves,
	    BaseLib::ConfigTree& config, int const tuple_size)
	{
		DBUG("Constructing UniformDisplacementBoundaryCondition from config.");
		config.checkConfParam("type", "UniformDisplacement");

		std::string const name_prefix = "component_value_";

		assert(0 <= tuple_size && tuple_size < 3);
		for (int i = 0; i < tuple_size; ++i)
		{
			auto const name = name_prefix + std::to_string(i);
			if (config.exists(name))
			{
				double const value = config.getConfParam<double>(name);
				_bcs[i].reset(
				    new UniformDirichletBoundaryCondition(geometry, value));
			}
		}
	}

	std::unique_ptr<UniformDirichletBoundaryCondition>
	get(int const component_id)
	{
		return std::move(_bcs[component_id]);
	}

private:
	std::array<std::unique_ptr<UniformDirichletBoundaryCondition>, 3> _bcs;
};

}  // namespace ProcessLib

#endif  // PROCESS_LIB_UNIFORM_DISPLACEMENT_BOUNDARY_CONDITION_H_
