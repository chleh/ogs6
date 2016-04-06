/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "InitialCondition.h"

#include <boost/optional.hpp>
#include <logog/include/logog.hpp>

#include "MathLib/Point3d.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"

#include "BaseLib/ConfigTree.h"

namespace ProcessLib
{
std::unique_ptr<InitialCondition> createUniformInitialCondition(
    BaseLib::ConfigTree const& config, int const n_components)
{
	config.checkConfParam("type", "Uniform");

	std::string const name_prefix = "component_value_";

	std::vector<double> values(n_components);  // parsed values for later ctor.
	for (int component_id = 0; component_id < n_components; ++component_id)
	{
		auto const name = name_prefix + std::to_string(component_id);
		double const value = config.getConfParam<double>(name);
		values[component_id] = value;
		DBUG("Using value %g for component %d.", value, component_id);
	}

	return std::unique_ptr<InitialCondition>{
	    new UniformInitialCondition{values}};
}

std::unique_ptr<InitialCondition> createMeshPropertyInitialCondition(
    BaseLib::ConfigTree const& config,
    MeshLib::Mesh const& mesh,
    int const n_components)
{
	config.checkConfParam("type", "MeshProperty");
	auto field_name = config.getConfParam<std::string>("field_name");
	DBUG("Using field_name %s", field_name.c_str());

	if (!mesh.getProperties().hasPropertyVector(field_name))
	{
		ERR("The required property %s does not exists in the mesh.",
		    field_name.c_str());
		std::abort();
	}
	auto const& property =
	    mesh.getProperties().template getPropertyVector<double>(field_name);
	if (!property)
	{
		ERR("The required property %s is not of the requested type.",
		    field_name.c_str());
		std::abort();
	}

	if (property->getNumberOfComponents() !=
	    static_cast<std::size_t>(n_components))
	{
		ERR("The required property %s has different number of components %d, "
		    "expected %d.",
		    field_name.c_str(), property->getNumberOfComponents(), n_components);
		std::abort();
	}
	return std::unique_ptr<InitialCondition>(
	    new MeshPropertyInitialCondition(*property));
}

}  // namespace ProcessLib
