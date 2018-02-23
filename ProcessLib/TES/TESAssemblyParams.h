/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Eigen>
#include <Eigen/Sparse>

#include "Material/DiffusionCoefficient.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/ReactionKinetics/ReactionRate.h"
#include "MaterialLib/ReactionKinetics/ReactiveSolidModel.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "VolumetricHeatLoss.h"

namespace ProcessLib
{
namespace TES
{
const unsigned NODAL_DOF = 3;
const unsigned COMPONENT_ID_PRESSURE = 0;
const unsigned COMPONENT_ID_TEMPERATURE = 1;
const unsigned COMPONENT_ID_MASS_FRACTION = 2;

struct AssemblyParams
{
    std::unique_ptr<MaterialLib::ReactiveSolidModel> reactive_solid;
    std::unique_ptr<MaterialLib::ReactionRate> reaction_rate;

    double fluid_specific_heat_source =
        std::numeric_limits<double>::quiet_NaN();
    double cpG = std::numeric_limits<
        double>::quiet_NaN();  // specific isobaric fluid heat capacity

    double solid_specific_heat_source =
        std::numeric_limits<double>::quiet_NaN();
    double solid_heat_cond = std::numeric_limits<double>::quiet_NaN();
    double cpS = std::numeric_limits<
        double>::quiet_NaN();  // specific isobaric solid heat capacity

    double tortuosity = std::numeric_limits<double>::quiet_NaN();
    std::unique_ptr<DiffusionCoefficient> diffusion_coefficient_component;

    Parameter<double>* poro = nullptr;
    Parameter<double>* permeability = nullptr;

    double rho_SR_dry = std::numeric_limits<double>::quiet_NaN();

    double rho_SR_lower = std::numeric_limits<double>::quiet_NaN();
    double rho_SR_upper = std::numeric_limits<double>::quiet_NaN();

    const double M_inert = MaterialLib::PhysicalConstant::MolarMass::N2;
    const double M_react = MaterialLib::PhysicalConstant::MolarMass::Water;

    // TODO unify variable names
    double initial_solid_density = std::numeric_limits<double>::quiet_NaN();
    std::string initial_solid_density_mesh_property;

    std::unique_ptr<VolumetricHeatLoss> volumetric_heat_loss;

    bool dielectric_heating_term_enabled = false;
    MathLib::PiecewiseLinearInterpolation heating_power_scaling =
        MathLib::PiecewiseLinearInterpolation(std::vector<double>{},
                                              std::vector<double>{}, true);

    double delta_t = std::numeric_limits<double>::quiet_NaN();

    bool output_element_matrices = false;

    double current_time = std::numeric_limits<double>::quiet_NaN();

    NumLib::LocalToGlobalIndexMap const* dof_table = nullptr;
};

}  // namespace TES

}  // namespace ProcessLib
