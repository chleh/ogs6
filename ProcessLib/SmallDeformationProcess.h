/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_SMALLDEFORMATIONPROCESS_H_
#define PROCESS_LIB_SMALLDEFORMATIONPROCESS_H_

#include <cassert>
#include <memory>

#include <boost/optional.hpp>
#include <boost/algorithm/string/erase.hpp>

#include "logog/include/logog.hpp"
#include "BaseLib/ConfigTree.h"

#ifdef USE_PETSC
#include "MeshLib/NodePartitionedMesh.h"
#include "MathLib/LinAlg/PETSc/PETScMatrixOption.h"
#endif

#include "AssemblerLib/LocalAssemblerBuilder.h"
#include "AssemblerLib/VectorMatrixAssembler.h"
#include "AssemblerLib/LocalDataInitializer.h"
#include "AssemblerLib/LocalToGlobalIndexMap.h"
#include "AssemblerLib/ComputeSparsityPattern.h"

#include "FileIO/VtkIO/VtuInterface.h"

#include "MathLib/LinAlg/ApplyKnownSolution.h"
#include "MathLib/LinAlg/SetMatrixSparsity.h"

#include "MeshLib/MeshSubset.h"
#include "MeshLib/MeshSubsets.h"
#include "MeshGeoToolsLib/MeshNodeSearcher.h"

#include "UniformDirichletBoundaryCondition.h"

#include "SmallDeformationFEM.h"
#include "NeumannBcAssembler.h"
#include "NeumannBc.h"
#include "Parameter.h"
#include "Process.h"
#include "ProcessVariable.h"

namespace MeshLib
{
class Element;
class Mesh;
template <typename PROP_VAL_TYPE>
class PropertyVector;
}

namespace ProcessLib
{
template <typename GlobalSetup, int DisplacementDim>
class SmallDeformationProcess : public Process<GlobalSetup>
{
public:
	using Base = Process<GlobalSetup>;
	using GlobalMatrix = typename GlobalSetup::MatrixType;
	using GlobalVector = typename GlobalSetup::VectorType;

public:
	SmallDeformationProcess(
	    MeshLib::Mesh& mesh,
	    typename Process<GlobalSetup>::NonlinearSolver& nonlinear_solver,
	    std::unique_ptr<typename Process<GlobalSetup>::TimeDiscretization>&&
	        time_discretization,
	    ProcessVariable& variable,
	    Parameter<double, MeshLib::Element const&> const& youngs_modulus,
	    Parameter<double, MeshLib::Element const&> const& poissons_ratio)
	    : Process<GlobalSetup>(mesh, nonlinear_solver,
	                           std::move(time_discretization)),
	      _youngs_modulus(youngs_modulus),
	      _poissons_ratio(poissons_ratio)
	{
		Base::_process_variables.emplace_back(variable);

		if (dynamic_cast<NumLib::ForwardEuler<GlobalVector>*>(
		        &Base::getTimeDiscretization()) != nullptr)
		{
			ERR(
			    "The SmallDeformationProcess can not be solved with the "
			    "ForwardEuler time discretization scheme. Aborting");
			// Because the M matrix is not assembled. Thus, the linearized
			// system would be singular. The same applies to CrankNicolson with
			// theta = 0.0, but this case is not checked here.  Anyway, the
			// SmallDeformationProcess shall be transferred to a simpler
			// ODESystemTag in the future.
			std::abort();
		}
	}

	template <unsigned GlobalDim>
	void createLocalAssemblers()
	{
		DBUG("Create local assemblers.");
		// Populate the vector of local assemblers.
		_local_assemblers.resize(Base::_mesh.getNElements());
		// Shape matrices initializer
		using LocalDataInitializer = AssemblerLib::LocalDataInitializer<
		    SmallDeformation::LocalAssemblerDataInterface,
		    SmallDeformation::LocalAssemblerData,
		    typename GlobalSetup::MatrixType, typename GlobalSetup::VectorType,
		    GlobalDim, Parameter<double, MeshLib::Element const&> const&,
		    Parameter<double, MeshLib::Element const&> const&>;

		LocalDataInitializer initializer;

		using LocalAssemblerBuilder =
		    AssemblerLib::LocalAssemblerBuilder<MeshLib::Element,
		                                        LocalDataInitializer>;

		LocalAssemblerBuilder local_asm_builder(
		    initializer, *Base::_local_to_global_index_map);

		DBUG("Calling local assembler builder for all mesh elements.");
		Base::_global_setup.transform(local_asm_builder,
		                              Base::_mesh.getElements(),
		                              _local_assemblers,
		                              Base::_integration_order,
		                              _youngs_modulus,
		                              _poissons_ratio);
	}

	void createLocalAssemblers() override
	{
		DBUG("Initialize SmallDeformationProcess.");

		if (Base::_mesh.getDimension() == 1)
			createLocalAssemblers<1>();
		else if (Base::_mesh.getDimension() == 2)
			createLocalAssemblers<2>();
		else if (Base::_mesh.getDimension() == 3)
			createLocalAssemblers<3>();
		else
			assert(false);
	}

	bool isLinear() const override { return false; }
private:
	void assembleConcreteProcess(const double t, GlobalVector const& x,
	                             GlobalMatrix& M, GlobalMatrix& K,
	                             GlobalVector& b) override
	{
		DBUG("Assemble SmallDeformationProcess.");

		// Call global assembler for each local assembly item.
		Base::_global_setup.executeDereferenced(
		    *Base::_global_assembler, _local_assemblers, t, x, M, K, b, _dt);
	}

	void assembleJacobianConcreteProcess(
	    const double t, GlobalVector const& x, GlobalVector const& xdot,
	    const double dxdot_dx, GlobalMatrix const& M, const double dx_dx,
	    GlobalMatrix const& K, GlobalMatrix& Jac) override
	{
		DBUG("AssembleJacobian SmallDeformationProcess.");

		(void) t;
		(void) x;
		(void) xdot;
		(void) dxdot_dx;
		(void) M;
		(void) dx_dx;
		(void) K;
		(void) Jac;


		// Call global assembler for each local assembly item.
		//Base::_global_setup.executeDereferenced(*Base::_global_assembler, _local_assemblers,
		 //                           t, x, M, K, b, _dt);
	}

	void preTimestep(GlobalVector const& x, double const /*t*/,
	                 double const dt) override
	{
		DBUG("PreTimestep SmallDeformationProcess.");

		_dt = dt;
		// Call global assembler for each local assembly item.
		Base::_global_setup.executeDereferenced(
		    [&](std::size_t const id, LocalAssembler& local_assembler,
		        GlobalVector const& x) -> void
			{
			    Base::_global_assembler->preTimestep(id, local_assembler, x);
			},
		    _local_assemblers, x);
	}

	void postTimestep(GlobalVector const& x, double const t) override
	{
		DBUG("PostTimestep SmallDeformationProcess.");

		// Call global assembler for each local assembly item.
		Base::_global_setup.executeDereferenced(
		    [&](std::size_t const id, LocalAssembler& local_assembler,
		        GlobalVector const& x, double const t) -> void
			{
			    Base::_global_assembler->postTimestep(id, local_assembler, x, t);
			},
		    _local_assemblers, x, t);
	}

private:
	struct FEMSharedData
	{
	};

	FEMSharedData _fem_shared_data;

	double _dt;
	Parameter<double, MeshLib::Element const&> const& _youngs_modulus;
	Parameter<double, MeshLib::Element const&> const& _poissons_ratio;

	using LocalAssembler = SmallDeformation::LocalAssemblerDataInterface<
	    typename GlobalSetup::MatrixType, typename GlobalSetup::VectorType>;

	std::vector<std::unique_ptr<LocalAssembler>> _local_assemblers;
};

template <typename GlobalSetup, int DisplacementDim>
std::unique_ptr<SmallDeformationProcess<GlobalSetup, DisplacementDim>>
createSmallDeformationProcess(
    MeshLib::Mesh& mesh,
    typename Process<GlobalSetup>::NonlinearSolver& nonlinear_solver,
    std::unique_ptr<typename Process<GlobalSetup>::TimeDiscretization>&&
        time_discretization,
    std::vector<ProcessVariable> const& variables,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config)
{
	config.checkConfParam("type", "SMALL_DEFORMATION");
	DBUG("Create SmallDeformationProcess.");

	// Process variable.
	ProcessVariable& process_variable =
	    findProcessVariable(config, "process_variable", variables);
	DBUG("Associate displacement with process variable \'%s\'.",
	     process_variable.getName().c_str());

	if (process_variable.getNumberOfComponents() != DisplacementDim)
	{
		ERR("Wrong number of components");
		std::abort();
	}

	// Young's modulus parameter.
	auto& youngs_modulus = findParameter<double, MeshLib::Element const&>(
	    config, "youngs_modulus", parameters);

	DBUG("Use \'%s\' as youngs modulus  parameter.",
	     youngs_modulus.name.c_str());

	// Poisson's ratio parameter.
	auto& poissons_ratio = findParameter<double, MeshLib::Element const&>(
	    config, "poissons_ratio", parameters);

	DBUG("Use \'%s\' as poissons ratio parameter.",
	     poissons_ratio.name.c_str());

	return std::unique_ptr<
	    SmallDeformationProcess<GlobalSetup, DisplacementDim>>{
	    new SmallDeformationProcess<GlobalSetup, DisplacementDim>{
	        mesh, nonlinear_solver, std::move(time_discretization),
	        process_variable, youngs_modulus, poissons_ratio}};
}
}  // namespace ProcessLib

#endif  // PROCESS_LIB_SMALLDEFORMATIONPROCESS_H_
