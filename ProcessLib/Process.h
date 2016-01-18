/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_PROCESS_H_
#define PROCESS_LIB_PROCESS_H_

#include <string>

#include "AssemblerLib/ComputeSparsityPattern.h"
#include "AssemblerLib/VectorMatrixAssembler.h"
#include "AssemblerLib/LocalToGlobalIndexMap.h"
#include "BaseLib/ConfigTreeNew.h"
#include "MathLib/LinAlg/SetMatrixSparsity.h"
#include "MeshLib/MeshSubsets.h"

#ifdef USE_PETSC
#include "MeshLib/NodePartitionedMesh.h"
#include "MathLib/LinAlg/PETSc/PETScMatrixOption.h"
#endif

#include "ProcessVariable.h"

namespace MeshLib
{
class Mesh;
}

namespace ProcessLib
{
template <typename GlobalSetup>
class Process
{
public:
	Process(MeshLib::Mesh& mesh) : _mesh(mesh) {}
	virtual ~Process()
	{
		for (auto p : _all_mesh_subsets)
			delete p;
	}

	/// Process specific initialization called by initialize().
	virtual void init() = 0;
	virtual bool assemble(const double delta_t) = 0;

	virtual std::string getLinearSolverName() const = 0;

	/// Postprocessing after solve().
	/// The file_name is indicating the name of possible output file.
	virtual void post(std::string const& file_name) = 0;
	virtual void postTimestep(std::string const& file_name,
	                          const unsigned timestep) = 0;

	/// Creates mesh subsets, i.e. components, for given mesh.
	virtual void initializeMeshSubsets(MeshLib::Mesh const& mesh) = 0;

	void initialize()
	{
		DBUG("Initialize process.");

		DBUG("Construct dof mappings.");
		initializeMeshSubsets(_mesh);

		_local_to_global_index_map.reset(
		    new AssemblerLib::LocalToGlobalIndexMap(
		        _all_mesh_subsets, AssemblerLib::ComponentOrder::BY_COMPONENT));

#ifndef USE_PETSC
		DBUG("Compute sparsity pattern");
		computeSparsityPattern();
#endif

		// create global vectors and linear solver
		createLinearSolver(getLinearSolverName());

		DBUG("Create global assembler.");
		_global_assembler.reset(
		    new GlobalAssembler(*_A, *_rhs, *_local_to_global_index_map));

		init();  // Execute proces specific initialization.
	}

	void initializeNeumannBcs(std::vector<NeumannBc<GlobalSetup>*> const& bcs)
	{
		for (auto bc : bcs)
			bc->initialize(_global_setup, *_A, *_rhs, _mesh.getDimension());
	}

	bool solve(const double delta_t)
	{
		_A->setZero();
		MathLib::setMatrixSparsity(*_A, _sparsity_pattern);

		bool const result = assemble(delta_t);

		_linear_solver->solve(*_rhs, *_x);
		return result;
	}

protected:
	/// Set linear solver options; called by the derived process which is
	/// parsing the configuration.
	void setLinearSolverOptions(BaseLib::ConfigTreeNew&& config)
	{
		_linear_solver_options.reset(new BaseLib::ConfigTreeNew(
			std::move(config)));
	}

	/// Sets the initial condition values in the solution vector x for a given
	/// process variable and component.
	void setInitialConditions(ProcessVariable const& variable,
	                          int const component_id)
	{
		std::size_t const n = _mesh.getNNodes();
		for (std::size_t i = 0; i < n; ++i)
		{
			MeshLib::Location const l(_mesh.getID(),
			                          MeshLib::MeshItemType::Node, i);
			auto const global_index = std::abs(
			    _local_to_global_index_map->getGlobalIndex(l, component_id));
			_x->set(global_index,
			        variable.getInitialConditionValue(*_mesh.getNode(i)));
		}
	}

private:
	/// Creates global matrix, rhs and solution vectors, and the linear solver.
	void createLinearSolver(std::string const& solver_name)
	{
		DBUG("Allocate global matrix, vectors, and linear solver.");
#ifdef USE_PETSC
		MathLib::PETScMatrixOption mat_opt;
		const MeshLib::NodePartitionedMesh& pmesh =
		    static_cast<const MeshLib::NodePartitionedMesh&>(_mesh);
		mat_opt.d_nz = pmesh.getMaximumNConnectedNodesToNode();
		mat_opt.o_nz = mat_opt.d_nz;
		const std::size_t num_unknowns =
		    _local_to_global_index_map->dofSizeGlobal();
		_A.reset(_global_setup.createMatrix(num_unknowns, mat_opt));
#else
		const std::size_t num_unknowns = _local_to_global_index_map->dofSize();
		_A.reset(_global_setup.createMatrix(num_unknowns));
#endif
		_x.reset(_global_setup.createVector(num_unknowns));
		_rhs.reset(_global_setup.createVector(num_unknowns));
		_linear_solver.reset(new typename GlobalSetup::LinearSolver(
		    *_A, solver_name, _linear_solver_options.get()));
		_linear_solver_options->checkAndInvalidate();
	}

	/// Computes and stores global matrix' sparsity pattern from given
	/// DOF-table.
	void computeSparsityPattern()
	{
		_sparsity_pattern = std::move(AssemblerLib::computeSparsityPattern(
		    *_local_to_global_index_map, _mesh));
	}

protected:
	MeshLib::Mesh& _mesh;
	std::vector<MeshLib::MeshSubsets*> _all_mesh_subsets;

	GlobalSetup _global_setup;

	using GlobalAssembler =
	    AssemblerLib::VectorMatrixAssembler<typename GlobalSetup::MatrixType,
	                                        typename GlobalSetup::VectorType>;

	std::unique_ptr<GlobalAssembler> _global_assembler;

	std::unique_ptr<AssemblerLib::LocalToGlobalIndexMap>
	    _local_to_global_index_map;

	std::unique_ptr<BaseLib::ConfigTreeNew> _linear_solver_options;
	std::unique_ptr<typename GlobalSetup::LinearSolver> _linear_solver;

	std::unique_ptr<typename GlobalSetup::MatrixType> _A;
	std::unique_ptr<typename GlobalSetup::VectorType> _rhs;
	std::unique_ptr<typename GlobalSetup::VectorType> _x;

	AssemblerLib::SparsityPattern _sparsity_pattern;
};

}  // namespace ProcessLib

#endif  // PROCESS_LIB_PROCESS_H_
