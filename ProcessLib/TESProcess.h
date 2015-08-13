/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_TESPROCESS_H_
#define PROCESS_LIB_TESPROCESS_H_

#include <memory>
#include <vector>

#include "AssemblerLib/LocalToGlobalIndexMap.h"
#include "AssemblerLib/VectorMatrixAssembler.h"

#include "FileIO/VtkIO/VtuInterface.h"

#include "MathLib/LinAlg/ApplyKnownSolution.h"
#include "MathLib/LinAlg/SetMatrixSparsity.h"

#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "MeshLib/NodeAdjacencyTable.h"

#include "ProcessVariable.h"
#include "Process.h"

#include "TESFEM.h"

namespace MeshLib
{
    class Element;
    class Mesh;
    template <typename PROP_VAL_TYPE> class PropertyVector;
}

namespace ProcessLib
{

namespace TES
{

const unsigned NODAL_DOF = 3;

template<typename GlobalSetup>
class TESProcess : public Process
{
    using ConfigTree = boost::property_tree::ptree;

    unsigned const _integration_order = 2;

public:
    TESProcess(MeshLib::Mesh& mesh,
            std::vector<ProcessVariable> const& variables,
            ConfigTree const& config);

    template <unsigned GlobalDim>
    void createLocalAssemblers();

    void initialize();

    void solve();

    void post(std::string const& file_name);

    ~TESProcess();

private:
    ProcessVariable* _process_vars[NODAL_DOF] = { nullptr, nullptr, nullptr };

    MeshLib::MeshSubset const* _mesh_subset_all_nodes = nullptr;
    std::vector<MeshLib::MeshSubsets*> _all_mesh_subsets;

    GlobalSetup _global_setup;
    std::unique_ptr<typename GlobalSetup::MatrixType> _A;
    std::unique_ptr<typename GlobalSetup::VectorType> _rhs;
    std::unique_ptr<typename GlobalSetup::VectorType> _x;

    using LocalAssembler = TES::LocalAssemblerDataInterface<
        typename GlobalSetup::MatrixType, typename GlobalSetup::VectorType>;

    std::vector<LocalAssembler*> _local_assemblers;

    using GlobalAssembler =
        AssemblerLib::VectorMatrixAssembler<
            typename GlobalSetup::MatrixType,
            typename GlobalSetup::VectorType>;


    std::unique_ptr<AssemblerLib::LocalToGlobalIndexMap> _local_to_global_index_map;

    std::unique_ptr<GlobalAssembler> _global_assembler;

    /// Global ids in the global matrix/vector where the dirichlet bc is
    /// imposed and their corresponding values.
    struct DirichletBC {
        std::vector<std::size_t> global_ids;
        std::vector<double> values;
    } _dirichlet_bc;

    std::vector<NeumannBc<GlobalSetup>*> _neumann_bcs;

    MeshLib::NodeAdjacencyTable _node_adjacency_table;
};

} // namespace TES

} // namespace ProcessLib

#include "TESProcess-impl.h"

#endif  // PROCESS_LIB_TESPROCESS_H_
