#include "EigenVector.h"

#define OGS_EIGENVECTOR_OUTPUT_PRE_POST_REFINE 0

#if OGS_EIGENVECTOR_OUTPUT_PRE_POST_REFINE
// TODO [DUNE] debug
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#endif
#include <sstream>

#include "BaseLib/Error.h"
#include "MeshLib/DUNEMesh.h"
#include "MeshLib/FEMMesh.h"
#include "NumLib/DOF/AbstractDOFTable.h"
#include "NumLib/DOF/refinementUtils.h"

/// Unregister pre and post refinement hooks with the mesh.
static void eigenVectorCleanupHooks(MeshLib::FEMMesh* mesh,
                                    std::size_t& key_pre, std::size_t& key_post)
{
    if (auto* mesh_ = dynamic_cast<MeshLib::DUNEMesh<2>*>(mesh))
    {
        mesh_->unPreRefine(key_pre);
        mesh_->unPostRefine(key_post);
    }
    else if (auto* mesh_ = dynamic_cast<MeshLib::DUNEMesh<3>*>(mesh))
    {
        mesh_->unPreRefine(key_pre);
        mesh_->unPostRefine(key_post);
    }

    key_pre = 0;
    key_post = 0;
}

/// Register pre and post refinement hooks with the mesh.
static void eigenVectorSetHooks(MathLib::EigenVector& vec,
                                NumLib::AbstractDOFTable const& dof_table,
                                std::size_t& key_pre,
                                std::size_t& key_post)
{
#if !OGS_EIGENVECTOR_OUTPUT_PRE_POST_REFINE
    (void)dof_table;
#endif

    auto* mesh = vec.getMesh();

    auto on_destroy = [&vec]() { vec.setMesh(nullptr); };

    auto post = [&vec
#if OGS_EIGENVECTOR_OUTPUT_PRE_POST_REFINE
                 ,
                 &dof_table
#endif
    ](auto& mesh, bool, MeshLib::DUNEIdToIdxMappings const& map_id_to_idx) {
        OGS_ALWAYS_ASSERT(&mesh == vec.getMesh());
        NumLib::adaptNodalFieldInplace(vec, map_id_to_idx);

#if OGS_EIGENVECTOR_OUTPUT_PRE_POST_REFINE
        auto const gridView = mesh.getMesh().leafGridView();
        Dune::VTKWriter<decltype(gridView)> vtkWriter(gridView);
        vtkWriter.addVertexData(vec.getRawVector(), "post",
                                dof_table.getNumberOfComponents());

        // TODO [DUNE] debug output
        std::stringstream name;
        name << "refine_" << ++eigenVectorDebugOutputCounterPost << '_' << &vec
             << "_post";
        vtkWriter.write(name.str(), Dune::VTK::OutputType::base64);
#endif
    };

    auto pre = [&vec
#if OGS_EIGENVECTOR_OUTPUT_PRE_POST_REFINE
                ,
                &dof_table
#endif
    ](auto& mesh, bool, MeshLib::DUNEIdToIdxMappings const& /*map_id_to_idx*/) {
#if OGS_EIGENVECTOR_OUTPUT_PRE_POST_REFINE
        auto const gridView = mesh.getMesh().leafGridView();
        Dune::VTKWriter<decltype(gridView)> vtkWriter(gridView);
        vtkWriter.addVertexData(vec.getRawVector(), "pre",
                                dof_table.getNumberOfComponents());

        // TODO [DUNE] debug output
        std::stringstream name;
        name << "refine_" << ++eigenVectorDebugOutputCounterPre << '_' << &vec
             << "_pre";
        vtkWriter.write(name.str(), Dune::VTK::OutputType::base64);
#else
        (void)mesh;
#endif
    };

    if (auto* mesh_ = dynamic_cast<MeshLib::DUNEMesh<2>*>(mesh))
    {
        key_pre = mesh_->onPreRefine(pre, on_destroy);
        key_post = mesh_->onPostRefine(post, on_destroy);
    }
    else if (auto* mesh_ = dynamic_cast<MeshLib::DUNEMesh<3>*>(mesh))
    {
        key_pre = mesh_->onPreRefine(pre, on_destroy);
        key_post = mesh_->onPostRefine(post, on_destroy);
    }
    else
    {
        key_pre = 0;
        key_post = 0;
    }
}

namespace MathLib
{
void EigenVector::setMesh(MeshLib::FEMMesh* mesh)
{
    eigenVectorCleanupHooks(_mesh, _key_pre, _key_post);

    _mesh = mesh;

    eigenVectorSetHooks(*this, *_dof_table, _key_pre, _key_post);
}

MeshLib::FEMMesh* EigenVector::getMesh() const
{
    return _mesh;
}

void EigenVector::setDOFTable(NumLib::AbstractDOFTable const* dof_table)
{
    _dof_table = dof_table;
}

NumLib::AbstractDOFTable const* EigenVector::getDOFTable() const
{
    return _dof_table;
}

EigenVector::~EigenVector()
{
    eigenVectorCleanupHooks(_mesh, _key_pre, _key_post);
}

EigenVector::EigenVector(EigenVector const& other)
    : _vec(other._vec), _mesh(other._mesh), _dof_table(other._dof_table)
{
    eigenVectorSetHooks(*this, *_dof_table, _key_pre, _key_post);
}

EigenVector& EigenVector::operator=(EigenVector const& other)
{
    eigenVectorCleanupHooks(_mesh, _key_pre, _key_post);

    _vec = other._vec;
    _mesh = other._mesh;
    _dof_table = other._dof_table;

    eigenVectorSetHooks(*this, *_dof_table, _key_pre, _key_post);

    return *this;
}

}  // namespace MathLib
