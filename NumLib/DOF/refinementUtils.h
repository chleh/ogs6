#pragma once

#include <unordered_map>

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/uggrid.hh>

#include "BaseLib/Error.h"
#include "MathLib/LinAlg/Eigen/EigenVector.h"
#include "MeshLib/DUNEMesh.h"

#include "DOFTableDUNE.h"

namespace NumLib
{
namespace detail
{
// TODO [DUNE] duplicated from SmallDeformationProcess-impl.h
template <int DisplacementDim>
decltype(auto) makeDisplacementBasis(
    typename BaseLib::DUNEGridType<DisplacementDim>::LeafGridView const&
        gridView)
{
    namespace DFB = Dune::Functions::BasisBuilder;
    return DFB::makeBasis(gridView,
                          DFB::power<DisplacementDim>(DFB::lagrange<1>(),
                                                      DFB::flatInterleaved()));
}

template <int DisplacementDim>
decltype(auto) makeScalarBasis(
    typename BaseLib::DUNEGridType<DisplacementDim>::LeafGridView const&
        gridView)
{
    namespace DFB = Dune::Functions::BasisBuilder;
    return DFB::makeBasis(gridView, DFB::lagrange<1>());
}

template <int GlobalDim, typename Basis>
void adaptNodalFieldInplace(MathLib::EigenVector& vec,
                            MeshLib::DUNEIdToIdxMappings const& map_id_to_idx,
                            MeshLib::DUNEMesh<GlobalDim> const& mesh,
                            DOFTableDUNE<Basis> const& dof_table)
{
    // TODO [DUNE] document vector are nodal d.o.f.
    auto const& grid = mesh.getMesh();
    auto const& gridView = grid.leafGridView();
    auto const& indexSet = gridView.indexSet();
    auto const& idSet = grid.localIdSet();

    auto& x_old = vec.getRawVector();
    auto x_new = Eigen::VectorXd(dof_table.dofSizeWithoutGhosts());

    std::vector<double> nodal_values;
    std::vector<Dune::FieldVector<double, 1>>
        shape_functions;  // temporary storage

    auto const& basis = dof_table.getBasis();
    auto localView = basis.localView();

    // TODO [DUNE] template param?
    auto const num_components = dof_table.getNumberOfComponents();

    for (auto const& element : Dune::elements(gridView))
    {
        if (element.isNew())
        {
            for (unsigned k = 0; k < element.subEntities(GlobalDim); k++)
            {
                auto father = element;
                auto positionInFather =
                    Dune::ReferenceElements<double, GlobalDim>::general(
                        element.type())
                        .position(k, GlobalDim);

                do
                {
                    positionInFather =
                        father.geometryInFather().global(positionInFather);
                    father = father.father();
                } while (father.isNew());

                localView.bind(father);

                for (int c = 0; c < num_components; ++c)
                {
                    // extract corner values
                    nodal_values.resize(father.subEntities(GlobalDim));
                    for (unsigned l = 0; l < father.subEntities(GlobalDim); l++)
                    {
                        auto const it = map_id_to_idx.map_nodes.find(
                            idSet.subId(father, l, GlobalDim));
                        assert(it != map_id_to_idx.map_nodes.end());
                        auto const old_idx = it->second;
                        auto const o = c + old_idx * num_components;

                        nodal_values[l] = x_old[o];
                    }

                    auto const& localFiniteElement =
                        localView.tree().child(c).finiteElement();
                    localFiniteElement.localBasis().evaluateFunction(
                        positionInFather, shape_functions);
                    assert(shape_functions.size() == nodal_values.size());

                    double interpolated_value = 0;
                    for (unsigned i = 0; i < shape_functions.size(); ++i)
                    {
                        interpolated_value +=
                            shape_functions[i] * nodal_values[i];
                    }

                    // TODO [DUNE] implies by-location DOF order
                    auto const new_idx =
                        indexSet.subIndex(element, k, GlobalDim);
                    auto const n = c + new_idx * num_components;

                    x_new[n] = interpolated_value;
                }
            }
        }
        else
        {
            for (unsigned k = 0; k < element.subEntities(GlobalDim); k++)
            {
                auto const it = map_id_to_idx.map_nodes.find(
                    idSet.subId(element, k, GlobalDim));
                assert(it != map_id_to_idx.map_nodes.end());
                auto const old_idx = it->second;
                auto const new_idx = indexSet.subIndex(element, k, GlobalDim);

                for (int c = 0; c < num_components; ++c)
                {
                    // TODO [DUNE] implies by-location DOF order
                    auto const n = c + new_idx * num_components;
                    auto const o = c + old_idx * num_components;
                    x_new[n] = x_old[o];
                }
            }
        }
    }

    x_old = std::move(x_new);
}

template <int GlobalDim>
void adaptNodalFieldInplace(MathLib::EigenVector& vec,
                            MeshLib::DUNEIdToIdxMappings const& map_id_to_idx,
                            MeshLib::DUNEMesh<GlobalDim> const& mesh)
{
    auto const& grid = mesh.getMesh();
    auto const& gridView = grid.leafGridView();
    auto const* dof_table = vec.getDOFTable();

    OGS_ALWAYS_ASSERT(dof_table != nullptr);

    // TODO [DUNE] very bad design
#if 0
    if (dof_table->getNumberOfComponents() == 1)
    {
        using Basis = decltype(makeScalarBasis<GlobalDim>(gridView));
        auto* dof_table_ = dynamic_cast<DOFTableDUNE<Basis> const*>(dof_table);
        OGS_ALWAYS_ASSERT(dof_table_ != nullptr);
        adaptNodalFieldInplace(vec, map_id_to_idx, mesh, *dof_table_);
        return;
    }
#endif
    if (dof_table->getNumberOfComponents() == GlobalDim)
    {
        using Basis = decltype(makeDisplacementBasis<GlobalDim>(gridView));
        auto* dof_table_ = dynamic_cast<DOFTableDUNE<Basis> const*>(dof_table);
        OGS_ALWAYS_ASSERT(dof_table_ != nullptr);
        adaptNodalFieldInplace(vec, map_id_to_idx, mesh, *dof_table_);
        return;
    }

    OGS_FATAL("unsupported finite element basis");
}

}  // namespace detail

/// Update a vector associated to a mesh after the mesh has been adapted.
/// \param vec The vector to be updated.
/// \param map_id_to_idx A mapping from DUNE's persistent IDs to DUNE's
/// non-persistent indices from right before the mesh adaptation.
void adaptNodalFieldInplace(MathLib::EigenVector& vec,
                            MeshLib::DUNEIdToIdxMappings const& map_id_to_idx)
{
    auto* mesh = vec.getMesh();

    OGS_ALWAYS_ASSERT(mesh != nullptr);

    if (auto* mesh_ = dynamic_cast<MeshLib::DUNEMesh<2>*>(mesh))
    {
        detail::adaptNodalFieldInplace(vec, map_id_to_idx, *mesh_);
        return;
    }
    if (auto* mesh_ = dynamic_cast<MeshLib::DUNEMesh<3>*>(mesh))
    {
        detail::adaptNodalFieldInplace(vec, map_id_to_idx, *mesh_);
        return;
    }

    OGS_FATAL("unsupported mesh type for refinement");
}

}  // namespace NumLib
