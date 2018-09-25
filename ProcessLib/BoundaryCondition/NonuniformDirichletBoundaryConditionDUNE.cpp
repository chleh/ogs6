/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "NonuniformDirichletBoundaryConditionDUNE.h"

#include <algorithm>
#include <set>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/grid/uggrid.hh>

#include "BoundaryConditionConfig.h"

#include "BaseLib/Functional.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace
{
template <std::size_t GlobalDim, typename A>
static double distanceSquared(A const& a,
                              std::array<double, GlobalDim> const& b)
{
    double dist2 = 0.0;
    for (std::size_t i = 0; i < GlobalDim; ++i)
    {
        dist2 += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return dist2;
}

struct CompareCoordinatesWithTolerance
{
    CompareCoordinatesWithTolerance(double tol) : tol2(tol * tol) {}

    template <std::size_t GlobalDim>
    bool operator()(std::array<double, GlobalDim> const& a,
                    std::array<double, GlobalDim> const& b) const
    {
        if (distanceSquared(a, b) < tol2)
            return false;
        return a < b;
    }

private:
    const double tol2;
};

struct CompareFieldVector
{
    template <int GlobalDim>
    bool operator()(Dune::FieldVector<double, GlobalDim> const& a,
                    Dune::FieldVector<double, GlobalDim> const& b) const
    {
        return std::lexicographical_compare(a.begin(), a.end(), b.begin(),
                                            b.end());
    }
};
}  // namespace

namespace ProcessLib
{
template <int GlobalDim>
NonuniformDirichletBoundaryConditionDUNE<GlobalDim>::
    NonuniformDirichletBoundaryConditionDUNE(
        const MeshLib::Mesh& boundary_mesh,
        MeshLib::PropertyVector<double> const& values,
        MeshLib::DUNEMesh<GlobalDim> const& bulk_mesh,
        int const variable_id_bulk,
        int const component_id_bulk)
    : _dirichlet_values(values), _component_id_bulk(component_id_bulk)
{
    DBUG("init dirichlet BC");

    OGS_ALWAYS_ASSERT(variable_id_bulk == 0);

    if (component_id_bulk < 0 || component_id_bulk >= GlobalDim)
    {
        OGS_FATAL(
            "Currently this BC is implemented only for a single vector-valued "
            "primary variable with %d components",
            GlobalDim);
    }

    // Initialization happens on the coarsest grid view.
    auto const gridView = bulk_mesh.getMesh().levelGridView(0);

    auto const& indexSet = gridView.indexSet();
    auto const& idSet = bulk_mesh.getMesh().localIdSet();

    double const coord_diff_tol = 1e-6;  // TODO [DUNE] input parameter

#if 0
    {
        int i = 0;
        for (auto& v : Dune::vertices(gridView))
        {
            std::cout << "Mesh vertex it " << i << ", idx " << indexSet.index(v)
                      << " -- " << v.geometry().corner(0) << '\n';
            i++;
        }

        i = 0;
        auto const* orig_ids =
            boundary_mesh.getProperties().getPropertyVector<unsigned long>(
                "OriginalSubsurfaceNodeIDs");
        for (auto* n : boundary_mesh.getNodes())
        {
            auto* cs = n->getCoords();
            std::cout << "BC mesh vertex it " << i << " orig "
                      << orig_ids->getComponent(i, 0) << " -- " << cs[0] << " "
                      << cs[1] << " " << cs[2] << '\n';
            i++;
        }
    }
#endif

    class VertexData
    {
        using Index = typename BaseLib::DUNEGridType<
            GlobalDim>::LeafGridView::IndexSet::IndexType;
        using Id =
            typename BaseLib::DUNEGridType<GlobalDim>::LocalIdSet::IdType;
        using BoundarySegments = std::set<std::size_t>;

    public:
        VertexData(Index idx_, Id id_) : idx(idx_), id(id_) {}

        Index idx;
        Id id;
        BoundarySegments boundary_segment_indices;
    };

    CompareCoordinatesWithTolerance compare(coord_diff_tol);

    // Will have an entry for each vertex at the domain boundary.
    std::map<std::array<double, GlobalDim>, VertexData,
             CompareCoordinatesWithTolerance>
        map_coords_to_vertex_data(compare);

    // Contains the number of corners of each boundary segment.
    std::map<std::size_t, unsigned> map_bdry_seg_idx_to_num_corners;

    // Iterate over all elements at the domain boundary.
    for (auto& element : Dune::elements(gridView))
    {
        // ??? Is there a more efficient way that only iterates over the
        // elements at the boundary?
        if (!element.hasBoundaryIntersections())
            continue;

        // Will be used (a) to check if a vertex is located at the domain
        // boundary and (b) to which boundary segments a vertex belongs.
        std::map<Dune::FieldVector<double, GlobalDim>, std::set<std::size_t>,
                 CompareFieldVector>
            map_element_vertex_coords_to_bdry_seg_idcs;

        // Collect boundary segments and their vertices.
        for (auto& intersection : Dune::intersections(gridView, element))
        {
            if (!intersection.boundary())
                continue;

            // Note: intersection does not provide direct access to its vertices
            // as grid entities, but only to their coordinates.
            // ??? Can we get access to the vertex id and index by using
            // intersection.inside() and the size(i, c, cc) and subEntities(i,
            // c, ii, cc) methods of the ReferenceElement? Or by using
            // intersection.indexInInside() and subIndex() subId()? That would
            // save the loop over the vertices of the element below.
            auto const seg = intersection.boundarySegmentIndex();
            auto const geom = intersection.geometry();
            auto const num_corners = geom.corners();

            {
                auto const it_and_succ =
                    map_bdry_seg_idx_to_num_corners.emplace(seg, num_corners);

                // Initially, each boundary segment must belong to one element
                // only. This Assumption is crucial when removing removing the
                // vertices of the intersection again later on.
                OGS_ALWAYS_ASSERT(it_and_succ.second);
            }

            for (int i = 0; i < num_corners; ++i)
            {
                auto const it_and_succ =
                    map_element_vertex_coords_to_bdry_seg_idcs.emplace(
                        std::piecewise_construct,
                        std::forward_as_tuple(geom.corner(i)),
                        std::forward_as_tuple());

                auto& bdry_seg_idcs = it_and_succ.first->second;
                bdry_seg_idcs.insert(seg);
            }
        }

        // Now map_element_vertex_coords_to_bdry_seg_idcs is established.

        // Save indices, ids and boundary segment indices  of all boundary
        // vertices.
        auto const num_vs = element.subEntities(GlobalDim);
        for (auto i = decltype(num_vs){0}; i < num_vs; ++i)
        {
            auto const vertex = element.template subEntity<GlobalDim>(i);
            auto const corner = vertex.geometry().corner(0);

            auto const it =
                map_element_vertex_coords_to_bdry_seg_idcs.find(corner);
            if (it == map_element_vertex_coords_to_bdry_seg_idcs.end())
            {
                // The current vertex is not located on the domain boundary.
                continue;
            }

            auto& bdry_segs = it->second;

            std::array<double, GlobalDim> vertex_coords;
            std::copy(corner.begin(), corner.end(), vertex_coords.begin());

            auto const id = idSet.id(vertex);
            auto const idx = indexSet.index(vertex);

            auto const it_and_succ = map_coords_to_vertex_data.emplace(
                std::piecewise_construct,
                std::forward_as_tuple(vertex_coords),
                std::forward_as_tuple(idx, id));
            auto& vertex_data = it_and_succ.first->second;

            if (it_and_succ.second)
            {
                vertex_data.boundary_segment_indices = std::move(bdry_segs);
            }
            else
            {
                // vertex already present
                assert(vertex_data.id == id);
                assert(vertex_data.idx == idx);

                // merge boundary segment indices
                vertex_data.boundary_segment_indices.insert(bdry_segs.begin(),
                                                            bdry_segs.end());
            }
        }
    }  // for each element

    // Now map_coords_to_vertex_data and map_bdry_seg_idx_to_num_corners are
    // established.

    _dirichlet_vertex_indices.reserve(boundary_mesh.getNumberOfNodes());

    {
        std::size_t dirichlet_dof_idx = 0;
        auto const& boundary_nodes = boundary_mesh.getNodes();

        // Iterate over Dirichlet DOFs and find the corresponding vertex indices
        // and ids.
        for (auto boundary_nodes_it = boundary_nodes.begin();
             boundary_nodes_it != boundary_nodes.end();
             ++boundary_nodes_it, ++dirichlet_dof_idx)
        {
            auto const* const boundary_coords =
                (*boundary_nodes_it)->getCoords();
            std::array<double, GlobalDim> vertex_coords;
            std::copy(boundary_coords, boundary_coords + GlobalDim,
                      vertex_coords.begin());

            auto const it = map_coords_to_vertex_data.find(vertex_coords);
            OGS_ALWAYS_ASSERT(it != map_coords_to_vertex_data.end());

#if 1
            {
                // Simple check for correspondence between bulk and boundary
                // mesh.
                auto const& bulk_coords = it->first;
                OGS_ALWAYS_ASSERT(
                    distanceSquared(boundary_coords, bulk_coords) <
                    coord_diff_tol * coord_diff_tol);

#if 0
                std::cout << "diri node idx " << it->second.idx << " --";
                for (int i = 0; i < GlobalDim; ++i)
                {
                    std::cout << " " << bulk_coords[i];
                }
                std::cout << " bdry_segs:";
                for (auto s : it->second.boundary_segment_indices)
                {
                    std::cout << " " << s;
                }
                std::cout << '\n';
#endif
            }
#endif

            auto const& vertex_data = it->second;
            _dirichlet_vertex_indices.push_back(vertex_data.idx);
            _map_boundary_vertex_id_to_value.emplace(
                vertex_data.id, _dirichlet_values[dirichlet_dof_idx]);

            for (auto seg : vertex_data.boundary_segment_indices)
            {
                if (--map_bdry_seg_idx_to_num_corners[seg] == 0)
                {
                    // All vertices of the boundary segment belong to this
                    // Dirichlet BC. Therefore the entire boundary segment
                    // belongs to this BC.
                    _boundary_segment_indices.insert(seg);
                }
            }

            map_coords_to_vertex_data.erase(it);
        }
    }

    // TODO [DUNE] fix const cast
    const_cast<MeshLib::DUNEMesh<GlobalDim>&>(bulk_mesh).onPostRefine(
        BaseLib::easyBind(
            &NonuniformDirichletBoundaryConditionDUNE<GlobalDim>::onPostRefine,
            *this),
        nullptr);

    DBUG(
        "Dirichlet BC will be set on %d vertices, which correspond to %d "
        "boundary segments",
        _dirichlet_values.size(), _boundary_segment_indices.size());
}

template <int GlobalDim>
void NonuniformDirichletBoundaryConditionDUNE<GlobalDim>::getEssentialBCValues(
    const double /*t*/, GlobalVector const& /*x*/,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    bc_values.ids.clear();
    bc_values.values.clear();  // test

    // Convert mesh node ids to global index for the given component.
    bc_values.ids.reserve(_dirichlet_values.size());
    bc_values.values.reserve(_dirichlet_values.size());

    // Map boundary dof indices to bulk dof indices and the corresponding
    // values.
    for (std::size_t i = 0; i < _dirichlet_values.size(); ++i)
    {
        auto const bulk_vertex_index = _dirichlet_vertex_indices[i];

        // TODO [DUNE] assumes one variable with GlobalDim components. And
        // global order by location
        auto const global_index =
            _component_id_bulk + bulk_vertex_index * GlobalDim;

        bc_values.ids.push_back(global_index);
        bc_values.values.push_back(_dirichlet_values[i]);
    }
}

template <int GlobalDim>
void NonuniformDirichletBoundaryConditionDUNE<GlobalDim>::onPostRefine(
    MeshLib::DUNEMesh<GlobalDim> const& mesh, bool globally_refined,
    MeshLib::DUNEIdToIdxMappings const&)
{
    auto const gridView = mesh.getMesh().leafGridView();
    auto const& idSet = mesh.getMesh().localIdSet();
    auto const& indexSet = gridView.indexSet();

    std::map<typename BaseLib::DUNEGridType<GlobalDim>::LocalIdSet::IdType,
             std::pair<typename BaseLib::DUNEGridType<
                           GlobalDim>::LeafGridView::IndexSet::IndexType,
                       double>>
        map_dirichlet_vertex_id_to_idx_and_value;

    // TODO [DUNE] here we have to use a FE basis that is compatible with the
    // one used in the process.
    auto const linear_lagrange_basis = Dune::Functions::BasisBuilder::makeBasis(
        gridView, Dune::Functions::BasisBuilder::lagrange<1>());
    auto localView = linear_lagrange_basis.localView();

    std::vector<double> nodal_values;  // temporary storage
    std::vector<Dune::FieldVector<double, 1>>
        shape_functions;  // temporary storage

    // Iterate over all elements at the domain boundary.
    for (auto& element_at_boundary : Dune::elements(gridView))
    {
        if (!element_at_boundary.hasBoundaryIntersections())
            continue;

        if (globally_refined || element_at_boundary.isNew())
        {
            std::set<Dune::FieldVector<double, GlobalDim>, CompareFieldVector>
                element_candidate_dirichlet_vertices;

            // Collect all vertices from all faces of element_at_boundary,
            // which belong to this Dirichlet BC.
            for (auto& intersection :
                 Dune::intersections(gridView, element_at_boundary))
            {
                if (!intersection.boundary() ||
                    _boundary_segment_indices.find(
                        intersection.boundarySegmentIndex()) ==
                        _boundary_segment_indices.end())
                {
                    // This intersection does not belong to this Dirichlet BC.
                    continue;
                }

                auto const geom = intersection.geometry();

                for (int i = 0; i < geom.corners(); ++i)
                {
                    element_candidate_dirichlet_vertices.insert(geom.corner(i));
                }
            }

            // Iterate over element_at_boundary's vertices.
            for (std::size_t k = 0;
                 k < element_at_boundary.subEntities(GlobalDim);
                 ++k)
            {
                auto const corner = element_at_boundary.geometry().corner(k);
                if (element_candidate_dirichlet_vertices.find(corner) ==
                    element_candidate_dirichlet_vertices.end())
                {
                    // This vertex does not belong to this Dirichlet BC.
                    continue;
                }

                auto father = element_at_boundary;
                auto positionInFather =
                    Dune::ReferenceElements<double, GlobalDim>::general(
                        element_at_boundary.type())
                        .position(k, GlobalDim);

                do
                {
                    positionInFather =
                        father.geometryInFather().global(positionInFather);
                    father = father.father();
                } while ((globally_refined && father.level() != 0) ||
                         father.isNew());  // TODO [DUNE] maybe always go back
                                           // to level 0?

                nodal_values.clear();

                // Extract corner values for subsequent interpolation.
                for (std::size_t l = 0; l < father.subEntities(GlobalDim); ++l)
                {
                    auto const it = _map_boundary_vertex_id_to_value.find(
                        idSet.subId(father, l, GlobalDim));
                    if (it != _map_boundary_vertex_id_to_value.end())
                    {
                        // There is a Dirichlet DOF for vertex l.
                        nodal_values.push_back(it->second);
                    }
                    else
                    {
                        // There is no Dirichlet DOF for vertex l.
                        // We use 0 as a dummy value. It is assumed that the
                        // actual value does not influence the interpolation
                        // procedure on father.
                        nodal_values.push_back(0);
                    }
                }

                // Interpolate using shape functions.
                localView.bind(father);
                auto const& localFiniteElement =
                    localView.tree().finiteElement();
                localFiniteElement.localBasis().evaluateFunction(
                    positionInFather, shape_functions);

                assert(shape_functions.size() == nodal_values.size());

                double dirichlet_value = 0;
                for (unsigned i = 0; i < shape_functions.size(); ++i)
                {
                    dirichlet_value += shape_functions[i] * nodal_values[i];
                }

                auto const id = idSet.subId(element_at_boundary, k, GlobalDim);
                auto const idx =
                    indexSet.subIndex(element_at_boundary, k, GlobalDim);

#if 0
                std::cout << "Dirichlet DOF id " << id << " idx " << idx
                          << " value " << dirichlet_value << " -- " << corner
                          << '\n';
#endif

                map_dirichlet_vertex_id_to_idx_and_value.emplace(
                    std::piecewise_construct,
                    std::forward_as_tuple(id),
                    std::forward_as_tuple(idx, dirichlet_value));
            }
        }
        else
        {
            // Not a refined element.
            for (unsigned k = 0; k < element_at_boundary.subEntities(GlobalDim);
                 ++k)
            {
                auto const id = idSet.subId(element_at_boundary, k, GlobalDim);
                auto const idx =
                    indexSet.subIndex(element_at_boundary, k, GlobalDim);
                auto const it = _map_boundary_vertex_id_to_value.find(id);

                if (it != _map_boundary_vertex_id_to_value.end())
                {
                    map_dirichlet_vertex_id_to_idx_and_value.emplace(
                        std::piecewise_construct,
                        std::forward_as_tuple(id),
                        std::forward_as_tuple(idx, it->second));
                }
                // TODO [DUNE] else fail?
            }
        }
    }  // for each element

    DBUG("Dirichlet BC pre refine %d d.o.f.", _dirichlet_vertex_indices.size());

    _map_boundary_vertex_id_to_value.clear();
    _dirichlet_vertex_indices.clear();
    _dirichlet_vertex_indices.reserve(
        map_dirichlet_vertex_id_to_idx_and_value.size());
    _dirichlet_values.clear();
    _dirichlet_values.reserve(map_dirichlet_vertex_id_to_idx_and_value.size());

    for (auto& id_and_idx_and_value : map_dirichlet_vertex_id_to_idx_and_value)
    {
        auto const id = id_and_idx_and_value.first;
        auto const idx = id_and_idx_and_value.second.first;
        auto const value = id_and_idx_and_value.second.second;

        _map_boundary_vertex_id_to_value.emplace(id, value);
        _dirichlet_vertex_indices.push_back(idx);
        _dirichlet_values.push_back(value);
    }

    DBUG("Dirichlet BC post refine %d d.o.f.",
         _dirichlet_vertex_indices.size());
}

template <int GlobalDim>
std::unique_ptr<NonuniformDirichletBoundaryConditionDUNE<GlobalDim>>
createNonuniformDirichletBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& boundary_mesh,
    int const variable_id, int const component_id,
    const MeshLib::DUNEMesh<GlobalDim>& bulk_mesh)
{
    DBUG("Constructing NonuniformDirichletDUNE BC from config.");
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "NonuniformDirichlet");

    // TODO [DUNE] finally use ProcessLib::Parameter here
    auto const field_name =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__NonuniformDirichletDUNE__field_name}
        config.getConfigParameter<std::string>("field_name");

    auto const* const property =
        boundary_mesh.getProperties().getPropertyVector<double>(field_name);

    if (!property)
    {
        OGS_FATAL("A property with name `%s' does not exist in `%s'.",
                  field_name.c_str(), boundary_mesh.getName().c_str());
    }

    if (property->getMeshItemType() != MeshLib::MeshItemType::Node)
    {
        OGS_FATAL(
            "Only nodal fields are supported for non-uniform fields. Field "
            "`%s' is not nodal.",
            field_name.c_str());
    }

    if (property->getNumberOfComponents() != 1)
    {
        OGS_FATAL("`%s' is not a one-component field.", field_name.c_str());
    }

    return std::make_unique<
        NonuniformDirichletBoundaryConditionDUNE<GlobalDim>>(
        boundary_mesh, *property, bulk_mesh, variable_id, component_id);
}

std::unique_ptr<BoundaryCondition> createNonuniformDirichletBoundaryCondition(
    const BoundaryConditionConfig& config,
    const NumLib::AbstractDOFTable& /*dof_table*/,
    const MeshLib::FEMMesh& bulk_mesh, const int variable_id)
{
    auto* boundary_mesh =
        dynamic_cast<MeshLib::Mesh const*>(&config.boundary_mesh);
    if (!boundary_mesh)
    {
        OGS_FATAL("Wrong type!");
    }

#if 0
    if (auto* mesh_ = dynamic_cast<MeshLib::Mesh const*>(&mesh))
    {
        return createNonuniformDirichletBoundaryCondition(
            config.config, dof_table, variable_id, *config.component_id,
            *mesh_);
    }
#endif
#if 0
    if (auto* mesh_ = dynamic_cast<MeshLib::DUNEMesh<1> const*>(&mesh))
    {
        return createNonuniformDirichletBoundaryCondition(
            config.config, variable_id, *config.component_id, *mesh_);
    }
#endif
    if (auto* bulk_mesh_ =
            dynamic_cast<MeshLib::DUNEMesh<2> const*>(&bulk_mesh))
    {
        return createNonuniformDirichletBoundaryCondition(
            config.config, *boundary_mesh, variable_id, *config.component_id,
            *bulk_mesh_);
    }
    if (auto* bulk_mesh_ =
            dynamic_cast<MeshLib::DUNEMesh<3> const*>(&bulk_mesh))
    {
        return createNonuniformDirichletBoundaryCondition(
            config.config, *boundary_mesh, variable_id, *config.component_id,
            *bulk_mesh_);
    }

    OGS_FATAL("unsupported mesh");
}

}  // namespace ProcessLib
