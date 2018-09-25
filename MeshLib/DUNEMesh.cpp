#include "DUNEMesh.h"

#include <logog/include/logog.hpp>

#include "BaseLib/Error.h"

namespace MeshLib
{
template <int GlobalDim>
DUNEMesh<GlobalDim>::DUNEMesh(
    const std::string& name,
    std::unique_ptr<BaseLib::DUNEGridType<GlobalDim>>&& mesh,
    unsigned global_refinement)
    : _mesh(std::move(mesh)), _global_refinement(global_refinement), _name(name)
{
    DBUG("DUNE mesh num boundary segments: %d.", _mesh->numBoundarySegments());
}

// only return const ref. grid shouldn't be changed directly
template <int GlobalDim>
BaseLib::DUNEGridType<GlobalDim> const& DUNEMesh<GlobalDim>::getMesh() const
{
    return *_mesh;
}

template <int GlobalDim>
unsigned DUNEMesh<GlobalDim>::getDimension() const
{
    return static_cast<unsigned>(BaseLib::DUNEGridType<GlobalDim>::dimension);
}

template <int GlobalDim>
void DUNEMesh<GlobalDim>::globalRefine()
{
    if (_global_refinement == 0)
        return;

    DBUG("DUNE mesh: globally refining grid");

    auto const gridView = _mesh->leafGridView();
    auto const& indexSet = gridView.indexSet();
    auto const& idSet = _mesh->localIdSet();

    DUNEIdToIdxMappings mappings;

    mappings.map_elements.reserve(gridView.size(0));
    for (auto& e : Dune::elements(gridView))
    {
        mappings.map_elements.emplace(idSet.id(e), indexSet.index(e));
    }

    mappings.map_nodes.reserve(gridView.size(GlobalDim));
    for (auto& v : Dune::vertices(gridView))
    {
        mappings.map_nodes.emplace(idSet.id(v), indexSet.index(v));
    }

    DBUG("DUNE mesh: %d pre refinement callbacks", _preRefine.size());
    for (auto& f : _preRefine)
        f.second.on_refine(*this, true, mappings);

    {
        auto const gridView = _mesh->leafGridView();
        INFO("DUNE mesh: pre refine %d elements, %d vertices",
             gridView.size(0),
             gridView.size(GlobalDim));
    }

    _mesh->globalRefine(static_cast<int>(_global_refinement));

    {
        auto const gridView = _mesh->leafGridView();
        INFO("DUNE mesh: post refine %d elements, %d vertices",
             gridView.size(0),
             gridView.size(GlobalDim));
    }

    DBUG("DUNE mesh: %d post refinement callbacks", _postRefine.size());
    for (auto& f : _postRefine)
        f.second.on_refine(*this, true, mappings);
}

template <int GlobalDim>
void DUNEMesh<GlobalDim>::mark(
    char indicator,
    typename BaseLib::DUNEGridType<GlobalDim>::Traits::template Codim<
        0>::Entity const& element)
{
    _mesh->mark(indicator, element);
}

template <int GlobalDim>
void DUNEMesh<GlobalDim>::adapt()
{
    DBUG("DUNE mesh: refining grid");

    auto const gridView = _mesh->leafGridView();
    auto const& indexSet = gridView.indexSet();
    auto const& idSet = _mesh->localIdSet();

    DUNEIdToIdxMappings mappings;

    mappings.map_elements.reserve(gridView.size(0));
    for (auto& e : Dune::elements(gridView))
    {
        mappings.map_elements.emplace(idSet.id(e), indexSet.index(e));
    }

    mappings.map_nodes.reserve(gridView.size(GlobalDim));
    for (auto& v : Dune::vertices(gridView))
    {
        mappings.map_nodes.emplace(idSet.id(v), indexSet.index(v));
    }

    DBUG("DUNE mesh: %d pre refinement callbacks", _preRefine.size());
    for (auto& f : _preRefine)
        f.second.on_refine(*this, false, mappings);

    {
        auto const gridView = _mesh->leafGridView();
        INFO("DUNE mesh: pre refine %d elements, %d vertices",
             gridView.size(0),
             gridView.size(GlobalDim));
    }

    _mesh->adapt();

    {
        auto const gridView = _mesh->leafGridView();
        INFO("DUNE mesh: post refine %d elements, %d vertices",
             gridView.size(0),
             gridView.size(GlobalDim));
    }

    DBUG("DUNE mesh: %d post refinement callbacks", _postRefine.size());
    for (auto& f : _postRefine)
        f.second.on_refine(*this, false, mappings);
}

template <int GlobalDim>
void DUNEMesh<GlobalDim>::postAdapt()
{
    _mesh->postAdapt();
}

template <int GlobalDim>
bool DUNEMesh<GlobalDim>::preAdapt()
{
    return _mesh->preAdapt();
}

template <int GlobalDim>
std::size_t DUNEMesh<GlobalDim>::onPreRefine(
    std::function<void(
        DUNEMesh<GlobalDim> const&, bool, MeshLib::DUNEIdToIdxMappings const&)>
        on_refine,
    std::function<void()>
        on_destroy)
{
    auto key = _key++;
    _preRefine.emplace(key,
                       Callbacks{std::move(on_refine), std::move(on_destroy)});
    return key;
}

template <int GlobalDim>
std::size_t DUNEMesh<GlobalDim>::onPostRefine(
    std::function<void(
        DUNEMesh<GlobalDim> const&, bool, MeshLib::DUNEIdToIdxMappings const&)>
        on_refine,
    std::function<void()>
        on_destroy)
{
    auto key = _key++;
    _postRefine.emplace(key,
                        Callbacks{std::move(on_refine), std::move(on_destroy)});
    return key;
}

template <int GlobalDim>
void DUNEMesh<GlobalDim>::unPreRefine(std::size_t key)
{
    if (key == 0)
        return;
    auto it = _preRefine.find(key);
    OGS_ALWAYS_ASSERT(it != _preRefine.end());
    _preRefine.erase(it);
}

template <int GlobalDim>
void DUNEMesh<GlobalDim>::unPostRefine(std::size_t key)
{
    if (key == 0)
        return;
    auto it = _postRefine.find(key);
    OGS_ALWAYS_ASSERT(it != _postRefine.end());
    _postRefine.erase(it);
}

template <int GlobalDim>
DUNEMesh<GlobalDim>::~DUNEMesh()
{
    // clean up hooks (otherwise references from the hooks to *this might
    // cause seg faults.
    while (!_preRefine.empty())
    {
        auto it = _preRefine.begin();
        if (it->second.on_destroy)
        {
            auto const old_size = _preRefine.size();
            it->second.on_destroy();
            if (_preRefine.size() >= old_size)
            {
                ERR("Cleaning up _preRefine went wrong: The on_destroy "
                    "hook did not clean up properly.");
                break;
            }
        }
        else
        {
            _preRefine.erase(it);
        }
    }

    while (!_postRefine.empty())
    {
        auto it = _postRefine.begin();
        if (it->second.on_destroy)
        {
            auto const old_size = _postRefine.size();
            it->second.on_destroy();
            if (_postRefine.size() >= old_size)
            {
                ERR("Cleaning up _postRefine went wrong: The on_destroy "
                    "hook did not clean up properly.");
                break;
            }
        }
        else
        {
            _postRefine.erase(it);
        }
    }
}

template class DUNEMesh<2>;
template class DUNEMesh<3>;

}  // namespace MeshLib
