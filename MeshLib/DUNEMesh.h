#pragma once
#include "BaseLib/DUNEConfig.h"

#if OGS_USE_DUNE

#include <functional>
#include <map>
#include <memory>
#include <unordered_map>

#include <dune/grid/uggrid.hh>

#include "FEMMesh.h"

namespace MeshLib
{
struct DUNEIdToIdxMappings
{
    std::unordered_map<std::size_t, std::size_t> map_elements;
    std::unordered_map<std::size_t, std::size_t> map_nodes;
};

template <int GlobalDim>
class DUNEMesh : public FEMMesh
{
public:
    explicit DUNEMesh(std::string const& name,
                      std::unique_ptr<BaseLib::DUNEGridType<GlobalDim>>&& mesh,
                      unsigned global_refinement);

    // only return const ref. grid shouldn't be changed directly
    BaseLib::DUNEGridType<GlobalDim> const& getMesh() const;

    std::size_t getID() const override { return 0 /* TODO [DUNE] fixme */; }

    std::string const& getName() const override { return _name; }

    unsigned getDimension() const override;

    /// Trigger the global mesh refinement.
    void globalRefine() override;

    /// Mark the \c element for refinement or coarsening.
    void mark(char indicator,
              typename BaseLib::DUNEGridType<GlobalDim>::Traits::template Codim<
                  0>::Entity const& element);

    /// Adapt the mesh locally.
    /// \pre preAdapt() must have been called before.
    void adapt();

    /// Declare the beginning of the mesh adaptation process.
    /// \pre All elements that should be adapted have to be mark()'ed before.
    /// \return Whether there are any elements marked for adaptation.
    bool preAdapt();

    /// Declare the end of the adaptation process. This methods cleans up the
    /// isNew status of the mesh elements. \pre adapt() must have been called
    /// before. \pre All post-adaptation data interpolation tasks must have been
    /// finished before.
    void postAdapt();

    /// Register a callback that is called prior to mesh adaptation.
    /// \param on_refine the callback function
    /// \param on_destroy a cleanup callback invoked if this mesh is destroyed.
    /// \return A key that can be used to unregister the callback again.
    std::size_t onPreRefine(
        std::function<void(DUNEMesh<GlobalDim> const&,
                           bool,
                           MeshLib::DUNEIdToIdxMappings const&)>
            on_refine,
        std::function<void()>
            on_destroy);

    /// Register a callback that is called after mesh adaptation.
    /// \param on_refine the callback function
    /// \param on_destroy a cleanup callback invoked if this mesh is destroyed.
    /// \return A key that can be used to unregister the callback again.
    std::size_t onPostRefine(
        std::function<void(DUNEMesh<GlobalDim> const&,
                           bool,
                           MeshLib::DUNEIdToIdxMappings const&)>
            on_refine,
        std::function<void()>
            on_destroy);

    /// Unregister the pre-adaptation callback identified by \c key.
    /// \see onPreRefine()
    void unPreRefine(std::size_t key);

    /// Unregister the post-adaptation callback identified by \c key.
    /// \see onPreRefine()
    void unPostRefine(std::size_t key);

    ~DUNEMesh();

private:
    std::unique_ptr<BaseLib::DUNEGridType<GlobalDim>> _mesh;

    /// The number of refinement levels initially applied globally to the mesh.
    unsigned const _global_refinement;

    std::string const _name;

    struct Callbacks
    {
        std::function<void(DUNEMesh<GlobalDim> const&,
                           bool,
                           MeshLib::DUNEIdToIdxMappings const&)>
            on_refine;

        // this function is expected to unregister itself from the DUNEMesh
        std::function<void()> on_destroy;
    };

    /// Pre-adaptation callbacks registered with this mesh.
    std::map<std::size_t, Callbacks> _preRefine;
    /// Pre-adaptation callbacks registered with this mesh.
    std::map<std::size_t, Callbacks> _postRefine;
    /// The next value for the key identifying pre- and post-adaptation
    /// callbacks.
    std::size_t _key = 1;
};

}  // namespace MeshLib
#endif  // OGS_USE_DUNE
