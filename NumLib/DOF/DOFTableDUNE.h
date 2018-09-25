#pragma once

#include <memory>
#include "AbstractDOFTable.h"

#include <iostream>

namespace NumLib
{
/// A glue class that translates between DUNE's and OGS's d.o.f. table
/// interfaces. Makes DUNE's d.o.f. table look like OGS.
template <typename Basis>
class DOFTableDUNE final : public AbstractDOFTable
{
    using Element = typename Basis::GridView::template Codim<0>::Entity;

public:
    explicit DOFTableDUNE(Basis&& basis) : _basis(std::move(basis))
    {
        std::cout << "basis size " << basis.size() << '\n';
        std::cout << "basis dim  " << basis.dimension() << '\n';
        std::cout << "size of an element: %d." << sizeof(Element) << '\n';

        cacheElements();
    }

    int getNumberOfComponents() const override
    {
        auto const& gridView = _basis.gridView();
        auto const dim = gridView.dimension;
        auto const num_vertices = gridView.size(dim);
        auto const num_dofs = _basis.size();
        assert(num_dofs % num_vertices == 0 &&
               "Currently DOFs must be given at vertivces only and each vertex "
               "must have the same DOFs.");
        return num_dofs / num_vertices;
    }
    int getNumberOfVariables() const override { return 1; }
    int getNumberOfVariableComponents(int /*variable_id*/) const override
    {
        // TODO [DUNE] implement
        return getNumberOfComponents();
    }
    std::size_t dofSizeWithoutGhosts() const override { return _basis.size(); }
    std::size_t size() const override { return _basis.gridView().size(0); }
    RowColumnIndices operator()(std::size_t const mesh_item_id,
                                const int component_id) const override
    {
        auto& e = _elements[mesh_item_id];

        auto localView = _basis.localView();
        auto localIndexSet = _basis.localIndexSet();
        localView.bind(e);
        localIndexSet.bind(localView);

        auto const num_comp = getNumberOfComponents();
        auto const num_local_dof = localView.size();
        assert(num_local_dof % num_comp == 0);
        auto const num_dof = num_local_dof / num_comp;

        // cf. VectorMatrixAssemblerDUNE
        _indices.clear();
        _indices.reserve(num_dof);

        for (auto i = decltype(num_dof){0}; i < num_dof; ++i)
        {
            auto const dof = num_dof * component_id +
                             i;  // local DOFs are ordered by component
            assert(localIndexSet.index(dof).size() == 1);
            _indices.push_back(
                static_cast<GlobalIndexType>(localIndexSet.index(dof)[0]));
        }

        return {_indices, _indices};
    }

    std::vector<GlobalIndexType> const& getGhostIndices() const override
    {
        // TODO [DUNE] add proper implementation.
        static std::vector<GlobalIndexType> ghosts;
        return ghosts;
    }

    /// Update the d.o.f. table, e.g., after mesh adaptation.
    void update(typename Basis::GridView const& gv)
    {
        _basis.update(gv);
        cacheElements();
    }

    Basis const& getBasis() const { return _basis; }

private:
    /// Cache elements for internal use.
    /// This is required to maintain the random access provided by operator().
    void cacheElements()
    {
        auto const gridView = _basis.gridView();
        _elements.resize(gridView.size(0));
        auto const& indexSet = gridView.indexSet();
        for (auto& element : Dune::elements(gridView))
        {
            _elements[indexSet.index(element)] = element;
        }
    }

    Basis _basis;

    /// Cached mesh elements.
    /// \todo This is probably a huge waste of memory.
    std::vector<Element> _elements;

    /// Local cache.
    /// \todo Not the best design.
    mutable std::vector<GlobalIndexType> _indices;
};

/// Helper function that enables template type deduction prior to C++17.
template <typename Basis>
std::unique_ptr<DOFTableDUNE<Basis>> makeDOFTableDUNE(Basis&& basis)
{
    return std::make_unique<DOFTableDUNE<Basis>>(std::move(basis));
}

}  // namespace NumLib
