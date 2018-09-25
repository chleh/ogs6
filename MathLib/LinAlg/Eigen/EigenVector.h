/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>
#ifndef NDEBUG
#include <fstream>
#include <string>
#endif

#include <Eigen/Eigen>
#include <Eigen/Sparse>

namespace MeshLib
{
class FEMMesh;
}

namespace NumLib
{
class AbstractDOFTable;
}

namespace MathLib
{

/// Global vector based on Eigen vector
class EigenVector final
{
public:
    using RawVectorType = Eigen::VectorXd;

    // The Index type of the Eigen::VectorXd class differs from the
    // Eigen::SparseMatrix<double> index type. Maybe an Eigen::SparseVector is a
    // more appropriate RawVectorType for the global vectors.
    using IndexType = Eigen::SparseMatrix<double>::Index;

    // TODO: preliminary
    EigenVector() = default;

    /// Constructor for initialization of the number of rows
    /// @param length number of rows
    explicit EigenVector(IndexType length) : _vec(length) {}

    /// copy constructor
    EigenVector(EigenVector const& other);
    EigenVector& operator=(EigenVector const& other);

    ~EigenVector();

    /// return a vector length
    IndexType size() const { return static_cast<IndexType>(_vec.size()); }

    /// return a start index of the active data range
    IndexType getRangeBegin() const { return 0;}

    /// return an end index of the active data range
    IndexType getRangeEnd() const { return size(); }

    // TODO preliminary
    void setZero() { _vec.setZero(); }

    /// access entry
    double const & operator[] (IndexType rowId) const { return _vec[rowId]; }
    double& operator[] (IndexType rowId) { return _vec[rowId]; }

    /// get entry
    double get(IndexType rowId) const
    {
        return _vec[rowId];
    }

    /// get entries
    std::vector<double> get(std::vector<IndexType> const& indices) const
    {
        std::vector<double> local_x;
        local_x.reserve(indices.size());

        for (auto i : indices) {
            local_x.emplace_back(_vec[i]);
        }

        return local_x;
    }

    /// set entry
    void set(IndexType rowId, double v)
    {
        _vec[rowId] = v;
    }

    /// add entry
    void add(IndexType rowId, double v)
    {
        _vec[rowId] += v;
    }

    /// add entries
    template<class T_SUBVEC>
    void add(const std::vector<IndexType> &pos, const T_SUBVEC &sub_vec)
    {
        auto const length = pos.size();
        for (std::size_t i=0; i<length; ++i) {
            add(pos[i], sub_vec[i]);
        }
    }

    /// Copy vector values.
    void copyValues(std::vector<double>& u) const
    {
        assert(u.size() == (std::size_t) _vec.size());
        copy_n(_vec.data(), _vec.size(), u.begin());
    }

#ifndef NDEBUG
    /// printout this equation for debugging
    void write (const std::string &filename) const { std::ofstream os(filename); os << _vec; }
#endif

    /// return a raw Eigen vector object
    RawVectorType& getRawVector() {return _vec; }

    /// return a raw Eigen vector object
    const RawVectorType& getRawVector() const {return _vec; }

    /// \defgroup DUNE
    /// The interpretation of the EigenVector has been changed. EigenVector is
    /// now taken to be the coefficient vector of some field that is defined on
    /// some domain.
    ///
    /// @{

    /// Set the mesh to which this EigenVector is associated.
    void setMesh(MeshLib::FEMMesh*);
    /// Return the mesh to which this EigenVector is associated.
    MeshLib::FEMMesh* getMesh() const;

    /// Set the d.o.f. table for accessing the d.o.f.s of thie EigenVector.
    void setDOFTable(NumLib::AbstractDOFTable const*);
    /// Return the d.o.f. table for accessing the d.o.f.s of thie EigenVector.
    NumLib::AbstractDOFTable const* getDOFTable() const;

    /// @}

private:
    RawVectorType _vec;

    /// \addtogroup DUNE
    /// @{

    /// The \c _mesh is the discretized domain that field is defined on. The \c
    /// _mesh is needed to update the EigenVector after the mesh has been
    /// adapted.
    MeshLib::FEMMesh* _mesh = nullptr;

    /// The \c _dof_table  maps mesh entities to indices in the coefficient
    /// vector. The \c _dof_table is needed to update the EigenVector after the
    /// mesh has been adapted.
    NumLib::AbstractDOFTable const* _dof_table = nullptr;

    /// Needed for cleaning up pre refinement hooks.
    std::size_t _key_pre = 0;

    /// Needed for cleaning up pre refinement hooks.
    std::size_t _key_post = 0;

    /// @}
};

} // MathLib
