/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Definition of the DenseVector class.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef DENSEVECTOR_H_
#define DENSEVECTOR_H_

#include <vector>
#include <valarray>
#include <fstream>
#include <iterator>

namespace MathLib
{

/**
 * Dense vector class
 */
template <typename T>
class DenseVector : public std::valarray<T>
{
public:
	typedef T FP_T;
	using IndexType = int;  // The type of valarray indices.

public:
	using std::valarray<T>::operator=;
	using std::valarray<T>::operator+=;
	using std::valarray<T>::operator-=;
	using std::valarray<T>::operator[];

	/**
	 * Constructor for initialization of the number of rows
	 * @param nrows number of rows
	 */
	explicit DenseVector(IndexType nrows=0)
	: std::valarray<T>(nrows)
	{}

	void resize(IndexType new_size)
	{
		std::valarray<T>::resize(new_size);
	}

	/// return a start index of the active data range
	IndexType getRangeBegin() const { return 0;}

	/// return an end index of the active data range
	IndexType getRangeEnd() const { return this->size(); }

	/// get entry
	T get(std::size_t i) const { return (*this)[i]; }

	/// set a value to entry
	void set(IndexType i, T v) { (*this)[i] = v; }

	/// add a value to entry
	void add(IndexType i, T v) { (*this)[i] += v; }

	/**
	 * add a sub vector
	 * @param pos       positions of each sub-vector entry in this vector
	 * @param sub_vec   sub-vector
	 */
	template<class T_SUBVEC>
	void add(const std::vector<IndexType> &pos, const T_SUBVEC &sub_vec)
	{
		for (std::size_t i=0; i<pos.size(); ++i) {
			this->add(pos[i], sub_vec[i]);
		}
	}

	void add(const std::vector<IndexType> &pos, const double v)
	{
		for (std::size_t i=0; i<pos.size(); ++i) {
			this->add(pos[i], v);
		}
	}

	/**
	 * writes the matrix entries into a file
	 * @param filename output file name
	 */
	void write (const std::string &filename) const
	{
		std::ofstream os(filename);
		write(os);
	}

	void write(std::ostream& os) const
	{
		os << static_cast<std::valarray<T> >(*this).size() << "\n";
		os << *this << "\n";
	}

    /// vector operation: add
    void operator+= (const DenseVector<T>& v)
    {
        *this += static_cast<std::valarray<T> >(v);
    }

    /// vector operation: subtract
    void operator-= (const DenseVector<T>& v)
    {
        *this -= static_cast<std::valarray<T> >(v);
    }

    void componentwiseDivide(DenseVector<T> const& v)
    {
        *this /= static_cast<std::valarray<T> >(v);
    }

};

/**
 * writes a vector content into the output stream
 * @param os the output stream
 */
template <typename T>
std::ostream& operator<<(std::ostream& os, DenseVector<T> const & v)
{
	std::copy(std::begin(v), std::end(v), std::ostream_iterator<T>(os, "\n"));
	return os;
}

}


#endif /* DENSEVECTOR_H_ */
