/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Definition of the LisVector class.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LISVECTOR_H_
#define LISVECTOR_H_

#include <string>
#include <vector>

#include <lis.h>

namespace MathLib
{
/**
 * \brief Lis vector wrapper class
 */
class LisVector
{
public:
	using IndexType = LIS_INT;
public:
	LisVector();

	/**
	 * Constructor for initialization of the number of rows
	 * @param length number of rows
	 */
	explicit LisVector(std::size_t length);

	/**
	 * Constructor using the given raw data
	 * @param length the length of the vector
	 * @param data   the raw data
	 */
	LisVector(std::size_t length, double* data);

	/// copy constructor
	LisVector(LisVector const& src);

	/**
	 *
	 */
	virtual ~LisVector();

	/// return a vector length
	std::size_t size() const;

	void resize(std::size_t length)
	{
		lis_vector_set_size(_vec, length, 0);
	}

	/// return a start index of the active data range
	std::size_t getRangeBegin() const { return 0; }
	/// return an end index of the active data range
	std::size_t getRangeEnd() const { return this->size(); }
	/// set all values in this vector
	LisVector& operator=(double v);

	/// access entry
	double operator[](IndexType rowId) const { return get(rowId); }
	/// get entry
	double get(IndexType rowId) const
	{
		double v = .0;
		lis_vector_get_value(_vec, rowId, &v);
		return v;
	}

	/// set entry
	void set(IndexType rowId, double v)
	{
		lis_vector_set_value(LIS_INS_VALUE, rowId, v, _vec);
	}

	/// add entry
	void add(IndexType rowId, double v)
	{
		lis_vector_set_value(LIS_ADD_VALUE, rowId, v, _vec);
	}

	/// printout this equation for debugging
	void write(const std::string& filename) const;

	/// return a raw Lis vector object
	LIS_VECTOR& getRawVector() { return _vec; }
	LIS_VECTOR const& getRawVector() const { return _vec; }
	/// vector operation: set data
	LisVector& operator=(const LisVector& src);

	/// vector operation: add
	void operator+=(const LisVector& v);

	/// vector operation: subtract
	void operator-=(const LisVector& v);

	///
	template <class T_SUBVEC>
	void add(const std::vector<IndexType>& pos, const T_SUBVEC& sub_vec)
	{
		for (std::size_t i = 0; i < pos.size(); ++i)
		{
			this->add(pos[i], sub_vec[i]);
		}
	}

	void add(const std::vector<IndexType>& pos, double x)
	{
		for (std::size_t i = 0; i < pos.size(); ++i)
		{
			this->add(pos[i], x);
		}
	}

	void componentwiseDivide(LisVector const& v)
	{
		// TODO: fragile: _vec passed twice
		lis_vector_pdiv(_vec, v.getRawVector(), _vec);
	}

private:
	LIS_VECTOR _vec;
};


inline double norm_1(const LisVector &v)
{
	double n;
	return lis_vector_nrm1(v.getRawVector(), &n);
	return n;
}

inline double norm_2(const LisVector &v)
{
	double n;
	return lis_vector_nrm2(v.getRawVector(), &n);
	return n;
}

inline double norm_max(const LisVector &v)
{
	double n;
	return lis_vector_nrmi(v.getRawVector(), &n);
	return n;
}


}  // MathLib

#endif  // LISVECTOR_H_
