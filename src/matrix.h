#pragma once

#include <cassert>
#include <vector>

namespace rlscore {
// these 4 lines are used 3 times verbatim - defining a temp macro to ease the pain
#define RL_MATRIX_DEFINE_OPERATORS \
	const T& operator()(Size_Type i) const { return m_data[i]; } \
	      T& operator()(Size_Type i)       { return m_data[i]; } \
	const T& operator()(Size_Type i, Size_Type j) const { return m_data[index(i, j)]; } \
	      T& operator()(Size_Type i, Size_Type j)       { return m_data[index(i, j)]; }

////////////////////////////////////////////////////////////////////////
template<typename T>
class Matrix {
	std::vector<T> m_data;
	Size_Type m_i, m_j;
  public:
	// column-major
	Matrix() : m_i(0), m_j(0) {}
	Matrix(Size_Type i, Size_Type j, const T& filler_val) : m_data(i*j, filler_val), m_i(i), m_j(j) {
		assert(m_j >=0);
		assert(m_i >=0);
	}
	Size_Type index(Size_Type i, Size_Type j) const {//return the index in m_data, column-major
		#ifdef MODE_DEBUG
		assert(j < m_j && j >=0);
		assert(i < m_i && i >=0);
		#endif
		return i + m_i*j; //column-major
	}
	void resize(Size_Type m, Size_Type n, const T& filler_val) {//resize to mxn, new sizes should be the same or greater than the old, preserves original data
		assert(m >= dim_1());//new sizes should be the same or greater than the old
		assert(n >= dim_2());
		if(m == dim_1() && n == dim_2()) return; // no-op
		std::vector<T> tmp(m*n, filler_val);
		for(Size_Type i = 0; i < m_i; i++) {
			for(Size_Type j = 0; j < m_j; j++) {
				tmp[i+m*j] = (*this)(i, j);
			}
		}
		m_data = tmp;
		m_i = m;
		m_j = n;
	}
	void append(const Matrix<T>& x, const T& filler_val) {//append Matrix x through the diagonal, it is a Matrix with two diagonal Matrix blocks
		Size_Type m = dim_1();
		Size_Type n = dim_2();
		resize(m + x.dim_1(), n + x.dim_2(), filler_val);
		for(Size_Type i = 0; i < x.dim_1(); i++) {
			for(Size_Type j = 0; j < x.dim_2(); j++) {
				(*this)(i+m, j+n) = x(i, j);
			}
		}
	}
	RL_MATRIX_DEFINE_OPERATORS // temp macro defined above, one for indexing vector, one for Matrix indexing
	Size_Type dim_1() const { return m_i; }
	Size_Type dim_2() const { return m_j; }
};
//////////////////////////////////////////////////////////////////////////




// inline Size_Type triangular_matrix_index(Size_Type n, Size_Type i, Size_Type j) {
// 	assert(j < n);
// 	assert(i <= j);
// 	return i + j*(j+1)/2;
// }
// //permissive
// inline Size_Type triangular_matrix_index_permissive(Size_Type n, Size_Type i, Size_Type j) {
// 	return (i <= j) ? triangular_matrix_index(n, i, j) : triangular_matrix_index(n, j, i);
// }
//only upper-right triangular{i,j}
// 0  1  3  6
//    2  4  7
//       5  8
//          9
template<typename T>
class Triangular_Matrix {
	std::vector<T> m_data;
	Size_Type m_dim;
  public:
	Triangular_Matrix() : m_dim(0) {}
	Triangular_Matrix(Size_Type n, const T& filler_val) : m_data(n*(n+1)/2, filler_val), m_dim(n) {
		assert(n > 0);
	}
	Size_Type index(Size_Type i, Size_Type j) const {//return corresponding index in vector m_data
		#ifdef MODE_DEBUG
		assert(i >= 0 && j >= 0);
		assert(j < m_dim);
		assert(i <= j);
		#endif
		return i + j*(j+1)/2;
	}
	Size_Type index_permissive(Size_Type i, Size_Type j) const {
		return (i < j) ? index(i, j) : index(j, i);
	}
	void resize(Size_Type n, const T& filler_val) {//resize m_data to has only the upper right Matrix elements, new sizes should be the same or greater than the old, preserves original data
		assert(n >= m_dim);
		if(n == m_dim) return; // no-op
		m_dim = n;
		m_data.resize(n*(n+1)/2, filler_val); // preserves original data
	}
	RL_MATRIX_DEFINE_OPERATORS // temp macro defined above, one for indexing vector, one for Matrix indexing
	Size_Type dim() const {//return m_dim, dimension of the Matrix
		return m_dim;
	}
};
///////////////////////////////////////////////////////////////////////////





//no diagonal elements
template<typename T>
class Strictly_Triangular_Matrix {
	std::vector<T> m_data;
	Size_Type m_dim;
  public:
	Strictly_Triangular_Matrix() : m_dim(0) {}
	Strictly_Triangular_Matrix(Size_Type n, const T& filler_val) : m_data(n*(n-1)/2, filler_val), m_dim(n) {
		assert(n > 0);
	}
	Size_Type index(Size_Type i, Size_Type j) const {//return corresponding index in vector m_data
		#ifdef MODE_DEBUG
		assert(i >= 0 && j > 0);
		assert(j < m_dim);
		assert(i < j);
		assert(j >= 1); // by implication, really
		#endif
		return i + j*(j-1)/2;
	}
	Size_Type index_permissive(Size_Type i, Size_Type j) const {
		return (i < j) ? index(i, j) : index(j, i);
	}
	void resize(Size_Type n, const T& filler_val) {//resize m_data to has only the upper right Matrix elements, new sizes should be the same or greater than the old, preserves original data
		assert(n >= m_dim);
		if(n == m_dim) return; // no-op
		m_dim = n;
		m_data.resize(n*(n-1)/2, filler_val); // preserves original data
	}
	void append(const Strictly_Triangular_Matrix<T>& m, const T& filler_val) {//append Matrix x through the diagonal, it is a Matrix with two diagonal Matrix blocks
		Size_Type n = dim();
		resize(n + m.dim(), filler_val);
		for(Size_Type i = 0; i < m.dim(); ++i) {
			for(Size_Type j = i+1; j < m.dim(); ++j) {
				(*this)(i+n, j+n) = m(i, j);
			}
		}
	}
	void append(const Matrix<T>& rectangular, const Strictly_Triangular_Matrix<T>& triangular) {//arguments rectangular and triangular must have the same dim
		//if(rectangular.dim_2() == 0) do nothing
		//if(rectangular.dim_1() == 0) assign argument triangular to this object
		//else append these three matrices like the following
		// i  i  i  r  r  r  r
		//    i  i  r  r  r  r
		//       i  r  r  r  r
		//          t  t  t  t
		//             t  t  t
		//                t  t
		//                   t
		//where i,r,t are this, rectangular and triangular, respectively.
		assert(dim() == rectangular.dim_1());
		assert(rectangular.dim_2() == triangular.dim());
		// a filler value is needed by append or resize
		// we will use a value from rectangular as the filler value
		// but it can not be obtained if dim_1 or dim_2 is 0
		// these cases have to be considered separately
		if(rectangular.dim_2() == 0) return;
		if(rectangular.dim_1() == 0) {
			(*this) = triangular;
			return;
		}
		const T& filler_val = rectangular(0, 0); // needed by 'append below'
		Size_Type n = dim();
		append(triangular, filler_val);
		for(Size_Type i = 0; i < rectangular.dim_1(); ++i) {
			for(Size_Type j = 0; j < rectangular.dim_2(); ++j) {
				(*this)(i, n + j) = rectangular(i, j);
			}
		}
	}
	RL_MATRIX_DEFINE_OPERATORS // temp macro defined above, one for indexing vector, one for Matrix indexing
	Size_Type dim() const {
		return m_dim;
	}
};
#undef RL_MATRIX_DEFINE_OPERATORS
////////////////////////////////////////////////////////////////////////////////////////////

}