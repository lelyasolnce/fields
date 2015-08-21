#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixMulKaratsuba_H
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixMulKaratsuba_H

#include <vector>
#include <stdlib.h>
#include <stdint.h>
#include "SlicedPolynomialMatrix.h"
#include "linbox/matrix/DenseMatrix/blas-matrix.h"
#include <givaro/givpoly1dense.h>
#include <givaro/givpoly1denseops.inl>

namespace LinBox
{ 
	template< class GField, class Operand1, class Operand2, class Operand3>
	class SlicedPolynomialMatrixMulKaratsuba
	{
	private:
		typedef typename Operand1::IntField IF;
		typedef typename Operand1::MatrixElement ME;
		typedef BlasMatrix<IF, std::vector<ME>> BM;
		typedef std::vector<BM> vec;
		typedef typename Operand1::polynomial polynomial;
		typedef Givaro::Poly1Dom<IF, Dense> PolyDom;		
		vec& karatsuba(IF& F, vec& A, vec& B, vec& C);
	public:
	   /* Matrix vectors are multiplied using Karatsuba algorithm.
	    */
		Operand1 &operator() (GField &GF, Operand1 &C, Operand2 &A, Operand3 &B);
	}; 
}

#include "SlicedPolynomialMatrixMulKaratsuba.inl"

#endif

