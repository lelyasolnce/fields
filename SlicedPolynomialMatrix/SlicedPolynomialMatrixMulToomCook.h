#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixMulToomCook_H
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixMulToomCook_H

#include "SlicedPolynomialMatrix.h"
#include "linbox/matrix/DenseMatrix/blas-matrix.h"
#include <givaro/givpoly1dense.h>
#include <givaro/givpoly1denseops.inl>
#include <stdlib.h>
#include <stdint.h>

namespace LinBox
{ 
	template< class GField, class Operand1, class Operand2, class Operand3>
	class SlicedPolynomialMatrixMulToomCook
	{
	private:
		typedef typename Operand1::IntField IntField;
		typedef typename IntField::Element Element;
		typedef Givaro::Poly1Dom<IntField, Dense> PolyDom;
		typedef typename Operand1::Rep Rep;
		typedef BlasMatrix<IntField, Rep> Matrix;
		typedef typename Operand1::polynomial polynomial;
		Matrix& EvaluationInterpolationMatrices (Matrix& TC, Matrix& iTC);
		Matrix& mul (IntField& F, Matrix& CMatBloc,  Matrix& AMatBloc,  Matrix& BMatBloc,
								    size_t m,  size_t k,  size_t n,  size_t e, polynomial irreducible);
	public:
	   /* All matrix classes should be SlicedPolynomialMatrices.
		* Matrices are multiplied in a classic way, while polynomial entries are multilpied using Toom-Cook algorithm.
		* For efficiency a vector of e BlasMatrices is transformed into a BlasMatrix of e rows.
		* Each row is a vector of all elements of one matrix coefficient.
		*/
		Operand1 &operator() ( GField &GF, Operand1 &C,  Operand2 &A,  Operand3 &B) ;
	}; 
}

#include "SlicedPolynomialMatrixMulToomCook.inl"

#endif


