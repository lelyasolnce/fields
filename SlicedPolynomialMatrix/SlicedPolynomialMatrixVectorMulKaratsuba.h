#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixVectorMulKaratsuba_H
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixVectorMulKaratsuba_H

#include "SlicedPolynomialMatrix.h"
#include "SlicedPolynomialVector.h"
#include <vector>
#include <stdlib.h>
#include <stdint.h>
#include "linbox/matrix/DenseMatrix/blas-matrix.h"
#include "linbox/vector/blas-vector.h"
#include <givaro/givpoly1dense.h>
#include <givaro/givpoly1denseops.inl>

namespace LinBox
{ 
	template< class Field, class Operand1, class Operand2, class Operand3>
	class SlicedPolynomialMatrixVectorMulKaratsuba
	{
	private:
		typedef typename Operand1::IntField IntField;
		typedef typename BlasMatrix<IntField> Matrix;
		typedef typename BlasVector<IntField> Vector;
		typedef std::vector<Matrix> vecM;
		typedef std::vector<Vector> vecV;
		typedef typename Operand1::polynomial polynomial;
		vec& karatsuba(IntField& F, vecV& C, vecM& A, vecV& B);
	public:
		Operand1 &operator() (const Field &GF, Operand1 &C, const Operand2 &A, const Operand3 &B) const;
	}; 
}

#include "SlicedPolynomialMatrixVectorMulKaratsuba.inl"

#endif