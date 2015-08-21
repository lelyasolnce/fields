#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixFastGauss_H
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixFastGauss_H

#include <stdlib.h>
#include <stdint.h>
#include "SlicedPolynomialMatrix.h"
#include "linbox/matrix/DenseMatrix/blas-matrix.h"
#include "linbox/matrix/DenseMatrix/blas-matrix-domain.h"
#include <givaro/gfq.h>
#include <givaro/modular.h>
#include "SlicedPolynomialMatrixConversion.h"
#include <fflas-ffpack/ffpack/ffpack.h>
#include <vector>

namespace LinBox
{ 
	template<class SPM, class Perm>
	class SlicedPolynomialMatrixFastGauss
	{
	private:
		typedef typename SPM::IntField Field;
		typedef typename Field::Element Element;
		typename BlasMatrix<Field, std::vector<Element>> Matrix;
		typedef typename Matrix::subMatrixType Submatrix;
		size_t PLS_decomposition(Submatrix &A, Perm &P, Perm &Q);
		size_t fast_PLS_decomposition(Submatrix &A, Perm &P, Perm &Q);
	public:
		SPM &operator() (SPM &A, Perm P, Perm Q);
	};
}

#endif