#ifndef __LINBOX_SlicedPolynomialMatrixMulParallel_H
#define __LINBOX_SlicedPolynomialMatrixMulParallel_H

#include "SlicedPolynomialMatrixParallelOperations.h"
#include "SlicedPolynomialMatrix.h"
#include "linbox/matrix/DenseMatrix/blas-matrix.h"
#include <vector>
#include <stdlib.h>
#include <stdint.h>
#include <omp.h>

namespace LinBox
{ 
	template< class GField, class Operand1, class Operand2, class Operand3>
	class SlicedPolynomialMatrixMulKaratsubaParallel
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
	    * BlasMatrix addition, subtraction and multiplication are parallelized using OpenMP.
	    */
		Operand1 &operator() (GField &GF, Operand1 &C, Operand2 &A, Operand3 &B);
	}; 

	template< class GField, class Operand1, class Operand2, class Operand3>
	class SlicedPolynomialMatrixMulToomCookParallel
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
		* BlasMatrix addition, subtraction and multiplication are parallelized using OpenMP.
		*/
		Operand1 &operator() ( GField &GF, Operand1 &C,  Operand2 &A,  Operand3 &B) ;
	}; 
}

#include "SlicedPolynomialMatrixMulParallel.inl"

#endif


