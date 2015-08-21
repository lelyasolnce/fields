#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixMulToomCook_INL
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixMulToomCook_INL

#include <fflas-ffpack/fflas/fflas.h>
#include <fflas-ffpack/ffpack/ffpack.h>
#include <givaro/givpoly1denseops.inl>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include "linbox/matrix/DenseMatrix/blas-matrix.h"

namespace LinBox
{
	template<class GField, class Operand1, class Operand2, class Operand3>
	BlasMatrix<typename Operand1::IntField, typename Operand1::Rep>&
	SlicedPolynomialMatrixMulToomCook<GField, Operand1, Operand2, Operand3 >::EvaluationInterpolationMatrices
				(Matrix& TC, Matrix& iTC)
	{
		size_t E = TC.rowdim();
		Element el;
		for (size_t i = 0 ; i < E ; ++i)
		{
			for (size_t j = 0 ; j < E ; ++j)
			{
				el = pow(i, j);
				TC.setEntry(i, j, el);
			}
		}
		int null;
		FFPACK::Invert<IntField>(TC.field(),E,TC.getPointer(),E,iTC.getWritePointer(),E,null);
		return TC;
	}

	template<class GField, class Operand1, class Operand2, class Operand3>
		BlasMatrix<typename Operand1::IntField, typename Operand1::Rep>&
		SlicedPolynomialMatrixMulToomCook<GField, Operand1, Operand2, Operand3 >::mul
			(IntField& F, Matrix& CMatBloc, Matrix& AMatBloc, Matrix& BMatBloc,
								    size_t m,
								    size_t k,
								    size_t n,  size_t e,
								   polynomial irreducible)
	{
		#if (__LINBOX_FFLAS_FFPACK_VERSION < 10501)F
				#warning "Invert is buggy in your fflas-ffpack version. please consider upgrading to >=1.5.1."
		#endif
		size_t E = 2*e - 1 ;

		Matrix TC    (F, E, E);
		Matrix iTC   (F, E, E);
		Matrix iEval (F, E, E);
		EvaluationInterpolationMatrices(TC, iTC);

		size_t s = m * n;
		Matrix TMatBloc( F, E, s);
		s = m * k;
		Matrix AEval( F , E, s);
		s = k * n;
		Matrix BEval( F , E, s);

		s = m * k;
		FFLAS::fgemm(F,
					 FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
					 E, s, e,
					 F.one,
					 TC.getPointer(),E,
					 AMatBloc.getPointer(), s,
					 F.zero,
					 AEval.getWritePointer(), s);


		s = n * k;
		FFLAS::fgemm(F,
					 FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
					 E ,s, e,
					 F.one,
					 TC.getPointer(),E,
					 BMatBloc.getPointer(), s,
					 F.zero,
					 BEval.getWritePointer(), s);


		for (size_t i = 0 ; i < E ; ++i)
		{
			FFLAS::fgemm(F,
						 FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
						 m, n, k,
						 F.one,
						 AEval.getPointer()+i*m*k, k,
						 BEval.getPointer()+i*n*k, n,
						 F.zero,
						 TMatBloc.getWritePointer()+i*m*n, n);
		}

		s = m * n;
		FFLAS::fgemm(F,
					 FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
				     	 E, s, E,
					 F.one,
					 iTC.getPointer(),E,
					 TMatBloc.getWritePointer(), s,
					 F.zero,
					 CMatBloc.getWritePointer(), s);

		return CMatBloc;
	}

	template<class GField, class Operand1, class Operand2, class Operand3>
	Operand1& SlicedPolynomialMatrixMulToomCook<GField, Operand1, Operand2, Operand3 >::operator()
									   ( GField& GF,
									   Operand1& C,
									    Operand2& A,
									    Operand3& B)
	{
		PolyDom polydom(C.fieldF());
		size_t e = C.length();
		size_t m = C.rowdim();
		size_t k = B.rowdim();
		size_t n = C.coldim();

		IntField F = A.fieldF();

		if (e == 1) {
			size_t n1 = A.rowdim();
			size_t n2 = A.coldim();
			size_t n3 = B.coldim();
			Matrix Am(F, n1, n2);
			Am = A.getMatrixCoefficient(0);
			Matrix Bm(F, n2, n3);
			Bm = B.getMatrixCoefficient(0);
			Matrix Cm(F, n1, n3);
			FFLAS::fgemm(F,
						 FFLAS::FflasNoTrans, FFLAS::FflasNoTrans,
						 m,n,k,
						 F.one,
						 Am.getPointer(), Am.getStride(),
						 Bm.getPointer(), Bm.getStride(),
						 F.zero,
						 Cm.getWritePointer(), Cm.getStride());
			C.setMatrixCoefficient(0, Cm);
			return C;
		}
		size_t s = m * n; size_t E = 2*e - 1;
		Matrix Cbloc(F,E,s);
		s = m * k;
		Matrix Abloc(F,e,s);
		s = k * n;
		Matrix Bbloc(F,e,s);

		Element el;
		for (size_t l = 0 ; l < e ; ++l)
		{
			for (size_t i = 0 ; i < m ; ++i)
			{
				for (size_t j = 0 ; j < k ; ++j)
				{
					s = i*k + j; el = A.getEntry(l, i, j);
					Abloc.setEntry(l, s, el);
				}
			}
		}

		for (size_t l = 0 ; l < e ; ++l)
		{
			for (size_t i = 0 ; i < k ; ++i)
			{
				for (size_t j = 0 ; j < n ; ++j)
				{
					s = i*n+j; el = B.getEntry(l, i, j);
					Bbloc.setEntry(l, s, el);
				}
			}
		}

		polynomial irred = C.irreducible();
		mul(C.fieldF(), Cbloc, Abloc, Bbloc, m, k, n, e, irred);
		
		polynomial x; for (int l = 0; l < E; l++) x.push_back(F.zero);
		for (size_t i = 0 ; i < m ; ++i)
		{
			for (size_t j = 0 ; j < n ; ++j)
			{
				for (size_t l = 0 ; l < E ; ++l)
				{
					s = i*n+j; x[l] = Cbloc.getEntry(l, s);
				}
				polynomial z; polydom.mod(z, x, irred);
				for (size_t l = 0 ; l < e ; ++l)
				{
					C.setEntry(l, i, j, z[l]);
				}
			}
		}

		return C;
	}
}

#endif