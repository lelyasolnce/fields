#ifndef __LINBOX_SlicedPolynomialMatrixMulParallel_INL
#define __LINBOX_SlicedPolynomialMatrixMulParallel_INL

#include <givaro/givpoly1dense.h>
#include <givaro/givpoly1denseops.inl>
#include "linbox/matrix/MatrixDomain/blas-matrix-domain.h"
#include <vector>
#include <stdlib.h>
#include <stdint.h>
#include <fflas-ffpack/fflas/fflas.h>
#include <fflas-ffpack/ffpack/ffpack.h>
#include <givaro/givpoly1dense.h>
#include <givaro/givpoly1denseops.inl>
#include <math.h>

namespace LinBox
{
	template<class GField, class Operand1, class Operand2, class Operand3>
	typename SlicedPolynomialMatrixMulKaratsubaParallel<GField, Operand1, Operand2, Operand3 >::vec&
		SlicedPolynomialMatrixMulKaratsubaParallel<GField, Operand1, Operand2, Operand3 >::karatsuba(IF& F, vec& C, vec& A, vec& B)
	{
	size_t n1 = A[0].rowdim();
	size_t n2 = A[0].coldim();
	size_t n3 = B[0].coldim();
		if (A.size() == 1)
		{
			for (int i = 0; i <B.size(); i++)
			{
				mulParallel(F, n1, n3, n2,
					A[0].getPointer(), n2,
					B[i].getPointer(), n3,
					C[i].getWritePointer(), n3);

			}
			return C;
		}
		if (B.size() == 1)
		{
			for (int i = 0; i < A.size(); i++)
			{
				mulParallel(F, n1, n3, n2,
					A[i].getPointer(), n2,
					B[0].getPointer(), n3,
					C[i].getWritePointer(), n3);
			}
			return C;
		}
		int m = (A.size() < B.size()) ? (B.size() / 2) : (A.size() / 2);
		if ((m < A.size()) && (m < B.size()))
		{
			vec A1(A.begin(), A.begin() + m);
			vec A2(A.begin() + m, A.end());
			vec B1(B.begin(), B.begin() + m);
			vec B2(B.begin() + m, B.end());
			vec A3;
			int minlength_a = (A1.size() < A2.size()) ? A1.size() : A2.size();
			int maxlength_a = (A1.size() < A2.size()) ? A2.size() : A1.size();
			for (int i = 0; i < minlength_a; i++)
			{
				BM AA(F, n1, n2);
				A3.push_back(BlasMatrixDomainAdd<IF, BM, BM, BM>()(F, AA, A1[i], A2[i]));
			}
			if (maxlength_a == A1.size())
			{
				for (int i = minlength_a; i < maxlength_a; i++)
				{
					A3.push_back(A1[i]);
				}
			}
			else
			{
				for (int i = minlength_a; i < maxlength_a; i++)
				{
					A3.push_back(A2[i]);
				}
			}
			vec B3;
			int minlength_b = (B1.size() < B2.size()) ? B1.size() : B2.size();
			int maxlength_b = (B1.size() < B2.size()) ? B2.size() : B1.size();
			for (int i = 0; i < minlength_b; i++)
			{
				BM BB(F, n2, n3);
				B3.push_back(BlasMatrixDomainAdd<IF, BM, BM, BM>()(F, BB, B1[i], B2[i]));
			}
			if (maxlength_b == B1.size())
			{
				for (int i = minlength_b; i < maxlength_b; i++)
				{
					B3.push_back(B1[i]);
				}
			}
			else
			{
				for (int i = minlength_b; i < maxlength_b; i++)
				{
					B3.push_back(B2[i]);
				}
			}
			vec C1;
			vec C2;
			vec C3;
			BM CC(F, n1, n3);
			int xx;
			xx = A1.size() + B1.size() - 1;
			for (int i = 0; i < xx; i++)
			{
				C1.push_back(CC);
			}
			xx = A2.size() + B2.size() - 1;
			for (int i = 0; i < xx; i++)
			{
				C2.push_back(CC);
			}
			xx = A3.size() + B3.size() - 1;
			for (int i = 0; i < xx; i++)
			{
				C3.push_back(CC);
			}
			karatsuba(F, C1, A1, B1);
			karatsuba(F, C2, A2, B2);
			karatsuba(F, C3, A3, B3);
			for (int i = 0; i < C1.size(); i++)
			{
				addParallel(F, n1, n3,
					C[i].getPointer(), n3,
					C1[i].getPointer(), n3,
					C[i].getWritePointer(), n3);
			}
			int mm = 2 * m;
			for (int i = 0; i < C2.size(); i++)
			{
				addParallel(F, n1, n3,
					C[mm + i].getPointer(), n3,
					C2[i].getPointer(), n3,
					C[mm + i].getWritePointer(), n3);
			}
			for (int i = 0; i < C3.size(); i++)
			{
				addParallel(F, n1, n3,
					C[m + i].getPointer(), n3,
					C3[i].getPointer(), n3,
					C[m + i].getWritePointer(), n3);
			}
			for (int i = 0; i < C1.size(); i++)
			{
				subParallel(F, n1, n3,
					C[m + i].getPointer(), n3,
					C1[i].getPointer(), n3,
					C[m + i].getWritePointer(), n3);
			}
			for (int i = 0; i < C2.size(); i++)
			{
				subParallel(F, n1, n3,
					C[m + i].getPointer(), n3,
					C2[i].getPointer(), n3,
					C[m + i].getWritePointer(), n3);
			}
			return C;
		}
		if (A.size() <= m)
		{
			vec B1(B.begin(), B.begin() + m);
			vec B2(B.begin() + m, B.end());
			vec C1;
			vec C2;
			BM CC(F, n1, n3);
			int xx;
			xx = A.size() + B1.size() - 1;
			for (int i = 0; i < xx; i++)
			{
				C1.push_back(CC);
			}
			xx = A.size() + B2.size() - 1;
			for (int i = 0; i < xx; i++)
			{
				C2.push_back(CC);
			}
			karatsuba(F, C1, A, B1);
			karatsuba(F, C2, A, B2);
			for (int i = 0; i < C1.size(); i++)
			{
				addParallel(F, n1, n3,
					C[i].getPointer(), n3,
					C1[i].getPointer(), n3,
					C[i].getWritePointer(), n3);
			}
			for (int i = 0; i < C2.size(); i++)
			{
				addParallel(F, n1, n3,
					C[m + i].getPointer(), n3,
					C2[i].getPointer(), n3,
					C[m + i].getWritePointer(), n3);
			}
			return C;
		}
		if (B.size() <= m)
		{
			vec A1(A.begin(), A.begin() + m);
			vec A2(A.begin() + m, A.end());
			vec C1;
			vec C2;
			BM CC(F, n1, n3);
			int xx;
			xx = A1.size() + B.size() - 1;
			for (int i = 0; i < xx; i++)
			{
				C1.push_back(CC);
			}
			xx = A2.size() + B.size() - 1;
			for (int i = 0; i < xx; i++)
			{
				C2.push_back(CC);
			}
			karatsuba(F, C1, A1, B);
			karatsuba(F, C2, A2, B);
			for (int i = 0; i < C1.size(); i++)
			{
				addParallel(F, n1, n3,
					C[i].getPointer(), n3,
					C1[i].getPointer(), n3,
					C[i].getWritePointer(), n3);
			}
			for (int i = 0; i < C2.size(); i++)
			{
				addParallel(F, n1, n3,
					C[m + i].getPointer(), n3,
					C2[i].getPointer(), n3,
					C[m + i].getWritePointer(), n3);
			}
			return C;
		}
		return C;
	}

	template<class GField, class Operand1, class Operand2, class Operand3>
	Operand1& SlicedPolynomialMatrixMulKaratsubaParallel<GField, Operand1, Operand2, Operand3 >::operator()(GField& GF, Operand1& C, Operand2& A, Operand3& B)
	{
		int xx;
		vec A1;
		vec B1;
		xx = A.length();
		for (int m = 0; m < xx; m++)
		{
			A1.push_back(A.getMatrixCoefficient(m));
		}
		xx = B.length();
		for (int m = 0; m < xx; m++)
		{
			B1.push_back(B.getMatrixCoefficient(m));
		}
		vec C1;
		xx = A1.size() + B1.size() - 1;
		size_t n1 = A.rowdim();
		size_t n2 = A.coldim();
		size_t n3 = B.coldim();
		BM CC(C.fieldF(), n1, n3);
		for (int i = 0; i < xx; i++)
		{
			C1.push_back(CC);
		}
		karatsuba(C.fieldF(), C1, A1, B1);
		
		PolyDom polydom(C.fieldF());
		polynomial _irred = C.irreducible();
		int mk = C1.size();
		int mi = C1[0].rowdim();
		int mj = C1[0].coldim();

		#pragma omp parallel for shared(C, C1, mi, mj, mk, _irred, polydom) private(entry, result, i, j, k)
		{
			for (int i = 0; i < mi; i++)
			{
				for (int j = 0; j < mj; j++)
				{
					polynomial entry;
					for (int k = 0; k < mk; k++)
					{
						entry.push_back(C1[k].getEntry(i, j));
					}
				
					polynomial result;
					polydom.mod(result, entry, _irred);
					for (size_t k = 0; k < C.length(); k++)
					{
						C.setEntry(k, i, j, result[k]);
					}
				}
			}
		}
		return C;
	}
}



namespace LinBox
{
	template<class GField, class Operand1, class Operand2, class Operand3>
	BlasMatrix<typename Operand1::IntField, typename Operand1::Rep>&
	SlicedPolynomialMatrixMulToomCookParallel<GField, Operand1, Operand2, Operand3 >::EvaluationInterpolationMatrices
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
		SlicedPolynomialMatrixMulToomCookParallel<GField, Operand1, Operand2, Operand3 >::mul
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
		mulParallel(F, E, s, e,
				TC.getPointer(), E,
				AMatBloc.getPointer(), s,
				AEval.getWritePointer(), s);
					 
		s = n * k;
		mulParallel(F, E ,s, e,
				 TC.getPointer(),E,
				 BMatBloc.getPointer(), s,
				 BEval.getWritePointer(), s);

		for (size_t i = 0 ; i < E ; ++i)
		{
			mulParallel(F, m, n, k,
					 AEval.getPointer()+i*m*k, k,
					 BEval.getPointer()+i*n*k, n,
					 TMatBloc.getWritePointer()+i*m*n, n);
		}

		s = m * n;
		mulParallel(F, E, s, E,
				 iTC.getPointer(),E,
				 TMatBloc.getWritePointer(), s, 
				 CMatBloc.getWritePointer(), s);

		return CMatBloc;
	}

	template<class GField, class Operand1, class Operand2, class Operand3>
	Operand1& SlicedPolynomialMatrixMulToomCookParallel<GField, Operand1, Operand2, Operand3 >::operator()
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
			mulParallel(F, m, n, k,
					Am.getPointer(), Am.getStride(),
					Bm.getPointer(), Bm.getStride(),
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
					C.setEntry(l, i, j, x[l]);
				}
			}
		}

		return C;
	}
}

#endif