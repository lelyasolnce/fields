#ifndef __LINBOX_BM_SlicedPolynomialBM_SlicedPolynomialMatrixMulKaratsuba_INL
#define __LINBOX_BM_SlicedPolynomialBM_SlicedPolynomialMatrixMulKaratsuba_INL

#include "SlicedPolynomialMatrixAddSub.h"
#include <vector>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>

namespace LinBox
{
	template<class GField, class Operand1, class Operand2, class Operand3>
	typename SlicedPolynomialMatrixMulKaratsuba<GField, Operand1, Operand2, Operand3 >::vec&
		SlicedPolynomialMatrixMulKaratsuba<GField, Operand1, Operand2, Operand3 >::karatsuba(IF& F, vec& C, vec& A, vec& B)
	{
		if (A.size() == 1)
		{
			for (int i = 0; i <B.size(); i++)
			{
				BlasMatrixDomainMul<IF, BM, BM, BM>()(F, C[i], A[0], B[i]);
			}
			return C;
		}
		if (B.size() == 1)
		{
			for (int i = 0; i < A.size(); i++)
			{
				BlasMatrixDomainMul<IF, BM, BM, BM>()(F, C[i], A[i], B[0]);
			}
			return C;
		}
		size_t n1 = A[0].rowdim();
		size_t n2 = A[0].coldim();
		size_t n3 = B[0].coldim();
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
				BlasMatrixDomainAddin<IF, BM, BM>()(F, C[i], C1[i]);
			}
			int mm = 2 * m;
			for (int i = 0; i < C2.size(); i++)
			{
				BlasMatrixDomainAddin<IF, BM, BM>()(F, C[mm + i], C2[i]);
			}
			for (int i = 0; i < C3.size(); i++)
			{
				BlasMatrixDomainAddin<IF, BM, BM>()(F, C[m + i], C3[i]);
			}
			for (int i = 0; i < C1.size(); i++)
			{
				BlasMatrixDomainSubin<IF, BM, BM>()(F, C[m + i], C1[i]);
			}
			for (int i = 0; i < C2.size(); i++)
			{
				BlasMatrixDomainSubin<IF, BM, BM>()(F, C[m + i], C2[i]);
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
				BlasMatrixDomainAddin<IF, BM, BM>()(F, C[i], C1[i]);
			}
			for (int i = 0; i < C2.size(); i++)
			{
				BlasMatrixDomainAddin<IF, BM, BM>()(F, C[m + i], C2[i]);
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
				BlasMatrixDomainAddin<IF, BM, BM>()(F, C[i], C1[i]);
			}
			for (int i = 0; i < C2.size(); i++)
			{
				BlasMatrixDomainAddin<IF, BM, BM>()(F, C[m + i], C2[i]);
			}
			return C;
		}
		return C;
	}

	template<class GField, class Operand1, class Operand2, class Operand3>
	Operand1& SlicedPolynomialMatrixMulKaratsuba<GField, Operand1, Operand2, Operand3 >::operator()(GField& GF, Operand1& C, Operand2& A, Operand3& B)
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
		return C;
	}
}

#endif