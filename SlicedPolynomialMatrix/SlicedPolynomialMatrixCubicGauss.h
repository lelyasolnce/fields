#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixCubicGauss_H
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixCubicGauss_H

#include <stdlib.h>
#include <stdint.h>
#include "SlicedPolynomialMatrix.h"
#include "linbox/matrix/DenseMatrix/blas-matrix.h"
#include <givaro/gfq.h>
#include <givaro/modular.h>
#include "SlicedPolynomialMatrixConversion.h"
#include <vector>

namespace LinBox
{ 
	template<class SPM>
	class SlicedPolynomialMatrixCubicGauss
	{
	private:
		typedef typename SPM::Field GField;
		typedef typename GField::Element Element;
	public:
		/* @param A1 is mxn SlicedPolynomialMatrix. A1 is reduced to row echelon form by row operations.
		 * @param C1 is mxm SlicedPolynomialMatrix. C1 is the corresponding permutation matrix.
		 */
		SPM &operator() (SPM &A1, SPM &C1)
		{
			size_t m = A1.rowdim(); size_t n = A1.coldim();
			BlasMatrix<SPM::GField, std::vector<Element>> A(A1.fieldGF(), m, n);
			BlasMatrix<SPM::GField, std::vector<Element>> C(A1.fieldGF(), m, m);
			SlicedPolynomialMatrixtoBlasMatrix<SPM, BM>()(A1, A);
			SlicedPolynomialMatrixtoBlasMatrix<SPM, BM>()(C1, C);
			Element el, el1, el2; GField GF = A.field();
			size_t M = A.rowdim();
			size_t N = A.coldim();
			for (size_t u = 0; u < M; u++)
			{
				for (size_t v = 0; v < N; v++)
				{
					el = (u == v) ? GF.one : GF.zero;
					C.setEntry(u, v, el);
				}
			}
			size_t i = 0;
			size_t j = 0;
			size_t s = 0;
			size_t t = 0;
			while ((i < M) && (j < N))
			{
				s = i;
				t = j;
				bool nonzero = false;
				while ((!nonzero) && (t < N))
				{
					if (A.getEntry(s, t) != F.zero)
					{
						nonzero = true;
					}
					else
					{
						if (s == M - 1)
						{
							s = i;
							t++;
						}
						else
						{
							s++;
						}
					}
				}
				if (!nonzero) return A;
				else
				{
					for (size_t v = 0; v < N; v++) 
					{
						el1 = A.getEntry(i, v);
						el2 = A.getEntry(s, v);
						A.setEntry(i, v, el2);
						A.setEntry(s, v, el1);
					}
					for (size_t v = 0; v < M; v++) 
					{
						el1 = C.getEntry(i, v);
						el2 = C.getEntry(s, v);
						C.setEntry(i, v, el2);
						C.setEntry(s, v, el1);

					}
					el1 = A.getEntry(i, t);
					GF.inv(el, el1);
					for (size_t v = 0; v < N; v++) 
					{
						GF.mul(el1, A.getEntry(i, v), el);
						A.setEntry(i, v, el1);
					}
					for (size_t v = 0; v < M; v++) 
					{
						GF.mul(el1, C.getEntry(i, v), el);
						C.setEntry(i, v, el1);
					}
					for (size_t u = 0; u < M; u++) 
					{
						if (u != i)
						{
							el = A.getEntry(u, t);
							for (size_t v = 0; v < N; v++) 
							{
								F.mul(el1, el, A.getEntry(i, v))
								F.sub(el2, A.getEntry(u, v), el1);
								A.setEntry(u, v, el2);
							}
							for (size_t v = 0; v < M; v++) 
							{
								F.mul(el1, el, C.getEntry(i, v))
								F.minus(el2, C.getEntry(u, v), el1);
								C.setEntry(u, v, el2);
							}
						}
					}
					i++;
					j++;
				}
			}
			BlasMatrixtoSlicedPolynomialMatrix<BM, SPM>()(A, A1);
			BlasMatrixtoSlicedPolynomialMatrix<BM, SPM>()(C, C1);
			return A1;
		}
	};
}

#endif
