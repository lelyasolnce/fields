#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixParallelOperations_H
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixParallelOperations_H

#include <omp.h>

namespace LinBox
{	
	template<class Field>
	void mulParallel(Field& F, size_t m, size_t n, size_t k,
		typename Field::Element_ptr A, size_t lda,
		typename Field::Element_ptr B, size_t ldb,
		typename Field::Element_ptr C, size_t ldc)
	{
		typename Field::Element p = F.characteristic();
		#pragma omp parallel for shared(A, B, C, m, k, n, p) private(i, l, j, el1, el2, el3)
		{
			typename Field::Element el1;
			typename Field::Element el2;
			typename Field::Element el3;
			for (size_t i = 0; i < m; i++)
			{
				for (size_t j = 0; j < n; j++)
				{
					el1 = 0;
					for (size_t l = 0; l < k; l++)
					{
						el2 = A[i*lda + l]; el3 = B[l*ldb + j];
						el1 += (el2 * el3); el1 = el1 % p;
					}					
					C[i*ldc + j] = el1;
				}
			}
		}
	}

	template<class Field>
	void addParallel(Field& F, size_t m, size_t n,
		typename Field::Element_ptr A, size_t lda,
		typename Field::Element_ptr B, size_t ldb,
		typename Field::Element_ptr C, size_t ldc)
	{
		typename Field::Element p = F.characteristic();
		#pragma omp parallel for shared(A, B, C, m, n, p) private(i, j, el)
		{
			typename Field::Element el;
			typename Field::Element el2;
			typename Field::Element el3;
			for (size_t i = 0; i < m; i++)
			{
				for (size_t j = 0; j < n; j++)
				{					
					C[i*ldc + j] = (A[i*lda + j] + B[i*ldb + j]) % p;
				}
			}
		}
	}

	template<class Field>
	void subParallel(Field& F, size_t m, size_t n,
		typename Field::Element_ptr A, size_t lda,
		typename Field::Element_ptr B, size_t ldb,
		typename Field::Element_ptr C, size_t ldc)
	{
		typename Field::Element p = F.characteristic();
		#pragma omp parallel for shared(A, B, C, m, n, p) private(i, j, el)
		{
			typename Field::Element el;
			typename Field::Element el2;
			typename Field::Element el3;
			for (size_t i = 0; i < m; i++)
			{
				for (size_t j = 0; j < n; j++)
				{					
					C[i*ldc + j] = (p + A[i*lda + j] - B[i*ldb + j]) % p;
				}
			}
		}
	}
}

#endif