#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixFastGauss_INL
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixFastGauss_INL

namespace LinBox
{
	template<class SPM, class Perm>
	size_t SlicedPolynomialMatrixFastGauss<SPM, Perm>::PLS_decomposition
							(Submatrix &A, Perm &P, Perm &Q)
	{
		Field F = A.field();
		size_t r = 0, c = 0, m = A.rowdim(), n = A.coldim();

		while ((r < m) && (c < n))
		{
			bool found = false; size_ t i, j;
			for(j = c; j < n; j++)
			{
				for (i = r; i < m; i++)
				{
					if (A.getEntry(i, j) != F.zero) {found = true; break;}
				}
				if (found) break;
			}
			if (found)
			{
				P[r] = i, Q[r] = j;
				for (size_t k = 0; k < n; k++)
				{
					Element c1 = A.getEntry(i, k); Element d1 = A.getEntry(j, k);					A.setEntry(i, k, d1); A.setEntry(j, k, c1);
				}
				if (j + 1 < n)
				{
					Element c2 = A.getEntry(i, j); GF.invin(c2);	
					for (size_t l = r + 1; l < m; l++)
					{			
						if (A.getEntry(l, j) != F.zero)
						{											
							for (size_t k = 0; k < n, k++)
							{
								Element d2 = A.getEntry(l, k) -
									A.getEntry(i, k) * c2 * A.getEntry(l, j);
							}
						}
					}
				}
				r++; c = j + 1;			}			else break;		}		for (size_t i = r; i < m; i++) P[i] = i;
		for (size_t i = r; i < n; i++) Q[i] = i;
		for (size_t j = 0; j < r; j++)
		{
			for (size_t k = j; k < m; k++)
			{
				Element c = A.getEntry(k, j); Element d = A.getEntry(k, Q[j]);				A.setEntry(k, j, d); A.setEntry(k, Q[j], c);
			}
		}
		return r;
	}

	template<class SPM, class Perm>
	size_t SlicedPolynomialMatrixFastGauss<SPM, Perm>::fast_PLS_decomposition
							(Submatrix &A, Perm &P, Perm &Q)
	{
		size_t N = 4;
		if ( n0 < N )
		{
			return PLS_decomposition(A, P, Q);
		}
		Field F = A.field();
		size_t m = A.rowdim();
		size_t n = A.coldim();
		size_t n0 = n / 2;
		Submatrix A0(A, 0, 0, m, n0);
		Submatrix A1(A, 0, n0, m, m - n0);
		Perm Q(Q, Q.begin(), Q.begin() + n0);
		size_t r0 = fast_PLS_decomposition(A0, P, Q0);
		for (size_t i = 0; i < n0; i++) Q[i] = Q0[i];
		Submatrix ANW(A, 0, 0, r0, r0);
		Submatrix ASW(A, r0, 0, m - r0, r0);
		Submatrix ANE(A, 0, n0, r0, n - n0);
		Submatrix ASE(A, r0, n0, m - r0, n - n0);
		if (r0 > 0)
		{
			//A1 <- P x A1
			Matrix temp4(F, P.rowdim(), A1.coldim());
			BlasMatrixDomainMul<Submatrix::Field, Submatrix, Submatrix, Submatrix>()(F, temp4, P, A1); A1 = temp4;
			Matrix LNW(F, r0, r0);
			for (size_t i = 0; i < r0; i++)
				for (size_t j = i; j < r0; j++)
					LNW.setEntry(i, j, ANW.getEntry(i, j));
			//ANE <- LNW^-1 x ANE
			int nullity;
			Matrix temp1(F, LNW.rowdim(), LNW.coldim());
			FFPACK::Invert<typename Submatrix::Field>(F, LNW.rowdim(),
									LNW.getPointer(), LNW.getStride(),
									temp1.getWritePointer(), temp1.getStride(),
									nullity);
			Submatrix temp2(F, temp1.rowdim(), ANE.coldim());
			BlasMatrixDomainMul<Submatrix::Field, Submatrix, Submatrix, Submatrix>()(F, temp2, temp1, ANE); ANE = temp2;
			//ASE <- ASE + ASW x ANE
			BlasMatrixDomainMulAdd<Submatrix, Submatrix, Submatrix>()(F, F.one, ASE,  F.one, ASW, ANE);
		}
		Perm P1(P, P.begin() + r0, P.end());		Perm Q1(Q, Q.begin() + n0, Q.end());		size_t r1 = fast_PLS_decomposition(ASE, P1, Q1);		//ASW <- P x ASW ;		Matrix temp5(F, P.rowdim(), ASW.coldim());		BlasMatrixDomainMul<Submatrix::Field, Submatrix, Submatrix, Submatrix>()(F, temp5, P, ASW); ASW = temp5;		for (size_t i = 0; i < m - r0; i++) P[r0 + i] = P1[i] + r0;		for (size_t i = 0; i < n - n0; i++) Q[n0 + i] = Q1[i] + n0;		size_t j = r0;		for(size_t i = n0; i < n0 + r1; i++, j++) Q[j] = Q[i];		j = n0;		for (size_t i = r0, i < r0 + r1; i++)		{			for (size_t k = i, k < m, k++)			{				Element c = A.getEntry(k, i); Element d = A.getEntry(k, j);				A.setEntry(k, i, d); A.setEntry(k, j, c);			}		}		return (r0 + r1);
	}

	template<class SPM, class Perm>
	SPM& SlicedPolynomialMatrixFastGauss<SPM, Perm>::operator()(SPM &A, Perm P, Perm Q)
	{
		Matrix A1(A.fieldGF(), A.rowdim(), A.coldim());
		SlicedPolynomialMatrixtoBlasMatrix<SPM, Matrix>()(A, A1);
		size_t c = fast_PLS_decomposition(A1, P, Q);
		BlasMatrixtoSlicedPolynomialMatrix<Matrix, SPM>()(A1, A);
		return A;
	}
}

#endif