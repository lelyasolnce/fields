#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_MATRIX_WINDOW_H
#define __LINBOX_matrix_SlicedPolynomialMatrix_MATRIX_WINDOW_H

namespace LinBox
{
	/*
	 *class MatrixWindow creates a Submatrix and secures that elements outside the window are not changed.
	 */

	template <class ParentMatrix>
	class SlicedPolynomialMatrixWindow
	{
		typedef ParentMatrix SPM;
		typedef typename SlicedPolynomialSubmatrix<ParentMatrix> SPSubM;
		SPSubM createMatrixWindow (SPM &M, size_t rowbeg, size_t colbeg, size_t Rowdim, size_t Coldim)
		{
			M.windowed = true;
			SPSubM x(M, rowbeg, colbeg, Rowdim, Coldim, true);
			return x;
		}
		SPSubM createMatrixWindow (SPM &M)
		{
			M.windowed = true;
			SPSubM x(M, true);
			return x;
		}
		SPSubM createMatrixWindow (SPSubM &SM, size_t rowbeg, size_t colbeg, size_t Rowdim, size_t Coldim)
		{
			SM.windowed = true;
			SPSubM x(SM, rowbeg, colbeg, Rowdim, Coldim, true);
			return x;
		}
		SPSubM createMatrixWindow (SPSubM &SM)
		{
			SM.windowed = true;
			SPSubM x(SM, true);
			return x;
		}
		void deleteMatrixWindow(SPSubM &x)
		{
			(x.pm).is_windowed = false;
			x.~SlicedPolynomialSubmatrix();
		}
	};
}

#endif