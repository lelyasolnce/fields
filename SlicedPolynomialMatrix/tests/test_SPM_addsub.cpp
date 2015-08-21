#include "linbox/matrix/SlicedPolynomialMatrix/SlicedPolynomialMatrix.h"
#include "linbox/matrix/SlicedPolynomialMatrix/SlicedPolynomialMatrixAddSub.h"
#include <givaro/gfq.h>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <fstream>

typedef int64_t TT;
typedef Givaro::GFqDom<TT> GField;
typedef LinBox::SlicedPolynomialMatrix<GField, TT> SPM;

/*
 * Tests SlicedPolynomialMatrices addition and subtraction
 */

int main()
{
	GField _gf(5, 3);

	size_t m = 2, n = 5;

	std::filebuf fb1; fb1.open("test_SPM_addsub_1.txt", std::ios::in); std::istream file1(&fb1);
	std::filebuf fb2; fb2.open("test_SPM_addsub_2.txt", std::ios::in); std::istream file2(&fb2);
	std::filebuf fb3; fb3.open("file3.txt", std::ios::out); std::ostream file3(&fb3);
	SPM spm1(_gf, m, n); SPM spm2(_gf, m, n); spm1.read(file1); spm2.read(file2);
	SPM spm3(_gf, m, n); 
	
	LinBox::SlicedPolynomialMatrixAdd<GField, SPM, SPM, SPM>()(_gf, spm3, spm1, spm2);
	file3<<"\nC = A + B\n"; spm3.write(file3);
	LinBox::SlicedPolynomialMatrixSub<GField, SPM, SPM, SPM>()(_gf, spm3, spm1, spm2);
	file3<<"\nC = A - B\n"; spm3.write(file3);
	LinBox::SlicedPolynomialMatrixAddin<GField, SPM, SPM>()(_gf, spm3, spm1);
	file3<<"\nC += A\n"; spm3.write(file3);
	LinBox::SlicedPolynomialMatrixSubin<GField, SPM, SPM>()(_gf, spm3, spm1);
	file3<<"\nC -= A\n"; spm3.write(file3);
}
