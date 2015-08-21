#include "linbox/matrix/SlicedPolynomialMatrix/SlicedPolynomialMatrix.h"
#include "linbox/matrix/SlicedPolynomialMatrix/SlicedPolynomialMatrixMulKaratsuba.h"
#include "linbox/matrix/SlicedPolynomialMatrix/SlicedPolynomialMatrixMulToomCook.h"
#include "linbox/matrix/SlicedPolynomialMatrix/SlicedPolynomialMatrixMulParallel.h"
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
	size_t m1 = 2, m2 = 3, m3 = 2;

	std::filebuf fb1; fb1.open("test_SPM_mul_1.txt", std::ios::in); std::istream file1(&fb1);
	std::filebuf fb2; fb2.open("test_SPM_mul_2.txt", std::ios::in); std::istream file2(&fb2);
	std::filebuf fb3; fb3.open("file3.txt", std::ios::out); std::ostream file3(&fb3);
	SPM spm1(_gf, m1, m2); SPM spm2(_gf, m2, m3); spm1.read(file1); spm2.read(file2);
	SPM spm3(_gf, m1, m3); SPM spm4(_gf, m1, m3); SPM spm5(_gf, m1, m3); SPM spm6(_gf, m1, m3); 
	
	LinBox::SlicedPolynomialMatrixMulKaratuba<GField, SPM, SPM, SPM>()(_gf, spm3, spm1, spm2);
	file3<<"\nKaratsuba\n"; spm3.write(file3);
	LinBox::SlicedPolynomialMatrixMulKaratsubaParallel<GField, SPM, SPM, SPM>()(_gf, spm4, spm1, spm2);
	file3<<"\nKaratsuba parallel version\n"; spm3.write(file4);
	LinBox::SlicedPolynomialMatrixMulToomCook<GField, SPM, SPM>()(_gf, spm5, spm3, spm1);
	file3<<"\nToom-Cook\n"; spm3.write(file5);
	LinBox::SlicedPolynomialMatrixMulToomCookParallel<GField, SPM, SPM>()(_gf, spm6, spm3, spm1);
	file3<<"\nToom-Cook parallel version\n"; spm3.write(file6);
}
