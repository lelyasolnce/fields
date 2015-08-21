#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixAddSub_H
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixAddSub_H

#include <vector>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include "linbox/matrix/SlicedPolynomialMatrix/SlicedPolynomialMatrix.h"

namespace LinBox
{
	/* 
	 * C = A + B
	 */
	template< class GField, class Operand1, class Operand2, class Operand3>
	class SlicedPolynomialMatrixAdd
	{
	private:
		typedef typename Operand1::MatrixElement ME1;
		typedef typename Operand1::IntField IF1;
		typedef BlasMatrix<IF1, std::vector<ME1>> BM1;
		typedef typename Operand2::MatrixElement ME2;
		typedef typename Operand2::IntField IF2;
		typedef BlasMatrix<IF2, std::vector<ME2>> BM2;
		typedef typename Operand3::MatrixElement ME3;
		typedef typename Operand3::IntField IF3;
		typedef BlasMatrix<IF3, std::vector<ME3>> BM3;
	public:
		Operand1 &operator() (GField &F, Operand1 &C, Operand2 &A, Operand3 &B);
	};

	/* 
	 * C = A - B
	 */
	template< class GField, class Operand1, class Operand2, class Operand3>
	class SlicedPolynomialMatrixSub
	{
	private:
		typedef typename Operand1::MatrixElement ME1;
		typedef typename Operand1::IntField IF1;
		typedef BlasMatrix<IF1, std::vector<ME1>> BM1;
		typedef typename Operand2::MatrixElement ME2;
		typedef typename Operand2::IntField IF2;
		typedef BlasMatrix<IF2, std::vector<ME2>> BM2;
		typedef typename Operand3::MatrixElement ME3;
		typedef typename Operand3::IntField IF3;
		typedef BlasMatrix<IF3, std::vector<ME3>> BM3;
	public:
		Operand1 &operator() (GField &F, Operand1 &C, Operand2 &A, Operand3 &B);
	};

	/*
	 * C += A
	 */
	template< class GField, class Operand1, class Operand2>
	class SlicedPolynomialMatrixAddin
	{
	private:
		typedef typename Operand1::MatrixElement ME1;
		typedef typename Operand1::IntField IF1;
		typedef BlasMatrix<IF1, std::vector<ME1>> BM1;
		typedef typename Operand2::MatrixElement ME2;
		typedef typename Operand2::IntField IF2;
		typedef BlasMatrix<IF2, std::vector<ME2>> BM2;
	public:
		Operand1 &operator() (GField &F, Operand1 &C, Operand2 &A);
	};

	/*
	 * C -= A
	 */
	template< class GField, class Operand1, class Operand2>
	class SlicedPolynomialMatrixSubin
	{
	private:
		typedef typename Operand1::MatrixElement ME1;
		typedef typename Operand1::IntField IF1;
		typedef BlasMatrix<IF1, std::vector<ME1>> BM1;
		typedef typename Operand2::MatrixElement ME2;
		typedef typename Operand2::IntField IF2;
		typedef BlasMatrix<IF2, std::vector<ME2>> BM2;
	public:
		Operand1 &operator() (GField &F, Operand1 &C, Operand2 &A);
	};
}

#include "SlicedPolynomialMatrixAddSub.inl"

#endif

