#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrix_H
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrix_H

#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdint.h>
#include "linbox/matrix/DenseMatrix/blas-matrix.h"
#include <givaro/modular.h>
#include <givaro/givpoly1dense.h>
#include <givaro/givpoly1denseops.inl>

namespace LinBox
{
	template <class _Field, class _MatrixElement = double>
	class SlicedPolynomialMatrix
	/*Requirement: variable class _Field must be a representation of a field of order p^e
	 *and have functions that return p (characteristic()) and e (exponent()).
	 *For example, Givaro::GfqDom.
	 */

    /*Description: this class implements representation of a matrix over GF
	 *as a vector of length n of matrices over F (BlasMatrices),
	 *which corresponds to representation of this matrix as a polynomial of degree e-1 with matrix-coefficients.
	 */
	{
	public:
		typedef _Field Field;
		typedef typename Field::Element Element;
		typedef _MatrixElement MatrixElement;
		typedef std::vector<MatrixElement> Rep;
		typedef Givaro::Modular<MatrixElement> IntField;
		typedef SlicedPolynomialMatrix<Field, MatrixElement> Self_t;
		typedef typename Givaro::Poly1Dom<IntField, Dense>::Rep polynomial;
	private:
		_Field *GF;
		IntField F;
	private:
		size_t e; // GF.cardinality == p^e
		std::vector<BlasMatrix<IntField, Rep>> V;
		typedef typename Field::Residu_t GFqDompolynomial;
	public:
		polynomial irreducible()
		{
			polynomial res;
			size_t p = GF->characteristic();
			GFqDompolynomial z = GF->irreducible();
			for (size_t i = 0; i <= e; i++)
			{
				res.push_back(z % p);
				z = z / p;
			}
			return res;
		}

		        			////////////////
		        			//Constructors//
		        			////////////////

	public:
		/*! Allocates a vector of new \f$ 0 \times 0\f$ matrices.
		 */
		SlicedPolynomialMatrix (_Field &BF);

		/*Allocates a vector of $ m1 \times m2\f$ zero matrices.
		 */
		SlicedPolynomialMatrix (_Field &BF, size_t & m1, size_t &m2);

		        			///////////////
		        			// Destructor//
		        			///////////////

	public:
		~SlicedPolynomialMatrix() {}

		        			////////////////////////
	          				//dimensions of vector//
		        			////////////////////////

	public:
		/*Get length of V.
		 * @returns length of V
		 */
		size_t length();
		
		/*Get the number of rows in a matrix.
		 * @returns Number of rows in a matrix
		 */
		size_t rowdim();

		/* Get the number of columns in a matrix.
		 * @returns Number of columns in a matrix
		 */
		size_t coldim();
		
		/* Get the stride of the matrix.
		 */
		size_t getStride()
		{
			return V[0].getStride();
		}
		        			/////////////////
		        			//return fields//
		        			/////////////////
	
	public:
		Field& fieldGF();

		IntField& fieldF();

		        			/////////////////////////
		        			//functions for entries//
		        			/////////////////////////

	public:
		/* Set the entry of the m-th matrix-coefficient at the (i, j) position to a_mij.
		 * @param m matrix-coefficient number, 0...length() - 1
		 * @param i Row number, 0...rowdim () - 1
		 * @param j Column number 0...coldim () - 1
		 * @param a_mij Element to set
		 */
		void setEntry (size_t m, size_t i, size_t j, _MatrixElement &a_mij);
		
	private:
		/* Get a writeable reference to the m-th matrix-coefficient at the (i, j) position.
		 * @param m matrix-coefficient number, 0...length() - 1
		 * @param i Row index of entry
		 * @param j Column index of entry
		 * @returns Reference to matrix entry
		 */
		_MatrixElement &refEntry (size_t m, size_t i, size_t j);

	public:
		/* Get a read-only reference to the m-th matrix-coefficient at the (i, j) position.
		 * @param m matrix-coefficient number, 0...length() - 1
		 * @param i Row index
		 * @param j Column index
		 * @returns Const reference to matrix entry
		 */
		_MatrixElement &getEntry (size_t m, size_t i, size_t j);
		
		        			/////////////////////////////////////
		        			//functions for matrix-coefficients//
		        			/////////////////////////////////////
						
	public:
		/* Set the m-th matrix-coefficient to V_m.
		 * @param m matrix-coefficient number, 0...length() - 1
		 * @param V_m matrix to set
		 */
		void setMatrixCoefficient (size_t m, BlasMatrix<IntField, Rep> &V_m) ;

	private:
		/* Get a writeable reference to the m-th matrix-coefficient.
		 * @param m matrix-coefficient number, 0...length() - 1
		 * @returns Reference to matrix-coefficent
		 */
		BlasMatrix<IntField, Rep> &refMatrixCoefficient (size_t m) ;

	public:
		/** Get a read-only reference to the m-th matrix-coefficient
		 * @param m matrix-coefficient number, 0...length() - 1
		 * @returns Const reference to matrix-coefficent
		 */
		BlasMatrix<IntField, Rep> &getMatrixCoefficient (size_t m);

		        			/////////
		        			//swaps//
		        			/////////

	public:
		/* Swap i1-th and i2-th rows of matrices.
		 * This is done inplace.
		 */
		void swapRows(size_t i1, size_t i2);

		/* Swap j1-th and j2-th columns of matrices.
		 * This is done inplace.
		 */
		void swapCols(size_t j1, size_t j2);
		
		        			/////////////
		        			//transpose//
		        			/////////////

	public:
		/* Creates a transposed polynomial matrix of this.
		 * @param[in] tV
		 * @return the transposed polynomial matrix of this.
		 */
		Self_t transpose(Self_t & tV);
		
		        			//////////////////
		        			//input / output//
		        			//////////////////
	
	public:
		std::istream &read (std::istream &file);

		std::ostream &write (std::ostream &os)
		{
			int K = this->length();
			int I = this->rowdim();
			int J = this->coldim();
			for (int k = 0; k < K; k++)
			{
				for (int i = 0; i < I; i++)
				{
					for (int j = 0; j < J; j++)
					{
						os << this->getEntry(k, i, j) << " ";
					}
					os << std::endl;
				}
				os << std::endl;
			}
			return os;
		}
	};
}

#include "SlicedPolynomialMatrix.inl"

#endif

