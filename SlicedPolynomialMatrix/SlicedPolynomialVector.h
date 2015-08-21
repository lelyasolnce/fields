#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialVector_H
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialVector_H

#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdint.h>
#include "linbox/vector/blas-vector.h"
#include <givaro/modular.h>
#include <givaro/givpoly1dense.h>
#include <givaro/givpoly1denseops.inl>

namespace LinBox
{
	template <class _Field, class _VectorElement = double>
	class SlicedPolynomialVector
	/*Requirement: variable class _Field must be a representation of a field of order p^e
	 *and have functions that return p (characteristic()) and e (exponent()).
	 *For example, Givaro::GfqDom.
	 */

    /*Description: this class implements representation of a vector over GF
	 *as a vector of length n of vectors over F (BlasVectors),
	 *which corresponds to representation of this vector as a polynomial of degree n-1 with vector-coefficients.
	 */
	{
	public:
		typedef _Field Field;
		typedef typename Field::Element Element;
		typedef std::vector<MatrixElement> Rep;
		typedef _VectorElement VectorElement;
		typedef typename Givaro::Modular<VectorElement> IntField;
		typedef typename SlicedPolynomialVector<Field, VectorElement> Self_t;
		typedef typename Givaro::Poly1Dom<IntField, Dense>::Rep polynomial;
	private:
		_Field *GF;
		IntField F;
	private:
		size_t e; // GF.cardinality == p^e
		std::vector<BlasVector<IntField>> V;
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
		/*! Allocates a vector of new zero vectors of size 0.
		 */
		SlicedPolynomialVector (const _Field &BF);

		/*Allocates a vector of new vectors of size m.
		 */
		SlicedPolynomialVector (const _Field &BF, const size_t &m);

		        			///////////////
		        			// Destructor//
		        			///////////////

	public:
		~SlicedPolynomialVector(){}

		        			////////////////////////
	          				//dimensions of vector//
		        			////////////////////////

	public:
		/*Get length of V.
		 * @returns length of V
		 */
		size_t length() const;
		
		/*Get the number of rows in a vector.
		 * @returns Number of rows in a vector
		 */
		size_t rowdim() const;
		
	                    			/////////////////
	                    			//return fields//
	                    			/////////////////
	
	public:
		const _Field& fieldGF() const;

		const IntField& fieldF() const;

	                    			/////////////////////////
	                    			//functions for entries//
	                    			/////////////////////////

	public:
		/* Set the entry of the m-th vector-coefficient at the (k) position to a_mk.
		 * @param m vector-coefficient number, 0...length() - 1
		 * @param k Row number 0...rowdim () - 1
		 * @param a_mk Element to set
		 */
		void setEntry (size_t m, size_t k, const _MatrixElement &a_mk);
		
	private:
		/* Get a writeable reference to the m-th vector-coefficient at the (k) position.
		 * @param m vector-coefficient number, 0...length() - 1
		 * @param k Row number 0...rowdim () - 1
		 * @returns Reference to vector entry
		 */
		_MatrixElement &refEntry (size_t m, size_t k);

	public:
		/* Get a read-only reference to the m-th vector-coefficient at the (k) position.
		 * @param m vector-coefficient number, 0...length() - 1
		 * @param k Row number 0...rowdim () - 1
		 * @returns Const reference to vector entry
		 */
		_MatrixElement &getEntry (size_t m, size_t k);
		
	                    			/////////////////////////////////////
	                    			//functions for matrix-coefficients//
	                    			/////////////////////////////////////
						
	public:
		/* Set the m-th vector-coefficient to V_m.
		 * @param m vector-coefficient number, 0...length() - 1
		 * @param V_m matrix to set
		 */
		void setVectorCoefficient (size_t m, const BlasVector<IntField> &V_m) ;

	private:
		/* Get a writeable reference to the m-th matrix-coefficient.
		 * @param m matrix-coefficient number, 0...length() - 1
		 * @returns Reference to matrix-coefficent
		 */
		BlasVector<IntField> &refVectorCoefficient (size_t m) ;

	public:
		/** Get a read-only reference to the m-th matrix-coefficient
		 * @param m matrix-coefficient number, 0...length() - 1
		 * @returns Const reference to matrix-coefficent
		 */
		const BlasVector<IntField> &getVectorCoefficient (size_t m) const ;

	                    			/////////
	                    			//swaps//
	                    			/////////

	public:
		/* Swap i1-th and i2-th rows of matrices.
		 * This is done inplace.
		 */
		void swapRows(size_t k1, size_t k2);
		
	                    			//////////////////
	                    			//input / output//
	                    			//////////////////
	
	public:
		std::istream &read (std::istream &file);

		std::ostream &write (std::ostream &os)
		{
			int K = this->length();
			int I = this->rowdim();
			for (int k = 0; k < K; k++)
			{
				for (int i = 0; i < I; i++)
				{
					os << this->getEntry(k, i, j) << " ";
				}
				os << std::endl;
			}
			return os;
		}
	};
}

#include "SlicedPolynomialVector.inl"

#endif