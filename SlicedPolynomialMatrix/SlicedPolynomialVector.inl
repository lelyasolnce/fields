#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialVector_INL
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialVector_INL

#include <linbox/vector/blas-vector.h>
#include <stdlib.h>
#include <stdint.h>

namespace Linbox
{
						////////////////
		        			//Constructors//
						////////////////

	template < class _Field, class _VectorElement >
	SlicedPolynomialVector< _Field, _VectorElement >::SlicedPolynomialVector (const _Field &BF)
	{
		GF = &BF;
		IntField F_temp(GF->characteristic());
		F = F_temp;
		e = GF->exponent(); //GF = GF(p^e)
		for (size_t m = 0; m < e; m++)
		{
			V.emplace_back(BlasVector<IntField>(F));
		}
	}

	template < class _Field, class _VectorElement >
	SlicedPolynomialVector< _Field, _VectorElement >::SlicedPolynomialVector (const _Field &BF, const size_t &m)
	{
		GF = &BF;
		IntField F_temp(GF->characteristic());
		F = F_temp;
		e = GF->exponent(); //GF = GF(p^e)
		for (size_t m = 0; m < e; m++)
		{
			V.emplace_back(BlasVector<IntField>(F, m));
		}
	}

						////////////////////////
		        			//dimensions of vector//
						////////////////////////

        template < class _Field, class _VectorElement >
	size_t SlicedPolynomialVector< _Field, _VectorElement >::length() const
	{
		return V.size();				
	}

	template < class _Field, class _VectorElement >
	size_t SlicedPolynomialVector< _Field, _VectorElement >::rowdim() const
	{
		return V[0].size();				
	}
	
	                    			/////////////////
	                    			//return fields//
	                    			/////////////////

	template < class _Field, class _VectorElement >
	const Field& SlicedPolynomialVector< _Field, _VectorElement >::fieldGF() const
	{
		return *GF;
	}

	template < class _Field, class _VectorElement >
	const IntField& SlicedPolynomialVector< _Field, _VectorElement >::fieldF() const
	{
		return F;
	}

						/////////////////////////
		        			//functions for entries//
						/////////////////////////
		
    template < class _Field, class _VectorElement >
	void SlicedPolynomialVector< _Field, _VectorElement >::setEntry (size_t m, size_t k, const _VectorElement &a_mk)
	{
		V[m].setEntry(k, a_mk);
	}

	template < class _Field, class _VectorElement >
	_VectorElement& SlicedPolynomialVector< _Field, _VectorElement >::refEntry (size_t m, size_t k)
	{
		return V[m].refEntry(k);

	}

	template < class _Field, class _VectorElement >
	_VectorElement & SlicedPolynomialVector< _Field, _VectorElement >::getEntry (size_t m, size_t k)
	{
		return V[m].getEntry(k);

	}
	
						/////////////////////////////////////
		                		//functions for matrix-coefficients//
						/////////////////////////////////////

	template < class _Field, class _VectorElement >
	void SlicedPolynomialVector< _Field, _VectorElement >::setMatrixCoefficient (size_t m, const BlasVector<IntField> &V_m)
	{
		V[m] = V_m;
	}

	template < class _Field, class _VectorElement >
	BlasVector<IntField> &SlicedPolynomialVector< _Field, _VectorElement >::refMatrixCoefficient (size_t m)
	{
		return V[m];
	}

	template < class _Field, class _VectorElement >
	const BlasVector<IntField> &SlicedPolynomialVector< _Field, _VectorElement >::getMatrixCoefficient (size_t m) const
	{
		return V[m];
	}

						/////////
		                		//swaps//
						/////////

	template < class _Field, class _VectorElement >
	void SlicedPolynomialVector< _Field, _VectorElement >::swapRows(size_t k1, size_t k2)
	{
		for (size_t m = 0; m < this->length(); m++)
		{
			VectorElement c = this->getEntry(m, k1);
			this->setEntry(m, k1, this->getEntry(m, k2));
			this->setEntry(m, k2, c);
		}
	}
	
						//////////////////
		                		//input / output//
						//////////////////	
	template < class _Field, class _VectorElement >
	std::istream& SlicedPolynomialVector< _Field, _VectorElement >::read (std::istream &file)
	{
		int M = this->length();
		int K = this->rowdim();
		VectorElement c;
		for (int m = 0; m < M; m++)
		{
			for (int k = 0; k < K; k++)
			{
				file >> c;
				this->setEntry(m, k, c);
			}
		}
		return file;
	}
}

#endif