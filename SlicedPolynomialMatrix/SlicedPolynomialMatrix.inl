#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrix_INL
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrix_INL

#include <stdlib.h>
#include <stdint.h>

namespace LinBox
{
		        			////////////////
		        			//Constructors//
		        			////////////////

	template < class _Field, class _MatrixElement >
	SlicedPolynomialMatrix< _Field, _MatrixElement >::SlicedPolynomialMatrix (_Field &BF)
	{
		GF = &BF;
		IntField F_temp(GF->characteristic());
		F = F_temp;
		e = GF->exponent(); //GF = GF(p^e)
		for (size_t m = 0; m < e; m++)
		{
			V.emplace_back(BlasMatrix<IntField, Rep>(F));
		}
	}

	template < class _Field, class _MatrixElement >
	SlicedPolynomialMatrix< _Field, _MatrixElement >::SlicedPolynomialMatrix (_Field &BF, size_t & m1, size_t &m2)
	{
		GF = &BF;
		IntField F_temp(GF->characteristic());
		F = F_temp;
		e = GF->exponent(); //GF = GF(p^e)
		for (size_t m = 0; m < e; m++)
		{
			V.emplace_back(BlasMatrix<IntField, Rep>(F, m1, m2));
		}
	}

		        			////////////////////////
		        			//dimensions of vector//
		        			////////////////////////

    template < class _Field, class _MatrixElement >
	size_t SlicedPolynomialMatrix< _Field, _MatrixElement >::length()
	{
		return V.size();				
	}

	template < class _Field, class _MatrixElement >
	size_t SlicedPolynomialMatrix< _Field, _MatrixElement >::rowdim()
	{
		return V[0].rowdim();				
	}

	template < class _Field, class _MatrixElement >
	size_t SlicedPolynomialMatrix< _Field, _MatrixElement >::coldim()
	{
		return V[0].coldim();				
	}
	
		        			/////////////////
		        			//return fields//
		        			/////////////////

	template < class _Field, class _MatrixElement >
	typename SlicedPolynomialMatrix< _Field, _MatrixElement >::Field&
	SlicedPolynomialMatrix< _Field, _MatrixElement >::fieldGF()
	{
		return *GF;
	}

	template < class _Field, class _MatrixElement >
	typename SlicedPolynomialMatrix< _Field, _MatrixElement >::IntField&
	SlicedPolynomialMatrix< _Field, _MatrixElement >::fieldF()
	{
		return F;
	}

		        			/////////////////////////
		        			//functions for entries//
		        			/////////////////////////
		
    template < class _Field, class _MatrixElement >
	void SlicedPolynomialMatrix< _Field,  _MatrixElement >::setEntry (size_t m, size_t i, size_t j, _MatrixElement &a_mij)
	{
		V[m].setEntry(i, j, a_mij);
	}

	template < class _Field, class _MatrixElement >
	_MatrixElement & SlicedPolynomialMatrix< _Field,  _MatrixElement >::refEntry (size_t m, size_t i, size_t j)
	{
		return V[m].refEntry(i, j);

	}

	template < class _Field, class _MatrixElement >
	_MatrixElement& SlicedPolynomialMatrix< _Field,  _MatrixElement >::getEntry (size_t m, size_t i, size_t j)
	{
		return V[m].getEntry(i, j);

	}
	
		        			/////////////////////////////////////
		        			//functions for matrix-coefficients//
		        			/////////////////////////////////////

	template < class _Field, class _MatrixElement >
	void SlicedPolynomialMatrix< _Field,  _MatrixElement >::setMatrixCoefficient (size_t m, BlasMatrix<IntField, Rep> &V_m)
	{
		V[m] = V_m;
	}

	template < class _Field, class _MatrixElement >
	BlasMatrix<Givaro::Modular<_MatrixElement>, std::vector<_MatrixElement>>&
	SlicedPolynomialMatrix< _Field,  _MatrixElement >::refMatrixCoefficient (size_t m)
	{
		return V[m];
	}

	template < class _Field, class _MatrixElement >
	BlasMatrix<Givaro::Modular<_MatrixElement>, std::vector<_MatrixElement>>&
	SlicedPolynomialMatrix< _Field,  _MatrixElement >::getMatrixCoefficient (size_t m)
	{
		return V[m];
	}

		        			/////////
		        			//swaps//
		        			/////////

	template < class _Field, class _MatrixElement >
	void SlicedPolynomialMatrix< _Field,  _MatrixElement >::swapRows(size_t i1, size_t i2)
	{
		for (size_t m = 0; m < this->length(); m++)
		{
			for (size_t j = 0; j < this->coldim(); j++)
			{
				MatrixElement c = this->getEntry(m, i1, j);
				this->setEntry(m, i1, j, this->getEntry(m, i2, j));
				this->setEntry(m, i2, j, c);
			}
		}
	}

	template < class _Field, class _MatrixElement >
	void SlicedPolynomialMatrix< _Field,  _MatrixElement >::swapCols(size_t j1, size_t j2)
	{
		for (size_t m = 0; m < this->length(); m++)
		{
			for (size_t i = 0; i < this->rowdim(); i++)
			{
				MatrixElement c = this->getEntry(m, i, j1);
				this->setEntry(m, i, j1, this->getEntry(m, i, j2));
				this->setEntry(m, i, j2, c);
			}
		}
	}
	
		        			/////////////
		        			//transpose//
		        			/////////////

	template < class _Field, class _MatrixElement >
	SlicedPolynomialMatrix< _Field,  _MatrixElement > SlicedPolynomialMatrix< _Field,  _MatrixElement >::transpose(SlicedPolynomialMatrix< _Field,  _MatrixElement > & tV)
	{
		for (size_t m = 0; m < this->length(); m++)
		{
			this->getMatrixCoefficent(m).transpose(tV.refMatrixCoefficent(m));
		} 
		return tV;
	}

		        			//////////////////
		        			//input / output//
		        			//////////////////	
	template < class _Field, class _MatrixElement >
	std::istream& SlicedPolynomialMatrix< _Field,  _MatrixElement >::read (std::istream &file)
	{
		int K = this->length();
		int I = this->rowdim();
		int J = this->coldim();
		MatrixElement c;
		for (int k = 0; k < K; k++)
		{
			for (int i = 0; i < I; i++)
			{
				for (int j = 0; j < J; j++)
				{
					file >> c;
					this->setEntry(k, i, j, c);
				}
			}
		}
		return file;
	}
	
}

#endif