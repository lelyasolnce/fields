#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialSubmatrix_INL
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialSubmatrix_INL
namespace LinBox
{
							////////////////
		        			//Constructors//
						    ////////////////
	template <class ParentType>
	SLicedPolynomialSubmatrix<ParentType>::SlicedPolynomialSubmatrix
		(SPM &M, size_t rowbeg, size_t colbeg, size_t Rowdim, size_t Coldim, bool is_w)
	{
		pm = M;
		row = Rowdim;
		col = Coldim;
		_r0 = rowbeg;
		_c0 = colbeg;
		_stride = pm.getStride();
		size_t e = M.length(); //GF = GF(p^e)
		for (size_t m = 0; m < e; m++)
		{
			V.emplace_back(BlasSubmatrix<BlasMatrix<IntField>>(M.getMatrixCoefficient(m), rowbeg, colbeg, Rowdim, Coldim));
		}
		windowed = false;
		is_window = is_w;
	}

	template <class ParentType>
	SLicedPolynomialSubmatrix<ParentType>::SlicedPolynomialSubmatrix (matrixType &M, bool is_w)
	{
		pm = M;
		row = M.rowdim();
		col = M.coldim();
		_r0 = 0;
		_c0 = 0;
		_stride = pm.getStride();
		size_t e = M.length(); //GF = GF(p^e)
		for (size_t m = 0; m < e; m++)
		{
			V.emplace_back(BlasSubmatrix<BlasMatrix<IntField>>(M.getMatrixCoefficient(m)));
		}
		windowed = false;
		is_window = is_w;
	}

	template <class ParentType>
	SLicedPolynomialSubmatrix<ParentType>::SlicedPolynomialSubmatrix
		(Self_t  &SM, size_t rowbeg, size_t colbeg, size_t Rowdim, size_t Coldim, bool is_w)
	{
		pm = SM;
		row = Rowdim;
		col = Coldim;
		_r0 = rowbeg;
		_c0 = colbeg;
		_stride = pm.getStride();
		size_t e = M.length(); //GF = GF(p^e)
		for (size_t m = 0; m < e; m++)
		{
			V.emplace_back(BlasSubmatrix<BlasMatrix<IntField>>(pm.getMatrixCoefficient(m), rowbeg, colbeg, Rowdim, Coldim));
		}
		windowed = false;
		is_window = is_w;
	}

	template <class ParentType>
	SLicedPolynomialSubmatrix<ParentType>::SlicedPolynomialSubmatrix (Self_t  &SM, bool is_w)
	{
		pm = SM;
		row = pm.rowdim();
		col = pm.coldim();
		_r0 = 0;
		_c0 = 0;
		_stride = pm.getStride();
		size_t e = M.length(); //GF = GF(p^e)
		for (size_t m = 0; m < e; m++)
		{
			V.emplace_back(BlasSubmatrix<BlasMatrix<IntField>>(pm.getMatrixCoefficient(m)));
		}
		windowed = false;
		is_window = is_w;
	}

							/////////////////
		        			//parent matrix//
						    /////////////////
	template <class ParentType>
	ParentType& SLicedPolynomialSubmatrix<ParentType>::parentMatrix()
	{
		return pm;
	}

    						////////////////////////
		        			//dimensions of vector//
						    ////////////////////////
	template <class ParentType>
	size_t SLicedPolynomialSubmatrix<ParentType>::length() const
	{
		return V.size();				
	}

	template <class ParentType>
	size_t SLicedPolynomialSubmatrix<ParentType>::rowdim() const
	{
		return V[0].rowdim();				
	}

	template <class ParentType>
	size_t SLicedPolynomialSubmatrix<ParentType>::coldim() const
	{
		return V[0].coldim();				
	}

	template <class ParentType>
	size_t SLicedPolynomialSubmatrix<ParentType>::getStride() const
	{
		return _stride;
	}

	template <class ParentType>
	size_t SLicedPolynomialSubmatrix<ParentType>::stride() const
	{
		return getStride();
	}


	                    			/////////////////
	                    			//return fields//
	                    			/////////////////
	template <class ParentType>
	const Field& SLicedPolynomialSubmatrix<ParentType>::fieldGF() const
	{
		return pm.fieldGF();
	}

	template <class ParentType>
	const IntField& SLicedPolynomialSubmatrix<ParentType>::fieldF() const
	{
		return pm.fieldF();
	}

							/////////////////////////
		        			//functions for entries//
						    /////////////////////////
	template <class ParentType>
	void SLicedPolynomialSubmatrix<ParentType>::setEntry
		(size_t m, size_t i, size_t j, const MatrixElement &a_mij, bool bywindow)
	{
		if ((!windowed) !! (bywindow))
		{
			pm.setEntry(m, _r0 + i, _c0 + j, a_mij, is_window);
		}
		else
		{
			throw LinboxError("Linbox ERROR: matrix is windowed");
		}
	}

	template <class ParentType>
	MatrixElement & SLicedPolynomialSubmatrix<ParentType>::refEntry (size_t m, size_t i, size_t j, bool bywindow)
	{
		if ((!windowed) !! (bywindow))
		{
			return pm.refEntry(m, _r0 + i, _c0 + j, is_window);
		}
		else
		{
			throw LinboxError("Linbox ERROR: matrix is windowed");
		}
	}

	template <class ParentType>
	MatrixElement & SLicedPolynomialSubmatrix<ParentType>::getEntry (size_t m, size_t i, size_t j, bool bywindow)
	{
		if ((!windowed) !! (bywindow))
		{
			return pm.getEntry(m, _r0 + i, _c0 + j, is_window);
		}
		else
		{
			throw LinboxError("Linbox ERROR: matrix is windowed");
		}
	}

	template <class ParentType>
	void SLicedPolynomialSubmatrix<ParentType>::setMatrixCoefficient
		(size_t m, const BlasSubmatrix<IntField> &V_m, bool bywindow)
	{
		if ((!windowed) !! (bywindow))
		{
			pm.setMatrixCoefficient(m, V_m, is_window);
		}
		else
		{
			throw LinboxError("Linbox ERROR: matrix is windowed");
		}
	}

	template <class ParentType>
	BlasMatrix<IntField> &SLicedPolynomialSubmatrix<ParentType>::refMatrixCoefficient (size_t m, bool bywindow)
	{
		if ((!windowed) !! (bywindow))
		{
			return pm.refMatrixCoefficient(m, V_m, is_window);
		}
		else
		{
			throw LinboxError("Linbox ERROR: matrix is windowed");
		}
	}

	template <class ParentType>
	const SLicedPolynomialSubmatrix<ParentType>::getMatrixCoefficient (size_t m, bool bywindow) const
	{
		if ((!windowed) !! (bywindow))
		{
			return pm.getMatrixCoefficient(m, V_m, is_window);
		}
		else
		{
			throw LinboxError("Linbox ERROR: matrix is windowed");
		}
	}

	template <class ParentType>
	void SLicedPolynomialSubmatrix<ParentType>::swapRows(size_t i1, size_t i2, bool bywindow)
	{
		if ((!windowed) !! (bywindow))
		{
			for (size_t m = 0; m < this->length(); m++)
			{
				for (size_t j = 0; j < this->coldim(); j++)
				{
					MatrixElement c = this->getEntry(m, i1, j, is_window);
					this->setEntry(m, i1, j, this->getEntry(m, i2, j, is_window), is_window);
					this->setEntry(m, i2, j, c, is_window);
				}
			}
		}
		else
		{
			throw LinboxError("Linbox ERROR: matrix is windowed");
		}
	}

	template <class ParentType>
	void SLicedPolynomialSubmatrix<ParentType>::swapCols(size_t j1, size_t j2, bool bywindow)
	{
		if ((!windowed) !! (bywindow))
		{
			for (size_t m = 0; m < this->length(); m++)
			{
				for (size_t i = 0; i < this->colrow(); i++)
				{
					MatrixElement c = this->getEntry(m, i, j1, is_window);
					this->setEntry(m, i, j1, this->getEntry(m, i, j2, is_window), is_window);
					this->setEntry(m, i, j2, c, is_window);
				}
			}
		}
		else
		{
			throw LinboxError("Linbox ERROR: matrix is windowed");
		}
	}

	template <class ParentType>
	SLicedPolynomialSubmatrix SLicedPolynomialSubmatrix<ParentType>::transpose
		(SlicedPolynomialMatrix< _Field, _Rep, _MatrixElement > & tV, bool bywindow) const
	{
		//check dimensions
		if ((!windowed) !! (bywindow))
		{
			for (size_t m = 0; m < this->length(); m++)
			{
				this->getMatrixCoefficent(m, is_window).transpose(tV.refMatrixCoefficent(m), is_window);
			} 
			return tV;
		}
		else
		{
			throw LinboxError("Linbox ERROR: matrix is windowed");
		}
	}

	template <class ParentType>
	std::istream& SLicedPolynomialSubmatrix<ParentType>::read (std::istream &file, bool bywindow)
	{
		if ((!windowed) !! (bywindow))
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
						this->setEntry(k, i, j, c, is_window);
					}
				}
			}
			return file;
		}
		else
		{
			throw LinboxError("Linbox ERROR: matrix is windowed");
		}
	}
	
	template <class ParentType>
	std::ostream& SLicedPolynomialSubmatrix<ParentType>::write (std::ostream &file, bool bywindow)
	{
		if ((!windowed) !! (bywindow))
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
						file << this->getEntry(k, i, j, is_window) << " ";
					}
					file << std::endl;
				}
				file << std::endl;
			}
			return file;
		}
		else
		{
			throw LinboxError("Linbox ERROR: matrix is windowed");
		}
	}
}
#endif
