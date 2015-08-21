#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialSubmatrix_H
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialSubmatrix_H
namespace Linbox
{
	template<class ParentType>
	class SlicedPolynomialSubmatrix
	{
	public:
		typedef typename SPM::Field           Field;
		typedef typename Field::Element         Element;    //!< Element type
		typedef typename SPM::Rep               Rep;    //!< Actually a <code>std::vector<Element></code> (or alike.)
		typedef SlicedPolynomialSubmatrix<typename SPM::Self_t>              Self_t;         //!< Self type
		typedef typename SPM::Element MatrixElement;
		typedef typename SPM::IntField IntField;
		typedef typename Givaro::Poly1Dom<IntField, Dense>::Rep polynomial;
		bool windowed;
		bool is_window;
	private:
		std::vector<BlasSubmatrix<IntField>> V;
		SPM &spm; // main SPM matrix
		ParentType &pm; // parent matrix    
		size_t row;                   //!< row dimension of Submatrix
		size_t col;                   //!< col dimension of Submatrix
		size_t _r0;                    //!< upper left corner row of Submatrix in \p _Mat
		size_t _c0;                    //!< upper left corner row of Submatrix in \p _Mat
		size_t _stride ;               //!< number of columns in \p _Mat (or stride of \p _Mat)
	
		        			////////////////
		        			//Constructors//
		        			////////////////

	public:
		/* Constructor from an existing SlicedPolynomialMatrix and dimensions.
		 * @ M is SlicedPolynomialMatrix of which to construct submatrix
		 * @ rowbeg is starting row
		 * @ colbeg is starting column
		 * @ Rowdim is row dimension
		 * @ Coldim is column dimension
		 * @ is_w is bool "submatrix created as a window"
		 * A pointer to the elements of the parent matrix are created, doesn't consume memory
		 */
		SlicedPolynomialSubmatrix (SPM &M, size_t rowbeg, size_t colbeg, size_t Rowdim, size_t Coldim, bool is_w = false);
		
		/* Constructor from an existing SlicedPolynomialMatrix
		 * @ M is SlicedPolynomialMatrix of which to construct submatrix
		 * @ is_w is bool "submatrix created as a window"
		 * A pointer to the elements of the parent matrix are created, doesn't consume memory
		 */
		SlicedPolynomialSubmatrix (SPM &M, bool is_w = false);

		/** Constructor from an existing submatrix and dimensions
		 * @ SM is SlicedPolynomialSubmatrix from which to construct submatrix
		 * @ rowbeg is starting row
		 * @ colbeg is starting column
		 * @ Rowdim is row dimension
		 * @ Coldim is column dimension
		 * @ is_w is bool "submatrix created as a window"
		 * A pointer to the elements of the parent matrix are created, doesn't consume memory
		 */
		SlicedPolynomialSubmatrix (Self_t  &SM, size_t rowbeg, size_t colbeg, size_t Rowdim, size_t Coldim, bool is_w = false);

		/* Copy constructor.
		 * @ SM is SlicedPolynomialSubmatrix to copy
		 * @ is_w is bool "submatrix created as a window"
		 * A pointer to the elements of the parent matrix are created, doesn't consume memory
		 */
		SlicedPolynomialSubmatrix (Self_t &SM, bool is_w = false);

	        			    ///////////////
		        			// Destructor//
		        			///////////////

	public:
		~SlicedPolynomialSubmatrix(){}

	        			    //////////////////
		        			// Parent matrix//
		        			//////////////////

	public:
		SPM& parentMatrix();
		
		        			/////////
		        			//swaps//
		        			/////////

	public:
		void swapRows(size_t i1, size_t i2, bool bywindow = false);
		void swapCols(size_t j1, size_t j2, bool bywindow = false);
		
		        			/////////////
		        			//transpose//
		        			/////////////

	public:
		Self_t transpose(Self_t & tV, bool bywindow = false) const;

							//////////////////
							//  DIMENSIONS  //
							//////////////////
	public:
		size_t length() const;
		size_t rowdim () const;
		size_t coldim () const ;
		size_t getStride() const;
		size_t stride() const;
		size_t offset() const;

							///////////////////
							//      I/O      //
							///////////////////

				                					//input / output//
	public:
		std::istream &read (std::istream &file);
		std::ostream &write (std::ostream &os) const
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


		        			/////////////////////////
		        			//functions for entries//
		        			/////////////////////////


	public:
		void setEntry (size_t m, size_t i, size_t j, const MatrixElement &a_mij, bool bywindow = false);
		
	private:
		MatrixElement &refEntry (size_t m, size_t i, size_t j, bool bywindow = false);

	public:
		MatrixElement &getEntry (size_t m, size_t i, size_t j, bool bywindow = false);
		
		        			/////////////////////////////////////
		        			//functions for matrix-coefficients//
		        			/////////////////////////////////////		
	public:
		void setMatrixCoefficient (size_t m, const BlasSubmatrix<IntField> &V_m, bool bywindow = false) ;

	private:
		BlasSubmatrix<IntField> &refMatrixCoefficient (size_t m, bool bywindow = false) ;

	public:
		const BlasSubmatrix<IntField> &getMatrixCoefficient (size_t m, bool bywindow = false) const ;

		        			/////////////////
		        			//return fields//
		        			/////////////////
	public:
		const Field& fieldGF() const;
		const IntField& fieldF() const;
	};
}
#endif