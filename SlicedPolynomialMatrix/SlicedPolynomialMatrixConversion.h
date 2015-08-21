#ifndef __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixConversion_H
#define __LINBOX_matrix_SlicedPolynomialMatrix_SlicedPolynomialMatrixConversion_H

#include "SlicedPolynomialMatrix.h"
#include "linbox/matrix/polynomial-matrix.h"
#include "linbox/matrix/DenseMatrix/blas-matrix.h"
#include <givaro/gfq.h>

namespace LinBox
{
	//class PM is supposed to be PolynomialMatrix<PMType::polfirst, PMStorage::plain, _FieldF>
	template <class SPM, class PM>
	class SlicedPolynomialMatrixtoPolynomialMatrix
	{
	public:
		PM& operator() (SPM &spm, PM &pm)
		{
			size_t rd = spm.rowdim();
			size_t cd = spm.coldim();
			size_t l = spm.length();
			for (size_t i = 0; i < rd; i++)
			{
				for (size_t j = 0; j < cd; j++)
				{
					for (size_t k = 0; k < l; k++)
					{
						pm.ref(i, j, k) = spm.getEntry(k, i, j);
					}
				}
			}
			return PM;
		}
	};

	template <class PM, class SPM>
	class PolynomialMatrixtoSlicedPolynomialMatrix
	{
	public:
		SPM &operator() (PM& pm, SPM& spm)
		{
			size_t rd = SPM.rowdim();
			size_t cd = SPM.coldim();
			size_t l = SPM.length();
			for (size_t k = 0; k < l; k++)
			{
				SPM.setMatrixCoefficient(k, PM[k]);
			}
			return SPM;
		}
	};

	template <class _FieldGF, class _Storage1, class _MatrixElement = double, class TT, class _Storage2>
	template <class SPM, class BM>
	class SlicedPolynomialMatrixtoBlasMatrix
	{
	private:
		typedef BM::Element Element;
		typedef typename Givaro::GFqDom<Element> GFqDom;
		
	public:
		BM &operator() (SPM &spm, BM &bm)
		{
			size_t rd = spm.rowdim();
			size_t cd = spm.coldim();
			size_t l = spm.length();
			Element el;

			typename GFqDom GF(spm.fieldGF().characteristic(), spm.fieldGF().exponent());
			PolynomialMatrix<typename PMType::polfirst, typename PMStorage::plain,
				typename SPM::IntField> pm(spm.fieldF(), rd, cd, l);
			SlicedPolynomialMatrixtoPolynomialMatrix<_FieldGF, _Storage1,
			_MatrixElement, spm::IntField>()(pm, spm);
			for (size_t i = 0; i < rd; i++)
			{
				for (size_t j = 0; j < cd; j++)
				{
					GF.init(el, pm(i, j));
					bm.setEntry(i, j, el);
				}
			}
			return BM;
		}
	};

	template <class BM, class SPM>
	class BlasMatrixtoSlicedPolynomialMatrix
	{
	private:
		typedef BM::Field Field;
		typedef Givaro::GFqDom<TT> GField;
	public:
		SPM& operator() (BM &bm, SPM &spm)
		{
			size_t rd = spm.rowdim();
			size_t cd = spm.coldim();
			size_t l = spm.length();
			if (l == 1)
			{
				for (size_t i = 0; i < rd; i++)
				{
					for (size_t j = 0; j < cd; j++)
					{
						spm.setEntry(0, i, j, BM.getEntry(i, j));
					}
				}
				return spm;
			}
			else
			{
				Field GF = spm.fieldDF();
				size_t p = GF.characteristic();				
				for (size_t i = 0; i < rd; i++)
				{
					for (size_t j = 0; j < cd; j++)
					{
						size_t el = GF.zech2padic(bm.getEntry(i, j));
						for (size_t k = 0; k < l; k++)
						{
							spm.setEntry(k, i, j, el % p);
							el /= p;
						}
					}
				}
				return SPM;
			}
		}
	};
}

#endif