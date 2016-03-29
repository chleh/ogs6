/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PHYSICALQUANITIES_SOLIDS_LINEARELASTICISOTROPIC_H_
#define PHYSICALQUANITIES_SOLIDS_LINEARELASTICISOTROPIC_H_

namespace Solids
{
namespace LinearElasticIsotropic
{
template <typename KelvinVectorType, typename KelvinMatrixType>
typename std::enable_if<KelvinVectorType::RowsAtCompileTime == 6 ||
                            KelvinVectorType::RowsAtCompileTime == 4,
                        void>::type
computeConstitutiveRelation(double const /*dt*/,
                            double const lambda,
                            double const mu,
                            KelvinVectorType const& eps_prev,
                            KelvinVectorType const& eps,
                            KelvinVectorType const& sigma_prev,
                            KelvinVectorType& sigma,
                            KelvinMatrixType& C)
{
	C.setZero();

	// C.topLeftCorner<3, 3>() = lambda;
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			C(i, j) = lambda;
	// C.diagonal() += 2 * mu;
	for (int i = 0; i < KelvinVectorType::RowsAtCompileTime; ++i)
		C(i, i) += 2 * mu;

	sigma = sigma_prev + C * (eps - eps_prev);
}

template <typename KelvinVectorType, typename KelvinMatrixType>
typename std::enable_if<KelvinVectorType::RowsAtCompileTime == 1, void>::type
computeConstitutiveRelation(double const /*dt*/,
                            double const lambda,
                            double const mu,
                            KelvinVectorType const& eps_prev,
                            KelvinVectorType const& eps,
                            KelvinVectorType const& sigma_prev,
                            KelvinVectorType& sigma,
                            KelvinMatrixType& C)
{
	C.setZero();

	C(0, 0) = mu * (3 * lambda + 2 * mu) / (lambda + mu);

	sigma = sigma_prev + C * (eps - eps_prev);
}
}
}

#endif  // PHYSICALQUANITIES_SOLIDS_LINEARELASTICISOTROPIC_H_
