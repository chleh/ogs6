/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_SMALLDEFORMATION_FEM_H_
#define PROCESS_LIB_SMALLDEFORMATION_FEM_H_

#include <memory>
#include <vector>

#include "Parameter.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Fem/Deformation/BMatrixPolicy.h"
#include "NumLib/Fem/Deformation/LinearBMatrix.h"

#include "PhysicalQuantities/Solids/LinearElasticIsotropic.h"

namespace ProcessLib
{
namespace SmallDeformation
{
template <typename GlobalMatrix, typename GlobalVector>
class LocalAssemblerDataInterface
{
public:
	virtual ~LocalAssemblerDataInterface() = default;

	virtual void assemble(double const t,
	                      std::vector<double> const& local_x) = 0;

	virtual void addToGlobal(
	    AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const&,
	    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) const = 0;
};

template <typename ShapeFunction_, typename IntegrationMethod_,
          typename GlobalMatrix, typename GlobalVector, unsigned GlobalDim>
class LocalAssemblerData final
    : public LocalAssemblerDataInterface<GlobalMatrix, GlobalVector>
{
public:
	using ShapeFunction = ShapeFunction_;
	using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
	using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
	using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
	using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

	using BMatricesType = BMatrixPolicyType<ShapeFunction, GlobalDim>;
	using BMatrixType = typename BMatricesType::BMatrixType;

	using StressVectorType = typename BMatricesType::StressVectorType;
	using StrainVectorType = typename BMatricesType::StressVectorType;

	using ModulusMatrixType = typename BMatricesType::ModulusMatrixType;
	using StiffnessMatrixType = typename BMatricesType::StiffnessMatrixType;
	using NodalForceVectorType = typename BMatricesType::NodalForceVectorType;
	using NodalDisplacementVectorType =
	    typename BMatricesType::NodalForceVectorType;

	/// The lambda and mu factors is directly integrated into the local
	/// element matrix.
	LocalAssemblerData(
	    MeshLib::Element const& e,
	    std::size_t const local_matrix_size,
	    unsigned const integration_order,
	    Parameter<double, MeshLib::Element const&> const& youngs_modulus,
	    Parameter<double, MeshLib::Element const&> const& poissons_ratio)
	{
		using FemType =
		    NumLib::TemplateIsoparametric<ShapeFunction, ShapeMatricesType>;

		FemType fe(
		    *static_cast<const typename ShapeFunction::MeshElement*>(&e));

		_integration_order = integration_order;
		IntegrationMethod_ integration_method(_integration_order);
		std::size_t const n_integration_points =
		    integration_method.getNPoints();

		//_shape_matrices.resize(n_integration_points); See note below for
        //resize.
		_b_matrices.resize(n_integration_points);
		_sigma.resize(n_integration_points);
		_sigma_prev.resize(n_integration_points);
		_eps.resize(n_integration_points);
		_eps_prev.resize(n_integration_points);
		_C.resize(n_integration_points);
		for (std::size_t ip(0); ip < n_integration_points; ip++)
		{
            // Previously it was necessary to resize the matrices, but now it is
            // handled by shape matrices' constructor. Isn't it?
			//_shape_matrices[ip].resize(ShapeFunction::DIM,
			//                           ShapeFunction::NPOINTS);
			_shape_matrices.emplace_back(ShapeFunction::DIM, GlobalDim,
			                             ShapeFunction::NPOINTS);

			fe.computeShapeFunctions(
			    integration_method.getWeightedPoint(ip).getCoords(),
			    _shape_matrices[ip]);

			_b_matrices[ip].resize(
			    BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::rows,
			    BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::columns);
			LinearBMatrix::computeBMatrix<
			    // Passing GlobalDim here is wrong. Should be DisplacementDim
			    // but this needs additional template argument and corresponding
			    // change in AssemblerLib::LocalDataInitializer
			    GlobalDim, ShapeFunction::NPOINTS,
			    typename ShapeMatricesType::GlobalDimNodalMatrixType,
			    BMatrixType>(_shape_matrices[ip].dNdx, _b_matrices[ip]);

			_sigma[ip].resize(
			    BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::rows);
			_sigma_prev[ip].resize(
			    BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::rows);
			_eps[ip].resize(
			    BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::rows);
			_eps_prev[ip].resize(
			    BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::rows);
			_C[ip].resize(
			    BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::rows,
			    BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::rows);
		}

		// Elementwise parameter checks.
		{

			if (youngs_modulus(e) <= 0.)
			{
				ERR(
				    "The Young's modulus parameter value for element %d is "
				    "negative.",
				    e.getID());
				std::abort();
			}

			if (poissons_ratio(e) <= -1 || poissons_ratio(e ) > 0.5)
			{
				ERR(
				    "The Poisson's ratio parameter value for element %d does "
				    "not lie in (-1, 0.5]",
				    e.getID());
				std::abort();
			}
		}
		_lambda = [&youngs_modulus, &poissons_ratio, &e]()
		{
			return youngs_modulus(e) * poissons_ratio(e) /
		       (1 + poissons_ratio(e)) / (1 - 2 * poissons_ratio(e));
		};

		_mu = [&youngs_modulus, &poissons_ratio, &e]()
		{
			return youngs_modulus(e)/(2*(1+poissons_ratio(e)));
		};

		_localA.reset(
		    new StiffnessMatrixType(local_matrix_size, local_matrix_size));
		_localRhs.reset(new NodalForceVectorType(local_matrix_size));

		for (int dim = 0; dim < 3; ++dim)
			for (unsigned i = 0; i < e.getNNodes(); ++i)
				_coords.push_back((*e.getNode(i))[dim]);
	}

	void assemble(double const /*t*/,
	              std::vector<double> const& local_x) override
	{
		_localA->setZero();
		_localRhs->setZero();

		IntegrationMethod_ integration_method(_integration_order);
		unsigned const n_integration_points = integration_method.getNPoints();

		NodalDisplacementVectorType local_displacement(
		    BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::columns);
		NodalDisplacementVectorType local_displacement_prev(
		    BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::columns);
		for (int i = 0;
		     i < BMatrixDimensions<ShapeFunction::NPOINTS, GlobalDim>::columns;
		     ++i)
		{
			local_displacement(i) = local_x[i];
			local_displacement_prev(i) = local_x_prev[i];
		}

		NodalVectorType coords_x(ShapeFunction::NPOINTS);
		NodalVectorType coords_y(ShapeFunction::NPOINTS);
		NodalVectorType coords_z(ShapeFunction::NPOINTS);

		for (int i = 0; i < ShapeFunction::NPOINTS; ++i)
		{
			coords_x[i] = _coords[0 * ShapeFunction::NPOINTS + i];
			coords_y[i] = _coords[1 * ShapeFunction::NPOINTS + i];
			coords_z[i] = _coords[2 * ShapeFunction::NPOINTS + i];
		}

		for (std::size_t ip(0); ip < n_integration_points; ip++)
		{
			auto const& sm = _shape_matrices[ip];
			auto const& B = _b_matrices[ip];
			auto const& wp = integration_method.getWeightedPoint(ip);

			LinearBMatrix::computeStrain(local_displacement, B, _eps[ip]);
			LinearBMatrix::computeStrain(local_displacement_prev, B,
			                             _eps_prev[ip]);

			Solids::LinearElasticIsotropic::computeConstitutiveRelation(
			    _lambda(), _mu(), _eps_prev[ip], _eps[ip], _sigma_prev[ip],
			    _sigma[ip], _C[ip]);
			std::cout << "ip " << ip << ":\t";
			std::cout << "C: " << _C[ip] << "\t";
			std::cout << "sigma: " << _sigma[ip] << "\t";
			std::cout << "sigma_prev: " << _sigma_prev[ip] << "\t";
			std::cout << "eps: " << _eps[ip] << "\t";
			std::cout << "eps_prev: " << _eps_prev[ip] << "\t";
			std::cout << std::endl;
			_localA->noalias() +=
			    B.transpose() * _C[ip] * B * sm.detJ * wp.getWeight();

			_localRhs->noalias() -=
				B.transpose() * _sigma[ip] * sm.detJ * wp.getWeight();
		}
	}

	void addToGlobal(
	    AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const& indices,
	    GlobalMatrix& /*M*/, GlobalMatrix& K, GlobalVector& b) const override
	{
		K.add(indices, *_localA);
		b.add(indices.rows, *_localRhs);
	}

private:
	std::vector<ShapeMatrices> _shape_matrices;
	std::vector<BMatrixType> _b_matrices;
	std::vector<StressVectorType> _sigma, _sigma_prev;
	std::vector<StrainVectorType> _eps, _eps_prev;
	std::vector<ModulusMatrixType> _C;

	std::vector<double> _coords;

	std::function<double(void)> _lambda;
	std::function<double(void)> _mu;

	std::unique_ptr<StiffnessMatrixType> _localA;
	std::unique_ptr<NodalForceVectorType> _localRhs;

	unsigned _integration_order = 2;
};

}  // namespace SmallDeformation
}  // namespace ProcessLib

#endif  // PROCESS_LIB_SMALLDEFORMATION_FEM_H_
