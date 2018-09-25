/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <iostream>
#include <memory>
#include <vector>

#include <dune/geometry/quadraturerules.hh>

#include "MaterialLib/SolidModels/LinearElasticIsotropic.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "LocalAssemblerInterface.h"
#include "SmallDeformationProcessData.h"

namespace ProcessLib
{
namespace SmallDeformationDUNE
{
template <typename BMatricesType, typename ShapeMatricesType,
          int DisplacementDim>
struct IntegrationPointData final
{
    explicit IntegrationPointData(
        MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material)
        : solid_material(solid_material),
          material_state_variables(
              solid_material.createMaterialStateVariables())
    {
    }

    typename BMatricesType::KelvinVectorType sigma, sigma_prev;
    typename BMatricesType::KelvinVectorType eps, eps_prev;

    MaterialLib::Solids::MechanicsBase<DisplacementDim>& solid_material;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;

    double integration_weight;
    typename ShapeMatricesType::NodalRowVectorType N;
    typename ShapeMatricesType::GlobalDimNodalMatrixType dNdx;

    void pushBackState()
    {
        eps_prev = eps;
        sigma_prev = sigma;
        material_state_variables->pushBackState();
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

using namespace MathLib::KelvinVector;

/// Used by for extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>> N;
};

// ////////////////////////////////////////////////////////////////////////////
// Note:
// The essential difference between this local assembler and the non-DUNE one is
// in the contstructor where the shape matrices and the quadrature scheme are
// evaluated. The assembly routine and other helper routines are basically the
// same for the DUNE and non-DUNE cases.
// ////////////////////////////////////////////////////////////////////////////
template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim, typename Basis>
class SmallDeformationLocalAssembler
    : public SmallDeformationLocalAssemblerInterface<DisplacementDim>
{
public:
    using ShapeMatricesType =
        ShapeMatrixPolicyType<ShapeFunction, DisplacementDim>;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;
    using BMatricesType = BMatrixPolicyType<ShapeFunction, DisplacementDim>;

    using BMatrixType = typename BMatricesType::BMatrixType;
    using StiffnessMatrixType = typename BMatricesType::StiffnessMatrixType;
    using NodalForceVectorType = typename BMatricesType::NodalForceVectorType;
    using NodalDisplacementVectorType =
        typename BMatricesType::NodalForceVectorType;

    SmallDeformationLocalAssembler(SmallDeformationLocalAssembler const&) =
        delete;
    SmallDeformationLocalAssembler(SmallDeformationLocalAssembler&&) = delete;

    SmallDeformationLocalAssembler(
        typename BaseLib::DUNEGridType<DisplacementDim>::Traits::template Codim<
            0>::Entity const& e,
        Basis const& basis,
        bool const is_axially_symmetric,
        unsigned const /*integration_order*/,
        SmallDeformationProcessData<DisplacementDim>& process_data)
        : _process_data(process_data),
          _is_axially_symmetric(is_axially_symmetric)
    {
        auto localView = basis.localView();
        localView.bind(e);
        const auto& localFiniteElement =
            localView.tree().child(0).finiteElement();

        auto order =
            2 * (DisplacementDim * localFiniteElement.localBasis().order() - 1);
        const auto& quad = Dune::QuadratureRules<double, DisplacementDim>::rule(
            e.type(), order);
        auto const n_integration_points = quad.size();
        // DBUG("n int pts: %d.", n_integration_points);

        _ip_data.reserve(n_integration_points);
        _secondary_data.N.resize(n_integration_points);

        auto geometry = e.geometry();  // Trafo from ref to actual element

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            // DBUG("ip %d.", ip);

            _ip_data.emplace_back(*_process_data.material);
            auto& ip_data = _ip_data[ip];

            // Position of the current quadrature point in the reference element
            const auto quadPos = quad[ip].position();

            // The transposed inverse Jacobian of the map from the reference
            // element to the element
            const auto jacobian = geometry.jacobianInverseTransposed(quadPos);

            // The multiplicative factor in the integral transformation formula
            const auto integrationElement =
                geometry.integrationElement(quadPos);

            // The gradients of the shape functions on the reference element
            std::vector<Dune::FieldMatrix<double, 1, DisplacementDim>>
                referenceGradients;
            localFiniteElement.localBasis().evaluateJacobian(
                quadPos, referenceGradients);

            // Compute the shape function gradients on the real element
            std::vector<Dune::FieldVector<double, DisplacementDim>> gradients(
                referenceGradients.size());
            for (std::size_t i = 0; i < gradients.size(); i++)
            {
                jacobian.mv(referenceGradients[i][0], gradients[i]);
                // std::cout << "grad: [" << i << "] " << gradients[i] << ".\n";

                _ip_data[ip].dNdx.resize(DisplacementDim, gradients.size());
                for (std::size_t d = 0; d < DisplacementDim; ++d)
                {
                    _ip_data[ip].dNdx(d, i) = gradients[i][d];
                }
            }

            // Compute the shape function
            std::vector<Dune::FieldVector<double, 1>> N;
            localFiniteElement.localBasis().evaluateFunction(quadPos, N);

            _ip_data[ip].N.resize(N.size());
            for (std::size_t i = 0; i < N.size(); i++)
            {
                // std::cout << "N: [" << i << "] " << N[i] << ".\n";
                _ip_data[ip].N[i] = N[i];
            }

            _ip_data[ip].integration_weight =
                quad[ip].weight() * integrationElement;

            // Initialize current time step values
            ip_data.sigma.setZero(
                KelvinVectorDimensions<DisplacementDim>::value);
            ip_data.eps.setZero(KelvinVectorDimensions<DisplacementDim>::value);

            // Previous time step values are not initialized and are set later.
            ip_data.sigma_prev.resize(
                KelvinVectorDimensions<DisplacementDim>::value);
            ip_data.eps_prev.resize(
                KelvinVectorDimensions<DisplacementDim>::value);

            // TODO [DUNE] remove. set only temporarily until ip data
            // interpolation works.
            ip_data.eps_prev.setZero();
            ip_data.sigma_prev.setZero();

            _secondary_data.N[ip] = _ip_data[ip].N;
        }
    }

    void assemble(double const /*t*/, std::vector<double> const& /*local_x*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_b_data*/) override
    {
        OGS_FATAL(
            "SmallDeformationLocalAssembler: assembly without jacobian is not "
            "implemented.");
    }

    void assembleWithJacobian(double const t,
                              std::vector<double> const& local_x,
                              std::vector<double> const& /*local_xdot*/,
                              const double /*dxdot_dx*/, const double /*dx_dx*/,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& local_b_data,
                              std::vector<double>& local_Jac_data) override
    {
        // DBUG("assembling SD process.");

        auto const local_matrix_size = local_x.size();

        auto local_Jac = MathLib::createZeroedMatrix<StiffnessMatrixType>(
            local_Jac_data, local_matrix_size, local_matrix_size);

        auto local_b = MathLib::createZeroedVector<NodalDisplacementVectorType>(
            local_b_data, local_matrix_size);

        unsigned const n_integration_points = _ip_data.size();

        SpatialPosition x_position;
        // TODO [DUNE] re-enable
        // x_position.setElementID(_element.getID());
        x_position.setElementID(0);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& w = _ip_data[ip].integration_weight;
            auto const& N = _ip_data[ip].N;
            auto const& dNdx = _ip_data[ip].dNdx;

            typename ShapeMatricesType::template MatrixType<DisplacementDim,
                                                            displacement_size>
                N_u_op = ShapeMatricesType::template MatrixType<
                    DisplacementDim,
                    displacement_size>::Zero(DisplacementDim,
                                             displacement_size);
            for (int i = 0; i < DisplacementDim; ++i)
                N_u_op
                    .template block<1, displacement_size / DisplacementDim>(
                        i, i * displacement_size / DisplacementDim)
                    .noalias() = N;

            auto const x_coord = 0.0;  // TODO [DUNE] use true value
#if 0
                interpolateXCoordinate<ShapeFunction, ShapeMatricesType>(
                    _element, N);
#endif
            auto const B = LinearBMatrix::computeBMatrix<
                DisplacementDim, ShapeFunction::NPOINTS,
                typename BMatricesType::BMatrixType>(dNdx, N, x_coord,
                                                     _is_axially_symmetric);

            auto const& eps_prev = _ip_data[ip].eps_prev;
            auto const& sigma_prev = _ip_data[ip].sigma_prev;

            auto& eps = _ip_data[ip].eps;
            auto& sigma = _ip_data[ip].sigma;
            auto& state = _ip_data[ip].material_state_variables;

            eps.noalias() =
                B *
                Eigen::Map<typename BMatricesType::NodalForceVectorType const>(
                    local_x.data(), ShapeFunction::NPOINTS * DisplacementDim);

            auto&& solution = _ip_data[ip].solid_material.integrateStress(
                t, x_position, _process_data.dt, eps_prev, eps, sigma_prev,
                *state, 0.0 /* T */);

            if (!solution)
                OGS_FATAL("Computation of local constitutive relation failed.");

            KelvinMatrixType<DisplacementDim> C;
            std::tie(sigma, state, C) = std::move(*solution);

            auto const rho = _process_data.solid_density(t, x_position)[0];
            auto const& b = _process_data.specific_body_force;
            local_b.noalias() -=
                (B.transpose() * sigma - N_u_op.transpose() * rho * b) * w;
            local_Jac.noalias() += B.transpose() * C * B * w;
        }
    }

    void preTimestepConcrete(std::vector<double> const& /*local_x*/,
                             double const /*t*/,
                             double const /*delta_t*/) override
    {
        preTimestep();
    }

    void preTimestep() override
    {
        for (auto& d : _ip_data)
        {
            d.pushBackState();
        }
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _secondary_data.N[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getIntPtSigma(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::AbstractDOFTable const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        using KelvinVectorType = typename BMatricesType::KelvinVectorType;
        auto const kelvin_vector_size =
            KelvinVectorDimensions<DisplacementDim>::value;
        auto const num_intpts = _ip_data.size();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, kelvin_vector_size, num_intpts);

        // TODO [DUNE] make a general implementation for converting
        // KelvinVectors back to symmetric rank-2 tensors.
        for (unsigned ip = 0; ip < num_intpts; ++ip)
        {
            auto const& sigma = _ip_data[ip].sigma;

            for (typename KelvinVectorType::Index component = 0;
                 component < kelvin_vector_size && component < 3;
                 ++component)
            {  // xx, yy, zz components
                cache_mat(component, ip) = sigma[component];
            }
            for (typename KelvinVectorType::Index component = 3;
                 component < kelvin_vector_size;
                 ++component)
            {  // mixed xy, yz, xz components
                cache_mat(component, ip) = sigma[component] / std::sqrt(2);
            }
        }

        return cache;
    }

    virtual std::vector<double> const& getIntPtEpsilon(
        const double /*t*/,
        GlobalVector const& /*current_solution*/,
        NumLib::AbstractDOFTable const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        using KelvinVectorType = typename BMatricesType::KelvinVectorType;
        auto const kelvin_vector_size =
            KelvinVectorDimensions<DisplacementDim>::value;
        auto const num_intpts = _ip_data.size();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, kelvin_vector_size, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, kelvin_vector_size, num_intpts);

        // TODO [DUNE] make a general implementation for converting
        // KelvinVectors back to symmetric rank-2 tensors.
        for (unsigned ip = 0; ip < num_intpts; ++ip)
        {
            auto const& eps = _ip_data[ip].eps;

            for (typename KelvinVectorType::Index component = 0;
                 component < kelvin_vector_size && component < 3;
                 ++component)
            {  // xx, yy, zz components
                cache_mat(component, ip) = eps[component];
            }
            for (typename KelvinVectorType::Index component = 3;
                 component < kelvin_vector_size;
                 ++component)
            {  // mixed xy, yz, xz components
                cache_mat(component, ip) = eps[component] / std::sqrt(2);
            }
        }

        return cache;
    }

    unsigned getNumberOfIntegrationPoints() const override
    {
        return _ip_data.size();
    }

    typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables const&
    getMaterialStateVariablesAt(unsigned integration_point) const override
    {
        return *_ip_data[integration_point].material_state_variables;
    }

    Eigen::Map<
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
    getIntPtSigma2(std::vector<double>& cache) const override
    {
        using KelvinVectorType = typename BMatricesType::KelvinVectorType;
        auto const kelvin_vector_size =
            KelvinVectorDimensions<DisplacementDim>::value;
        auto const num_intpts = _ip_data.size();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, kelvin_vector_size, num_intpts);

        for (unsigned ip = 0; ip < num_intpts; ++ip)
        {
            auto const& sigma = _ip_data[ip].sigma;

            for (typename KelvinVectorType::Index component = 0;
                 component < kelvin_vector_size;
                 ++component)
            {
                cache_mat(component, ip) = sigma[component];
            }
        }

        return cache_mat;
    }

    Eigen::Map<
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
    getIntPtSigmaPrev2(std::vector<double>& cache) const override
    {
        using KelvinVectorType = typename BMatricesType::KelvinVectorType;
        auto const kelvin_vector_size =
            KelvinVectorDimensions<DisplacementDim>::value;
        auto const num_intpts = _ip_data.size();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, kelvin_vector_size, num_intpts);

        for (unsigned ip = 0; ip < num_intpts; ++ip)
        {
            auto const& sigma = _ip_data[ip].sigma_prev;

            for (typename KelvinVectorType::Index component = 0;
                 component < kelvin_vector_size;
                 ++component)
            {
                cache_mat(component, ip) = sigma[component];
            }
        }

        return cache_mat;
    }

    Eigen::Map<
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
    getIntPtEpsilon2(std::vector<double>& cache) const override
    {
        using KelvinVectorType = typename BMatricesType::KelvinVectorType;
        auto const kelvin_vector_size =
            KelvinVectorDimensions<DisplacementDim>::value;
        auto const num_intpts = _ip_data.size();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, kelvin_vector_size, num_intpts);

        // TODO [DUNE] make a general implementation for converting
        // KelvinVectors back to symmetric rank-2 tensors.
        for (unsigned ip = 0; ip < num_intpts; ++ip)
        {
            auto const& eps = _ip_data[ip].eps;

            for (typename KelvinVectorType::Index component = 0;
                 component < kelvin_vector_size;
                 ++component)
            {
                cache_mat(component, ip) = eps[component];
            }
        }

        return cache_mat;
    }

    Eigen::Map<
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>
    getIntPtEpsilonPrev2(std::vector<double>& cache) const override
    {
        using KelvinVectorType = typename BMatricesType::KelvinVectorType;
        auto const kelvin_vector_size =
            KelvinVectorDimensions<DisplacementDim>::value;
        auto const num_intpts = _ip_data.size();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
            double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, kelvin_vector_size, num_intpts);

        // TODO [DUNE] make a general implementation for converting
        // KelvinVectors back to symmetric rank-2 tensors.
        for (unsigned ip = 0; ip < num_intpts; ++ip)
        {
            auto const& eps = _ip_data[ip].eps_prev;

            for (typename KelvinVectorType::Index component = 0;
                 component < kelvin_vector_size;
                 ++component)
            {
                cache_mat(component, ip) = eps[component];
            }
        }

        return cache_mat;
    }

    void setIntPtSigma(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                     Eigen::RowMajor> const& values) override
    {
        using KelvinVectorType = typename BMatricesType::KelvinVectorType;
        auto const kelvin_vector_size =
            KelvinVectorDimensions<DisplacementDim>::value;
        auto const num_intpts = _ip_data.size();

        for (unsigned ip = 0; ip < num_intpts; ++ip)
        {
            auto& sigma = _ip_data[ip].sigma;

            for (typename KelvinVectorType::Index component = 0;
                 component < kelvin_vector_size;
                 ++component)
            {
                sigma[component] = values(component, ip);
            }
        }
    }
    void setIntPtSigmaPrev(
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                      Eigen::RowMajor> const& values) override
    {
        using KelvinVectorType = typename BMatricesType::KelvinVectorType;
        auto const kelvin_vector_size =
            KelvinVectorDimensions<DisplacementDim>::value;
        auto const num_intpts = _ip_data.size();

        for (unsigned ip = 0; ip < num_intpts; ++ip)
        {
            auto& sigma = _ip_data[ip].sigma_prev;

            for (typename KelvinVectorType::Index component = 0;
                 component < kelvin_vector_size;
                 ++component)
            {
                sigma[component] = values(component, ip);
            }
        }
    }

    void setIntPtEpsilon(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                       Eigen::RowMajor> const& values) override
    {
        using KelvinVectorType = typename BMatricesType::KelvinVectorType;
        auto const kelvin_vector_size =
            KelvinVectorDimensions<DisplacementDim>::value;
        auto const num_intpts = _ip_data.size();

        for (unsigned ip = 0; ip < num_intpts; ++ip)
        {
            auto& eps = _ip_data[ip].eps;

            for (typename KelvinVectorType::Index component = 0;
                 component < kelvin_vector_size;
                 ++component)
            {
                eps[component] = values(component, ip);
            }
        }
    }

    void setIntPtEpsilonPrev(
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                      Eigen::RowMajor> const& values) override
    {
        using KelvinVectorType = typename BMatricesType::KelvinVectorType;
        auto const kelvin_vector_size =
            KelvinVectorDimensions<DisplacementDim>::value;
        auto const num_intpts = _ip_data.size();

        for (unsigned ip = 0; ip < num_intpts; ++ip)
        {
            auto& eps = _ip_data[ip].eps_prev;

            for (typename KelvinVectorType::Index component = 0;
                 component < kelvin_vector_size;
                 ++component)
            {
                eps[component] = values(component, ip);
            }
        }
    }

private:
    SmallDeformationProcessData<DisplacementDim>& _process_data;

    std::vector<
        IntegrationPointData<BMatricesType, ShapeMatricesType, DisplacementDim>,
        Eigen::aligned_allocator<IntegrationPointData<
            BMatricesType, ShapeMatricesType, DisplacementDim>>>
        _ip_data;

    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;
    bool const _is_axially_symmetric;

    static const int displacement_size =
        ShapeFunction::NPOINTS * DisplacementDim;
};

}  // namespace SmallDeformationDUNE
}  // namespace ProcessLib
