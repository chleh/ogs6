/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cassert>
#include <vector>

#include "BaseLib/DUNEConfig.h"

#if HAVE_UG
#pragma message "HAVE_UG set #####"
#else
#pragma message "HAVE_UG not set #####"
#endif

#include <dune/functions/functionspacebases/defaultglobalbasis.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>

#include "BaseLib/Functional.h"
#include "NumLib/DOF/DOFTableDUNE.h"
#include "NumLib/Extrapolation/ExtrapolateSingleElement.h"
#include "NumLib/Extrapolation/LocalLinearLeastSquaresExtrapolator.h"
#include "NumLib/Extrapolation/ZienkiewiczZhuGradientRecoveryEstimator.h"
#include "ProcessLib/Process.h"

#include "CreateLocalAssemblers.h"
#include "SmallDeformationFEM.h"
#include "SmallDeformationProcess.h"

static std::size_t smallDeformationDUNERefinementOutputCounter = 0;

namespace
{
template <int DisplacementDim, typename Element, typename LocalView>
decltype(auto) get_quadrature(Element const& element,
                              LocalView const& localView)
{
    const auto& localFiniteElement = localView.tree().child(0).finiteElement();
    auto const integration_order =
        2 * (DisplacementDim * localFiniteElement.localBasis().order() - 1);

    const auto& quadrature =
        Dune::QuadratureRules<double, DisplacementDim>::rule(element.type(),
                                                             integration_order);
    return quadrature;
}

template <typename Position, typename Element, typename Predicate>
void map_local_coordinates_to_father_recursive(Position& position,
                                               Element& element,
                                               Predicate&& loop_while)
{
    do
    {
        position = element.geometryInFather().global(position);
        element = element.father();
    } while (loop_while(element));
}

template <typename Matrix, typename Quadrature, typename LocalFiniteElement,
          typename ShapeFunctions>
decltype(auto) interpolate_nodal_data_to_integration_points(
    Matrix const& nodal_data,
    Quadrature const& quadrature,
    LocalFiniteElement const& localFiniteElement,
    ShapeFunctions& shape_functions)
{
    auto const num_components = nodal_data.rows();
    auto const num_int_pts = quadrature.size();

    // (num_components x num_element_int_pts) as wanted by local
    // assembler setters
    auto int_pt_data = Eigen::MatrixXd(num_components, num_int_pts);

    for (unsigned ip = 0; ip < num_int_pts; ++ip)
    {
        // interpolate nodal data to integration points
        for (unsigned c = 0; c < num_components; ++c)
        {
            localFiniteElement.localBasis().evaluateFunction(
                quadrature[ip].position(), shape_functions);

            double interpolated_value = 0;
            for (unsigned i = 0; i < shape_functions.size(); ++i)
            {
                interpolated_value += shape_functions[i] * nodal_data(c, i);
            }

            int_pt_data(c, ip) = interpolated_value;
        }
    }  // for each int pt

    return int_pt_data;
}

template <typename Matrix, typename Element, typename Quadrature,
          typename LocalFiniteElement, typename ShapeFunctions>
decltype(auto) interpolate_father_nodal_data_to_new_childs_integration_points(
    Matrix const& father_nodal_data, Element const& element,
    Quadrature const& quadrature_element,
    LocalFiniteElement const& localFiniteElement_father,
    ShapeFunctions& shape_functions)
{
    auto const num_components = father_nodal_data.rows();
    auto const num_element_int_pts = quadrature_element.size();

    // (num_components x num_element_int_pts) as wanted by local
    // assembler setters
    auto element_int_pt_data =
        Eigen::MatrixXd(num_components, num_element_int_pts);

    for (unsigned ip = 0; ip < num_element_int_pts; ++ip)
    {
        auto e = element;
        auto int_pt_position = quadrature_element[ip].position();

        map_local_coordinates_to_father_recursive(
            int_pt_position, e, [](auto const& e) { return e.isNew(); });

        // interpolate father's nodal data to new element's
        // integration points
        for (unsigned c = 0; c < num_components; ++c)
        {
            localFiniteElement_father.localBasis().evaluateFunction(
                int_pt_position, shape_functions);
            assert(shape_functions.size() ==
                   static_cast<std::size_t>(father_nodal_data.cols()));

            double interpolated_value = 0;
            for (unsigned i = 0; i < shape_functions.size(); ++i)
            {
                interpolated_value +=
                    shape_functions[i] * father_nodal_data(c, i);
            }

            element_int_pt_data(c, ip) = interpolated_value;
        }
    }  // for each int pt

    return element_int_pt_data;
}

template <typename ChildData, typename ShapeFunctionsCache,
          typename LocalFiniteElement>
decltype(auto) get_father_shape_functions_at_integration_points_of_all_children(
    ChildData const& child_data, ShapeFunctionsCache& shape_functions,
    std::vector<double>& all_shape_functions,
    LocalFiniteElement const& localFiniteElement_father)
{
    all_shape_functions.clear();
    std::size_t int_pt_cumulative = 0;
    std::size_t num_father_nodes = 0;
    for (auto const& child_idx_and_int_pts : child_data)
    {
        auto const& child_int_pts =
            child_idx_and_int_pts.int_pt_coords_in_father;
        auto const num_child_int_pts = child_int_pts.size();

        // get father's shape functions at the integration points of the
        // children
        for (std::size_t int_pt = 0; int_pt < num_child_int_pts;
             ++int_pt, ++int_pt_cumulative)
        {
            auto const& pt = child_int_pts[int_pt];

            // evaluate father's shape function at child's
            // integration point
            localFiniteElement_father.localBasis().evaluateFunction(
                pt, shape_functions);

            // now actually collect all shape functions into one
            // matrix
            num_father_nodes = shape_functions.size();
            if (int_pt == 0)
                all_shape_functions.resize(all_shape_functions.size() +
                                           num_father_nodes *
                                               num_child_int_pts);
            auto const num_int_pts_so_far =
                all_shape_functions.size() / num_father_nodes;
            auto all_shape_functions_mat = MathLib::toMatrix(
                all_shape_functions, num_int_pts_so_far, num_father_nodes);

            for (std::size_t n = 0; n < num_father_nodes; ++n)
                all_shape_functions_mat(int_pt_cumulative, n) =
                    shape_functions[n];
        }
    }

    auto const all_shape_functions_mat = MathLib::toMatrix(
        all_shape_functions, int_pt_cumulative, num_father_nodes);
    return all_shape_functions_mat;
}

template <typename ChildData, typename LocalAssemblers, typename Getter>
decltype(auto) get_integration_point_data_of_all_children(
    ChildData const& child_data,
    LocalAssemblers const& local_assemblers,
    Getter&& getter,
    std::vector<double>& all_int_pt_data,
    std::vector<double>& single_element_int_pt_data)
{
    all_int_pt_data.clear();

    std::size_t num_components = 0;

    // collect integration point data from all children
    for (auto const& child_idx_and_int_pts : child_data)
    {
        auto const& child_int_pts =
            child_idx_and_int_pts.int_pt_coords_in_father;
        auto const num_child_int_pts = child_int_pts.size();

        // get child's integration point data
        auto child_idx = child_idx_and_int_pts.index;
        assert(local_assemblers[child_idx]);
        auto const& child_asm = *local_assemblers[child_idx];
        auto const child_int_pt_data =
            (child_asm.*getter)(single_element_int_pt_data);
        assert(num_components == 0 ||
               num_components ==
                   static_cast<std::size_t>(child_int_pt_data.rows()));
        num_components = static_cast<std::size_t>(child_int_pt_data.rows());
        assert(static_cast<std::size_t>(child_int_pt_data.cols()) ==
               num_child_int_pts);

        // collect all integration point data into one matrix
        all_int_pt_data.resize(all_int_pt_data.size() +
                               num_components * num_child_int_pts);
        auto const num_int_pts_so_far = all_int_pt_data.size() / num_components;
        auto all_int_pt_data_mat = MathLib::toMatrix(
            all_int_pt_data, num_int_pts_so_far, num_components);
        all_int_pt_data_mat.bottomRows(num_child_int_pts).noalias() =
            child_int_pt_data.transpose();
    }

    auto const all_int_pt_data_mat =
        MathLib::toMatrix(all_int_pt_data,
                          all_int_pt_data.size() / num_components,
                          num_components);
    return all_int_pt_data_mat;
}

}  // anonymous namespace

namespace ProcessLib
{
namespace SmallDeformationDUNE
{
// DUNE docs says there exists dynamic polymorphism for that, but how?
template <int DisplacementDim>
decltype(auto) makeDisplacementBasis(
    typename BaseLib::DUNEGridType<DisplacementDim>::LeafGridView const&
        gridView)
{
    namespace DFB = Dune::Functions::BasisBuilder;
    return DFB::makeBasis(gridView,
                          DFB::power<DisplacementDim>(DFB::lagrange<1>(),
                                                      DFB::flatInterleaved()));
}

template <int DisplacementDim>
decltype(auto) makeScalarBasis(
    typename BaseLib::DUNEGridType<DisplacementDim>::LeafGridView const&
        gridView)
{
    namespace DFB = Dune::Functions::BasisBuilder;
    return DFB::makeBasis(gridView, DFB::lagrange<1>());
}

template <int DisplacementDim>
SmallDeformationProcess<DisplacementDim>::SmallDeformationProcess(
    MeshLib::FEMMesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    SmallDeformationProcessData<DisplacementDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller)
    : Process(mesh, nullptr /* global assembler */, parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller)),
      _process_data(std::move(process_data)),
      _global_assembler_dune(std::move(jacobian_assembler))
{
// TODO [DUNE] re-enable
#if 0
    _nodal_forces = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "NodalForces", MeshLib::MeshItemType::Node, DisplacementDim);
#endif
}

template <int DisplacementDim>
bool SmallDeformationProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::setInitialConditions(
    const int process_id, const double t, GlobalVector& x)
{
    DBUG("Set initial conditions.");

    SpatialPosition pos;

    auto const& pvs = getProcessVariables(process_id);
    assert(pvs.size() == 1);
    ProcessVariable& pv = pvs.front();
    auto const& ic = pv.getInitialCondition();

    auto const num_comp = pv.getNumberOfComponents();

    // TODO [DUNE] assumes by location ordering
    auto gridView = _process_data.grid.getMesh().leafGridView();
    auto const num_nodes = gridView.size(DisplacementDim);

    for (int component_id = 0; component_id < num_comp; ++component_id)
    {
        for (GlobalIndexType node = 0; node < num_nodes; ++node)
        {
            auto const global_index = component_id * num_nodes + node;

            pos.setNodeID(node);
            auto const& ic_value = ic(t, pos);
            x.set(global_index, ic_value[component_id]);
        }
    }
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::AbstractDOFTable const& /*dof_table*/,
    MeshLib::FEMMesh const& /*mesh*/,
    unsigned const /*integration_order*/)
{
    DBUG("init SmallDefDUNE.");

    DBUG("DUNE Mesh %p", &_process_data.grid.getMesh());

    createLocalAssemblers(_process_data.grid, true,
                          MeshLib::DUNEIdToIdxMappings{});

    initializeExtrapolator();

// TODO [DUNE] re-enable
#if 0
    // TODO [DUNE] move the two data members somewhere else.
    // for extrapolation of secondary variables
    std::vector<MeshLib::MeshSubsets> all_mesh_subsets_single_component;
    all_mesh_subsets_single_component.emplace_back(
        _mesh_subset_all_nodes.get());
    _local_to_global_index_map_single_component =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);
    _nodal_forces->resize(DisplacementDim * mesh.getNumberOfNodes());
#endif

    Base::_secondary_variables.addSecondaryVariable(
        "sigma",
        makeExtrapolator(MathLib::KelvinVector::KelvinVectorType<
                             DisplacementDim>::RowsAtCompileTime,
                         getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigma));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon",
        makeExtrapolator(MathLib::KelvinVector::KelvinVectorType<
                             DisplacementDim>::RowsAtCompileTime,
                         getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsilon));

    // enable output of internal variables defined by material models
    auto const internal_variables =
        _process_data.material->getInternalVariables();
    for (auto const& internal_variable : internal_variables)
    {
        auto const& name = internal_variable.name;
        auto const& fct = internal_variable.getter;
        auto const num_components = internal_variable.num_components;
        DBUG("Registering internal variable %s.", name.c_str());

        auto getIntPtValues = BaseLib::easyBind(
            [fct, num_components](
                LocalAssemblerInterface const& loc_asm,
                const double /*t*/,
                GlobalVector const& /*current_solution*/,
                NumLib::AbstractDOFTable const& /*dof_table*/,
                std::vector<double>& cache) -> std::vector<double> const& {
                const unsigned num_int_pts =
                    loc_asm.getNumberOfIntegrationPoints();

                cache.clear();
                auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
                    double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
                    cache, num_components, num_int_pts);

                // TODO [DUNE] avoid the heap allocation (one per finite
                // element)
                std::vector<double> cache_column(num_int_pts);

                for (unsigned i = 0; i < num_int_pts; ++i)
                {
                    auto const& state = loc_asm.getMaterialStateVariablesAt(i);

                    auto const& int_pt_values = fct(state, cache_column);
                    assert(int_pt_values.size() == num_components);
                    auto const int_pt_values_vec =
                        MathLib::toVector(int_pt_values);

                    cache_mat.col(i).noalias() = int_pt_values_vec;
                }

                return cache;
            });

        Base::_secondary_variables.addSecondaryVariable(
            name,
            makeExtrapolator(num_components, getExtrapolator(),
                             _local_assemblers, std::move(getIntPtValues)));
    }

    // TODO [DUNE] might do duplicate work in case of initial global refinement
    computeSparsityPattern(_process_data.grid, false,
                           MeshLib::DUNEIdToIdxMappings{});

    _process_data.grid.onPostRefine(
        BaseLib::easyBind(
            &SmallDeformationProcess<DisplacementDim>::computeSparsityPattern,
            *this),
        nullptr);
    _process_data.grid.onPostRefine(
        BaseLib::easyBind(
            &SmallDeformationProcess<DisplacementDim>::createLocalAssemblers,
            *this),
        nullptr);
    _process_data.grid.onPostRefine(
        BaseLib::easyBind(
            &SmallDeformationProcess<DisplacementDim>::updateDOFtables, *this),
        nullptr);
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::computeSparsityPattern(
    MeshLib::DUNEMesh<DisplacementDim> const& mesh, bool,
    MeshLib::DUNEIdToIdxMappings const&)
{
    DBUG("Computing Sparsity Pattern.");
    OGS_ALWAYS_ASSERT(&mesh == &_process_data.grid);
    auto gridView = _process_data.grid.getMesh().leafGridView();

    auto basis = makeDisplacementBasis<DisplacementDim>(gridView);
    DBUG("Basis type %s\nof size %d.", Dune::className(basis).c_str(),
         sizeof(basis));

    // compute sparsity pattern
    std::vector<std::vector<GlobalIndexType>> global_idcs(basis.size());
    auto localView = basis.localView();
    auto localIndexSet = basis.localIndexSet();

    for (auto const& e : Dune::elements(gridView))
    {
        localView.bind(e);
        localIndexSet.bind(localView);
        for (size_t i = 0; i < localIndexSet.size(); i++)
        {
            auto const row = localIndexSet.index(i);
            assert(row.size() == 1);

            for (size_t j = 0; j < localIndexSet.size(); j++)
            {
                auto const col = localIndexSet.index(j);
                assert(col.size() == 1);

                global_idcs[row[0]].push_back(col[0]);
            }
        }
    }

    _sparsity_pattern.clear();
    _sparsity_pattern.reserve(basis.size());
    for (auto const& row : global_idcs)
    {
        _sparsity_pattern.push_back(row.size());
    }
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::initializeExtrapolator()
{
    auto gridView = _process_data.grid.getMesh().leafGridView();
    auto basis = makeDisplacementBasis<DisplacementDim>(gridView);

    _local_to_global_index_map = NumLib::makeDOFTableDUNE(std::move(basis));

    auto scalar_basis = makeScalarBasis<DisplacementDim>(gridView);
    auto* dof_table_single_component =
        NumLib::makeDOFTableDUNE(std::move(scalar_basis)).release();

    auto extrapolator =
        std::make_unique<NumLib::LocalLinearLeastSquaresExtrapolator>(
            *dof_table_single_component);

    // TODO [DUNE] Later on the DOF table can change during the simulation!
    _extrapolator_data =
        ExtrapolatorData(std::move(extrapolator), dof_table_single_component,
                         /*manage_storage*/ true);

    auto error_estimator =
        std::make_unique<NumLib::ZienkiewiczZhuGradientRecoveryEstimator>();

    // Init error estimator
    _error_estimator_data = ErrorEstimatorData{std::move(error_estimator),
                                               dof_table_single_component};
}

template <int DisplacementDim>
const GlobalVector* SmallDeformationProcess<DisplacementDim>::estimateError(
    GlobalVector const& /*x*/, double& global_relative_error)
{
    _error_estimator_data.estimator_algorithm->estimate(
        _process_data.grid.getMesh(), _local_assemblers,
        &LocalAssemblerInterface::getIntPtSigma2);
    return &_error_estimator_data.estimator_algorithm->getErrorEstimate(
        global_relative_error);
}

template <int DisplacementDim>
bool SmallDeformationProcess<DisplacementDim>::refine(
    std::vector<char> const& elements_for_refinement)
{
    auto const& grid = _process_data.grid.getMesh();
    auto const gridView = grid.leafGridView();

    auto const& indexSet = gridView.indexSet();
    auto const& idSet = grid.localIdSet();

    // mark and adapt grid
    for (auto& e : Dune::elements(gridView))
    {
        auto const indicator = elements_for_refinement[indexSet.index(e)];
        assert(indicator == -1 || indicator == 0 || indicator == 1);
        _process_data.grid.mark(indicator, e);
    }

    if (!_process_data.grid.preAdapt())
    {
        INFO("mesh does not need to be adapted");
        return false;
    }

    // TODO [DUNE] debug output
    Dune::VTKWriter<decltype(gridView)> vtkWriter(gridView);
    vtkWriter.addCellData(elements_for_refinement, "refinement_indicator", 1);
    vtkWriter.write(
        "refinement_indicator_" +
            std::to_string(++smallDeformationDUNERefinementOutputCounter),
        Dune::VTK::OutputType::base64);

    auto const getters_setters = {
        std::make_pair(&LocalAssemblerInterface::getIntPtEpsilon2,
                       &LocalAssemblerInterface::setIntPtEpsilon),
        std::make_pair(&LocalAssemblerInterface::getIntPtEpsilonPrev2,
                       &LocalAssemblerInterface::setIntPtEpsilonPrev),
        std::make_pair(&LocalAssemblerInterface::getIntPtSigma2,
                       &LocalAssemblerInterface::setIntPtSigma),
        std::make_pair(&LocalAssemblerInterface::getIntPtSigmaPrev2,
                       &LocalAssemblerInterface::setIntPtSigmaPrev)};

    auto basis = makeDisplacementBasis<DisplacementDim>(gridView);
    auto localView_element = basis.localView();

    // establish mapping from ids to indices
    std::unordered_map<std::size_t, std::size_t> map_id_to_idx;
    map_id_to_idx.reserve(gridView.size(0));

    // using GlobalCoordinate = typename
    // std::decay<decltype(gridView)>::type::GlobalCoordinate;
    using GlobalCoordinate = Dune::FieldVector<double, DisplacementDim>;

    struct ChildData
    {
        ChildData(std::size_t index_, std::vector<GlobalCoordinate>&& coords)
            : index(index_), int_pt_coords_in_father(std::move(coords))
        {
        }
        std::size_t index;  //!< index of the (child) element

        //! coordinates of the (child) element's integration points in the
        //! mapped to the reference element of the father
        std::vector<GlobalCoordinate> int_pt_coords_in_father;
    };

    std::unordered_map<std::size_t, std::vector<ChildData>>
        map_father_id_to_vanishing_child_data;

    // establish id to index mapping and save some extra data for elements that
    // might vanish
    for (auto const& element : Dune::elements(gridView))
    {
        assert(!element.isNew());
        map_id_to_idx.emplace(idSet.id(element), indexSet.index(element));
        assert(indexSet.index(element) < _local_assemblers.size());

        if (element.mightVanish() && element.hasFather())
        {
            // save quadrature point coords in father's reference element
            localView_element.bind(element);
            const auto& quadrature =
                get_quadrature<DisplacementDim>(element, localView_element);
            auto const num_int_pts = quadrature.size();

            std::vector<GlobalCoordinate> int_pt_coords_in_father;
            int_pt_coords_in_father.reserve(num_int_pts);

            typename std::decay<decltype(element)>::type father;
            for (unsigned ip = 0; ip < num_int_pts; ++ip)
            {
                father = element;
                auto int_pt_position = quadrature[ip].position();
                map_local_coordinates_to_father_recursive(
                    int_pt_position, father,
                    [](auto const& e) { return e.mightVanish(); });

                int_pt_coords_in_father.push_back(int_pt_position);
            }

            auto it_succ = map_father_id_to_vanishing_child_data.emplace(
                std::piecewise_construct,
                std::forward_as_tuple(idSet.id(father)),
                std::forward_as_tuple());
            it_succ.first->second.emplace_back(
                indexSet.index(element), std::move(int_pt_coords_in_father));
        }
    }

    _process_data.grid.adapt();

    basis.update(gridView);
    auto localView_father_ = basis.localView();

    using LDI =
        LocalDataInitializer<LocalAssemblerInterface,
                             SmallDeformationLocalAssembler, DisplacementDim,
                             decltype(basis), bool const, unsigned const,
                             SmallDeformationProcessData<DisplacementDim>&>;
    LDI initializer;

    auto const num_elements = gridView.size(0);

    decltype(_local_assemblers) local_assemblers(num_elements);

    // temporary storage
    std::vector<double> single_element_int_pt_data;
    std::vector<Dune::FieldVector<double, 1>> shape_functions;
    std::vector<double> all_shape_functions;
    std::vector<double> all_int_pt_data;

    // counters for statistical output
    std::size_t num_new_elements = 0;
    std::size_t num_coarsened_elements = 0;
    std::size_t num_kept_elements = 0;

    // Treat new elements
    for (auto& element : Dune::elements(gridView))
    {
        if (!element.isNew())
        {
            // "old" elements will be covered subsequently
            continue;
        }

        ++num_new_elements;

        auto const element_idx = indexSet.index(element);
        assert(!local_assemblers[element_idx]);

        local_assemblers[element_idx] =
            initializer(element, basis, false /*is_axially_symmetric*/,
                        0 /*integration_order*/, _process_data);
        auto& element_asm = *local_assemblers[element_idx];
        auto const num_element_int_pts =
            element_asm.getNumberOfIntegrationPoints();

        // quadrature variables
        localView_element.bind(element);
        const auto& quadrature_element =
            get_quadrature<DisplacementDim>(element, localView_element);
        assert(quadrature_element.size() == num_element_int_pts);

        auto father = element;
        do
        {
            father = father.father();
        } while (father.isNew());

        auto& localView_father =
            localView_father_;  // this way this loop and the "old" element loop
                                // use the same variable names
        localView_father.bind(father);
        // simplification: father's nodal data will all be interpolated with
        // child(0)
        auto const& localFiniteElement_father =
            localView_father.tree().child(0).finiteElement();

        auto const father_id = idSet.id(father);
        auto const father_id_idx_it = map_id_to_idx.find(father_id);

        /* Case 1: father easily found (father was a leaf element in the grid
         * before the current adaptation).
         *
         * Strategy:
         * 1. Extrapolate father's integration point data to father's nodes.
         * 2. Interpolate those nodal data to the integration points of the new
         *    (child) element.
         */
        if (father_id_idx_it != map_id_to_idx.end())
        {
            auto const father_idx = father_id_idx_it->second;
            auto const& father_asm = *_local_assemblers[father_idx];

            for (auto const& get_set : getters_setters)
            {
                auto const getter = get_set.first;
                auto const setter = get_set.second;

                auto const father_int_pt_data =
                    (father_asm.*getter)(single_element_int_pt_data);
#ifndef NDEBUG
                auto const father_num_int_pts =
                    father_asm.getNumberOfIntegrationPoints();
                auto const num_components =
                    KelvinVectorDimensions<DisplacementDim>::value;

                assert(father_int_pt_data.rows() == num_components);
                assert(father_int_pt_data.cols() == father_num_int_pts);
#endif
                // Extrapolate father's integration point data to nodes.
                // Note: The extrapolation minimizes the error in the l^2 norm,
                // i.e., every integration point has weight 1. This is
                // different, e.g., from the L^2 norm.
                // TODO [DUNE] don't use father's assembler
                // (extrapolateSingleElement) here
                auto const father_nodal_data = NumLib::extrapolateSingleElement(
                    father_asm, father_int_pt_data);

                // interpolate the nodal data just obtained
                auto element_int_pt_data =
                    interpolate_father_nodal_data_to_new_childs_integration_points(
                        father_nodal_data, element, quadrature_element,
                        localFiniteElement_father, shape_functions);

                (element_asm.*setter)(element_int_pt_data);
            }  // for each getter/setter
        }
        /* Case 2: father not easily found (father was a leaf element in the
         * grid before the current adaptation).
         *
         * Strategy:
         * 1. The father did not exist as a leaf before this adaptation. He must
         *    be part of a coarser level of the grid.
         * 2. Get all former children of the father that were leafs before the
         *    adaptation.
         * 3. Extrapolate the integration point data of those children to the
         *    nodes of the father.
         * 4. Interpolate those nodal data to the integration points of the new
         *    (child) element.
         */
        else
        {
            auto const father_id_and_child_idcs_it =
                map_father_id_to_vanishing_child_data.find(father_id);
            assert(father_id_and_child_idcs_it !=
                   map_father_id_to_vanishing_child_data.end());

            auto const& child_data = father_id_and_child_idcs_it->second;

            // collect all shape functions into one matrix
            auto const all_shape_functions_mat =
                get_father_shape_functions_at_integration_points_of_all_children(
                    child_data, shape_functions, all_shape_functions,
                    localFiniteElement_father);

            auto const decomposition =
                (all_shape_functions_mat.transpose() * all_shape_functions_mat)
                    .ldlt();

            for (auto const& get_set : getters_setters)
            {
                auto const getter = get_set.first;
                auto const setter = get_set.second;

                auto const all_int_pt_data_mat =
                    get_integration_point_data_of_all_children(
                        child_data, _local_assemblers, getter, all_int_pt_data,
                        single_element_int_pt_data);

                // Compute the extrapolated father nodal data from all
                // children's integration point data
                // Note: The extrapolation minimizes the error in the l^2 norm,
                // i.e., every integration point has weight 1. This is
                // different, e.g., from the L^2 norm.
                Eigen::MatrixXd const father_nodal_data =
                    decomposition
                        .solve(all_shape_functions_mat.transpose() *
                               all_int_pt_data_mat)
                        .transpose();
#ifndef NDEBUG
                auto const num_components =
                    KelvinVectorDimensions<DisplacementDim>::value;
                assert(father_nodal_data.rows() == num_components);
                // the following assertion holds only for linear shape
                // fcts.
                assert(father_nodal_data.cols() ==
                       father.subEntities(DisplacementDim));
#endif

                // interpolate the nodal data just obtained
                auto const element_int_pt_data =
                    interpolate_father_nodal_data_to_new_childs_integration_points(
                        father_nodal_data, element, quadrature_element,
                        localFiniteElement_father, shape_functions);

                (element_asm.*setter)(element_int_pt_data);
            }  // for each getter/setter
        }      // if father easily found
    }          // for each element

    // Treat coarsened elements
    for (auto& element : Dune::elements(gridView))
    {
        if (element.isNew())
        {
            // new elements already covered by the preceding loop
            continue;
        }

        auto const element_id = idSet.id(element);

        if (map_id_to_idx.find(element_id) != map_id_to_idx.end())
        {
            // "normal" old element, will be treated in the next loop
            continue;
        }

        ++num_coarsened_elements;

        /* Strategy:
         * 1. The element is neither new nor was it a leaf before this
         *    adaptation. Therefore it must be part of coarser grid level.
         * 2. Get all former children of the element that were leafs before the
         *    adaptation.
         * 3. Extrapolate the integration point data of those children to the
         *    nodes of the element.
         * 4. Interpolate those nodal data to the integration points of the
         *    element.
         */

        // TODO [DUNE] the following is mostly a 1:1 copy of the second part of
        // the adaptation loop above.

        auto const element_idx = indexSet.index(element);
        assert(!local_assemblers[element_idx]);

        local_assemblers[element_idx] =
            initializer(element, basis, false /*is_axially_symmetric*/,
                        0 /*integration_order*/, _process_data);
        auto& element_asm = *local_assemblers[element_idx];
        auto const num_element_int_pts =
            element_asm.getNumberOfIntegrationPoints();

        // quadrature variables
        localView_element.bind(element);
        const auto& quadrature_element =
            get_quadrature<DisplacementDim>(element, localView_element);
        assert(quadrature_element.size() == num_element_int_pts);

        auto const& father = element;

        auto const& localView_father = localView_element;
        // father's nodal data will all be interpolated with child(0)
        auto const& localFiniteElement_father =
            localView_father.tree().child(0).finiteElement();

        auto const father_id = element_id;

        auto const father_id_and_child_idcs_it =
            map_father_id_to_vanishing_child_data.find(father_id);
        assert(father_id_and_child_idcs_it !=
               map_father_id_to_vanishing_child_data.end());

        auto const& child_data = father_id_and_child_idcs_it->second;

        // collect all shape functions into one matrix
        auto const all_shape_functions_mat =
            get_father_shape_functions_at_integration_points_of_all_children(
                child_data, shape_functions, all_shape_functions,
                localFiniteElement_father);

        auto const decomposition =
            (all_shape_functions_mat.transpose() * all_shape_functions_mat)
                .ldlt();

        for (auto const& get_set : getters_setters)
        {
            auto const getter = get_set.first;
            auto const setter = get_set.second;

            auto const all_int_pt_data_mat =
                get_integration_point_data_of_all_children(
                    child_data, _local_assemblers, getter, all_int_pt_data,
                    single_element_int_pt_data);

            // compute the extrapolated father nodal data from all
            // children's integration point data
            // Note: The extrapolation minimizes the error in the l^2 norm,
            // i.e., every integration point has weight 1. This is
            // different, e.g., from the L^2 norm.
            Eigen::MatrixXd const father_nodal_data =
                decomposition
                    .solve(all_shape_functions_mat.transpose() *
                           all_int_pt_data_mat)
                    .transpose();
#ifndef NDEBUG
            auto const num_components =
                KelvinVectorDimensions<DisplacementDim>::value;
            assert(father_nodal_data.rows() == num_components);
            // the following assertion holds only for linear shape fcts.
            assert(father_nodal_data.cols() ==
                   father.subEntities(DisplacementDim));
#endif

            // interpolate the nodal data just obtained
            auto const element_int_pt_data =
                /* note: this function is essential and different from the
                   procedure above */
                interpolate_nodal_data_to_integration_points(
                    father_nodal_data, quadrature_element,
                    localFiniteElement_father, shape_functions);

            (element_asm.*setter)(element_int_pt_data);
        }  // for each getter/setter
    }

    // Treat old elements. This must be done after the new/coarsened elements,
    // because otherwise access to old local assemblers that have been moved
    // away will lead to segmentation faults.
    for (auto& element : Dune::elements(gridView))
    {
        if (element.isNew())
        {
            // new elements already covered by the pre-preceding loop
            continue;
        }

        auto const element_id = idSet.id(element);
        auto const element_idx = indexSet.index(element);

        auto const it = map_id_to_idx.find(element_id);
        if (it == map_id_to_idx.end())
        {
            // coarsened element, was treated in the preceding loop
            continue;
        }

        ++num_kept_elements;

        local_assemblers[element_idx] =
            std::move(_local_assemblers[it->second]);
    }

    DBUG("refinement element counts: %d total, %d new, %d coarsened, %d kept.",
         num_new_elements + num_coarsened_elements + num_kept_elements,
         num_new_elements, num_coarsened_elements, num_kept_elements);

    _local_assemblers = std::move(local_assemblers);

    _process_data.grid.postAdapt();

    for (auto& e : Dune::elements(gridView))
    {
        assert(!e.isNew());
    }

    return true;
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::updateDOFtables(
    MeshLib::DUNEMesh<DisplacementDim> const&, bool,
    MeshLib::DUNEIdToIdxMappings const&)
{
    auto gridView = _process_data.grid.getMesh().leafGridView();

    using DT = decltype(makeDisplacementBasis<DisplacementDim>(gridView));
    auto dt = dynamic_cast<NumLib::DOFTableDUNE<DT>*>(
        _local_to_global_index_map.get());
    OGS_ALWAYS_ASSERT(dt != nullptr);
    dt->update(gridView);

    using DT1 = decltype(makeScalarBasis<DisplacementDim>(gridView));
    auto* dt1 = dynamic_cast<NumLib::DOFTableDUNE<DT1>*>(
        &_extrapolator_data.getDOFTable());
    OGS_ALWAYS_ASSERT(dt1 != nullptr);
    dt1->update(gridView);
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::createLocalAssemblers(
    MeshLib::DUNEMesh<DisplacementDim> const& mesh, bool global_refine,
    MeshLib::DUNEIdToIdxMappings const&)
{
    if (!global_refine)
        return;

    OGS_ALWAYS_ASSERT(&mesh == &_process_data.grid);
    auto gridView = mesh.getMesh().leafGridView();

    namespace DFB = Dune::Functions::BasisBuilder;
    auto basis = makeDisplacementBasis<DisplacementDim>(gridView);
    DBUG("Basis type %s\nof size %d.", Dune::className(basis).c_str(),
         sizeof(basis));

    DBUG("Create local assemblers.");

    ProcessLib::SmallDeformationDUNE::createLocalAssemblers<
        DisplacementDim, decltype(basis), SmallDeformationLocalAssembler>(
        basis, _local_assemblers, false /* axially symmetric */,
        static_cast<unsigned>(2) /* integration order */, _process_data);
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::assembleConcreteProcess(
    const double /*t*/, GlobalVector const& /*x*/, GlobalMatrix& /*M*/,
    GlobalMatrix& /*K*/, GlobalVector& /*b*/)
{
    DBUG("Assemble SmallDeformationProcess.");

    OGS_FATAL("not implemented.");
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(const double t, GlobalVector const& x,
                                        GlobalVector const& xdot,
                                        const double dxdot_dx,
                                        const double dx_dx, GlobalMatrix& M,
                                        GlobalMatrix& K, GlobalVector& b,
                                        GlobalMatrix& Jac)
{
    DBUG("AssembleWithJacobian SmallDeformationProcess.");

    auto gridView = _process_data.grid.getMesh().leafGridView();
    auto basis = makeDisplacementBasis<DisplacementDim>(gridView);

    // Shortcut: except CrankNicolson, matrices are never stored in
    // intermediate stages. Therefore we don't need to add pre/post
    // refinement hooks to EigenMatrix. Instead we adapt the matrices here.
    if (static_cast<std::size_t>(Jac.getRawMatrix().rows()) != basis.size())
    {
        Jac = std::move(*MathLib::MatrixVectorTraits<GlobalMatrix>::newInstance(
            getMatrixSpecifications(0 /* TODO [DUNE] fixme */)));
        Jac.setZero();
        M = Jac;
        K = Jac;
    }

    auto localView = basis.localView();
    auto localIndexSet = basis.localIndexSet();

    auto e_it = Dune::elements(gridView).begin();
    auto const e_end = Dune::elements(gridView).end();
    auto l_it = std::begin(_local_assemblers);
    auto const l_end = std::end(_local_assemblers);

    for (; e_it != e_end && l_it != l_end; ++e_it, ++l_it)
    {
        localView.bind(*e_it);
        localIndexSet.bind(localView);
        _global_assembler_dune.assembleWithJacobian(
            *e_it, **l_it, basis, localView, localIndexSet, t, x, xdot,
            dxdot_dx, dx_dx, M, K, b, Jac, _coupled_solutions);
    }

    // TODO [DUNE] re-enable
#if 0
    b.copyValues(*_nodal_forces);
    std::transform(_nodal_forces->begin(), _nodal_forces->end(),
                   _nodal_forces->begin(), [](double val) { return -val; });
#endif
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::preTimestepConcreteProcess(
    GlobalVector const& /*x*/, double const t, double const dt,
    const int /*process_id*/)
{
    DBUG("PreTimestep SmallDeformationProcess.");

    _process_data.dt = dt;
    _process_data.t = t;

    for (auto& l : _local_assemblers)
    {
        l->preTimestep();
    }
}

}  // namespace SmallDeformationDUNE
}  // namespace ProcessLib
