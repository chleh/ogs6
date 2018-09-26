/** This file is taken from Oliver Sander's book draft! It must not be published
 * as is!
 */

#if HAVE_CONFIG_H
#include <config.h>
#endif  // HAVE_CONFIG_H

#include "dune-config.h"

#include <vector>

#include <dune/common/function.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/matrix.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/matrixmarket.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/pqknodalbasis.hh>

using namespace Dune;

// Compute the stiffness matrix for a single element
template <class LocalView, class MatrixType>
void assembleElementStiffnessMatrix(const LocalView& localView,
                                    MatrixType& elementMatrix)
{
    using Element = typename LocalView::Element;
    auto element = localView.element();
    const int dim = Element::dimension;
    auto geometry = element.geometry();  // Trafo from ref to actual element

    // Get set of shape functions for this element
    const auto& localFiniteElement = localView.tree().finiteElement();

    // Set all matrix entries to zero
    elementMatrix.setSize(localFiniteElement.size(), localFiniteElement.size());
    elementMatrix = 0;  // fills the entire matrix with zeros

    // Get a quadrature rule
    int order = 2 * (dim * localFiniteElement.localBasis().order() - 1);
    const auto& quad =
        QuadratureRules<double, dim>::rule(element.type(), order);

    // Loop over all quadrature points
    for (size_t pt = 0; pt < quad.size(); pt++)
    {
        // Position of the current quadrature point in the reference element
        const auto quadPos = quad[pt].position();

        // The transposed inverse Jacobian of the map from the reference element
        // to the element
        const auto jacobian = geometry.jacobianInverseTransposed(quadPos);

        // The multiplicative factor in the integral transformation formula
        const auto integrationElement = geometry.integrationElement(quadPos);

        // The gradients of the shape functions on the reference element
        std::vector<FieldMatrix<double, 1, dim>> referenceGradients;
        localFiniteElement.localBasis().evaluateJacobian(quadPos,
                                                         referenceGradients);

        // Compute the shape function gradients on the real element
        std::vector<FieldVector<double, dim>> gradients(
            referenceGradients.size());
        for (size_t i = 0; i < gradients.size(); i++)
            jacobian.mv(referenceGradients[i][0], gradients[i]);

        // Compute the actual matrix entries
        for (size_t i = 0; i < elementMatrix.N(); i++)
            for (size_t j = 0; j < elementMatrix.M(); j++)
                elementMatrix[localView.tree().localIndex(i)]
                             [localView.tree().localIndex(j)] +=
                    (gradients[i] * gradients[j]) * quad[pt].weight() *
                    integrationElement;
    }
}

// Compute the source term for a single element
template <class LocalView>
void getVolumeTerm(const LocalView& localView,
                   BlockVector<FieldVector<double, 1>>& localRhs,
                   const std::function<double(
                       FieldVector<double, LocalView::Element::dimension>)>
                       volumeTerm)
{
    using Element = typename LocalView::Element;
    auto element = localView.element();
    const int dim = Element::dimension;  // TODO compile-time constant!

    // Set of shape functions for a single element
    const auto& localFiniteElement = localView.tree().finiteElement();

    // Set all entries to zero
    localRhs.resize(localFiniteElement.size());
    localRhs = 0;

    // A quadrature rule
    int order = dim;
    const QuadratureRule<double, dim>& quad =
        QuadratureRules<double, dim>::rule(element.type(), order);

    // Loop over all quadrature points
    for (size_t pt = 0; pt < quad.size(); pt++)
    {
        // Position of the current quadrature point in the reference element
        const FieldVector<double, dim>& quadPos = quad[pt].position();

        // The multiplicative factor in the integral transformation formula
        const double integrationElement =
            element.geometry().integrationElement(quadPos);

        double functionValue = volumeTerm(element.geometry().global(quadPos));

        // Evaluate all shape function values at this point
        std::vector<FieldVector<double, 1>> shapeFunctionValues;
        localFiniteElement.localBasis().evaluateFunction(quadPos,
                                                         shapeFunctionValues);

        // Actually compute the vector entries
        for (size_t i = 0; i < localRhs.size(); i++)
            localRhs[i] += shapeFunctionValues[i] * functionValue *
                           quad[pt].weight() * integrationElement;
    }
}

// Get the occupation pattern of the stiffness matrix
template <class Basis>
void getOccupationPattern(const Basis& basis, MatrixIndexSet& nb)
{
    nb.resize(basis.size(), basis.size());

    auto gridView = basis.gridView();

    // A loop over all elements of the grid
    auto localView = basis.localView();
    auto localIndexSet = basis.localIndexSet();

    for (const auto& element : elements(gridView))
    {
        localView.bind(element);
        localIndexSet.bind(localView);

        for (size_t i = 0; i < localIndexSet.size(); i++)
        {
            // The global index of the i-th vertex of the element
            auto row = localIndexSet.index(i);

            for (size_t j = 0; j < localIndexSet.size(); j++)
            {
                // The global index of the j-th vertex of the element
                auto col = localIndexSet.index(j);
                nb.add(row, col);
            }
        }
    }
}

/** \brief Assemble the Laplace stiffness matrix on the given grid view */
template <class Basis>
void assemblePoissonProblem(
    const Basis& basis,
    BCRSMatrix<FieldMatrix<double, 1, 1>>& matrix,
    BlockVector<FieldVector<double, 1>>& rhs,
    const std::function<double(FieldVector<double, Basis::GridView::dimension>)>
        volumeTerm)
{
    auto gridView = basis.gridView();

    // MatrixIndexSets store the occupation pattern of a sparse matrix.
    // They are not particularly efficient, but simple to use.
    MatrixIndexSet occupationPattern;
    getOccupationPattern(basis, occupationPattern);

    // ... and give it the occupation pattern we want.
    occupationPattern.exportIdx(matrix);

    // set rhs to correct length
    rhs.resize(basis.size());

    // Set all entries to zero
    matrix = 0;
    rhs = 0;

    // A loop over all elements of the grid //////////////////////////
    // TODO ??? How to store element data?
    auto localView = basis.localView();
    auto localIndexSet = basis.localIndexSet();

    // TODO Note: element has the same type for all elements in the gridView!
    for (const auto& element : elements(gridView))
    {
        // Now let's get the element stiffness matrix
        // A dense matrix is used for the ***element*** stiffness matrix
        localView.bind(element);
        localIndexSet.bind(localView);

        Matrix<FieldMatrix<double, 1, 1>> elementMatrix;
        assembleElementStiffnessMatrix(localView, elementMatrix);

        // add element matrix to global matrix
        for (size_t i = 0; i < elementMatrix.N(); i++)
        {
            // The global index of the i-th vertex of the element
            auto row = localIndexSet.index(i);

            for (size_t j = 0; j < elementMatrix.M(); j++)
            {
                // The global index of the j-th vertex of the element
                auto col = localIndexSet.index(j);
                matrix[row][col] += elementMatrix[i][j];
            }
        }

        // Now get the local contribution to the right-hand side vector
        BlockVector<FieldVector<double, 1>> localRhs;
        getVolumeTerm(localView, localRhs, volumeTerm);

        for (size_t i = 0; i < localRhs.size(); i++)
        {
            // The global index of the i-th vertex of the element
            auto row = localIndexSet.index(i);
            rhs[row] += localRhs[i];
        }
    }
}

int main(int argc, char* argv[]) try
{
    // Set up MPI, if available
    MPIHelper::instance(argc, argv);

    // ////////////////////////////////
    //   Generate the grid
    // ////////////////////////////////

    if (argc != 2)
    {
        std::cerr << "exactly one argument taken\n";
        return 1;
    }
    const int dim = 2;
    typedef UGGrid<dim> GridType;
    std::shared_ptr<GridType> grid(GmshReader<GridType>::read(argv[1]));

    grid->globalRefine(2);

    typedef GridType::LeafGridView GridView;
    GridView gridView = grid->leafGridView();

    // ///////////////////////////////////////////////////////
    //   Stiffness matrix and right hand side vector
    // ///////////////////////////////////////////////////////

    typedef BlockVector<FieldVector<double, 1>> VectorType;
    typedef BCRSMatrix<FieldMatrix<double, 1, 1>> MatrixType;

    VectorType rhs;
    MatrixType stiffnessMatrix;

    // ///////////////////////////////////////////////////////
    //   Assemble the system
    // ///////////////////////////////////////////////////////

    // 1st order Lagrange elements
    Functions::PQkNodalBasis<GridView, 1> basis(gridView);

    auto rightHandSide = [](const FieldVector<double, dim>& /*x*/) {
        return -5.0;
    };
    assemblePoissonProblem(basis, stiffnessMatrix, rhs, rightHandSide);

    // Determine Dirichlet dofs by marking all degrees of freedom whose Lagrange
    // nodes comply with a given predicate. Since interpolating into a
    // vector<bool> is currently not supported, we use a vector<char> which, in
    // contrast to vector<bool> is a real container.
    auto predicate = [](auto p) {
        return p[0] < 1e-8 or p[1] < 1e-8 or (p[0] > 0.4999 and p[1] > 0.4999);
    };

    // Interpolating the predicate will mark all desired Dirichlet degrees of
    // freedom
    std::vector<char> dirichletNodes;
    Functions::interpolate(basis, dirichletNodes, predicate);

    // /////////////////////////////////////////
    //   Modify Dirichlet rows
    // /////////////////////////////////////////
    // loop over the matrix rows
    for (size_t i = 0; i < stiffnessMatrix.N(); i++)
    {
        // set diagonal of dirichlet rows to one
        if (dirichletNodes[i])
        {
            auto cIt = stiffnessMatrix[i].begin();
            auto cEndIt = stiffnessMatrix[i].end();
            // loop over nonzero matrix entries in current row
            for (; cIt != cEndIt; ++cIt)
                *cIt = (i == cIt.index()) ? 1.0 : 0.0;
        }
    }

    // Set Dirichlet values
    auto dirichletValues = [](auto p) {
        return (p[0] < 1e-8 or p[1] < 1e-8) ? 0 : 0.5;
    };
    Functions::interpolate(basis, rhs, dirichletValues, dirichletNodes);

    /////////////////////////////////////////////////////////////////////////////
    // Write matrix and load vector to files, to be used in later examples
    /////////////////////////////////////////////////////////////////////////////
    storeMatrixMarket(stiffnessMatrix, "poisson-matrix.mm");
    storeMatrixMarket(rhs, "poisson-rhs.mm");

    // /////////////////////////
    //   Compute solution
    // /////////////////////////

    // Choose an initial iterate
    VectorType x(basis.size());
    x = 0;

    // Technicality:  turn the matrix into a linear operator
    MatrixAdapter<MatrixType, VectorType, VectorType> linearOperator(
        stiffnessMatrix);

    // Sequential incomplete LU decomposition as the preconditioner
    SeqILU0<MatrixType, VectorType, VectorType> preconditioner(stiffnessMatrix,
                                                               1.0);

    // Preconditioned conjugate-gradient solver
    CGSolver<VectorType> cg(linearOperator,
                            preconditioner,
                            1e-4,  // desired residual reduction factor
                            50,    // maximum number of iterations
                            2);    // verbosity of the solver

    // Object storing some statistics about the solving process
    InverseOperatorResult statistics;

    // Solve!
    cg.apply(x, rhs, statistics);

    // Output result
    VTKWriter<GridView> vtkWriter(gridView);
    vtkWriter.addVertexData(x, "solution");
    vtkWriter.write("poissonequation_result");

    return 0;
}
// Error handling
catch (std::exception& e)
{
    std::cout << e.what() << std::endl;
    return 1;
}
