/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "logog/include/logog.hpp"

#include "PardisoLinearSolver.h"

#include "MathLib/LinAlg/Sparse/LOLCSRTools.h"

extern "C"
{
#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_spblas.h"
}


namespace MathLib
{


class PardisoLinearSolverImpl
{
public:
    void solve(LOLMatrix &A, DenseVector<double> &rhs, DenseVector<double> &result);

    ~PardisoLinearSolverImpl();

private:
    void initParams();

    void*   pt[64]    = { 0 };
    MKL_INT iparm[64] = { 0 };

    MKL_INT mtype;		/* Real unsymmetric matrix */

    MKL_INT maxfct;			/* Maximum number of numerical factorizations. */
    MKL_INT mnum;			/* Which factorization to use. */
    MKL_INT msglvl;			/* Print statistical information in file */
    MKL_INT error;			/* Initialize error flag */

    MKL_INT nrhs;
    MKL_INT nrows; // number of rows

    /* Auxiliary variables. */
    double  ddum;			/* Double dummy */
    MKL_INT idum;			/* Integer dummy. */

};

void PardisoLinearSolverImpl::initParams()
{
    /* -------------------------------------------------------------------- */
    /* .. Setup Pardiso control parameters. */
    /* -------------------------------------------------------------------- */
    iparm[0] = 1;			/* No solver default */
    iparm[1] = 2;			/* Fill-in reordering from METIS */
    /* Numbers of processors, value of OMP_NUM_THREADS */
    iparm[2] = 1;
    iparm[3] = 0;			/* No iterative-direct algorithm */
    iparm[4] = 0;			/* No user fill-in reducing permutation */
    iparm[5] = 0;			/* Write solution into x */
    iparm[6] = 0;			/* Not in use */
    iparm[7] = 2;			/* Max numbers of iterative refinement steps */
    iparm[8] = 0;			/* Not in use */
    iparm[9] = 13;		/* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;		/* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;		/* Conjugate transposed/transpose solve */
    iparm[12] = 1;		/* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
    iparm[13] = 0;		/* Output: Number of perturbed pivots */
    iparm[14] = 0;		/* Not in use */
    iparm[15] = 0;		/* Not in use */
    iparm[16] = 0;		/* Not in use */
    iparm[17] = -1;		/* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1;		/* Output: Mflops for LU factorization */
    iparm[19] = 0;		/* Output: Numbers of CG Iterations */
    iparm[34] = 1;      /* zero based indexing */

    mtype = 11;		/* Real unsymmetric matrix */

    maxfct = 1;			/* Maximum number of numerical factorizations. */
    mnum = 1;			/* Which factorization to use. */
    msglvl = 1;			/* Print statistical information in file */
    error = 0;			/* Initialize error flag */

    nrhs = 1;
}


void
PardisoLinearSolverImpl::
solve(LOLMatrix &A, DenseVector<double> &rhs, DenseVector<double> &result)
{
    initParams();

    auto A_csr = toCSR(A); // TODO: do not always allocate
    double* a   = const_cast<double*>(A_csr.getValues().data());
    MKL_INT* ia = const_cast<int*   >(A_csr.getRowIndices().data());
    MKL_INT* ja = const_cast<int*   >(A_csr.getColumnIndices().data());

    nrows = A_csr.getNRows();
    assert(nrows == static_cast<MKL_INT>(rhs.size()));

    // TODO: needs only be done once
    /* -------------------------------------------------------------------- */
    /* .. Reordering and Symbolic Factorization. This step also allocates */
    /* all memory that is necessary for the factorization. */
    /* -------------------------------------------------------------------- */
    MKL_INT phase = 11;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
            &nrows, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if (error != 0)
    {
        ERR("ERROR during symbolic factorization: %d", error);
        std::abort();
    }
    DBUG("Reordering completed ... ");
    DBUG("Number of nonzeros in factors = %d", iparm[17]);
    DBUG("Number of factorization MFLOPS = %d\n", iparm[18]);


    /* -------------------------------------------------------------------- */
    /* .. Numerical factorization. */
    /* -------------------------------------------------------------------- */
    phase = 22;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
            &nrows, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if (error != 0)
    {
        ERR("ERROR during numerical factorization: %d", error);
        std::abort();
    }
    DBUG("Factorization completed ... ");


    /* -------------------------------------------------------------------- */
    /* .. Back substitution and iterative refinement. */
    /* -------------------------------------------------------------------- */
    phase = 33;

    double* b = const_cast<double*>(&rhs[0]);

    double* x = const_cast<double*>(&result[0]);

    iparm[11] = 0;
    DBUG("\nSolving system...");
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
            &nrows, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
    if (error != 0)
    {
        ERR("ERROR during solution: %d", error);
        std::abort();
    }

#ifndef NDEBUG
    // Compute residual
    std::vector<double> bsv(rhs.size());
    double* bs = bsv.data();

    char* uplo = (char*) "non-transposed";
    mkl_dcsrgemv (uplo, &nrows, a, ia, ja, x, bs);

    double res = 0.0;
    double res0 = 0.0;
    for (int j = 1; j <= nrows; j++)
    {
        res += (bs[j - 1] - b[j - 1]) * (bs[j - 1] - b[j - 1]);
        res0 += b[j - 1] * b[j - 1];
    }
    res = sqrt (res) / sqrt (res0);

    DBUG("Relative residual = %e", res);
    // Check residual
    /*
    if (res > 1e-10)
    {
        printf ("Error: residual is too high!\n");
        exit (10);
    }
    */
#endif
}

PardisoLinearSolverImpl::~PardisoLinearSolverImpl()
{
    /* -------------------------------------------------------------------- */
    /* .. Termination and release of memory. */
    /* -------------------------------------------------------------------- */
    MKL_INT phase = -1;			/* Release internal memory. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
            &nrows, &ddum, &idum, &idum, &idum, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);
    // TODO: hopefully there won't be leaks here, since I have omitted row/column arrays
    //       ia and ja.
}


PardisoLinearSolver::PardisoLinearSolver(
        LOLMatrix &A, const std::string /*solver_name*/,
        BaseLib::ConfigTree const*const /*option*/)
    : _A{A}
    , _data{new PardisoLinearSolverImpl}
{
    /*
    boost::optional<std::string> solver_type
            = option.get_optional<std::string>("solver_type");
    if (solver_type)
    {

    }
    else
    {
        ERR("option <solver_type> not given");
    }
    */
}

void PardisoLinearSolver
::solve(DenseVector<double> &rhs, DenseVector<double> &result)
{
    _data->solve(_A, rhs, result);
}


PardisoLinearSolver::~PardisoLinearSolver()
{}

}
