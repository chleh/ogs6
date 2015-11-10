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

PardisoLinearSolver::PardisoLinearSolver(const BaseLib::ConfigTree &/*option*/)
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
::solve(LOLMatrix &A, DenseVector<double> &rhs, DenseVector<double> &result)
{
    auto A_csr = toCSR(A);
    double* a   = const_cast<double*>(A_csr.getValues().data());
    MKL_INT* ia = const_cast<int*   >(A_csr.getRowIndices().data());
    MKL_INT* ja = const_cast<int*   >(A_csr.getColumnIndices().data());
    MKL_INT mtype = 11;		/* Real unsymmetric matrix */

    void* pt[64]; for (int i=0; i<64; ++i) { pt[i] = 0; }

    /* -------------------------------------------------------------------- */
    /* .. Setup Pardiso control parameters. */
    /* -------------------------------------------------------------------- */
    MKL_INT iparm[64] = { 0 };
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

    MKL_INT maxfct = 1;			/* Maximum number of numerical factorizations. */
    MKL_INT mnum = 1;			/* Which factorization to use. */
    MKL_INT msglvl = 1;			/* Print statistical information in file */
    MKL_INT error = 0;			/* Initialize error flag */

    MKL_INT nrhs = 1;
    MKL_INT nrows = A_csr.getNRows(); // number of rows
    assert(nrows == static_cast<MKL_INT>(rhs.size()));

    /* Auxiliary variables. */
    double ddum;			/* Double dummy */
    MKL_INT idum;			/* Integer dummy. */

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

    std::vector<double> bsv(rhs.size());
    double* bs = bsv.data();

    double* x = const_cast<double*>(&result[0]);

    double res, res0;

    iparm[11] = 0;
    DBUG("\nSolving system...");
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
            &nrows, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
    if (error != 0)
    {
        ERR("ERROR during solution: %d", error);
        std::abort();
    }

    // Compute residual
    char* uplo = (char*) "non-transposed";
    mkl_dcsrgemv (uplo, &nrows, a, ia, ja, x, bs);
    res = 0.0;
    res0 = 0.0;
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

    /* -------------------------------------------------------------------- */
    /* .. Termination and release of memory. */
    /* -------------------------------------------------------------------- */
    phase = -1;			/* Release internal memory. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
            &nrows, &ddum, ia, ja, &idum, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);

}

}
