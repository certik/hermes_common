// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "matrix.h"
#include "solvers.h"

#ifdef COMMON_WITH_SUPERLU
#include <superlu/slu_ddefs.h>

void solve_linear_system_superlu(Matrix *mat, double *res)
{
    CooMatrix *mcoo = dynamic_cast<CooMatrix*>(mat);
    CSCMatrix mcsc(mcoo);

    int nnz = mcsc.get_nnz();
    int *Ap = mcsc.get_Ap();
    int *Ai = mcsc.get_Ai();
    double *Ax = mcsc.get_A();

    SuperMatrix A;
    SuperMatrix B;
    SuperMatrix L;      /* factor L */
    SuperMatrix U;      /* factor U */

    int nrhs, info;

    superlu_options_t options;
    SuperLUStat_t stat;

    // Set the default input options:
    /*
    options.Fact = DOFACT;
    options.Equil = YES;
    options.ColPerm = COLAMD;
    options.DiagPivotThresh = 1.0;
    options.Trans = NOTRANS;
    options.IterRefine = NOREFINE;
    options.SymmetricMode = NO;
    options.PivotGrowth = NO;
    options.ConditionNumber = NO;
    options.PrintStat = YES;
    */

    set_default_options(&options);

    // create csc matrix
    dCreate_CompCol_Matrix(&A, mcsc.get_size(), mcsc.get_size(), nnz, Ax, Ai, Ap, SLU_NC, SLU_D, SLU_GE);
    // dPrint_CompCol_Matrix("A", &A);

    nrhs = 1;
    // create rhs matrix
    dCreate_Dense_Matrix(&B, mcsc.get_size(), nrhs, res, mcsc.get_size(), SLU_DN, SLU_D, SLU_GE);
    // dPrint_Dense_Matrix("B", &B);

    // column permutation vector
    int *perm_c = intMalloc(mcsc.get_size());
    // row permutations from partial pivoting
    int *perm_r = intMalloc(mcsc.get_size());
    if (!perm_c) ABORT("Malloc fails for perm_c[].");
    if (!perm_r) ABORT("Malloc fails for perm_r[].");

    // initialize the statistics variables
    StatInit(&stat);
    dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

    mem_usage_t mem_usage;
    if ( info == 0 )
    {
        // solution
        double *x = (double*) ((DNformat*) B.Store)->nzval;

        // copy result
        memcpy(res, x, mcsc.get_size()*sizeof(double));

        // print_vector("X", x, mcsc.get_size());
        /*
        SCformat *Lstore = (SCformat *) L.Store;
        NCformat *Ustore = (NCformat *) U.Store;
        printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
        printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
        printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - mcsc.get_size());
        printf("FILL ratio = %.1f\n", (float)(Lstore->nnz + Ustore->nnz - mcsc.get_size())/nnz);

        mem_usage_t mem_usage;
        dQuerySpace(&L, &U, &mem_usage);
        printf("L\\U MB %.3f\ttotal MB needed %.3f\n", mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
        */
    }
    else
    {
        printf("dgssv() error returns INFO = %d\n", info);
        if ( info <= mcsc.get_size() )
        {
            // factorization completes
            dQuerySpace(&L, &U, &mem_usage);
            printf("L\\U MB %.3f\ttotal MB needed %.3f\n", mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
        }
    }

    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    /*
    FIXME
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    */
}

#else

void solve_linear_system_superlu(Matrix *mat, double *res)
{
    _error("hermes_common: solve_linear_system_superlu - SUPERLU is not available.");
}
#endif
