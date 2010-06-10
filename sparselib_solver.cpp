// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "matrix.h"
#include "solvers.h"

// #ifdef COMMON_WITH_SPARSELIB
#include <coord_double.h>
#include <compcol_double.h>
#include <mvvd.h>
#include <ilupre_double.h>
#include <bicg.h>
#include <cg.h>
#include <cgs.h>
#include <bicgstab.h>
#include <cheby.h>
#include <gmres.h>
#include <ir.h>
#include <qmr.h>

void solve_linear_system_sparselib(Matrix *mat, double *res)
{
    CooMatrix *mcoo = dynamic_cast<CooMatrix*>(mat);

    int nnz = mcoo->triplets_len();
    int *row = new int[nnz];
    int *col = new int[nnz];
    double *data = new double[nnz];

    mcoo->get_row_col_data(row, col, data);

    // matrix
    Coord_Mat_double Acoo(mcoo->get_size(), mcoo->get_size(), nnz, data, row, col);
    CompCol_Mat_double Acsc(Acoo);

    // rhs
    VECTOR_double rhs(res, mcoo->get_size());

    // preconditioner
    CompCol_ILUPreconditioner_double ILU(Acsc);

    // solve
    VECTOR_double xv = ILU.solve(rhs);

    int maxiter = 1000;
    double tol = 1e-6;

    int result = CGS(Acsc, xv, rhs, ILU, maxiter, tol);
    if (result == 0)
        printf("SparseLib CGS: maxiter: %i, tol: %f\n", maxiter, tol);
    else
        _error("SparseLib error.");

    delete[] data;
    delete[] row;
    delete[] col;

    double *x;
    x = (double*) malloc(mcoo->get_size() * sizeof(double));

    for (int i = 0 ; i < xv.size() ; i++)
    {
        // printf("(%i): %f\n", i, xv(i));
        x[i] = xv(i);
    }

    memcpy(res, x, mcoo->get_size()*sizeof(double));
}

/*
#else

void solve_linear_system_sparselib(Matrix *mat, double *res)
{
    _error("hermes_common: solve_linear_system_sparselib - SPARSELIB is not available.");
}
#endif
*/

