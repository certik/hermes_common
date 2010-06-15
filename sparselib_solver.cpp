// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "matrix.h"
#include "solvers.h"

#ifdef COMMON_WITH_SPARSELIB
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

bool CommonSolverSparseLib::solve(Matrix *mat, double *res)
{    
    CooMatrix *mcoo = dynamic_cast<CooMatrix*>(mat);

    int nnz = mcoo->get_nnz();
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
    VECTOR_double xv = ILU.solve(rhs);

    // method
    int result = -1;
    switch (method)
    {
    case CommonSolverSparseLibSolver_ConjugateGradientSquared:
        result = CGS(Acsc, xv, rhs, ILU, maxiter, tolerance);
        break;
    case CommonSolverSparseLibSolver_RichardsonIterativeRefinement:
        result = IR(Acsc, xv, rhs, ILU, maxiter, tolerance);
        break;
    default:
        _error("SparseLib error. Method is not defined.");
    }

    if (result == 0)
        printf("SparseLib solver: maxiter: %i, tol: %e\n", maxiter, tolerance);
    else
        _error("SparseLib error.");

    delete[] data;
    delete[] row;
    delete[] col;

    double *x;
    x = (double*) malloc(mcoo->get_size() * sizeof(double));

    for (int i = 0 ; i < xv.size() ; i++)
        x[i] = xv(i);

    memcpy(res, x, mcoo->get_size()*sizeof(double));
}

bool CommonSolverSparseLib::solve(Matrix *mat, cplx *res)
{
    _error("CommonSolverSparseLib::solve(Matrix *mat, cplx *res) not implemented.");
}

#else

bool CommonSolverSparseLib::solve(Matrix *mat, double *res)
{
    _error("CommonSolverSparseLib::solve(Matrix *mat, double *res) not implemented.");
}

bool CommonSolverSparseLib::solve(Matrix *mat, cplx *res)
{
    _error("CommonSolverSparseLib::solve(Matrix *mat, cplx *res) not implemented.");
}

#endif
