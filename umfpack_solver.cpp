// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "matrix.h"
#include "solvers.h"

#ifdef COMMON_WITH_UMFPACK
#include <umfpack.h>

double control_array[UMFPACK_CONTROL];
double info_array[UMFPACK_INFO];

void print_status(int status)
{
    switch (status)
    {
    case UMFPACK_OK:
        break;
    case UMFPACK_WARNING_singular_matrix:
        _error("UMFPACK: singular stiffness matrix!");
        break;

    case UMFPACK_ERROR_out_of_memory:           _error("UMFPACK: out of memory!");
    case UMFPACK_ERROR_argument_missing:        _error("UMFPACK: argument missing");
    case UMFPACK_ERROR_invalid_Symbolic_object: _error("UMFPACK: invalid Symbolic object");
    case UMFPACK_ERROR_invalid_Numeric_object:  _error("UMFPACK: invalid Numeric object");
    case UMFPACK_ERROR_different_pattern:       _error("UMFPACK: different pattern");
    case UMFPACK_ERROR_invalid_system:          _error("UMFPACK: invalid system");
    case UMFPACK_ERROR_n_nonpositive:           _error("UMFPACK: n nonpositive");
    case UMFPACK_ERROR_invalid_matrix:          _error("UMFPACK: invalid matrix");
    case UMFPACK_ERROR_internal_error:          _error("UMFPACK: internal error");
    default:                                    _error("UMFPACK: unknown error");
    }
}

void solve_linear_system_umfpack(Matrix *mat, double *res)
{
    CooMatrix *mcoo = dynamic_cast<CooMatrix*>(mat);

    int nnz = mcoo->triplets_len();
    int *row = new int[nnz];
    int *col = new int[nnz];
    double *data = new double[nnz];

    mcoo->get_row_col_data(row, col, data);

    int *Ap = new int[mcoo->get_size() + 1];
    int *Ai = new int[nnz];
    double *Ax = new double[nnz];

    umfpack_di_defaults(control_array);

    /* convert matrix from triplet form to compressed-column form */
    int status_triplet_to_col = umfpack_di_triplet_to_col(mcoo->get_size(), mat->get_size(),
                                                          nnz, row, col, data, Ap, Ai, Ax, NULL);
    print_status(status_triplet_to_col);

    /* symbolic analysis */
    void *symbolic, *numeric;
    int status_symbolic = umfpack_di_symbolic(mcoo->get_size(), mcoo->get_size(),
                                              Ap, Ai, NULL, &symbolic,
                                              control_array, info_array);
    print_status(status_symbolic);

    /* LU factorization */
    int status_numeric = umfpack_di_numeric(Ap, Ai, Ax, symbolic, &numeric,
                                            control_array, info_array);
    print_status(status_numeric);

    umfpack_di_free_symbolic(&symbolic);

    double *x;
    x = (double*) malloc(mcoo->get_size() * sizeof(double));

    /* solve system */
    int status_solve = umfpack_di_solve(UMFPACK_A,
                                        Ap, Ai, Ax, x, res, numeric,
                                        control_array, info_array);

    print_status(status_solve);

    umfpack_di_free_numeric(&numeric);

    delete[] data;
    delete[] row;
    delete[] col;
    delete[] Ap;
    delete[] Ai;
    delete[] Ax;

    if (symbolic) umfpack_di_free_symbolic(&symbolic);
    if (numeric) umfpack_di_free_numeric(&numeric);

    memcpy(res, x, mcoo->get_size()*sizeof(double));
}

void solve_linear_system_umfpack(Matrix *mat, cplx *res)
{
    _error("hermes_common: solve_linear_system_umfpack - not implemented.");
}

#else

void solve_linear_system_umfpack(Matrix *mat, double *res)
{
    _error("hermes_common: solve_linear_system_umfpack - UMFPACK is not available.");
}
void solve_linear_system_umfpack(Matrix *mat, cplx *res)
{
    _error("hermes_common: solve_linear_system_umfpack - UMFPACK is not available.");
}
#endif
