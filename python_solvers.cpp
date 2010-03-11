// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "matrix.h"

#include "_hermes_common_api.h"

void solve_linear_system_numpy(CooMatrix *mat, double *res)
{
    if (import__hermes_common())
        throw std::runtime_error("hermes_common failed to import.");
    CSRMatrix M(mat);
    insert_object("m", c2py_CSRMatrix(&M));
    insert_object("rhs", c2numpy_double_inplace(res, mat->get_size()));
    cmd("A = m.to_scipy_csr().todense()");
    cmd("from numpy.linalg import solve");
    cmd("x = solve(A, rhs)");
    double *x;
    int n;
    numpy2c_double_inplace(get_object("x"), &x, &n);
    memcpy(res, x, n*sizeof(double));
}

void solve_linear_system_scipy_umfpack(CooMatrix *mat, double *res)
{
    if (import__hermes_common())
        throw std::runtime_error("hermes_common failed to import.");
    CSCMatrix M(mat);
    insert_object("m", c2py_CSCMatrix(&M));
    insert_object("rhs", c2numpy_double_inplace(res, mat->get_size()));
    cmd("A = m.to_scipy_csc()");
    cmd("from scipy.sparse.linalg import spsolve");
    cmd("x = spsolve(A, rhs)");
    double *x;
    int n;
    numpy2c_double_inplace(get_object("x"), &x, &n);
    memcpy(res, x, n*sizeof(double));
}

void solve_linear_system_scipy_cg(CooMatrix *mat, double *res)
{
    if (import__hermes_common())
        throw std::runtime_error("hermes_common failed to import.");
    CSRMatrix M(mat);
    insert_object("m", c2py_CSRMatrix(&M));
    insert_object("rhs", c2numpy_double_inplace(res, mat->get_size()));
    cmd("A = m.to_scipy_csr()");
    cmd("from scipy.sparse.linalg import cg");
    cmd("x, res = cg(A, rhs)");
    double *x;
    int n;
    numpy2c_double_inplace(get_object("x"), &x, &n);
    memcpy(res, x, n*sizeof(double));
}

void solve_linear_system_scipy_gmres(CooMatrix *mat, double *res)
{
    if (import__hermes_common())
        throw std::runtime_error("hermes_common failed to import.");
    CSRMatrix M(mat);
    insert_object("m", c2py_CSRMatrix(&M));
    insert_object("rhs", c2numpy_double_inplace(res, mat->get_size()));
    cmd("A = m.to_scipy_csr()");
    cmd("from scipy.sparse.linalg import gmres");
    cmd("x, res = gmres(A, rhs)");
    double *x;
    int n;
    numpy2c_double_inplace(get_object("x"), &x, &n);
    memcpy(res, x, n*sizeof(double));
}
