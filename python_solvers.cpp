// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "matrix.h"

#include "python_api.h"

void solve_linear_system_numpy(CooMatrix *mat, double *res)
{
    CSRMatrix M(mat);
    Python *p = new Python();
    p->push("m", c2py_CSRMatrix(&M));
    p->push("rhs", c2numpy_double_inplace(res, mat->get_size()));
    p->exec("A = m.to_scipy_csr().todense()");
    p->exec("from numpy.linalg import solve");
    p->exec("x = solve(A, rhs)");
    double *x;
    int n;
    numpy2c_double_inplace(p->pull("x"), &x, &n);
    memcpy(res, x, n*sizeof(double));
    delete p;
}

void solve_linear_system_scipy_umfpack(CooMatrix *mat, double *res)
{
    CSCMatrix M(mat);
    Python *p = new Python();
    p->push("m", c2py_CSCMatrix(&M));
    p->push("rhs", c2numpy_double_inplace(res, mat->get_size()));
    p->exec("A = m.to_scipy_csc()");
    p->exec("from scipy.sparse.linalg import spsolve");
    p->exec("x = spsolve(A, rhs)");
    double *x;
    int n;
    numpy2c_double_inplace(p->pull("x"), &x, &n);
    memcpy(res, x, n*sizeof(double));
    delete p;
}

void solve_linear_system_scipy_cg(CooMatrix *mat, double *res)
{
    CSRMatrix M(mat);
    Python *p = new Python();
    p->push("m", c2py_CSRMatrix(&M));
    p->push("rhs", c2numpy_double_inplace(res, mat->get_size()));
    p->exec("A = m.to_scipy_csr()");
    p->exec("from scipy.sparse.linalg import cg");
    p->exec("x, res = cg(A, rhs)");
    double *x;
    int n;
    numpy2c_double_inplace(p->pull("x"), &x, &n);
    memcpy(res, x, n*sizeof(double));
    delete p;
}

void solve_linear_system_scipy_gmres(CooMatrix *mat, double *res)
{
    CSRMatrix M(mat);
    Python *p = new Python();
    p->push("m", c2py_CSRMatrix(&M));
    p->push("rhs", c2numpy_double_inplace(res, mat->get_size()));
    p->exec("A = m.to_scipy_csr()");
    p->exec("from scipy.sparse.linalg import gmres");
    p->exec("x, res = gmres(A, rhs)");
    double *x;
    int n;
    numpy2c_double_inplace(p->pull("x"), &x, &n);
    memcpy(res, x, n*sizeof(double));
    delete p;
}
