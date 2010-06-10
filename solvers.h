// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef __HERMES_COMMON_SOLVERS_H
#define __HERMES_COMMON_SOLVERS_H

// default solvers
void solve_linear_system_dense_lu(Matrix *mat, double *res);
int solve_linear_system_cg(Matrix* A, double *x,
                           double matrix_solver_tol = 1e-6,
                           int matrix_solver_maxiter = 30);

// python numpy - optional
void solve_linear_system_numpy(Matrix *mat, double *res);
void solve_linear_system_numpy(Matrix *mat, cplx *res);

// python scipy - optional
void solve_linear_system_scipy_umfpack(Matrix *mat, double *res);
void solve_linear_system_scipy_umfpack(Matrix *mat, cplx *res);
void solve_linear_system_scipy_cg(Matrix *mat, double *res);
void solve_linear_system_scipy_gmres(Matrix *mat, double *res);

// c++ umfpack - optional
void solve_linear_system_umfpack(Matrix *mat, double *res);
void solve_linear_system_umfpack(Matrix *mat, cplx *res);

// c++ superlu - optional
void solve_linear_system_superlu(Matrix *mat, double *res);

// c++ sparselib - optional
void solve_linear_system_sparselib(Matrix *mat, double *res);

#endif
