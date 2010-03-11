// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "matrix.h"

#include "_hermes_common_api.h"

#define TINY 1e-20

/// Transposes an m by n matrix. If m != n, the array matrix in fact has to be
/// a square matrix of the size max(m, n) in order for the transpose to fit inside it.
template<typename T>
void transpose(T** matrix, int m, int n)
{
  int min = std::min(m, n);
  for (int i = 0; i < min; i++)
    for (int j = i+1; j < min; j++)
      std::swap(matrix[i][j], matrix[j][i]);

  if (m < n)
    for (int i = 0; i < m; i++)
      for (int j = m; j < n; j++)
        matrix[j][i] = matrix[i][j];
  else if (n < m)
    for (int i = n; i < m; i++)
      for (int j = 0; j < n; j++)
        matrix[j][i] = matrix[i][j];
}


/// Changes the sign of a matrix
template<typename T>
void chsgn(T** matrix, int m, int n)
{
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      matrix[i][j] = -matrix[i][j];
}


/// Given a matrix a[n][n], this routine replaces it by the LU decomposition of a rowwise
/// permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;
/// indx[n] is an output vector that records the row permutation effected by the partial
/// pivoting; d is output as +-1 depending on whether the number of row interchanges was even
/// or odd, respectively. This routine is used in combination with lubksb to solve linear equations
/// or invert a matrix.
void ludcmp(double** a, int n, int* indx, double* d);

/// Solves the set of n linear equations AX = B. Here a[n][n] is input, not as the matrix
/// A but rather as its LU decomposition, determined by the routine ludcmp. indx[n] is input
/// as the permutation vector returned by ludcmp. b[n] is input as the right-hand side vector
/// B, and returns with the solution vector X. a, n, and indx are not modified by this routine
/// and can be left in place for successive calls with different right-hand sides b. This routine takes
/// into account the possibility that b will begin with many zero elements, so it is efficient for use
/// in matrix inversion.
template<typename T>
void lubksb(double** a, int n, int* indx, T* b)
{
  int i, ip, j;
  T sum;

  for (i = 0; i < n; i++)
  {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    for (j = 0; j < i; j++) sum -= a[i][j]*b[j];
    b[i] = sum;
  }
  for (i = n-1; i >= 0; i--)
  {
    sum = b[i];
    for (j = i+1; j < n; j++) sum -= a[i][j]*b[j];
    b[i] = sum / a[i][i];
  }
}

/// Given a positive-definite symmetric matrix a[n][n], this routine constructs its Cholesky
/// decomposition, A = L*L^T . On input, only the upper triangle of a need be given; it is not
/// modified. The Cholesky factor L is returned in the lower triangle of a, except for its diagonal
/// elements which are returned in p[n].
void choldc(double **a, int n, double p[]);

/// Solves the set of n linear equations A*x = b, where a is a positive-definite symmetric matrix.
/// a[n][n] and p[n] are input as the output of the routine choldc. Only the lower
/// subdiagonal portion of a is accessed. b[n] is input as the right-hand side vector. The
/// solution vector is returned in x[n]. a, n, and p are not modified and can be left in place
/// for successive calls with different right-hand sides b. b is not modified unless you identify b and
/// x in the calling sequence, which is allowed. The right-hand side b can be complex, in which case
/// the solution x is also complex.
template<typename T>
void cholsl(double **a, int n, double p[], T b[], T x[])
{
  int i, k;
  T sum;

  for (i = 0; i < n; i++)
  {
    sum = b[i];
    k = i;
    while (--k >= 0)
      sum -= a[i][k] * x[k];
    x[i] = sum / p[i];
  }

  for (i = n-1; i >= 0; i--)
  {
    sum = x[i];
    k = i;
    while (++k < n)
      sum -= a[k][i] * x[k];
    x[i] = sum / p[i];
  }
}


void ludcmp(double** a, int n, int* indx, double* d)
{
  int i, imax = 0, j, k;
  double big, dum, sum, temp;
  double* vv = new double[n];

  *d = 1.0;
  for (i = 0; i < n; i++)
  {
    big=0.0;
    for (j = 0; j < n; j++) 
      if ((temp = fabs(a[i][j])) > big) 
        big = temp;
    if (big == 0.0) _error("Singular matrix!");
    vv[i] = 1.0 / big;
  }
  for (j = 0; j < n; j++)
  {
    for (i = 0; i < j; i++)
    {
      sum = a[i][j];
      for (k = 0; k < i; k++) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
    }
    big = 0.0;
    for (i = j; i < n; i++)
    {
      sum = a[i][j];
      for (k = 0; k < j; k++) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
      if ((dum = vv[i]*fabs(sum)) >= big)
      {
        big = dum;
        imax = i;
      }
    }
    if (j != imax)
    {
      for (k = 0; k < n; k++)
      {
        dum = a[imax][k];
        a[imax][k] = a[j][k];
        a[j][k] = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if (a[j][j] == 0.0) a[j][j] = TINY;
    if (j != n-1) 
    {
      dum = 1.0 / (a[j][j]);
      for (i = j+1; i < n; i++) a[i][j] *= dum;
    }
  }
  delete [] vv;
}


void choldc(double **a, int n, double p[])
{
  int i, j, k;
  for (i = 0; i < n; i++)
  {
    for (j = i; j < n; j++)
    {
      double sum = a[i][j];
      k = i;
      while (--k >= 0)
        sum -= a[i][k] * a[j][k];
      
      if (i == j)
      {
        if (sum <= 0.0)
          _error("CHOLDC failed!");
        else
          p[i] = sqrt(sum);
      }
      else
        a[j][i] = sum / p[i];
    }
  }
}

int initialized=0;

Matrix::Matrix()
{
    if (!initialized) {
        // this is necessary, so that we can use Python from matrix.cpp:
        if (import__hermes_common())
            throw std::runtime_error("hermes_common failed to import.");
        initialized=1;
    }
}

CSRMatrix::CSRMatrix(CooMatrix *m):Matrix()
{
    insert_object("m", c2py_CooMatrix(m));
    cmd("n = m.to_scipy_coo().tocsr()");
    cmd("A = n.data");
    cmd("IA = n.indptr");
    cmd("JA = n.indices");
    //XXX: this should *not* be inplace, we need to fix it.
    numpy2c_double_inplace(get_object("A"), &(this->A), &(this->nnz));
    numpy2c_int_inplace(get_object("IA"), &(this->IA), &(this->size));
    numpy2c_int_inplace(get_object("JA"), &(this->JA), &(this->nnz));
    this->size--;

    // Original C++ implementation using a DenseMatrix:
    //DenseMatrix *dmat = new DenseMatrix(m);
    //this->add_from_dense_matrix(dmat);
    //delete dmat;
}

void CSRMatrix::print()
{
    insert_object("m", c2py_CSRMatrix(this));
    cmd("S = str(m.to_scipy_csr())");
    printf("%s\n", py2c_str(get_object("S")));
}

CSCMatrix::CSCMatrix(CooMatrix *m):Matrix()
{
    insert_object("m", c2py_CooMatrix(m));
    cmd("n = m.to_scipy_coo().tocsc()");
    cmd("A = n.data");
    cmd("IA = n.indices");
    cmd("JA = n.indptr");
    //XXX: this should *not* be inplace, we need to fix it.
    numpy2c_double_inplace(get_object("A"), &(this->A), &(this->nnz));
    numpy2c_int_inplace(get_object("IA"), &(this->IA), &(this->nnz));
    numpy2c_int_inplace(get_object("JA"), &(this->JA), &(this->size));
    this->size--;
}

void CSCMatrix::print()
{
    insert_object("m", c2py_CSCMatrix(this));
    cmd("S = str(m.to_scipy_csc())");
    printf("%s\n", py2c_str(get_object("S")));
}

void CooMatrix::print()
{
    insert_object("m", c2py_CooMatrix(this));
    cmd("S = str(m.to_scipy_coo())");
    printf("%s\n", py2c_str(get_object("S")));
}

void solve_linear_system_dense(DenseMatrix *mat, double *res)
{
    int n = mat->get_size();
    int *indx = new int[n];
    double **_mat = mat->get_mat();
    double d;
    ludcmp(_mat, n, indx, &d);
    lubksb(_mat, n, indx, res);
}


void solve_linear_system(Matrix *mat, double *res)
{
    DenseMatrix *dmat = new DenseMatrix(mat);
    solve_linear_system_dense(dmat, res);
    delete dmat;
}
