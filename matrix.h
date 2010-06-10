// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef __HERMES_COMMON_MATRIX_H
#define __HERMES_COMMON_MATRIX_H

#include <math.h>
#include <stdlib.h>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <typeinfo>
#include <math.h>
#include <string.h>
#include <complex>

typedef std::complex<double> cplx;
class Matrix;

#include "solvers.h"

// printf debug information about the stiffness/Jacobi matrix
#define _error(x) throw std::runtime_error(x)
#define DEBUG_MATRIX 0


/// Creates a new (full) matrix with m rows and n columns with entries of the type T.
/// The entries can be accessed by matrix[i][j]. To delete the matrix, just
/// do "delete matrix".
template<typename T>
T** _new_matrix(int m, int n = 0)
{
    if (!n) n = m;
    T** vec = (T**) new char[sizeof(T*)*m + sizeof(T)*m*n];
    if (vec == NULL) _error("Out of memory.");
    T* row = (T*) (vec + m);
    for (int i = 0; i < m; i++, row += n)
        vec[i] = row;
    return vec;
}

template <typename SCALAR>
        class Triple {
        public:
    int i;
    int j;
    SCALAR v;
    Triple<SCALAR> *next;

    Triple(int i, int j, SCALAR v) {
        this->i = i;
        this->j = j;
        this->v = v;
        this->next = NULL;
    }
};

class Matrix {
public:
    Matrix();
    virtual ~Matrix();
    virtual int get_size() = 0;
    virtual void add(int m, int n, double v) = 0;
    virtual void add(int m, int n, cplx v) {
	    _error("internal error: add(int, int, cplx) not implemented.");
    }
    virtual void add_block(int *iidx, int ilen, int *jidx, int jlen,
                           double** mat) {
        for (int i = 0; i < ilen; i++)
            for (int j=0; j < jlen; j++)
                if (iidx[i] >= 0 && jidx[j] >= 0)
                    this->add(iidx[i], jidx[j], mat[i][j]);
    }
    virtual void add_block(int *iidx, int ilen, int *jidx, int jlen,
                           cplx** mat) {
        for (int i = 0; i < ilen; i++)
            for (int j=0; j < jlen; j++)
                if (iidx[i] >= 0 && jidx[j] >= 0)
                    this->add(iidx[i], jidx[j], mat[i][j]);
    }
    virtual void set_zero() = 0;
    virtual double get(int m, int n) = 0;
    virtual cplx get_cplx(int m, int n) {
	    _error("internal error: get_cplx() not implemented.");
    }
    virtual void copy_into(Matrix *m) = 0;
    virtual void print() = 0;
    virtual bool is_complex() {
        return this->_is_complex;
    }
    virtual void times_vector(double* vec, double* result, int rank) {
	    _error("internal error: times_vector() not implemented.");
    }
protected:
    bool _is_complex;
};


class CSRMatrix;
class CSCMatrix;

class CooMatrix : public Matrix {
public:
    CooMatrix(bool is_complex=false):Matrix() {
        init();

        this->_is_complex = is_complex;
    }
    CooMatrix(int size, bool is_complex=false):Matrix() {
        init();

        this->_is_complex = is_complex;
        this->size = size;
    }
    ~CooMatrix() {
        this->set_zero();
    }

    virtual void init() {
        this->_is_complex = false;
        this->size = 0;
        this->list = NULL;
        this->list_last = NULL;
        this->list_cplx = NULL;
        this->list_last_cplx = NULL;
    }

    template <typename SCALAR>
    void free_data() {
        Triple<SCALAR> *t = this->get_list<SCALAR>();
        while (t != NULL) {
            Triple<SCALAR> *t_old = t;
            t = t->next;
            delete t_old;
        }
    }

    virtual void set_zero();

    virtual void add(int m, int n, double v) {
        if (this->_is_complex)
            _error("can't use add(int, int, double) for complex matrix");

        // adjusting size if necessary
        if (m+1 > this->size) this->size = m+1;
        if (n+1 > this->size) this->size = n+1;

        // debug
        if (DEBUG_MATRIX) {
            printf("Matrix_add %d %d %g -> size = %d\n", m, n, v, this->size);
        }

        /*
        // add into existing coo
        Triple<double> *l = this->list;
        while (l) {
            if (m == l->i && n == l->j)
            {
                l->v += v;
                return;
            }
            l = l->next;
        }
        */

        // add new triple
        Triple<double> *t = new Triple<double>(m, n, v);
        if (this->list == NULL) {
            this->list = t;
            this->list_last = t;
        } else {
            this->list_last->next = t;
            this->list_last = this->list_last->next;
        }
    }

    virtual void add(int m, int n, cplx v) {
        if (!(this->_is_complex))
            _error("can't use add(int, int, cplx) for real matrix");

        // adjusting size if necessary
        if (m+1 > this->size) this->size = m+1;
        if (n+1 > this->size) this->size = n+1;

        // add into existing coo
        Triple<cplx> *l = this->list_cplx;
        while (l) {
            if (m == l->i && n == l->j)
            {
                l->v += v;
                return;
            }
            l = l->next;
        }

        Triple<cplx> *t = new Triple<cplx>(m, n, v);
        if (this->list_cplx == NULL) {
            this->list_cplx = t;
            this->list_last_cplx = t;
        } else {
            this->list_last_cplx->next = t;
            this->list_last_cplx = this->list_last_cplx->next;
        }
    }

    virtual void copy_into(Matrix *m) {
        m->set_zero();
        if (this->_is_complex) {
            Triple<cplx> *t = this->list_cplx;
            while (t != NULL) {
                m->add(t->i, t->j, t->v);
                t = t->next;
            }
        } else {
            Triple<double> *t = this->list;
            while (t != NULL) {
                m->add(t->i, t->j, t->v);
                t = t->next;
            }
        }
    }

    // Returns the number of triplets
    int triplets_len() {
        Triple<double> *t = this->list;
        int len = 0;
        while (t != NULL) {
            len++;
            t = t->next;
        }
        return len;
    }

    int triplets_len_cplx() {
        Triple<cplx> *t = this->list_cplx;
        int len = 0;
        while (t != NULL) {
            len++;
            t = t->next;
        }
        return len;
    }

    // Returns the row/col indices, together with the data. All row, col
    // and data arrays have to be initialized (use this->triplets_len).
    void get_row_col_data(int *row, int *col, double *data) {
        Triple<double> *t = this->list;
        int i = 0;
        while (t) {
            row[i] = t->i;
            col[i] = t->j;
            data[i] = t->v;
            i++;
            t = t->next;
        }
    }

    void get_row_col_data(int *row, int *col, cplx *data) {
        Triple<cplx> *t = this->list_cplx;
        int i = 0;
        while (t) {
            row[i] = t->i;
            col[i] = t->j;
            data[i] = t->v;
            i++;
            t = t->next;
        }
    }

    virtual double get(int m, int n) {
        Triple<double> *t = this->list;
        while (t) {
            if (m == t->i && n == t->j)
                return t->v;
            t = t->next;
        }
        return 0.0;
    }

    virtual int get_size() {
        return this->size;
    }

    virtual void print();

    virtual void times_vector(double* vec, double* result, int rank) {
        for (int i=0; i < rank; i++) result[i] = 0;
        Triple<double> *t = this->list;
        while (t != NULL) {
            result[t->i] += t->v * vec[t->j];
            t = t->next;
        }
    }

    template <typename SCALAR>
            Triple<SCALAR> *get_list() {
        // this works for SCALAR=double, and in matrix.cpp, we specialize
        // it for SCALAR=cplx
        return this->list;
    }


private:
    int size;
    /*
           We represent the COO matrix as a list of Triples (i, j, v), where
           (i, j) can be redundant (then the corresponding "v" have to be
           summed). We remember the pointers to the first (*list) and last
           (*list_last) element.
           */
    Triple<double> *list;
    Triple<double> *list_last;

    Triple<cplx> *list_cplx;
    Triple<cplx> *list_last_cplx;
};

class DenseMatrix : public Matrix {
public:
    DenseMatrix(int size, bool is_complex=false) {
        this->_is_complex = is_complex;

        if (is_complex)
            this->mat_cplx = _new_matrix<cplx>(size, size);
        else
            this->mat = _new_matrix<double>(size, size);
        this->size = size;

        if (is_complex) {
            for (int i = 0; i<size; i++)
                for (int j = 0; j<size; j++)
                    this->mat_cplx[i][j] = 0;
        } else {
            for (int i = 0; i<size; i++)
                for (int j = 0; j<size; j++)
                    this->mat[i][j] = 0;
        }
    }
    DenseMatrix(Matrix *m, bool is_complex=false) {
        this->_is_complex = is_complex;
        if (is_complex)
            this->mat_cplx =
                    _new_matrix<cplx>(m->get_size(), m->get_size());
        else
            this->mat = _new_matrix<double>(m->get_size(), m->get_size());
        this->size = m->get_size();
        m->copy_into(this);
    }
    ~DenseMatrix() {
        free_data();
    }

    virtual void free_data() {
        if (this->mat) { delete[] this->mat; this->mat = NULL; };
        if (this->mat_cplx) { delete[] this->mat_cplx; this->mat_cplx = NULL; };
    }

    virtual void add(int m, int n, double v) {
        this->mat[m][n] += v;
    }
    virtual void add(int m, int n, cplx v) {
        this->mat_cplx[m][n] += v;
    }
    virtual void set_zero() {
        if (this->_is_complex) {
            for (int i = 0; i<size; i++)
                for (int j = 0; j<size; j++)
                    this->mat_cplx[i][j] = 0;
        } else {
            for (int i = 0; i<size; i++)
                for (int j = 0; j<size; j++)
                    this->mat[i][j] = 0;
        }
    }
    virtual double get(int m, int n) {
        return this->mat[m][n];
    }

    virtual int get_size() {
        return this->size;
    }
    virtual void copy_into(Matrix *m) {
	    m->set_zero();
        for (int i = 0; i < this->size; i++)
            for (int j = 0; j < this->size; j++) {
            double v = this->get(i, j);
            if (fabs(v) > 1e-12)
                m->add(i, j, v);
        }
    }

    virtual void print() {
        for (int i = 0; i < this->size; i++) {
            for (int j = 0; j < this->size; j++) {
                double v = this->get(i, j);
                printf("%f ", v);
                /*
                    if (fabs(v) > 1e-12)
                        printf("%f ", v);
                    else
                        printf("0 ");
                        */
            }
            printf("\n");
        }
    }

    // Return the internal matrix.
    double **get_mat() {
        return this->mat;
    }


private:
    int size;

    double **mat;
    cplx **mat_cplx;
};

class CSRMatrix : public Matrix {
public:
    CSRMatrix(int size) : Matrix() {
        init();

        this->size = size;
    }

    CSRMatrix(Matrix *m);
    CSRMatrix(CooMatrix *m):Matrix() {
        init();

        this->add_from_CooMatrix(m);
    }
    CSRMatrix(CSCMatrix *m):Matrix() {
        init();

        this->add_from_CSCMatrix(m);
    }
    CSRMatrix(DenseMatrix *m):Matrix() {
        init();

        this->add_from_DenseMatrix(m);
    }

    ~CSRMatrix() {
        free_data();
    }

    virtual void init() {
        this->_is_complex = false;
        this->size = 0;
        this->nnz = 0;
        this->A = NULL;
        this->A_cplx = NULL;
        this->Ap = NULL;
        this->Ai = NULL;
    }

    virtual void free_data() {
        if (this->Ap) { delete[] this->Ap; this->Ap = NULL; }
        if (this->Ai) { delete[] this->Ai; this->Ai = NULL; }
        if (this->A) { delete[] this->A; this->A = NULL; }
        if (this->A_cplx) { delete[] this->A_cplx; this->A_cplx = NULL; }
    }

    void add_from_DenseMatrix(DenseMatrix *m);
    void add_from_CooMatrix(CooMatrix *m);
    void add_from_CSCMatrix(CSCMatrix *m);

    virtual void add(int m, int n, double v) {
        _error("CSR matrix add() not implemented.");
    }
    virtual void set_zero() {
        _error("CSR matrix set_zero() not implemented.");
    }
    virtual double get(int m, int n) {
        _error("CSR matrix get() not implemented.");
    }

    virtual int get_size() {
        return this->size;
    }
    int get_nnz() {
        return this->nnz;
    }
    virtual void copy_into(Matrix *m) {
        _error("CSR matrix copy_into() not implemented.");
    }

    virtual void print();

    int *get_Ap() {
        return this->Ap;
    }
    int *get_Ai() {
        return this->Ai;
    }
    double *get_A() {
        return this->A;
    }
    cplx *get_A_cplx() {
        return this->A_cplx;
    }

private:
    // matrix size
    int size;
    // number of non-zeros
    int nnz;

    double *A;
    cplx *A_cplx;

    int *Ap;
    int *Ai;
};

class CSCMatrix : public Matrix {
public:
    CSCMatrix(int size):Matrix() {
        init();

        this->size = size;
    }
    CSCMatrix(Matrix *m):Matrix() {
        init();

        if (dynamic_cast<CooMatrix*>(m))
            this->add_from_CooMatrix((CooMatrix*)m);
        else if (dynamic_cast<CSRMatrix*>(m))
            this->add_from_CSRMatrix((CSRMatrix*)m);
        else
            _error("Matrix type not supported.");
    }
    CSCMatrix(CooMatrix *m):Matrix() {
        init();

        this->add_from_CooMatrix(m);
    }
    CSCMatrix(CSRMatrix *m):Matrix() {
        init();

        this->add_from_CSRMatrix(m);
    }
    ~CSCMatrix() {
        free_data();
    }

    virtual void init() {
        this->_is_complex = false;
        this->size = 0;
        this->nnz = 0;
        this->A = NULL;
        this->A_cplx = NULL;
        this->Ap = NULL;
        this->Ai = NULL;
    }

    virtual void free_data() {
        if (this->Ap) { delete[] this->Ap; this->Ap = NULL; }
        if (this->Ai) { delete[] this->Ai; this->Ai = NULL; }
        if (this->A) { delete[] this->A; this->A = NULL; }
        if (this->A_cplx) { delete[] this->A_cplx; this->A_cplx = NULL; }
    }

    void add_from_CooMatrix(CooMatrix *m);
    void add_from_CSRMatrix(CSRMatrix *m);

    virtual void add(int m, int n, double v) {
        _error("CSC matrix add() not implemented.");
    }
    virtual void set_zero() {
        _error("CSC matrix set_zero() not implemented.");
    }
    virtual double get(int m, int n) {
        _error("CSC matrix get() not implemented.");
    }

    virtual int get_size() {
        return this->size;
    }
    int get_nnz() {
        return this->nnz;
    }
    virtual void copy_into(Matrix *m) {
        _error("CSC matrix copy_into() not implemented.");
    }

    virtual void print();

    int *get_Ap() {
        return this->Ap;
    }
    int *get_Ai() {
        return this->Ai;
    }
    double *get_A() {
        return this->A;
    }
    cplx *get_A_cplx() {
        return this->A_cplx;
    }

private:
    // matrix size
    int size;
    // number of non-zeros
    int nnz;

    double *A;
    cplx *A_cplx;

    int *Ap;
    int *Ai;
};

// print vector - int
void print_vector(const char *label, int *value, int size);
// print vector - double
void print_vector(const char *label, double *value, int size);
// print vector - cplx
void print_vector(const char *label, cplx *value, int size);

#endif
