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

// printf debug information about the stiffness/Jacobi matrix
#define _error(x) throw std::runtime_error(x)
#define DEBUG_MATRIX 0

class Matrix {
public:
    Matrix();
    virtual ~Matrix() { }
    virtual int get_size() = 0;
    virtual void add(int m, int n, double v) = 0;
    virtual void add_block(int *iidx, int ilen, int *jidx, int jlen,
            double** mat) {
        for (int i = 0; i < ilen; i++)
            for (int j=0; j < jlen; j++)
                if (iidx[i] >= 0 && jidx[j] >= 0)
                    this->add(iidx[i], jidx[j], mat[i][j]);
    }
    virtual void set_zero() = 0;
    virtual double get(int m, int n) = 0;
    virtual void copy_into(Matrix *m) = 0;
    virtual void print() = 0;
    virtual void times_vector(double* vec, double* result, int rank) = 0;
};

class Triple {
    public:
        int i;
        int j;
        double v;
        Triple *next;

        Triple(int i, int j, double v) {
            this->i = i;
            this->j = j;
            this->v = v;
            this->next = NULL;
        }
};

class CooMatrix : public Matrix {
    public:
        CooMatrix():Matrix() {
            this->size = 0;
            this->list = NULL;
            this->list_last = NULL;
        }
        CooMatrix(int size):Matrix() {
            this->size = size;
            this->list = NULL;
            this->list_last = NULL;
        }
        ~CooMatrix() {
            this->free_data();
            this->list = NULL;
            this->list_last = NULL;
            this->size = 0;
        }
        void free_data() {
            Triple *t = this->list;
            while (t != NULL) {
                Triple *t_old = t;
                t = t->next;
                delete t_old;
            }
        }
        virtual void set_zero() {
            this->free_data();
            this->list = NULL;
            this->list_last = NULL;
        }
        virtual void add(int m, int n, double v) {
            // adjusting size if necessary
            if (m+1 > this->size) this->size = m+1;
            if (n+1 > this->size) this->size = n+1;
            // debug
            if (DEBUG_MATRIX) {
                printf("Matrix_add %d %d %g -> size = %d\n", m, n, v, this->size);
            }
            Triple *t = new Triple(m, n, v);
            if (this->list == NULL) {
                this->list = t;
                this->list_last = t;
            } else {
                this->list_last->next = t;
                this->list_last = this->list_last->next;
            }
        }
        virtual void copy_into(Matrix *m) {
            m->set_zero();
            Triple *t = this->list;
            while (t != NULL) {
                m->add(t->i, t->j, t->v);
                t = t->next;
            }
        }

        // Returns the number of triplets
        int triplets_len() {
            Triple *t = this->list;
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
            Triple *t = this->list;
            int i = 0;
            while (t != NULL) {
                row[i] = t->i;
                col[i] = t->j;
                data[i] = t->v;
                i++;
                t = t->next;
            }
        }

        virtual double get(int m, int n) {
            double v=0;
            Triple *t = this->list;
            if (t != NULL) {
                while (t->next != NULL) {
                    if (m == t->i && n == t->j)
                        v += t->v;
                    t = t->next;
                }
            }
            return v;
        }

        virtual int get_size() {
            return this->size;
        }

        virtual void print();

        virtual void times_vector(double* vec, double* result, int rank) {
            for (int i=0; i < rank; i++) result[i] = 0;
            Triple *t = this->list;
            while (t != NULL) {
                result[t->i] += t->v * vec[t->j];
                t = t->next;
            }
        }

    private:
        int size;
        /*
           We represent the COO matrix as a list of Triples (i, j, v), where
           (i, j) can be redundant (then the corresponding "v" have to be
           summed). We remember the pointers to the first (*list) and last
           (*list_last) element.
           */
        Triple *list;
        Triple *list_last;
};


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

class DenseMatrix : public Matrix {
    public:
        DenseMatrix(int size) {
            this->mat = _new_matrix<double>(size, size);
            this->size = size;
            for (int i = 0; i<size; i++)
              for (int j = 0; j<size; j++) this->mat[i][j] = 0;
        }
        DenseMatrix(Matrix *m) {
            this->mat = _new_matrix<double>(m->get_size(), m->get_size());
            //printf("%d %d", this->size, m->get_size());
            //exit(1);
            this->size = m->get_size();
            //this->size = size;
            m->copy_into(this);

        }
        ~DenseMatrix() {
            delete[] this->mat;
        }
        virtual void add(int m, int n, double v) {
            this->mat[m][n] += v;
            //printf("calling add: %d %d %f\n", m, n, v);
        }
        virtual void set_zero() {
            for (int i = 0; i<size; i++)
              for (int j = 0; j<size; j++) this->mat[i][j] = 0;
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

        virtual void times_vector(double* vec, double* result, int rank) {
	    _error("times_vector() in dense matrix not implemented yet.");
        }

        // Return the internal matrix.
        double **get_mat() {
            return this->mat;
        }
        double **mat;

    private:
        int size;

};

class CSRMatrix : public Matrix {
    public:
        CSRMatrix(int size):Matrix() {
            this->size = size;
            this->A = NULL;
            this->IA = NULL;
            this->JA = NULL;
        }
        CSRMatrix(CooMatrix *m);
        CSRMatrix(DenseMatrix *m):Matrix() {
            this->add_from_dense_matrix(m);
        }
        ~CSRMatrix() {
            delete[] this->A;
            delete[] this->IA;
            delete[] this->JA;
        }

        void add_from_dense_matrix(DenseMatrix *m) {
            this->size = m->get_size();
            this->nnz = 0;
            for(int i = 0; i < this->size; i++)
                for(int j = 0; j < this->size; j++) {
                    double v = m->get(i, j);
                    if (fabs(v) > 1e-12)
                        this->nnz++;
                }
            this->A = new double[this->nnz];
            this->IA = new int[this->size+1];
            this->JA = new int[this->nnz];
            int count = 0;
            this->IA[0] = 0;
            for(int i = 0; i < this->size; i++) {
                for(int j = 0; j < this->size; j++) {
                    double v = m->get(i, j);
                    if (fabs(v) > 1e-12) {
                        this->A[count] = v;
                        this->JA[count] = j;
                        count++;
                    }
                }
                this->IA[i+1] = count;
            }
        }

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

        virtual void times_vector(double* vec, double* result, int rank) {
	    _error("times_vector() in CSR matrix not implemented yet.");
        }

        int *get_IA() {
            return this->IA;
        }
        int *get_JA() {
            return this->JA;
        }
        double *get_A() {
            return this->A;
        }

    private:
        int size;
        int nnz;
        double *A;
        int *IA;
        int *JA;

};

// solve linear system
void solve_linear_system(Matrix *mat, double *res);
void solve_linear_system_dense(DenseMatrix *mat, double *res);

#endif
