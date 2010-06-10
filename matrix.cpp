// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "matrix.h"

// print vector - int
void print_vector(const char *label, int *value, int size) {
    printf("%s [", label);
    for (int i = 0; i < size; i++) {
        if (i < size-1)
            printf("%i, ", value[i]);
        else
            printf("%i", value[i]);
    }
    printf("]\n");
}

// print vector - double
void print_vector(const char *label, double *value, int size) {
    printf("%s [", label);
    for (int i = 0; i < size; i++) {
        if (i < size-1)
            printf("%f, ", value[i]);
        else
            printf("%f", value[i]);
    }
    printf("]\n");
}

// print vector - cplx
void print_vector(const char *label, cplx *value, int size) {
    printf("%s [", label);
    for (int i = 0; i < size; i++) {
        if (i < size-1)
            printf("(%f, %f), ", value[i].real(), value[i].imag());
        else
            printf("(%f, %f)", value[i].real(), value[i].imag());
    }
    printf("]\n");
}

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




Matrix::Matrix()
{
}

Matrix::~Matrix()
{
}

// *********************************************************************************************************************

void CooMatrix::print()
{
    printf("\nCoo Matrix:\n");
    if (is_complex()) {
        Triple<cplx> *t = this->list_cplx;
        while (t != NULL) {
            printf("(%i, %i): (%f, %f) \n", t->i, t->j, t->v.real(), t->v.imag());
            t = t->next;
        }
    } else {
        Triple<double> *t = this->list;
        while (t != NULL) {
            printf("(%i, %i): %f\n", t->i, t->j, t->v);
            t = t->next;
        }
    }
}

// *********************************************************************************************************************

CSRMatrix::CSRMatrix(Matrix *m):Matrix() {
    init();

    if (dynamic_cast<CooMatrix*>(m))
        this->add_from_CooMatrix((CooMatrix*)m);
    else if (dynamic_cast<CSCMatrix*>(m))
        this->add_from_CSCMatrix((CSCMatrix*)m);
    else if (dynamic_cast<DenseMatrix*>(m))
        this->add_from_DenseMatrix((DenseMatrix*)m);
    else
        _error("Matrix type not supported.");
}

void CSRMatrix::add_from_DenseMatrix(DenseMatrix *m) {
    this->size = m->get_size();

    // count nnz
    this->nnz = 0;
    for(int i = 0; i < this->size; i++) {
        for(int j = 0; j < this->size; j++) {
            double v = m->get(i, j);
            if (fabs(v) > 1e-12)
                this->nnz++;
        }
    }

    // allocate arrays
    this->A = new double[this->nnz];
    this->Ap = new int[this->size+1];
    this->Ai = new int[this->nnz];

    int count = 0;
    this->Ap[0] = 0;
    for(int i = 0; i < this->size; i++) {
        for(int j = 0; j < this->size; j++) {
            double v = m->get(i, j);
            if (fabs(v) > 1e-12) {
                this->A[count] = v;
                this->Ai[count] = j;
                count++;
            }
        }
        this->Ap[i+1] = count;
    }
}

void CSRMatrix::add_from_CooMatrix(CooMatrix *m)
{
    this->size = m->get_size();
    this->nnz = m->is_complex() ? m->triplets_len_cplx() : m->triplets_len();
    this->_is_complex = m->is_complex();

    // allocate data
    free_data();

    this->Ap = new int[this->size+1];
    this->Ai = new int[this->nnz];
    if (is_complex())
        this->A_cplx = new cplx[this->nnz];
    else
        this->A = new double[this->nnz];

    // get data
    int *row = new int[nnz];
    int *col = new int[nnz];
    double *data = NULL;
    cplx *data_cplx = NULL;

    if (is_complex()) {
        data_cplx = new cplx[nnz];
        m->get_row_col_data(row, col, data_cplx);
    } else {
        data = new double[nnz];
        m->get_row_col_data(row, col, data);
    }

    // compute number of non-zero entries per row of A
    std::fill(this->Ap, this->Ap + this->size, 0);

    for (int n = 0; n < nnz; n++) {
        this->Ap[row[n]]++;
    }

    // cumsum the nnz per row to get this->row_ptr[]
    for(int i = 0, cumsum = 0; i < this->size; i++){
        int temp = this->Ap[i];
        this->Ap[i] = cumsum;
        cumsum += temp;
    }
    this->Ap[this->size] = nnz;

    // write Aj,Ax into Bj,Bx
    for(int n = 0; n < nnz; n++){
        int index  = row[n];
        int dest = this->Ap[index];

        this->Ai[dest] = col[n];
        if (is_complex())
            this->A_cplx[dest] = data_cplx[n];
        else
            this->A[dest] = data[n];

        this->Ap[index]++;
    }

    for(int i = 0, last = 0; i <= this->size; i++){
        int temp = this->Ap[i];
        this->Ap[i]  = last;
        last   = temp;
    }

    // free data
    delete[] row;
    delete[] col;
    if (data) delete[] data;
    if (data_cplx) delete[] data_cplx;
}

void CSRMatrix::add_from_CSCMatrix(CSCMatrix *m)
{
    _error("internal error: CSRMatrix::add_from_CSCMatrix(CSCMatrix *m) not implemented.");
}

void CSRMatrix::print()
{
    printf("\nCSR Matrix:\n");
    printf("size: %i\n", this->size);
    printf("nzz: %i\n", this->nnz);

    print_vector("row_ptr", this->Ap, this->size+1);
    print_vector("col_ind", this->Ai, this->nnz);
    if (is_complex())
        print_vector("data", this->A_cplx, this->nnz);
    else
        print_vector("data", this->A, this->nnz);
}

// *********************************************************************************************************************

void CSCMatrix::add_from_CooMatrix(CooMatrix *m)
{
    this->size = m->get_size();
    this->nnz = m->is_complex() ? m->triplets_len_cplx() : m->triplets_len();
    this->_is_complex = m->is_complex();

    // allocate data
    free_data();

    this->Ap = new int[this->size+1];
    this->Ai = new int[this->nnz];
    if (is_complex())
        this->A_cplx = new cplx[this->nnz];
    else
        this->A = new double[this->nnz];

    // get data
    int *row = new int[nnz];
    int *col = new int[nnz];
    double *data = NULL;
    cplx *data_cplx = NULL;

    if (is_complex()) {
        data_cplx = new cplx[nnz];
        m->get_row_col_data(col, row, data_cplx);
    } else {
        data = new double[nnz];
        m->get_row_col_data(col, row, data);
    }

    // compute number of non-zero entries per row of A
    std::fill(this->Ap, this->Ap + this->size, 0);

    for (int n = 0; n < nnz; n++) {
        this->Ap[row[n]]++;
    }

    // cumsum the nnz per row to get this->row_ptr[]
    for(int i = 0, cumsum = 0; i < this->size; i++){
        int temp = this->Ap[i];
        this->Ap[i] = cumsum;
        cumsum += temp;
    }
    this->Ap[this->size] = nnz;

    // write Aj, Ax into Bj, Bx
    for(int n = 0; n < nnz; n++){
        int index  = row[n];
        int dest = this->Ap[index];

        this->Ai[dest] = col[n];
        if (is_complex())
            this->A_cplx[dest] = data_cplx[n];
        else
            this->A[dest] = data[n];

        this->Ap[index]++;
    }

    for(int i = 0, last = 0; i <= this->size; i++){
        int temp = this->Ap[i];
        this->Ap[i]  = last;
        last   = temp;
    }

    // free data
    delete[] row;
    delete[] col;
    if (data) delete[] data;
    if (data_cplx) delete[] data_cplx;
}

void CSCMatrix::add_from_CSRMatrix(CSRMatrix *m)
{
    _error("internal error: CSRMatrix::add_from_CSCMatrix(CSCMatrix *m) not implemented.");
    /*
    this->size = m->get_size();
    this->nnz = m->get_nnz();
    this->_is_complex = m->is_complex();

    // allocate data
    free_data();

    this->Ap = new int[this->size+1];
    this->Ai = new int[this->nnz];
    if (is_complex())
        this->A_cplx = new cplx[this->nnz];
    else
        this->A = new double[this->nnz];

    // get data
    int *row = m->get_Ap();
    int *col = m->get_Ai();
    double *data = NULL;
    cplx *data_cplx = NULL;

    if (is_complex()) {
        data_cplx = m->get_A_cplx();
    } else {
        data = m->get_A();
    }

    //compute number of non-zero entries per column of A
    std::fill(this->Ap, this->Ap + this->size, 0);

    for (int n = 0; n < nnz; n++){
        this->Ap[this->Ai[n]]++;
    }

    //cumsum the nnz per column to get Bp[]
    for(int col = 0, cumsum = 0; col < this->size; col++){
        int temp  = this->Ap[col];
        this->Ap[col] = cumsum;
        cumsum += temp;
    }
    this->Ap[this->size] = nnz;

    for(int l = 0; l < this->size; l++){
        for(int jj = row[l]; jj < row[l+1]; jj++){
            int temp  = col[jj];
            int dest = this->Ap[temp];

            this->Ai[dest] = l;
            if (is_complex())
                this->A_cplx[dest] = data_cplx[jj];
            else
                this->A[dest] = data[jj];

            this->Ap[temp]++;
        }
    }

    for(int l = 0, last = 0; l <= this->size; l++){
        int temp = this->Ap[l];
        this->Ap[l] = last;
        last = temp;
    }
    */
}

void CSCMatrix::print()
{
    printf("\nCSC Matrix:\n");
    printf("size: %i\n", this->size);
    printf("nzz: %i\n", this->nnz);

    print_vector("col_ptr", this->Ap, this->size+1);
    print_vector("row_ind", this->Ai, this->nnz);
    if (is_complex())
        print_vector("data", this->A_cplx, this->nnz);
    else
        print_vector("data", this->A, this->nnz);
}

template <> Triple<cplx> *CooMatrix::get_list<cplx>() {
    return this->list_cplx;
}

void CooMatrix::set_zero() {
    if (this->_is_complex) {
        this->free_data<cplx>();
        this->list_cplx = NULL;
        this->list_last_cplx = NULL;
    } else {
        this->free_data<double>();
        this->list = NULL;
        this->list_last = NULL;
    }
    this->size = 0;
}
