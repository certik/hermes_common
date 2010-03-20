#include <iostream>
#include <stdexcept>

#include "_hermes_common_api_new.h"

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                              -1

void _assert(bool a)
{
    if (!a) throw std::runtime_error("Assertion failed.");
}

void test_basic()
{
    insert_object("i", c2py_int(5));
    cmd("i = i*2");
    int i = py2c_int(get_object("i"));
    _assert(i == 10);
}

void test_numpy()
{
    // double arrays
    cmd("from numpy import array");
    cmd("a = array(range(20), dtype='double')");
    cmd("assert a.strides == (8,)");
    cmd("b = a[::5]");
    cmd("assert b.strides == (40,)");
    double *A;
    int n;
    numpy2c_double_inplace(get_object("a"), &A, &n);
    _assert(n == 20);
    _assert(
            (fabs(A[0] - 0.)  < 1e-10) &&
            (fabs(A[1] - 1.)  < 1e-10) &&
            (fabs(A[2] - 2.) < 1e-10) &&
            (fabs(A[3] - 3.) < 1e-10)
           );
    numpy2c_double_inplace(get_object("b"), &A, &n);
    _assert(n == 4);
    _assert(
            (fabs(A[0] - 0.)  < 1e-10) &&
            (fabs(A[1] - 5.)  < 1e-10) &&
            (fabs(A[2] - 10.) < 1e-10) &&
            (fabs(A[3] - 15.) < 1e-10)
           );

    double a[3] = {1., 5., 3.};
    insert_object("A", c2numpy_double(a, 3));
    cmd("assert (A == array([1., 5., 3.])).all()");

    // integer arrays
    cmd("a = array(range(20), dtype='int32')");
    cmd("assert a.strides == (4,)");
    cmd("b = a[::5]");
    cmd("assert b.strides == (20,)");
    int *B;
    numpy2c_int_inplace(get_object("a"), &B, &n);
    _assert(n == 20);
    _assert(
            (B[0] == 0) &&
            (B[1] == 1) &&
            (B[2] == 2) &&
            (B[3] == 3)
           );
    numpy2c_int_inplace(get_object("b"), &B, &n);
    _assert(n == 4);
    _assert(
            (B[0] == 0) &&
            (B[1] == 5) &&
            (B[2] == 10) &&
            (B[3] == 15)
           );

    int b[3] = {1, 5, 3};
    insert_object("B", c2numpy_int(b, 3));
    cmd("assert (B == array([1, 5, 3])).all()");
}

int main(int argc, char* argv[])
{
    try {
        // This is a hack, this should be set somewhere else:
        putenv((char *)"PYTHONPATH=../..");
        Py_Initialize();
        PySys_SetArgv(argc, argv);
        if (import__hermes_common())
            throw std::runtime_error("hermes_common failed to import.");

        test_basic();
        test_numpy();

        return ERROR_SUCCESS;
    } catch(std::exception const &ex) {
        std::cout << "Exception raised: " << ex.what() << "\n";
        return ERROR_FAILURE;
    } catch(...) {
        std::cout << "Exception raised." << "\n";
        return ERROR_FAILURE;
    }
}
