#include <iostream>
#include <stdexcept>

#include "_hermes_common_api.h"

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                              -1

void _assert(bool a)
{
    if (!a) throw std::runtime_error("Assertion failed.");
}

void test_basic()
{
    cmd("print 'ok'");

    insert_object("i", c2py_int(5));
    cmd("print i");
    cmd("i = i*2");
    int i = py2c_int(get_object("i"));
    _assert(i == 10);
}

void test_numpy()
{
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
