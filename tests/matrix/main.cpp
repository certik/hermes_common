#include <iostream>
#include <stdexcept>

#include "matrix.h"

#include "_hermes_common_api.h"

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                              -1

void _assert(bool a)
{
    if (!a) throw std::runtime_error("Assertion failed.");
}

void test_matrix1()
{
    CooMatrix m(5);
    m.add(1, 3, 3.5);
    m.add(2, 3, 4.5);
    m.add(3, 4, 1.5);
    m.add(4, 2, 1.5);
    m.add(2, 3, 1);
    m.print();

    printf("----\n");
    CSRMatrix n(&m);
    n.print();
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

        test_matrix1();

        return ERROR_SUCCESS;
    } catch(std::exception const &ex) {
        std::cout << "Exception raised: " << ex.what() << "\n";
        return ERROR_FAILURE;
    } catch(...) {
        std::cout << "Exception raised." << "\n";
        return ERROR_FAILURE;
    }
}
