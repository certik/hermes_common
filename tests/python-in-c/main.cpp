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

        return ERROR_SUCCESS;
    } catch(std::exception const &ex) {
        std::cout << "Exception raised: " << ex.what() << "\n";
        return ERROR_FAILURE;
    } catch(...) {
        std::cout << "Exception raised." << "\n";
        return ERROR_FAILURE;
    }
}
