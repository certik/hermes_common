#include <iostream>
#include <stdexcept>

#include "python_api.h"

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                              -1

Python *python;

void _assert(bool a)
{
    if (!a) throw std::runtime_error("Assertion failed.");
}

void test_basic()
{
    python->insert_object("i", c2py_int(5));
    python->eval("i = i*2");
    int i = py2c_int(python->get_object("i"));
    _assert(i == 10);
}

int main(int argc, char* argv[])
{
    try {
        python = new Python();

        test_basic();

        delete python;

        return ERROR_SUCCESS;
    } catch(std::exception const &ex) {
        std::cout << "Exception raised: " << ex.what() << "\n";
        return ERROR_FAILURE;
    } catch(...) {
        std::cout << "Exception raised." << "\n";
        return ERROR_FAILURE;
    }
}
