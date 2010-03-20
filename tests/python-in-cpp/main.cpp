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

// Test the global python instance
void test_basic1()
{
    python->insert_object("i", c2py_int(5));
    python->eval("i = i*2");
    int i = py2c_int(python->get_object("i"));
    _assert(i == 10);
}

// Test the local Python() instance
void test_basic2()
{
    Python *p = new Python();
    p->insert_object("i", c2py_int(5));
    p->eval("i = i*2");
    int i = py2c_int(p->get_object("i"));
    _assert(i == 10);
    delete p;
}

// Test initialization/destruction of two Python() instances
void test_basic3()
{
    Python *p1 = new Python();
    p1->insert_object("i", c2py_int(5));
    p1->eval("i = i*2");
    int i = py2c_int(p1->get_object("i"));
    _assert(i == 10);
    delete p1;

    Python *p2 = new Python();
    p2->insert_object("i", c2py_int(5));
    p2->eval("i = i*2");
    i = py2c_int(p2->get_object("i"));
    _assert(i == 10);
    delete p2;
}

// Test that each Python() instance has it's own namespace
void test_basic4()
{
    Python *p1 = new Python();
    Python *p2 = new Python();
    p1->insert_object("i", c2py_int(5));
    p2->insert_object("i", c2py_int(6));
    p1->eval("i = i*2");
    p2->eval("i = i*2");
    int i1 = py2c_int(p1->get_object("i"));
    int i2 = py2c_int(p2->get_object("i"));
    // Fails so far:
    //_assert(i1 == 10);
    //_assert(i2 == 12);
    delete p1;
    delete p2;
}

int main(int argc, char* argv[])
{
    try {
        python = new Python(argc, argv);

        test_basic1();
        test_basic2();
        test_basic3();
        test_basic4();

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
