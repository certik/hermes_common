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

// Test the creation of a big numpy array
void use_big_numpy_array1(int size)
{
    Python *p = new Python();
    p->push("i", c2py_int(size));
    p->exec("import numpy");
    p->exec("A = numpy.zeros((i,), dtype='double')");
    p->exec("print A.shape, 'size (MB):', A.nbytes/1024.**2");
    delete p;
}

void test_leaks1()
{
    // If the matrices are not deallocated, this should take ~760GB:
    for (int i=0; i<5; i++)
        // Roughly 0.76GB:
        use_big_numpy_array1(5*pow(10,8));
}

void test_dealloc()
{
    Python *p = new Python();
    p->exec("from test_dealloc import A");
    p->exec("a = A()");
    p->exec("print 'A() was created'");
    delete p;
}

// Test gettting the data from the array back
// ...

int main(int argc, char* argv[])
{
    try {
        test_leaks1();
        test_dealloc();

        return ERROR_SUCCESS;
    } catch(std::exception const &ex) {
        std::cout << "Exception raised: " << ex.what() << "\n";
        return ERROR_FAILURE;
    } catch(...) {
        std::cout << "Exception raised." << "\n";
        return ERROR_FAILURE;
    }
}
