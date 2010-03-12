#include "python_api.h"

static int python_count=0;

Python::Python()
{
    python_count++;
    if (python_count == 1) {
        Py_Initialize();
    }
}

Python::~Python()
{
    python_count--;
    if (python_count == 0) {
        Py_Finalize();
    }
}

void Python::eval(const char *text)
{
    cmd(text);
}

void Python::insert_object(const char *name, PyObject *o)
{
    insert_object(name, o);
}

PyObject *Python::get_object(const char *name)
{
    get_object(name);
}
