#include <stdexcept>

#include "python_api.h"

static int python_count=0;

Python::Python()
{
    Python::Python(-1, NULL);
}

Python::Python(int argc, char* argv[])
{
    python_count++;
    if (python_count == 1) {
        // This is a hack:
        putenv((char *)"PYTHONPATH=../..");
        Py_Initialize();
        if (argc >= 0)
            PySys_SetArgv(argc, argv);
        if (import__hermes_common())
            throw std::runtime_error("hermes_common failed to import.");
    }
}

Python::~Python()
{
    // Don't finalize the interpreter, because for some reason,
    // import__hermes_common() doesn't work when called again in the second
    // interpreter (segfaults). So for now we just keep one interpreter and
    // that's it.
    /*
    python_count--;
    if (python_count == 0) {
        Py_Finalize();
    }
    */
}

void Python::eval(const char *text)
{
    cmd(text);
}

void Python::insert_object(const char *name, PyObject *o)
{
    ::insert_object(name, o);
}

PyObject *Python::get_object(const char *name)
{
    ::get_object(name);
}
