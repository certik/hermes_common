#include <stdexcept>

#include "python_api.h"

static int python_count=0;

Python::Python()
{
    this->_init(-1, NULL);
}

Python::Python(int argc, char* argv[])
{
    this->_init(argc, argv);
}

void Python::_init(int argc, char* argv[])
{
    python_count++;
    if (python_count == 1) {
        // This is a hack, so that we can load hermes_common below. Some better
        // mechanism should be used instead, once we figure out how to do it:
        putenv((char *)"PYTHONPATH=../..");
        Py_Initialize();
        if (argc >= 0)
            PySys_SetArgv(argc, argv);
        if (import__hermes_common())
            throw std::runtime_error("hermes_common failed to import.");
    }
    this->_namespace = namespace_create();
}

Python::~Python()
{
    // free the namespace:
    Py_DECREF(this->_namespace);

    // free the interpreter if this was the last instance using it:
    python_count--;
    if (python_count == 0) {
        // don't finalize python, because the numpy package segfaults when
        // imported again:
        //Py_Finalize();
    }
}

void Python::print()
{
    namespace_print(_namespace);
}

void Python::exec(const char *text)
{
    run_cmd(text, this->_namespace);
}

void Python::push(const char *name, PyObject *o)
{
    namespace_push(this->_namespace, name, o);
    Py_DECREF(o);
}

PyObject *Python::pull(const char *name)
{
    return namespace_pull(this->_namespace, name);
}
