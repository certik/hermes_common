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
        putenv((char *)"PYTHONPATH=.:../..:../../../python");
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
    // Free the namespace. This frees all the dictionary items, so if there
    // are some numpy arrays (or your classes) in the namespace, they will be
    // deallocated at this time.
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
    // namespace_push() is a regular Cython function and
    // as such, it increfs the object "o" before storing it in the namespace,
    // but we want to steal the reference, so we decref it here (there is still
    // at least one reference stored in the dictionary this->_namespace, so
    // it's safe). This is so that
    //     this->push("i", c2py_int(5));
    // doesn't leak (c2py_int() creates a python reference and push() destroys
    // this python reference)
    Py_DECREF(o);
}

PyObject *Python::pull(const char *name)
{
    PyObject *tmp = namespace_pull(this->_namespace, name);
    // namespace_pull() is a regular Cython function and
    // as such, it increfs the result before returning it, but we only want to
    // borrow a reference, so we decref it here (there is still at least one
    // reference stored in the dictionary this->_namespace, so it's safe)
    // This is so that
    //     int i = py2c_int(this->pull("i"));
    // doesn't leak (pull() borrows the reference, py2c_int() doesn't do
    // anything with the reference, so no leak nor segfault happens)
    Py_DECREF(tmp);
    return tmp;
}
