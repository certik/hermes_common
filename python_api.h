// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef __HERMES_COMMON_PYTHON_API_H
#define __HERMES_COMMON_PYTHON_API_H

#include "_hermes_common_api_new.h"

/*
   This is a nice C++ Python API and the only header file that you should
   include in your code.

   You have to create an instance of the Python() class first, which will
   initialize the pointers to the conversion methods (like py2c_int, c2py_int,
   ...), defined in _hermes_common_api_new.h. Once you instantiated one
   Python() instance, then you don't have to worry about this at all. If you
   call the c2py_int (and similar methods) without instantiating Python()
   first, it will segfault (as they point to NULL).
*/

class Python {
public:
    Python();
    Python(int argc, char* argv[]);
    ~Python();
    void print();
    void eval(const char *text);
    void insert_object(const char *name, PyObject *o);
    PyObject *get_object(const char *name);
private:
    PyObject *_namespace;
    void _init(int argc, char* argv[]);
};

#endif
