// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#ifndef __HERMES_COMMON_MATRIX_H
#define __HERMES_COMMON_MATRIX_H

#include "_hermes_common_api_new.h"

class Python {
public:
    Python();
    Python(int argc, char* argv[]);
    ~Python();
    void eval(const char *text);
    void insert_object(const char *name, PyObject *o);
    PyObject *get_object(const char *name);
};

#endif
