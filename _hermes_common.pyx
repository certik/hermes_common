# Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
# Distributed under the terms of the BSD license (see the LICENSE
# file for the exact terms).
# Email: hermes1d@googlegroups.com, home page: http://hpfem.org/


#-----------------------------------------------------------------------
# Common C++ <-> Python+NumPy conversion tools:

import sys
import traceback
# this is important to be called here, otherwise we can't use the NumPy C/API:
import_array()

global_namespace = {"verbose": False}

cdef api void cmd(char *text):
    """
    Runs the command "text" in the Python namespace.
    """
    n = run_cmd(text, global_namespace)
    global_namespace.update(n)

cdef api void set_verbose_cmd(int verbose):
    global_namespace["verbose"] = verbose

cdef api void insert_object(char *name, object o):
    """
    Inserts an object into the global namespace.

    Example 1:

    insert_object("a", c2py_int(3));
    cmd("print a");

    This prints "3".

    Example 2:

    int a[3] = {1, 5, 3};
    insert_object("A", c2numpy_int(a, 3));
    cmd("print A");

    This prints "[1  5  3]" (this is how the NumPy array is printed).

    Example 3:

    double a[3] = {1, 5, 3};
    insert_object("A", c2numpy_double(a, 3));
    cmd("print A");

    This prints "[ 1.  5.  3.]" (this is how the NumPy array is printed).
    """
    global_namespace.update({name: o})

cdef api object get_object(char *name):
    """
    Retrieves an object from the Python namespace.

    Example:

    // put into python:
    int a[3] = {1, 5, 3};
    insert_object("A", c2numpy_int(a, 3));

    // retrieve from python:
    double *A;
    int n;
    numpy2c_double_inplace(get_object("A"), &A, &n);
    """
    return global_namespace.get(name)

cdef api object c2py_int(int i):
    return i

cdef api int py2c_int(object i):
    return i

cdef api double py2c_double(object i):
    return i

cdef api object c2numpy_int(int *A, int len):
    """
    Construct the integer NumPy array by copying the data.
    """
    from numpy import empty
    cdef ndarray vec = empty([len], dtype="int32")
    cdef int *pvec = <int *>vec.data
    memcpy(pvec, A, len*sizeof(int))
    return vec

cdef api object c2numpy_int_inplace(int *A, int len):
    """
    Construct the integer NumPy array inplace (don't copy any data).
    """
    cdef npy_intp dim = len
    return PyArray_SimpleNewFromData(1, &dim, NPY_INT, A)

cdef api object c2numpy_double(double *A, int len):
    """
    Construct the double NumPy array by copying the data.
    """
    from numpy import empty
    cdef ndarray vec = empty([len], dtype="double")
    cdef double *pvec = <double *>vec.data
    memcpy(pvec, A, len*sizeof(double))
    return vec

cdef api object c2numpy_double_inplace(double *A, int len):
    """
    Construct the double NumPy array inplace (don't copy any data).
    """
    cdef npy_intp dim = len
    return PyArray_SimpleNewFromData(1, &dim, NPY_DOUBLE, A)

_AA = None

cdef api void numpy2c_int_inplace(object A_n, int **A_c, int *n):
    """
    Returns the C array, that points to the numpy array (inplace).

    Only if strides != sizeof(int), the data get copied first.
    """
    cdef ndarray A = A_n
    if not (A.nd == 1 and A.strides[0] == sizeof(int)):
        from numpy import array
        A = array(A.flat, dtype="int32")
        # this is needed so that numpy doesn't dealocate the arrays
        global _AA
        _AA = A
    n[0] = len(A)
    A_c[0] = <int *>(A.data)

cdef api void numpy2c_double_inplace(object A_n, double **A_c, int *n):
    """
    Returns the C array, that points to the numpy array (inplace).

    Only if strides != sizeof(double), the data get copied first.
    """
    cdef ndarray A = A_n
    if not (A.nd == 1 and A.strides[0] == sizeof(double)):
        from numpy import array
        A = array(A.flat, dtype="double")
        # this is needed so that numpy doesn't dealocate the arrays
        global _AA
        _AA = A
    n[0] = len(A)
    A_c[0] = <double *>(A.data)

cdef api object run_cmd(char *text, object namespace):
    try:
        verbose = namespace.get("verbose")
        if verbose:
            print "got a text:", text
        if verbose:
            print "evaluting in the namespace:"
            print namespace
        code = compile(text, "", "exec")
        eval(code, {}, namespace)
        if verbose:
            print "new namespace:"
            print namespace
        return namespace
    except SystemExit, e:
        try:
            exit_code = int(e)
        except:
            exit_code = -1
        exit(exit_code)
    except:
        etype, value, tb = sys.exc_info()
        s = "".join(traceback.format_exception(etype, value, tb))
        s = "Exception raised in the Python code:\n" + s
        throw_exception(s)
