#ifndef STDCYTHON_H
#define STDCYTHON_H

// To avoid compilation warnings:
#undef _XOPEN_SOURCE
#undef _POSIX_C_SOURCE

#include "Python.h"

#ifdef __cplusplus
extern "C" {
#endif

void throw_exception(char *text);

#ifdef __cplusplus
}
#endif

#endif

