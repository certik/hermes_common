#! /usr/bin/env python

from time import strftime, localtime

def convert_h():
    f = open("_hermes_common_api_new.h", "w")
    f.write("/* Generated by convert_api.py on %s */\n\n" % \
            strftime("%a %b %d %H:%M:%S %Y", localtime()));
    lines = open("_hermes_common_api.h").readlines()
    line = lines[0]
    while not line.startswith("static"):
        f.write(line)
        del lines[0]
        line = lines[0]

    while line.startswith("static"):
        line = line.replace("static", "extern")
        f.write(line)
        del lines[0]
        line = lines[0]

    f.write("""
extern int import__hermes_common(void);

#endif\n""")

def convert_cpp():
    f = open("_hermes_common_api_new.cpp", "w")
    f.write("/* Generated by convert_api.py on %s */\n\n" % \
            strftime("%a %b %d %H:%M:%S %Y", localtime()));
    lines = open("_hermes_common_api.h").readlines()
    line = lines[0]
    while not line.startswith("static"):
        del lines[0]
        line = lines[0]

    f.write("""\
#include "_hermes_common_api_new.h"

""")

    while line.startswith("static"):
        line = line.replace("static ", "")
        f.write(line)
        del lines[0]
        line = lines[0]

    f.write(line)
    line_old = line
    for line in lines:
        if line.startswith("#ifndef"):
            continue
        if line.startswith("#define"):
            continue
        if line.startswith("#endif"):
            if not line_old.startswith("    "):
                continue
        line = line.replace("static ", "")
        f.write(line)
        line_old = line

convert_h()
convert_cpp()
