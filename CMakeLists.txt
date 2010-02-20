project(hermes_common)
cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(PythonLibs REQUIRED)
find_package(NumPy REQUIRED)

include(UseCython)
# if the user specified PYTHON_INCLUDE_PATH, let's use that, otherwise let's
# try to guess using cmake/FindPythonLibs.cmake
if (NOT PYTHON_INCLUDE_PATH)
    find_package(PythonLibs REQUIRED)
endif (NOT PYTHON_INCLUDE_PATH)
find_package(NumPy REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH} ${NUMPY_INCLUDE_PATH})

enable_testing()

# this is done by cmake anyway, no matter what actual COMMAND we put here (!):
add_custom_target(test
    COMMAND ctest
)

#PYTHONPATH=${hermes_common_SOURCE_DIR}

add_subdirectory(tests)



set(HERMES_COMMON hermes_common)

set(SRC
    common.cpp
    matrix.cpp
    )

add_library(${HERMES_COMMON} SHARED ${SRC})

cython_add_module(_hermes_common stdcython.cpp)
target_link_libraries(_hermes_common ${HERMES_COMMON})
set_property(TARGET _hermes_common PROPERTY COMPILE_FLAGS "-g")
