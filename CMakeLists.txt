cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(HERMES_COMMON hermes_common)

set(SRC
    common.cpp
    matrix.cpp
    )

add_library(${HERMES_COMMON} SHARED ${SRC})
