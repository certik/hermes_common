project(python-in-c)
add_executable(${PROJECT_NAME} main.cpp)

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(python-in-c ${BIN})
