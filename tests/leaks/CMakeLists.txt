include_directories(${hermes_common_SOURCE_DIR})

project(leaks)
add_executable(${PROJECT_NAME} main.cpp cpp_api.cpp)
target_link_libraries(${PROJECT_NAME} ${HERMES_COMMON})
add_dependencies(${PROJECT_NAME} ${HERMES_COMMON})

cython_add_module(my_api2 cpp_api.cpp)

# tests:
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(leaks ${BIN})
