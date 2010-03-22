include_directories(${hermes_common_SOURCE_DIR})

project(cpp-callbacks)
add_executable(${PROJECT_NAME} main.cpp cpp_api.cpp)
target_link_libraries(${PROJECT_NAME} ${HERMES_COMMON})
add_dependencies(${PROJECT_NAME} ${HERMES_COMMON})

cython_add_module(my_api cpp_api.cpp)


# tests:
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(cpp-callbacks ${BIN})
