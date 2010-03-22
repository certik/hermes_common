include_directories(${hermes_common_SOURCE_DIR})

project(leaks)
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} ${HERMES_COMMON})
add_dependencies(${PROJECT_NAME} ${HERMES_COMMON})




# tests:
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
# Disable this test --- run it manually, as it requires some memory
#add_test(leaks ${BIN})
