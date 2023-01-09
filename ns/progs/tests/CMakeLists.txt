# ============= test_jacoby_solver
add_executable(test_jacoby_solver test_jacoby_solver.cpp)
target_link_libraries(test_jacoby_solver PUBLIC ${CFDLIB_TARGET})