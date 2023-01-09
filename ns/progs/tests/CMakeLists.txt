# ============= test_jacoby_solver
add_executable(test_jacoby_solver test_jacoby_solver.cpp)
target_link_libraries(test_jacoby_solver PUBLIC ${CFDLIB_TARGET})

# ============= test_seidel_solver
add_executable(test_seidel_solver test_seidel_solver.cpp)
target_link_libraries(test_seidel_solver PUBLIC ${CFDLIB_TARGET})