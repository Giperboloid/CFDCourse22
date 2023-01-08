# ============= sandbox
add_executable(sandbox sandbox.cpp)
target_link_libraries(sandbox PUBLIC ${CFDLIB_TARGET})

# ============= grids_test
add_executable(grids_test grids_test.cpp)
target_link_libraries(grids_test PUBLIC ${CFDLIB_TARGET})

# ============= poisson
add_executable(poisson poisson.cpp)
target_link_libraries(poisson PUBLIC ${CFDLIB_TARGET})

# ============= psi_omega_semi_imp
add_executable(psi_omega_semi_imp psi_omega_semi_imp.cpp)
target_link_libraries(psi_omega_semi_imp PUBLIC ${CFDLIB_TARGET})

# ============= test_jacoby_solver
add_executable(test_jacoby_solver test_jacoby_solver.cpp)
target_link_libraries(test_jacoby_solver PUBLIC ${CFDLIB_TARGET})
