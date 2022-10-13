# ============= grids_test
add_executable(grids_test grids_test.cpp)
target_link_libraries(grids_test PUBLIC ${CFDLIB_TARGET})

# ============= poisson
add_executable(poisson poisson.cpp)
target_link_libraries(poisson PUBLIC ${CFDLIB_TARGET})

# ============= psi_omega_semi_imp
add_executable(psi_omega_semi_imp psi_omega_semi_imp.cpp)
target_link_libraries(psi_omega_semi_imp PUBLIC ${CFDLIB_TARGET})
