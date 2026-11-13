cmake_minimum_required (VERSION 2.6)
project (fragPLUSrecomb Fortran)

add_subdirectory (src)

# install dependencies
install (FILES parameters_hd.dat max_wigner_c.dat max_wigner_b.dat recomb_c_raw.dat recomb_b_raw.dat DESTINATION ${PROJECT_NAME})
