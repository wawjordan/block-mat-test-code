add_subdirectory( constants )
add_subdirectory( functions )
add_subdirectory( types )
# Propagate lib_srcs list to parent
add_lib_sources( set_precision.f90 )
# add_exe_sources( test_program main_program.f90 )
# add_exe_sources( test_program dget_test.f90 )
# add_exe_sources( test_program band_test.f90 )
# add_exe_sources( test_program lu_test.f90 )
# add_exe_sources( test_program blk_band_test.f90 )
add_exe_sources( test_program tridiag_test.f90 )
