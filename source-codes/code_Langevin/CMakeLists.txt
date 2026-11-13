cmake_minimum_required (VERSION 2.6)
project (JETSCAPE CXX Fortran)

set (CMAKE_CXX_FLAGS_OSG "${CMAKE_CXX_FLAGS_RELEASE} -static")
set (CMAKE_Fortran_FLAGS_OSG "${CMAKE_Fortran_FLAGS_RELEASE} -static")

# additional flags for Intel compilers
# the Intel compilers can generate multiple binaries for different CPUs and store them in the same file
# the correct binary will then be chosen at runtime
# here, the -ax option generates two binaries for Intel CPUs with SSE4.2 and 4.1 support
# the -m option generates a fallback for generic CPUs with SSE2 support
# this appears to work with every CPU on the OSG
# the remaining flags are suggested by the Intel docs for performance
set (INTEL_FLAGS "-axsse4.2,sse4.1 -msse2 -ipo -ipo-jobs8 -no-prec-div")

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
  set (CMAKE_CXX_FLAGS_OSG "${CMAKE_CXX_FLAGS_OSG} ${INTEL_FLAGS}")
endif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")

if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  set (CMAKE_Fortran_FLAGS_OSG "${CMAKE_Fortran_FLAGS_OSG} ${INTEL_FLAGS} -heap-arrays")
endif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")



# add subprojects

# add_subdirectory (LBT)
add_subdirectory (diffusion)
add_subdirectory (fragPLUSrecomb)
#add_subdirectory (B2D)
#add_subdirectory (B2Jpsi)
#add_subdirectory (HMdecay)
#add_subdirectory (bmsap)
#add_subdirectory (RAA)

#
# installation
#

# install run script and dependency
install (PROGRAMS run-patch DESTINATION .)
install (FILES msg DESTINATION .)
#install (DIRECTORY hydroProfile DESTINATION .)

# create empty results directory
install (DIRECTORY DESTINATION results1)
install (DIRECTORY DESTINATION results2)
# install (DIRECTORY DESTINATION hydroProfile)

