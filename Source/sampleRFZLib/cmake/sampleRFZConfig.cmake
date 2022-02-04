include("${CMAKE_CURRENT_LIST_DIR}/sampleRFZLibTargets.cmake")


set(sampleRFZ_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/../../../include")
set(sampleRFZ_LIB_DIRS "${CMAKE_CURRENT_LIST_DIR}/../../../lib;${CMAKE_CURRENT_LIST_DIR}/../../../bin")

set(sampleRFZ_Fortran_COMPILER_NAME @Fortran_COMPILER_NAME@)
set(sampleRFZ_BUILD_SHARED_LIBS "@BUILD_SHARED_LIBS@")

if (sampleRFZ_Fortran_COMPILER_NAME MATCHES "gfortran.*")
  set(sampleRFZ_Fortran_RT_Libs gfortran @sampleRFZ_FORTRAN_SUPPORT_LIBS@)
  set(sampleRFZ_Fortran_Lib_Dir @GFortran_LIB_DIR@)
endif()
