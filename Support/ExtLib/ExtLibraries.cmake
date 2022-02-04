
#------------------------------------------------------------------------------
# Find the Intel Math Kernel Library (MKL) which has FFT functions
# On mac systems, we will also need to build up the RPATH
if (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # Define the interface layers and link type for MKL
  set(MKL_Link_Type Static)
  set(MKL_Interface_Layer 32)
  set(MKL_ThreadingLayer Sequential)
  set(MKL_OpenMP_Library iomp5)
  set(MKL_F95Interface BLAS95 LAPACK95)
  find_package(MKL REQUIRED)
  if(NOT MKL_FOUND)
    message(FATAL_ERROR "MKL is Required when using the Intel Fortran Compiler")
  endif()
  set(sampleRFZ_OpenMP_LIBRARY ${MKL_${MKL_OpenMP_Library}_LIBRARY})
  if(sampleRFZ_OpenMP_LIBRARY)
    get_filename_component(sampleRFZ_OpenMP_LIB_DIR ${sampleRFZ_OpenMP_LIBRARY} DIRECTORY)
    message(STATUS "sampleRFZ_OpenMP_LIB_DIR: ${sampleRFZ_OpenMP_LIB_DIR}")
    set(sampleRFZ_OpenMP_LIB_DIR ${sampleRFZ_OpenMP_LIB_DIR} CACHE PATH "")
    get_property(sampleRFZSearchDirs GLOBAL PROPERTY sampleRFZSearchDirs)
    file(APPEND "${sampleRFZSearchDirs}" "${sampleRFZ_OpenMP_LIB_DIR};")
  endif()


endif()


#------------------------------------------------------------------------------
# Find the GFotran Specific or matched libraries
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")

  set(sampleRFZ_FORTRAN_SUPPORT_LIBS ${sampleRFZ_FORTRAN_SUPPORT_LIBS} gcc_eh gomp)
  # Find an OpenMP Library
  set(sampleRFZ_OpenMP_LIBRARY gomp)

  # Find BLAS/LAPACK Library
  if(APPLE)
    find_library(sampleRFZ_BLAS_LAPACK_LIBS Accelerate)
  else()
    find_package(LAPACK REQUIRED)
    set(sampleRFZ_BLAS_LAPACK_LIBS ${LAPACK_LIBRARIES})
  endif()
endif()

#include_directories(${JSONFORTRAN_INCLUDE_DIR} ${FFTW3_INCLUDE_DIR} ${CLFortran_INCLUDE_DIR})


