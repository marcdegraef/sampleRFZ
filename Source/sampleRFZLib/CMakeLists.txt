project(sampleRFZLib)

include (${CMP_SOURCE_DIR}/cmpCMakeMacros.cmake)


#------------------------------------------------------------------------------
# Setup all the source files
set(sampleRFZLib_SRCS

  ${sampleRFZLib_SOURCE_DIR}/mod_kinds.f90
  ${sampleRFZLib_SOURCE_DIR}/mod_global.f90
  ${sampleRFZLib_SOURCE_DIR}/mod_io.f90
  ${sampleRFZLib_SOURCE_DIR}/mod_math.f90
  ${sampleRFZLib_SOURCE_DIR}/mod_rotations.f90
  ${sampleRFZLib_SOURCE_DIR}/mod_Lambert.f90
  ${sampleRFZLib_SOURCE_DIR}/mod_quaternions.f90
  ${sampleRFZLib_SOURCE_DIR}/mod_symmetry.f90
  ${sampleRFZLib_SOURCE_DIR}/mod_crystallography.f90
  ${sampleRFZLib_SOURCE_DIR}/mod_so3.f90
  ${sampleRFZLib_SOURCE_DIR}/mod_OMPsupport.f90

  ${sampleRFZLib_SOURCE_DIR}/program_mods/mod_sampleRFZ.f90

  ${sampleRFZLib_Additional_SRCS}
)

#-------------------------------------------------------------------------------
# Put just the *.c files into their own library so Visual Studio will compile
# everything correctly. VS can NOT have both fortran and C files in the same target
# add_library(sampleRFZLib_C STATIC ${sampleRFZLib_C_SRCS})
# set_target_properties (sampleRFZLib_C PROPERTIES
#   LINKER_LANGUAGE Fortran
#    # Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
#   FOLDER sampleRFZPublic
# )
# add_library(sampleRFZLib_Cpp STATIC ${sampleRFZLib_Cpp_SRCS})
# set_target_properties (sampleRFZLib_Cpp PROPERTIES
#   LINKER_LANGUAGE Fortran
#    # Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
#   FOLDER sampleRFZPublic
# )
# target_include_directories(sampleRFZLib_Cpp PRIVATE ${SHTfile_BINARY_DIR})


add_library(sampleRFZLib ${LIB_TYPE} ${sampleRFZLib_SRCS})

target_link_libraries(sampleRFZLib)

if (Fortran_COMPILER_NAME MATCHES "ifort.*")
  target_link_libraries(sampleRFZLib
                  ${MKL_LIBRARIES}
  )
endif()

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  target_link_libraries(sampleRFZLib
                  ${sampleRFZ_BLAS_LAPACK_LIBS}
                  ${sampleRFZ_FORTRAN_SUPPORT_LIBS}
                  gomp
  )
endif()

target_include_directories(sampleRFZLib
                          PRIVATE
                            "${sampleRFZ_SOURCE_DIR}/Source"
                          )

set_target_properties (sampleRFZLib PROPERTIES
                  LINKER_LANGUAGE Fortran
                  FOLDER sampleRFZPublic
                  )

# On Windows we need a few extra compile and linker flags. By making these PUBLIC they should
# flow down to every dependent target
if(WIN32)
  target_compile_options(sampleRFZLib PUBLIC
                          $<$<COMPILE_LANGUAGE:Fortran>:/assume:byterecl>
                          $<$<COMPILE_LANGUAGE:Fortran>:/assume:nounderscore>
                          $<$<COMPILE_LANGUAGE:Fortran>:/Qopenmp>
                          $<$<COMPILE_LANGUAGE:Fortran>:/Qdiag-disable:11082>
                          $<$<COMPILE_LANGUAGE:Fortran>:/Qdiag-disable:7025>
                          $<$<COMPILE_LANGUAGE:Fortran>:/Qip>
                          $<$<COMPILE_LANGUAGE:Fortran>:/traceback>
                          #$<$<COMPILE_LANGUAGE:Fortran>:/warn:all>
                        #  $<$<COMPILE_LANGUAGE:Fortran>:/heap-arrays>
                          )
  target_link_options(sampleRFZLib PUBLIC /STACK:2000000000 /machine:x64)
else()
  if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    # First, Flags that should be used under any build condition
    target_compile_options(sampleRFZLib PUBLIC $<$<COMPILE_LANGUAGE:Fortran>:-fopenmp>
        $<$<COMPILE_LANGUAGE:Fortran>:-fno-sign-zero>
        $<$<COMPILE_LANGUAGE:Fortran>:-fno-f2c>
        $<$<COMPILE_LANGUAGE:Fortran>:-std=gnu>
        $<$<COMPILE_LANGUAGE:Fortran>:-fbacktrace>
    )
    # Next, DEBUG compiler options.
    target_compile_options(sampleRFZLib PUBLIC $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-g>
                          $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-fcheck=bounds>
                          $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-O0>
    )
    # Next, RELEASE compiler options.
    target_compile_options(sampleRFZLib PUBLIC $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Release>>:-funroll-all-loops>
                          $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Release>>:-O3>
    )

  elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # First, Flags that should be used under any build condition
    target_compile_options(sampleRFZLib PUBLIC
                                            #$<$<COMPILE_LANGUAGE:Fortran>:-heap-arrays>
                                            $<$<COMPILE_LANGUAGE:Fortran>:-qopenmp>
                                            $<$<COMPILE_LANGUAGE:Fortran>:-assume byterecl>
                                            $<$<COMPILE_LANGUAGE:Fortran>:/assume:nounderscore>
                                            $<$<COMPILE_LANGUAGE:Fortran>:-stand f08 -diag-disable=7025>
                                            $<$<COMPILE_LANGUAGE:Fortran>:/traceback>
    )
    # Next, DEBUG compiler options.
    target_compile_options(sampleRFZLib PUBLIC $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-debug all>
                                            $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-O0>
                                            $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Debug>>:-g>
    )
    # Next, RELEASE compiler options.
    target_compile_options(sampleRFZLib PUBLIC $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CONFIG:Release>>:-O3>)
  else()
    message(WARNING "Fortran compiler not recognized. ${Fortran_COMPILER_NAME}")
  endif()
endif()


#-------------------------------------------------------------------------------
# This sets up the two variables install_dir and lib_install_dir
sampleRFZ_SetupInstallDirs()

#if(BUILD_SHARED_LIBS)
  INSTALL(TARGETS sampleRFZLib # sampleRFZLib_C sampleRFZLib_Cpp
    COMPONENT Applications
    EXPORT sampleRFZLibTargets
    RUNTIME DESTINATION ${install_dir}
    LIBRARY DESTINATION ${lib_install_dir}
    ARCHIVE DESTINATION ${lib_install_dir}
    BUNDLE DESTINATION "."
  )
#endif()

if(APPLE AND BUILD_SHARED_LIBS)
  set(QAB_INSTALL_DEST "lib")
  set(osx_app_name "libsampleRFZLib.dylib")
  get_property(sampleRFZBinaryDir GLOBAL PROPERTY sampleRFZBinaryDir)
  set(OSX_MAKE_STANDALONE_BUNDLE_CMAKE_SCRIPT "${sampleRFZ_BINARY_DIR}/sampleRFZLib/sampleRFZLib_CompleteTool.cmake")
  set(OPTIMIZE_BUNDLE_SHELL_SCRIPT            "${sampleRFZ_BINARY_DIR}/sampleRFZLib/sampleRFZLib_OptimizeTool.sh")
  set(PROJECT_INSTALL_DIR "lib")
  configure_file("${CMP_OSX_TOOLS_SOURCE_DIR}/CompleteTool.cmake.in"
                "${OSX_MAKE_STANDALONE_BUNDLE_CMAKE_SCRIPT}" @ONLY IMMEDIATE)
  configure_file("${CMP_OSX_TOOLS_SOURCE_DIR}/CompleteTool.sh.in"
                "${OPTIMIZE_BUNDLE_SHELL_SCRIPT}" @ONLY IMMEDIATE)
  #install(SCRIPT "${OSX_MAKE_STANDALONE_BUNDLE_CMAKE_SCRIPT}" COMPONENT ${QAB_COMPONENT})
endif()


# --------------------------------------------------------------------
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/sampleRFZ/sampleRFZLibTargetsConfigVersion.cmake"
  VERSION ${sampleRFZ_VERSION}
  COMPATIBILITY AnyNewerVersion
)
#if(BUILD_SHARED_LIBS)
  export(EXPORT sampleRFZLibTargets
    FILE "${CMAKE_CURRENT_BINARY_DIR}/sampleRFZ/sampleRFZLibTargets.cmake"
    NAMESPACE sampleRFZ::
  )
#endif()

configure_file(cmake/sampleRFZConfig.cmake
  "${CMAKE_CURRENT_BINARY_DIR}/sampleRFZ/sampleRFZConfig.cmake"
  @ONLY
)

set(ConfigPackageLocation ${lib_install_dir}/cmake/sampleRFZ)

#if(BUILD_SHARED_LIBS)
  install(EXPORT sampleRFZLibTargets
    FILE
      sampleRFZLibTargets.cmake
    NAMESPACE
      sampleRFZ::
    DESTINATION
      ${ConfigPackageLocation}
  )

#endif()

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/sampleRFZ/sampleRFZConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/sampleRFZ/sampleRFZLibTargetsConfigVersion.cmake"
  DESTINATION
    ${ConfigPackageLocation}
  COMPONENT
    Devel
)
