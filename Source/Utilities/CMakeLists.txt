set(APP_DIR "${sampleRFZ_SOURCE_DIR}/Source/Utilities")
set(TMPLT_DIR "${sampleRFZ_SOURCE_DIR}/NamelistTemplates")
set(LIB_SEARCH_DIRS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})


Add_sampleRFZ_Executable(TARGET sampleRFZ
                       SOURCES ${APP_DIR}/sampleRFZ.f90
                       LINK_LIBRARIES sampleRFZLib 
                       SOLUTION_FOLDER sampleRFZPublic/Utilities
                       INSTALL_PROGRAM TRUE
                       INCLUDE_DIRS ${sampleRFZLib_BINARY_DIR}
                       )


