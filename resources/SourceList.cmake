
#---------------------------------------------------------------------
# Set some variables to shorten up the call to the function below
set(APP_DIR ${sampleRFZ_SOURCE_DIR}/resources)

#---------------------------------------------------------------------
# Aggregate all the files that are needed
set(sampleRFZ_RESOURCE_FILES
  ${APP_DIR}/templatecodes.txt
)

if(NOT EXISTS "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/resources")
  file(MAKE_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/resources")
endif()

foreach(file ${sampleRFZ_RESOURCE_FILES})
  file(COPY "${file}" DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/resources")
endforeach()


#---------------------------------------------------------------------
# This sets up the two variables install_dir and lib_install_dir
sampleRFZ_SetupInstallDirs()

#---------------------------------------------------------------------
# Create the Installation Rules
INSTALL(FILES ${sampleRFZ_RESOURCE_FILES}
  COMPONENT Applications
  DESTINATION ${extra_install_dir}/resources
)

