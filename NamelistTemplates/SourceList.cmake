
#---------------------------------------------------------------------
# Set some variables to shorten up the call to the function below
set(APP_DIR ${sampleRFZ_SOURCE_DIR}/NamelistTemplates)

#---------------------------------------------------------------------
# Aggregate all the OpenCL files that are needed
set(EMSoft_RESOURCE_FILES
  ${APP_DIR}/sampleRFZ.template
)

if(NOT EXISTS "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/NamelistTemplates")
  file(MAKE_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/NamelistTemplates")
endif()

foreach(file ${sampleRFZ_RESOURCE_FILES})
  file(COPY "${file}" DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/NamelistTemplates")
endforeach()

#---------------------------------------------------------------------
# This sets up the two variables install_dir and lib_install_dir
sampleRFZ_SetupInstallDirs()

#---------------------------------------------------------------------
# Create the Installation Rules
INSTALL(FILES ${sampleRFZ_RESOURCE_FILES}
  COMPONENT Applications
  DESTINATION ${extra_install_dir}/NamelistTemplates
)
