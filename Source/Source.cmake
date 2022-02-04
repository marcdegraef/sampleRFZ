

set_property(GLOBAL PROPERTY sampleRFZ_PACKAGE_DEST_PREFIX ".")
# -----------------------------------------------------------------------
#
# -----------------------------------------------------------------------

get_property(sampleRFZ_PACKAGE_DEST_PREFIX GLOBAL PROPERTY sampleRFZ_PACKAGE_DEST_PREFIX)

include("${sampleRFZ_SOURCE_DIR}/Source/sampleRFZ_Functions.cmake")

add_subdirectory(${PROJECT_SOURCE_DIR}/Source/sampleRFZLib ${PROJECT_BINARY_DIR}/sampleRFZLib)


set(MODALITY_DIRS
    Utilities
)
# -----------------------------------------------------------------------
# Establish which modalities are going to be compiled
# -----------------------------------------------------------------------
foreach(MODALITY ${MODALITY_DIRS})
  option(sampleRFZ_ENABLE_${MODALITY} "Build sources and programs related to ${MODALITY}" ON)
endforeach()


# -----------------------------------------------------------------------
# Add a wrapper lib thats uses the enabled modality options to compile itself
# -----------------------------------------------------------------------
# add_subdirectory(${PROJECT_SOURCE_DIR}/Source/sampleRFZWrapperLib ${PROJECT_BINARY_DIR}/sampleRFZWrapperLib)

# -----------------------------------------------------------------------
# Add the executables
# -----------------------------------------------------------------------
foreach(MODALITY ${MODALITY_DIRS})
  if( "${sampleRFZ_ENABLE_${MODALITY}}" STREQUAL "ON" )
    message(STATUS "sampleRFZ: Enabling public ${MODALITY} Modality")
    add_subdirectory( ${PROJECT_SOURCE_DIR}/Source/${MODALITY} ${PROJECT_BINARY_DIR}/${MODALITY})
  endif()
endforeach()



# -----------------------------------------------------------------------
# Does the developer want to compile the GUI for sampleRFZ?
# -----------------------------------------------------------------------
# if( sampleRFZ_ENABLE_sampleRFZWorkbench )

#   INCLUDE (${sampleRFZ_SOURCE_DIR}/Support/cmp/cmpCMakeMacros.cmake )
#   # --------------------------------------------------------------------
#   # Find and Use the Qt5 Libraries
#   include(${sampleRFZ_SOURCE_DIR}/Support/cmp/ExtLib/Qt5Support.cmake)
#   set(sampleRFZWorkbench_Qt5_Components Core Widgets Network Gui Concurrent Svg Xml OpenGL PrintSupport )
#   CMP_AddQt5Support( "${sampleRFZWorkbench_Qt5_Components}"
#                     "FALSE"
#                     "${sampleRFZ_BINARY_DIR}"
#                     "sampleRFZWorkbench")

#   include(${PROJECT_SOURCE_DIR}/Source/sampleRFZWorkbench/SourceList.cmake)
# endif()
