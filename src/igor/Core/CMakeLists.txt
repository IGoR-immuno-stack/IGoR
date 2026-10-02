### CMakeLists.txt ---

project(Core
VERSION
  ${${CMAKE_PROJECT_NAME}_VERSION}
LANGUAGES
  CXX)

## #################################################################
## Core Target Definition (add_library, export header, version, alias)
## #################################################################

add_library(${PROJECT_NAME} SHARED)

# Step 1a keeps Core's flat export header (<igorCoreExport.h>); aligning it on the other
# layers' <igor/<Layer>/Export.h> scheme is step 1b of doc/LAYER_REFACTORING_PROPOSAL.md.
generate_export_header(${PROJECT_NAME}
  EXPORT_FILE_NAME "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}${PROJECT_NAME}Export.h")

set_target_properties(${PROJECT_NAME}
  PROPERTIES
    OUTPUT_NAME ${CMAKE_PROJECT_NAME}${PROJECT_NAME}
    VERSION   ${${PROJECT_NAME}_VERSION}
    SOVERSION ${${PROJECT_NAME}_VERSION_MAJOR})

add_library(${CMAKE_PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

## #################################################################
## Sources, Headers, Includes, Linking, Install
## #################################################################

target_sources(${PROJECT_NAME}
  PRIVATE
    Legacy/GeneticCode.cpp
    Legacy/IntStr.cpp
    Legacy/Utils.cpp

  PUBLIC FILE_SET HEADERS
    BASE_DIRS
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${CMAKE_CURRENT_BINARY_DIR}
    FILES
      Legacy/CoreEnums.h
      Legacy/DynamicSequenceMap.h
      Legacy/GeneticCode.h
      Legacy/IntStr.h
      Legacy/LayeredArray.h
      Legacy/SegmentSpan.h
      Legacy/SeqOffsetsMap.h
      Legacy/SeqTypeRegistry.h
      Legacy/StdTypedefs.h
      Legacy/Typedef.h
      Legacy/Utils.h
      ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}${PROJECT_NAME}Export.h
      ${CMAKE_CURRENT_BINARY_DIR}/Config
      ${CMAKE_CURRENT_BINARY_DIR}/Config.h
)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Config.h")
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Config")

target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/src>
)

if(WIN32)
  target_link_libraries(${PROJECT_NAME} PRIVATE ws2_32)
endif()

## ###################################################################
## Install rules
## ###################################################################

install(
  TARGETS
    ${PROJECT_NAME}
  EXPORT
    ${CMAKE_PROJECT_NAME}${PROJECT_NAME}-targets
  LIBRARY
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME
    DESTINATION ${CMAKE_INSTALL_BINDIR}
  FILE_SET HEADERS
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}
  INCLUDES
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
  EXPORT      ${CMAKE_PROJECT_NAME}${PROJECT_NAME}-targets
  FILE        ${CMAKE_PROJECT_NAME}${PROJECT_NAME}Targets.cmake
  NAMESPACE   ${CMAKE_PROJECT_NAME}::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${CMAKE_PROJECT_NAME}
)

export(
  EXPORT    ${CMAKE_PROJECT_NAME}${PROJECT_NAME}-targets
  NAMESPACE ${CMAKE_PROJECT_NAME}::
  FILE      ${CMAKE_BINARY_DIR}/${CMAKE_PROJECT_NAME}${PROJECT_NAME}Targets.cmake
)

######################################################################
### CMakeLists.txt ends here
