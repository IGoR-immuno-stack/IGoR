### CMakeLists.txt ---

project(Core
VERSION
  ${${CMAKE_PROJECT_NAME}_VERSION}
LANGUAGES
  CXX)

## #################################################################
## Target Definition (add_library, export header, version, alias)
## #################################################################

add_library(${PROJECT_NAME} SHARED)

generate_export_header(${PROJECT_NAME}
  EXPORT_FILE_NAME "${CMAKE_CURRENT_BINARY_DIR}/Export.h")

set_target_properties(${PROJECT_NAME}
  PROPERTIES
    OUTPUT_NAME ${CMAKE_PROJECT_NAME}${PROJECT_NAME}
    VERSION   ${${PROJECT_NAME}_VERSION}
    SOVERSION ${${PROJECT_NAME}_VERSION_MAJOR})

add_library(${CMAKE_PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

## #################################################################
## Sources, Headers, Includes, Linking, Install
## #################################################################

set(LAYER_DIR ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME})
set(LAYER_BINARY_DIR ${CMAKE_BINARY_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME})

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Config.h")

target_sources(${PROJECT_NAME}
  PRIVATE
    AlignmentData.cpp
    Types.cpp
    GeneticCode.cpp
    IntStr.cpp
  PUBLIC FILE_SET HEADERS
    BASE_DIRS
      ${CMAKE_SOURCE_DIR}/src
      ${CMAKE_BINARY_DIR}/src
    FILES
      ${LAYER_BINARY_DIR}/Export.h
      ${LAYER_BINARY_DIR}/Config.h
      ${LAYER_DIR}/AlignmentData.h
      ${LAYER_DIR}/Platform.h
      ${LAYER_DIR}/Types.h
      ${LAYER_DIR}/DynamicSequenceMap.h
      ${LAYER_DIR}/GeneticCode.h
      ${LAYER_DIR}/IntStr.h
      ${LAYER_DIR}/LayeredArray.h
      ${LAYER_DIR}/SegmentSpan.h
      ${LAYER_DIR}/SeqOffsetsMap.h
      ${LAYER_DIR}/SeqTypeRegistry.h
)

# Legacy files: everything under Legacy/, the stubs of promoted headers included.
target_sources(${PROJECT_NAME}
  PRIVATE
    Legacy/Utils.cpp
  PUBLIC FILE_SET HEADERS
    BASE_DIRS
      ${CMAKE_SOURCE_DIR}/src
      ${CMAKE_BINARY_DIR}/src
    FILES
      ${LAYER_DIR}/Legacy/AlignmentData.h
      ${LAYER_DIR}/Legacy/CoreEnums.h
      ${LAYER_DIR}/Legacy/DynamicSequenceMap.h
      ${LAYER_DIR}/Legacy/GeneticCode.h
      ${LAYER_DIR}/Legacy/IntStr.h
      ${LAYER_DIR}/Legacy/LayeredArray.h
      ${LAYER_DIR}/Legacy/SegmentSpan.h
      ${LAYER_DIR}/Legacy/SeqOffsetsMap.h
      ${LAYER_DIR}/Legacy/SeqTypeRegistry.h
      ${LAYER_DIR}/Legacy/StdTypedefs.h
      ${LAYER_DIR}/Legacy/Typedef.h
      ${LAYER_DIR}/Legacy/Utils.h
)

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
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
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
