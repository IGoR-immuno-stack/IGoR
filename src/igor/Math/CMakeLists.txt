### CMakeLists.txt ---

project(Math
VERSION
  ${${CMAKE_PROJECT_NAME}_VERSION}
LANGUAGES
  CXX)

## #################################################################
## Core Target Definition (add_library, export header, version, alias)
## #################################################################

add_library(${PROJECT_NAME} INTERFACE)

set_target_properties(${PROJECT_NAME}
  PROPERTIES
    VERSION   ${${PROJECT_NAME}_VERSION}
    SOVERSION ${${PROJECT_NAME}_VERSION_MAJOR})

add_library(${CMAKE_PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

## #################################################################
## Sources, Headers, Includes, Linking, Install
## #################################################################

target_sources(${PROJECT_NAME}
  PUBLIC FILE_SET HEADERS
    BASE_DIRS
      ${CMAKE_SOURCE_DIR}/src
    FILES
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}.h
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/HybridBuffer.h
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/HybridBuffer.tpp
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/Linalg.h
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/Linalg.tpp
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/MdspanCompat.h
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/Tensor.h
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/Tensor.tpp
)

target_include_directories(${PROJECT_NAME} INTERFACE
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)

## #################################################################
## Compile Definitions
## #################################################################

# Define IGOR_MATH_DEBUG for Debug builds
target_compile_definitions(${PROJECT_NAME} INTERFACE
  $<$<CONFIG:Debug>:IGOR_MATH_DEBUG>
)

## ###################################################################
## Install rules
## ###################################################################

install(
  TARGETS
    ${PROJECT_NAME}
  EXPORT
    ${CMAKE_PROJECT_NAME}${PROJECT_NAME}-targets
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
