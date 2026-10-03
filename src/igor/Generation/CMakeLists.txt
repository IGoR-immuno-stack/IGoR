### CMakeLists.txt ---

project(Generation
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

target_sources(${PROJECT_NAME}
  PRIVATE
    SamplingHandlerFactory.cpp
    Legacy/FastGenerator.cpp
    Legacy/FastSampling.cpp
    Legacy/SequenceGenerator.cpp
  PUBLIC FILE_SET HEADERS
    BASE_DIRS
      ${CMAKE_SOURCE_DIR}/src
      ${CMAKE_BINARY_DIR}/src
    FILES
      ${CMAKE_BINARY_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/Export.h
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}.h
      ${LAYER_DIR}/CategoricalSamplingHandler.h
      ${LAYER_DIR}/CategoricalSamplingHandler.tpp
      ${LAYER_DIR}/MarkovSamplingHandler.h
      ${LAYER_DIR}/MarkovSamplingHandler.tpp
      ${LAYER_DIR}/SamplingEngine.h
      ${LAYER_DIR}/SamplingEngine.tpp
      ${LAYER_DIR}/SamplingHandler.h
      ${LAYER_DIR}/SamplingHandler.tpp
      ${LAYER_DIR}/SamplingHandlerFactory.h
      ${LAYER_DIR}/SamplingHandlerFactory.tpp
      ${LAYER_DIR}/Legacy/FastGenerator.h
      ${LAYER_DIR}/Legacy/FastSampling.h
      ${LAYER_DIR}/Legacy/SequenceGenerator.h
)

target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/src>
)

# OpenMP is PUBLIC because Legacy/FastGenerator.h includes <omp.h>.
target_link_libraries(${PROJECT_NAME} PUBLIC
  ${CMAKE_PROJECT_NAME}::Model
  OpenMP::OpenMP_CXX
)

## #################################################################
## Install rules
## #################################################################

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
