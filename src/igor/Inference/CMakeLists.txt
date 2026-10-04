### CMakeLists.txt ---

project(Inference
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
    InferenceHandlerFactory.cpp
  PUBLIC FILE_SET HEADERS
    BASE_DIRS
      ${CMAKE_SOURCE_DIR}/src
      ${CMAKE_BINARY_DIR}/src
    FILES
      ${CMAKE_BINARY_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/Export.h
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}.h
      ${LAYER_DIR}/CategoricalInferenceHandler.h
      ${LAYER_DIR}/CategoricalInferenceHandler.tpp
      ${LAYER_DIR}/InferenceEngine.h
      ${LAYER_DIR}/InferenceEngine.tpp
      ${LAYER_DIR}/InferenceHandler.h
      ${LAYER_DIR}/InferenceHandlerFactory.h
      ${LAYER_DIR}/InferenceHandlerFactory.tpp
      ${LAYER_DIR}/MarkovInferenceHandler.h
      ${LAYER_DIR}/MarkovInferenceHandler.tpp
)

# Legacy files: everything under Legacy/, the stubs of promoted headers included.
target_sources(${PROJECT_NAME}
  PRIVATE
    Legacy/Bestscenarioscounter.cpp
    Legacy/Coverageerrcounter.cpp
    Legacy/Errorscounter.cpp
    Legacy/GenModel.cpp
    Legacy/Pgencounter.cpp
  PUBLIC FILE_SET HEADERS
    BASE_DIRS
      ${CMAKE_SOURCE_DIR}/src
      ${CMAKE_BINARY_DIR}/src
    FILES
      ${LAYER_DIR}/Legacy/Bestscenarioscounter.h
      ${LAYER_DIR}/Legacy/Coverageerrcounter.h
      ${LAYER_DIR}/Legacy/Errorscounter.h
      ${LAYER_DIR}/Legacy/GenModel.h
      ${LAYER_DIR}/Legacy/Pgencounter.h
)

target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/src>
)

# OpenMP is PUBLIC because Legacy/GenModel.h includes <omp.h>.
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
