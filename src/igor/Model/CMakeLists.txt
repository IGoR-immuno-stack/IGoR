### CMakeLists.txt ---

project(Model
VERSION
  ${${CMAKE_PROJECT_NAME}_VERSION}
LANGUAGES
  CXX)

## #################################################################
## Core Target Definition (add_library, export header, version, alias)
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

target_sources(${PROJECT_NAME}
  PUBLIC FILE_SET HEADERS
    BASE_DIRS
      ${CMAKE_BINARY_DIR}/src
    FILES
      ${CMAKE_BINARY_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/Export.h
)

target_sources(${PROJECT_NAME}
  PRIVATE
    EventFactory.cpp
    InferenceHandlerFactory.cpp
    SamplingHandlerFactory.cpp
    Topology.cpp
  PUBLIC FILE_SET HEADERS
    BASE_DIRS
      ${CMAKE_SOURCE_DIR}/src
    FILES
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}.h
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/CategoricalInferenceHandler.h
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/CategoricalInferenceHandler.tpp
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/CategoricalSamplingHandler.h
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/CategoricalSamplingHandler.tpp
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/EventFactory.h
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/InferenceEngine.h
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/InferenceEngine.tpp
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/InferenceHandler.h
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/InferenceHandlerFactory.h
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/InferenceHandlerFactory.tpp
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/LegacyBridge.h
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/LegacyBridge.tpp
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/MarkovInferenceHandler.h
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/MarkovInferenceHandler.tpp
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/MarkovSamplingHandler.h
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/MarkovSamplingHandler.tpp
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/Navigator.h
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/RecombinationModel.h
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/RecombinationModel.tpp
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/SamplingEngine.h
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/SamplingEngine.tpp
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/SamplingHandler.h
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/SamplingHandler.tpp
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/SamplingHandlerFactory.h
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/SamplingHandlerFactory.tpp
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/Scenario.h
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/Topology.h
      )

target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/src>
)

target_link_libraries(${PROJECT_NAME} PUBLIC
  ${CMAKE_PROJECT_NAME}::Math
  ${CMAKE_PROJECT_NAME}::Core
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
