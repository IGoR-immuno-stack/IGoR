### CMakeLists.txt ---

project(Model
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
    EventFactory.cpp
    Topology.cpp
    Legacy/Counter.cpp
    Legacy/Deletion.cpp
    Legacy/Dinuclmarkov.cpp
    Legacy/Errorrate.cpp
    Legacy/EventUtils.cpp
    Legacy/Genechoice.cpp
    Legacy/gene_to_seqtype_migr.cpp
    Legacy/HypermutationfullNmererrorrate.cpp
    Legacy/Hypermutationglobalerrorrate.cpp
    Legacy/Insertion.cpp
    Legacy/ModelJson.cpp
    Legacy/Model_Parms.cpp
    Legacy/Model_marginals.cpp
    Legacy/Rec_Event.cpp
    Legacy/Singleerrorrate.cpp
  PUBLIC FILE_SET HEADERS
    BASE_DIRS
      ${CMAKE_SOURCE_DIR}/src
      ${CMAKE_BINARY_DIR}/src
    FILES
      ${CMAKE_BINARY_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}/Export.h
      ${CMAKE_SOURCE_DIR}/src/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}.h
      ${LAYER_DIR}/EventFactory.h
      ${LAYER_DIR}/Forward.h
      ${LAYER_DIR}/LegacyBridge.h
      ${LAYER_DIR}/LegacyBridge.tpp
      ${LAYER_DIR}/Navigator.h
      ${LAYER_DIR}/RecombinationModel.h
      ${LAYER_DIR}/RecombinationModel.tpp
      ${LAYER_DIR}/Scenario.h
      ${LAYER_DIR}/Topology.h
      ${LAYER_DIR}/Legacy/AccumulationContext.h
      ${LAYER_DIR}/BoundTightness.h
      ${LAYER_DIR}/Legacy/BoundTightness.h
      ${LAYER_DIR}/Legacy/Counter.h
      ${LAYER_DIR}/Legacy/Deletion.h
      ${LAYER_DIR}/Legacy/Dinuclmarkov.h
      ${LAYER_DIR}/Legacy/Errorrate.h
      ${LAYER_DIR}/Legacy/EventUtils.h
      ${LAYER_DIR}/Legacy/EventTypedefs.h
      ${LAYER_DIR}/Legacy/ExplorationContext.h
      ${LAYER_DIR}/Legacy/Genechoice.h
      ${LAYER_DIR}/Legacy/gene_to_seqtype_migr.h
      ${LAYER_DIR}/Legacy/HypermutationfullNmererrorrate.h
      ${LAYER_DIR}/Legacy/Hypermutationglobalerrorrate.h
      ${LAYER_DIR}/Legacy/Insertion.h
      ${LAYER_DIR}/Legacy/JsonDetail.h
      ${LAYER_DIR}/JunctionGeometry.h
      ${LAYER_DIR}/Legacy/JunctionGeometry.h
      ${LAYER_DIR}/Legacy/ModelContext.h
      ${LAYER_DIR}/Legacy/ModelJson.h
      ${LAYER_DIR}/Legacy/Model_Parms.h
      ${LAYER_DIR}/Legacy/Model_marginals.h
      ${LAYER_DIR}/Legacy/QuerySequenceContext.h
      ${LAYER_DIR}/Legacy/Rec_Event.h
      ${LAYER_DIR}/SafetyMatrix.h
      ${LAYER_DIR}/Legacy/SafetyMatrix.h
      ${LAYER_DIR}/Legacy/Scenario.h
      ${LAYER_DIR}/Legacy/ScenarioContext.h
      ${LAYER_DIR}/Legacy/Singleerrorrate.h
      ${LAYER_DIR}/SpanProfile.h
      ${LAYER_DIR}/Legacy/SpanProfile.h
      ${LAYER_DIR}/UnfilledSegmentLengths.h
      ${LAYER_DIR}/Legacy/UnfilledSegmentLengths.h
)

target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/src>
)

# GSL is PUBLIC because Hypermutationglobalerrorrate.h includes gsl headers; nlohmann_json
# because Rec_Event.h names nlohmann::json in to_json()'s signature (headers only include
# <nlohmann/json_fwd.hpp>; the full json.hpp stays in the .cpp files).
# igor::Alignment is a transitional link: Legacy/QuerySequenceContext.h holds a JournaledQuery
# by value, and the events call the Aligner's nucleotide helpers. It goes away when the contexts
# leave Model (stage 2).
target_link_libraries(${PROJECT_NAME} PUBLIC
  ${CMAKE_PROJECT_NAME}::Math
  ${CMAKE_PROJECT_NAME}::Core
  ${CMAKE_PROJECT_NAME}::Alignment
  GSL::gsl
  nlohmann_json::nlohmann_json
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
