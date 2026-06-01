### CMakeLists.txt ---

project(igor
VERSION
  ${igor_VERSION}
LANGUAGES
  CXX)

## #################################################################
## Build rules
## #################################################################

add_executable(${PROJECT_NAME}
  main.cpp
  legacy_main.cpp)

find_package(Git QUIET)
set(IGOR_GIT_COMMIT "unknown")
if(Git_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE IGOR_GIT_COMMIT
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET)
endif()
target_compile_definitions(${PROJECT_NAME} PRIVATE IGOR_GIT_COMMIT="${IGOR_GIT_COMMIT}")

target_link_libraries(${PROJECT_NAME} PRIVATE igor::Core)
target_link_libraries(${PROJECT_NAME} PRIVATE GSL::gsl)
target_link_libraries(${PROJECT_NAME} PRIVATE CLI11::CLI11)
target_link_libraries(${PROJECT_NAME} PRIVATE tomlplusplus::tomlplusplus)

if(APPLE)
  target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP)
else()
  target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
endif()

## #################################################################
## Install rules
## #################################################################

install(TARGETS ${PROJECT_NAME}
   BUNDLE DESTINATION bin
  RUNTIME DESTINATION bin)

######################################################################
### CMakeLists.txt ends here
