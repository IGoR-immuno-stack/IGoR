### CMakeLists.txt ---

# A consumer of the installed IGoR package, outside the IGoR build. It is not part of the main
# build: scripts/tests/test_package.sh installs IGoR into a scratch prefix and configures this
# project against it, once per scenario. What it guards is cmake/igorConfig.cmake.in, which
# nothing inside the repository consumes and which was broken for months without a failing test.

cmake_minimum_required(VERSION 3.28)

project(igor_package_test LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SCENARIO "all" CACHE STRING "all, model, core, streaming, umbrella or unknown")

if(SCENARIO STREQUAL "all")
  # No component list: every layer, each with its dependencies.
  find_package(igor REQUIRED)
  foreach(layer Core Math Alignment Model Inference Generation Streaming)
    if(NOT TARGET igor::${layer})
      message(FATAL_ERROR "find_package(igor) did not import igor::${layer}")
    endif()
  endforeach()
  add_executable(consumer model.cpp)
  target_link_libraries(consumer PRIVATE igor::Inference igor::Generation igor::Streaming)

elseif(SCENARIO STREQUAL "model")
  # One component: it comes with what it depends on, and with nothing else.
  find_package(igor REQUIRED COMPONENTS Model)
  foreach(layer Core Math Alignment Model)
    if(NOT TARGET igor::${layer})
      message(FATAL_ERROR "COMPONENTS Model did not import igor::${layer}")
    endif()
  endforeach()
  foreach(layer Inference Generation Streaming)
    if(TARGET igor::${layer})
      message(FATAL_ERROR "COMPONENTS Model imported igor::${layer}, which it does not depend on")
    endif()
  endforeach()
  add_executable(consumer model.cpp)
  target_link_libraries(consumer PRIVATE igor::Model)

elseif(SCENARIO STREQUAL "core")
  find_package(igor REQUIRED COMPONENTS Core)
  if(TARGET igor::Model OR TARGET igor::Alignment)
    message(FATAL_ERROR "COMPONENTS Core imported more than Core")
  endif()
  add_executable(consumer core.cpp)
  target_link_libraries(consumer PRIVATE igor::Core)

elseif(SCENARIO STREQUAL "streaming")
  # The one layer with an external dependency of its own (sparrow).
  find_package(igor REQUIRED COMPONENTS Streaming)
  add_executable(consumer streaming.cpp)
  target_link_libraries(consumer PRIVATE igor::Streaming)

elseif(SCENARIO STREQUAL "umbrella")
  find_package(igor REQUIRED QUIET)
  add_executable(consumer model.cpp)
  target_link_libraries(consumer PRIVATE igor::igor)

elseif(SCENARIO STREQUAL "unknown")
  # An unknown component is reported: not found, and no half-imported package.
  find_package(igor QUIET COMPONENTS NoSuchLayer)
  if(igor_FOUND)
    message(FATAL_ERROR "find_package(igor COMPONENTS NoSuchLayer) reported igor as found")
  endif()
  if(TARGET igor::Core)
    message(FATAL_ERROR "an unknown component still imported targets")
  endif()

else()
  message(FATAL_ERROR "unknown SCENARIO '${SCENARIO}'")
endif()

######################################################################
### CMakeLists.txt ends here
