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
  main.cpp)

target_link_libraries(${PROJECT_NAME} PRIVATE igor::Core)
target_link_libraries(${PROJECT_NAME} PRIVATE GSL::gsl)
target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)

## #################################################################
## Install rules
## #################################################################

install(TARGETS ${PROJECT_NAME}
   BUNDLE DESTINATION bin
  RUNTIME DESTINATION bin)

######################################################################
### CMakeLists.txt ends here
