### CMakeLists.txt for igor-convert ---

project(igor-convert
  VERSION ${${CMAKE_PROJECT_NAME}_VERSION}
  LANGUAGES CXX)

add_executable(${PROJECT_NAME}
  igor-convert.cpp)

target_link_libraries(${PROJECT_NAME}
  PRIVATE
    igor::Streaming
    igor::Core)

install(
  TARGETS ${PROJECT_NAME}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

######################################################################
### CMakeLists.txt ends here
