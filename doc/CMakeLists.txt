# configured documentation tools and intermediate build results
set(BINARY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_build")
set_directory_properties(PROPERTIES
  ADDITIONAL_MAKE_CLEAN_FILES "${BINARY_BUILD_DIR};${CMAKE_CURRENT_BINARY_DIR}/public"
  )

if(BUILD_DOC_DOXYGEN)
  find_package(Doxygen REQUIRED)

  configure_file(Doxyfile.in ${BINARY_BUILD_DIR}/Doxyfile IMMEDIATE @ONLY)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/doxygen-theme/custom.css
    ${BINARY_BUILD_DIR}/doxygen-theme/custom.css
    COPYONLY)
  add_custom_target(doc
      WORKING_DIRECTORY ${BINARY_BUILD_DIR}
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/public/doxygen
      COMMAND ${DOXYGEN_EXECUTABLE} ${BINARY_BUILD_DIR}/Doxyfile
      COMMENT "Building Doxygen doc")
endif()
