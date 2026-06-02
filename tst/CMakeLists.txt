### CMakeLists.txt ---

enable_testing()

find_package(Catch2 3 REQUIRED)

if(WIN32 AND DEFINED ENV{CONDA_PREFIX})
    message(STATUS "Copying Conda DLLs to ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} to avoid Windows DLL shadowing...")
    file(GLOB ALL_CONDA_DLLS "$ENV{CONDA_PREFIX}/Library/bin/*.dll")
    if(ALL_CONDA_DLLS)
        file(COPY ${ALL_CONDA_DLLS} DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
    endif()
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/igor)

add_subdirectory(igor/Core)
add_subdirectory(igor/Streaming)
add_subdirectory(igor/Math)
add_subdirectory(igor/Model)

######################################################################
### CMakeLists.txt ends here
