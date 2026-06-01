### CMakeLists.txt ---

enable_testing()

find_package(Catch2 3 REQUIRED)

# On Windows, copy runtime DLLs (project shared libs + imported DLLs) next to
# test executables so the loader finds them in the application directory.
# This avoids 0xc0000142 (STATUS_DLL_NOT_FOUND) errors during Catch2 test
# discovery, which can occur when DL_PATHS / PATH manipulation is unreliable.
function(igor_copy_runtime_dlls target)
    if(WIN32)
        add_custom_command(TARGET ${target} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_RUNTIME_DLLS:${target}> $<TARGET_FILE_DIR:${target}>
            COMMAND_EXPAND_LISTS
        )
    endif()
endfunction()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/igor)

add_subdirectory(igor/Core)
add_subdirectory(igor/Streaming)
add_subdirectory(igor/Math)
add_subdirectory(igor/Model)

######################################################################
### CMakeLists.txt ends here
