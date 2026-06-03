### CMakeLists.txt ---

enable_testing()

find_package(Catch2 3 REQUIRED)

function(igor_copy_runtime_dlls target)
    if(WIN32)
        add_custom_command(TARGET ${target} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_RUNTIME_DLLS:${target}>
                $<TARGET_FILE_DIR:${target}>
            COMMAND_EXPAND_LISTS
            VERBATIM
        )

        if(DEFINED ENV{CONDA_PREFIX})
            set(_conda_bin "$ENV{CONDA_PREFIX}/Library/bin")
            file(GLOB _igor_runtime_extra
                "${_conda_bin}/libomp140.x86_64.dll"
                "${_conda_bin}/gsl-*.dll"
            )
            if(_igor_runtime_extra)
                add_custom_command(TARGET ${target} POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E copy_if_different
                        ${_igor_runtime_extra}
                        $<TARGET_FILE_DIR:${target}>
                    COMMAND_EXPAND_LISTS
                    VERBATIM
                )
            endif()
        endif()
    endif()
endfunction()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/igor)

add_subdirectory(igor/Core)
add_subdirectory(igor/Math)
add_subdirectory(igor/Model)
add_subdirectory(igor/Streaming)

######################################################################
### CMakeLists.txt ends here
