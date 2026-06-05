### igorSetupOpenMP.cmake ---

if(APPLE AND NOT DEFINED ENV{CONDA_PREFIX})
  execute_process(
    COMMAND brew --prefix libomp
    OUTPUT_VARIABLE _brew_libomp_prefix
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  if(_brew_libomp_prefix)
    set(OpenMP_libomp_LIBRARY "${_brew_libomp_prefix}/lib/libomp.dylib"
        CACHE FILEPATH "Path to libomp for OpenMP")
  endif()
endif()

if(WIN32 AND MSVC)
  set(OpenMP_RUNTIME_MSVC "llvm")
endif()

find_package(OpenMP REQUIRED COMPONENTS CXX)

######################################################################
### igorSetupOpenMP.cmake ends here
