### igorSetupOpenMP.cmake ---

if (APPLE)
  find_library(OpenMP_LIBRARY
    NAMES omp
  )
  find_path(OpenMP_INCLUDE_DIR
    omp.h
    PATHS
    $ENV{CONDA_PREFIX}/lib/clang/4.0.1/include
  )
  mark_as_advanced(OpenMP_LIBRARY OpenMP_INCLUDE_DIR)
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(OpenMP DEFAULT_MSG
    OpenMP_LIBRARY OpenMP_INCLUDE_DIR)

  if (OpenMP_FOUND)
    set(OpenMP_LIBRARIES ${OpenMP_LIBRARY})
    set(OpenMP_INCLUDE_DIRS ${OpenMP_INCLUDE_DIR})
    set(OpenMP_COMPILE_OPTIONS -Xpreprocessor -fopenmp)

    add_library(OpenMP::OpenMP SHARED IMPORTED)
    set_target_properties(OpenMP::OpenMP PROPERTIES
      IMPORTED_LOCATION ${OpenMP_LIBRARIES}
      INTERFACE_INCLUDE_DIRECTORIES "${OpenMP_INCLUDE_DIRS}"
      INTERFACE_COMPILE_OPTIONS "${OpenMP_COMPILE_OPTIONS}"
    )
  endif()
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lpthread")
elseif(WIN32)
  find_library(OpenMP_LIBRARY
    NAMES libomp
    PATHS "$ENV{CONDA_PREFIX}/Library/lib"
    NO_DEFAULT_PATH
  )
  find_file(OpenMP_DLL
    NAMES libomp.dll
    PATHS "$ENV{CONDA_PREFIX}/Library/bin"
    NO_DEFAULT_PATH
  )
  find_path(OpenMP_INCLUDE_DIR
    omp.h
    PATHS "$ENV{CONDA_PREFIX}/Library/include"
    NO_DEFAULT_PATH
  )
  mark_as_advanced(OpenMP_LIBRARY OpenMP_DLL OpenMP_INCLUDE_DIR)
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(OpenMP DEFAULT_MSG
    OpenMP_LIBRARY OpenMP_DLL OpenMP_INCLUDE_DIR)

  if(OpenMP_FOUND)
    add_library(OpenMP::OpenMP_CXX SHARED IMPORTED GLOBAL)
    set_target_properties(OpenMP::OpenMP_CXX PROPERTIES
      IMPORTED_LOCATION "${OpenMP_DLL}"
      IMPORTED_IMPLIB   "${OpenMP_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${OpenMP_INCLUDE_DIR}"
      INTERFACE_COMPILE_OPTIONS "/openmp:llvm"
    )
  endif()
else()
  find_package(OpenMP REQUIRED)
endif()

######################################################################
### igorSetupOpenMP.cmake ends here
