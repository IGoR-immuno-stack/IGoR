include(FetchContent)

FetchContent_Declare(
  mdspan
  GIT_REPOSITORY https://github.com/kokkos/mdspan.git
  GIT_TAG        mdspan-0.6.0
)

FetchContent_MakeAvailable(mdspan)
