// MdspanCompat.h ---

#pragma once

/**
 * @file MdspanCompat.h
 * @brief Compatibility header for std::mdspan across different platforms.
 *
 * - macOS/Windows: Uses standard <mdspan> (libc++ / MSVC STL)
 * - Linux (GCC): Uses Kokkos reference implementation
 *
 * After including this header, std::mdspan and std::dextents are available.
 */

// Try standard <mdspan> first (Clang libc++, MSVC)
#if __has_include(<mdspan>) && !defined(ADJOINTCHECK_FORCE_KOKKOS_MDSPAN)
    #include <mdspan>

// Fallback to Kokkos mdspan (for GCC on Linux)
#elif __has_include(<experimental/mdspan>)
    #include <experimental/mdspan>
    // Bring Kokkos mdspan into std namespace for compatibility
    namespace std {
        using ::std::experimental::mdspan;
        using ::std::experimental::dextents;
        using ::std::experimental::extents;
        using ::std::experimental::layout_right;
        using ::std::experimental::layout_left;
        using ::std::experimental::default_accessor;
    }

#else
    #error "No mdspan implementation found. Install 'mdspan' package (conda-forge) or use a compiler with C++23 mdspan support."
#endif

//
// MdspanCompat.h ends here
