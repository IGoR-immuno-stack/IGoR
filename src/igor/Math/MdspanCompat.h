// MdspanCompat.h ---

#pragma once

/**
 * @file MdspanCompat.h
 * @brief Compatibility header for std::mdspan.
 *
 * Configured to prefer Kokkos reference implementation (experimental/mdspan)
 * to ensure access to submdspan for slicing operations.
 */

// Prefer Kokkos mdspan (experimental) to get submdspan support
#if __has_include(<experimental/mdspan>)
    #include <experimental/mdspan>
    
    // Bring Kokkos mdspan into std namespace for compatibility
    namespace std {
        using ::std::experimental::mdspan;
        using ::std::experimental::dextents;
        using ::std::experimental::extents;
        using ::std::experimental::layout_right;
        using ::std::experimental::layout_left;
        using ::std::experimental::layout_stride;
        using ::std::experimental::default_accessor;

        // Slicing support (submdspan)
        using ::std::experimental::submdspan;
        using ::std::experimental::full_extent;
        using ::std::experimental::full_extent_t;
    }

// Fallback to standard <mdspan> (C++23) if Kokkos is missing
// Note: Slicing (submdspan) might be missing here until C++26
#elif __has_include(<mdspan>)
    #include <mdspan>
    
    #ifndef IGOR_NO_SUBMDSPAN
        #define IGOR_NO_SUBMDSPAN
        #ifdef _MSC_VER
            #pragma message("Using standard <mdspan> which may lack submdspan. Slicing disabled.")
        #else
            #warning "Using standard <mdspan> which may lack submdspan. Slicing disabled."
        #endif
    #endif

#else
    #error "No mdspan implementation found. Install 'mdspan' package (conda-forge)."
#endif

// MdspanCompat.h ends here
