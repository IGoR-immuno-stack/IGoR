// MdspanCompat.h ---

#pragma once

/**
 * @file MdspanCompat.h
 * @brief Compatibility header for std::mdspan.
 *
 * Configured to prefer Kokkos reference implementation (experimental/mdspan)
 * to ensure access to submdspan for slicing operations.
 */

// Prefer Kokkos mdspan (experimental) to get submdspan support.
// NOTE: Requires NOMINMAX on Windows to prevent the max macro from breaking
// std::numeric_limits<size_t>::max() inside the Kokkos dynamic_extent header.
// To avoid clashes with GCC 13/14's incomplete <experimental/mdspan> which
// lacks dextents and submdspan, we check for Kokkos's internal headers directly
// and include them.
#if __has_include(<experimental/__p0009_bits/mdspan.hpp>)
    #include <experimental/__p0009_bits/default_accessor.hpp>
    #include <experimental/__p0009_bits/full_extent_t.hpp>
    #include <experimental/__p0009_bits/mdspan.hpp>
    #include <experimental/__p0009_bits/dynamic_extent.hpp>
    #include <experimental/__p0009_bits/extents.hpp>
    #include <experimental/__p0009_bits/layout_stride.hpp>
    #include <experimental/__p0009_bits/layout_left.hpp>
    #include <experimental/__p0009_bits/layout_right.hpp>
    
    #if __has_include(<experimental/__p2630_bits/submdspan.hpp>)
        #include <experimental/__p2630_bits/submdspan.hpp>
    #endif
    
    // Bring Kokkos mdspan into std namespace for compatibility
    namespace std {
#if defined(MDSPAN_IMPL_STANDARD_NAMESPACE)
        using ::MDSPAN_IMPL_STANDARD_NAMESPACE::mdspan;
        using ::MDSPAN_IMPL_STANDARD_NAMESPACE::dextents;
        using ::MDSPAN_IMPL_STANDARD_NAMESPACE::extents;
        using ::MDSPAN_IMPL_STANDARD_NAMESPACE::layout_right;
        using ::MDSPAN_IMPL_STANDARD_NAMESPACE::layout_left;
        using ::MDSPAN_IMPL_STANDARD_NAMESPACE::layout_stride;
        using ::MDSPAN_IMPL_STANDARD_NAMESPACE::default_accessor;
        
        #if MDSPAN_HAS_CXX_17
        using ::MDSPAN_IMPL_STANDARD_NAMESPACE::submdspan;
        #endif
        using ::MDSPAN_IMPL_STANDARD_NAMESPACE::full_extent;
        using ::MDSPAN_IMPL_STANDARD_NAMESPACE::full_extent_t;
#else
        using ::std::experimental::mdspan;
        using ::std::experimental::dextents;
        using ::std::experimental::extents;
        using ::std::experimental::layout_right;
        using ::std::experimental::layout_left;
        using ::std::experimental::layout_stride;
        using ::std::experimental::default_accessor;

        // Slicing support (submdspan)
        #if MDSPAN_HAS_CXX_17
        using ::std::experimental::submdspan;
        #endif
        using ::std::experimental::full_extent;
        using ::std::experimental::full_extent_t;
#endif
    }

// Fallback to standard <experimental/mdspan> if Kokkos internal bits are missing
#elif __has_include(<experimental/mdspan>)
    #include <experimental/mdspan>
    
    namespace std {
#if defined(MDSPAN_IMPL_STANDARD_NAMESPACE)
        using ::MDSPAN_IMPL_STANDARD_NAMESPACE::mdspan;
        using ::MDSPAN_IMPL_STANDARD_NAMESPACE::dextents;
        using ::MDSPAN_IMPL_STANDARD_NAMESPACE::extents;
        using ::MDSPAN_IMPL_STANDARD_NAMESPACE::layout_right;
        using ::MDSPAN_IMPL_STANDARD_NAMESPACE::layout_left;
        using ::MDSPAN_IMPL_STANDARD_NAMESPACE::layout_stride;
        using ::MDSPAN_IMPL_STANDARD_NAMESPACE::default_accessor;

        #if MDSPAN_HAS_CXX_17
        using ::MDSPAN_IMPL_STANDARD_NAMESPACE::submdspan;
        #endif
        using ::MDSPAN_IMPL_STANDARD_NAMESPACE::full_extent;
        using ::MDSPAN_IMPL_STANDARD_NAMESPACE::full_extent_t;
#else
        using ::std::experimental::mdspan;
        using ::std::experimental::dextents;
        using ::std::experimental::extents;
        using ::std::experimental::layout_right;
        using ::std::experimental::layout_left;
        using ::std::experimental::layout_stride;
        using ::std::experimental::default_accessor;

        // Slicing support (submdspan)
        #if MDSPAN_HAS_CXX_17
        using ::std::experimental::submdspan;
        #endif
        using ::std::experimental::full_extent;
        using ::std::experimental::full_extent_t;
#endif
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
