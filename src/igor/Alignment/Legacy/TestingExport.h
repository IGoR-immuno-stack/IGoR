#pragma once

#include <igor/Alignment/Export.h>

// Marks a symbol that is internal to Alignment (not part of its installed public API) but still
// needs to cross the shared library boundary for whitebox tests to link against it directly
// (see AlignerInternal.h). Resolves to a real export/import only when ALIGNMENT_TESTING_ENABLED
// is defined (see tst/igor/Alignment/CMakeLists.txt), so production builds keep these symbols
// hidden.
#ifdef ALIGNMENT_TESTING_ENABLED
#  define ALIGNMENT_TESTING_EXPORT ALIGNMENT_EXPORT
#else
#  define ALIGNMENT_TESTING_EXPORT ALIGNMENT_NO_EXPORT
#endif
