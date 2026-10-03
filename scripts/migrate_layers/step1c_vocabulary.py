#!/usr/bin/env python3
"""Step 1c of doc/LAYER_REFACTORING_PROPOSAL.md, item 4: the Core vocabulary.

  1. Core/Legacy/{CoreEnums,StdTypedefs,Typedef}.h merge into Core/Types.h, in igor::core, with
     PascalCase type names (EventType, SeqSide, SeqType, SeqTypeString, EventName,
     MarginalArrayPtr, SeqOffset, CodonTable). Enumerators and index_type keep their names. A stub
     stays at each old path and aliases the legacy names.
  2. The getpid/gethostid shims leave Utils.h and the global namespace for Core/Platform.h.
  3. `promote.py core-vocabulary` then promotes IntStr and GeneticCode.
  4. Promoted and new code (everything under src/igor outside Legacy/) drops the legacy spelling
     of these names.

Run from the repository root; idempotent.
"""

import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import promote  # noqa: E402  (also changes directory to the repository root)

LICENSE_FROM = "src/igor/Core/SegmentSpan.h"

TYPES_BODY = """#pragma once

#include <cstdint>
#include <iosfwd>
#include <memory>
#include <string>
#include <unordered_map>

#include <igor/Core/Export.h>

/**
 * \\file Types.h
 * \\brief The dependency-free vocabulary of IGoR: the enums describing events and sequence
 *        segments, and the small typedefs every layer names.
 *
 * Merged from the legacy CoreEnums.h, StdTypedefs.h and Typedef.h (step 1c of
 * doc/LAYER_REFACTORING_PROPOSAL.md). Nothing here depends on any other IGoR header. The
 * Gene_class enums stay in Legacy/Utils.h for now because they carry export annotations and a
 * cluster of conversion functions.
 */

namespace igor::core {

enum EventType { GeneChoice_t, Deletion_t, Insertion_t, Dinuclmarkov_t, Undefined_t };

/// Which end of a constructed sequence segment an offset or a deletion refers to.
enum SeqSide { Five_prime = 0, Three_prime = 1, Undefined_side = 2 };

/**
 * The six sequence types of the legacy VDJ topology. Their values are the SeqTypeIds the
 * legacy registry assigns, which Model_Parms::read_model_parms() relies on, so enum-keyed code
 * still addresses the right slots.
 */
enum SeqType { V_gene_seq = 0, VD_ins_seq = 1, D_gene_seq = 2, DJ_ins_seq = 3, J_gene_seq = 4, VJ_ins_seq = 5 };

/// Name of a sequence type, as written in model files.
using SeqTypeString = std::string;

/// Type used as key for unordered maps, since an event cannot be instantiated.
using EventName = std::string;

/// Array of long doubles holding the marginal values.
using MarginalArrayPtr = std::unique_ptr<long double[]>;

/// Offset of an aligned sequence in the sequence_offsets maps. Characterizes the beginning
/// and the end of a sequence piece on the data sequence.
using SeqOffset = int;

using CodonTable = std::unordered_map<std::string, std::string>;

/// Index of a node, an event or a tensor dimension.
using index_type = std::int64_t;

/// Textual forms of the vocabulary enums, as written in model files and messages.
CORE_EXPORT SeqTypeString to_string(const SeqType);
CORE_EXPORT std::string to_string(const SeqSide);
CORE_EXPORT std::ostream &operator<<(std::ostream &, SeqSide);
CORE_EXPORT std::string operator+(const std::string &, SeqSide);
CORE_EXPORT std::string operator+(const std::string &, EventType);

} // namespace igor::core
"""

STUBS = {
    "CoreEnums": """using Event_type = igor::core::EventType;
using igor::core::GeneChoice_t;
using igor::core::Deletion_t;
using igor::core::Insertion_t;
using igor::core::Dinuclmarkov_t;
using igor::core::Undefined_t;

using Seq_side = igor::core::SeqSide;
using igor::core::Five_prime;
using igor::core::Three_prime;
using igor::core::Undefined_side;

using Seq_type = igor::core::SeqType;
using igor::core::V_gene_seq;
using igor::core::VD_ins_seq;
using igor::core::D_gene_seq;
using igor::core::DJ_ins_seq;
using igor::core::J_gene_seq;
using igor::core::VJ_ins_seq;

// The textual forms of these enums moved with them; legacy code calls them unqualified.
using igor::core::to_string;
using igor::core::operator<<;
using igor::core::operator+;
""",
    "StdTypedefs": """using Seq_type_String = igor::core::SeqTypeString;
using Rec_Event_name = igor::core::EventName;
using Marginal_array_p = igor::core::MarginalArrayPtr;
using Seq_Offset = igor::core::SeqOffset;
using UMCodonTable = igor::core::CodonTable;
""",
    "Typedef": """using igor::core::index_type;
""",
}

PLATFORM = """#pragma once

/**
 * \\file Platform.h
 * \\brief The two operating-system calls IGoR needs under one name on every platform.
 *
 * Out of Legacy/Utils.h since step 1c: the system headers are included here, at global scope,
 * and the shims live in igor::core instead of the global namespace.
 */

#include <cstdint>

#if defined(_WIN32)

#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif

#  include <process.h>
#  include <winsock2.h>
#  include <windows.h>
#  include <ws2tcpip.h>

#else

#  include <unistd.h>

#endif

namespace igor::core {

#if defined(_WIN32)

inline int portable_getpid()
{
    return _getpid();
}

inline uint32_t portable_gethostid()
{
    WSADATA wsaData;
    WSAStartup(MAKEWORD(2, 2), &wsaData);

    char hostname[256];
    gethostname(hostname, sizeof(hostname));

    struct addrinfo hints{};
    hints.ai_family = AF_INET;

    struct addrinfo *info;
    if (getaddrinfo(hostname, nullptr, &hints, &info) != 0)
        return 0;

    uint32_t res = ((struct sockaddr_in *)info->ai_addr)->sin_addr.S_un.S_addr;

    freeaddrinfo(info);
    WSACleanup();
    return res;
}

#else

inline int portable_getpid()
{
    return getpid();
}

inline uint32_t portable_gethostid()
{
    return gethostid();
}

#endif

} // namespace igor::core
"""


def stub(name, body, license_text):
    return (license_text + "#pragma once\n\n#include <igor/Core/Types.h>\n\n"
            f"// {name}.h was merged into igor/Core/Types.h (step 1c of doc/LAYER_REFACTORING_PROPOSAL.md).\n"
            "// This stub keeps the legacy include path and the legacy names for the code that has not\n"
            "// been promoted yet; it goes when its last consumer switches.\n"
            "namespace igor::core::legacy {\n\n" + body + "\n} // namespace igor::core::legacy\n")


def make_types():
    if os.path.exists("src/igor/Core/Types.h"):
        return
    license_text = promote.license_of(promote.read(LICENSE_FROM)).replace("SegmentSpan.h", "Types.h")
    # keep the history of the largest of the three under the new name
    subprocess.check_call(["git", "mv", "src/igor/Core/Legacy/CoreEnums.h", "src/igor/Core/Types.h"])
    promote.write("src/igor/Core/Types.h", license_text + TYPES_BODY)
    for name, body in STUBS.items():
        path = f"src/igor/Core/Legacy/{name}.h"
        promote.write(path, stub(name, body, license_text.replace("Types.h", f"{name}.h")))
        subprocess.check_call(["git", "add", path])
    cmake = "src/igor/Core/CMakeLists.txt"
    text = promote.read(cmake)
    anchor = "      ${LAYER_DIR}/Legacy/CoreEnums.h\n"
    if "${LAYER_DIR}/Types.h" not in text:
        text = text.replace(anchor, "      ${LAYER_DIR}/Platform.h\n      ${LAYER_DIR}/Types.h\n" + anchor, 1)
    promote.write(cmake, text)


def make_platform():
    path = "src/igor/Core/Platform.h"
    if os.path.exists(path):
        return
    license_text = promote.license_of(promote.read(LICENSE_FROM)).replace("SegmentSpan.h", "Platform.h")
    promote.write(path, license_text + PLATFORM)
    subprocess.check_call(["git", "add", path])
    utils = "src/igor/Core/Legacy/Utils.h"
    text = promote.read(utils)
    start = text.index("#if defined(_WIN32)\n\n#  ifndef WIN32_LEAN_AND_MEAN")
    end = text.index("#if defined(_MSC_VER)\n#  include <intrin.h>")
    text = text[:start] + "#include <igor/Core/Platform.h>\n\n" + text[end:]
    promote.write(utils, text)
    source = "src/igor/Core/Legacy/Utils.cpp"
    text = promote.read(source)
    text = text.replace("::portable_getpid()", "igor::core::portable_getpid()")
    text = text.replace("::portable_gethostid()", "igor::core::portable_gethostid()")
    text = text.replace("igor::core::igor::core::", "igor::core::")
    promote.write(source, text)


def modernize_new_code():
    for dirpath, _, files in os.walk("src/igor"):
        if "/Legacy" in dirpath:
            continue
        for name in files:
            if not name.endswith((".h", ".cpp", ".tpp")):
                continue
            path = os.path.join(dirpath, name)
            text = promote.read(path)
            new = promote.modernize_core_names(text, in_core=path.startswith("src/igor/Core/"))
            if new != text:
                promote.write(path, new)


def main():
    make_types()
    make_platform()
    for stem, symbols, renames in promote.PROMOTIONS["core-vocabulary"][1]:
        done = promote.promote("Core", stem, symbols, renames)
        print(f"{'promoted' if done else 'skipped '} Core/{stem}")
    promote.restub()
    modernize_new_code()
    print("step1c_vocabulary: done")


if __name__ == "__main__":
    main()
