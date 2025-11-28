/*
 * ArrowCompatibility.h
 *
 *  Created on: Nov 28, 2025
 *      Author: IGoR Development Team
 *
 *  This source code is distributed as part of the IGoR software.
 *
 *  Compatibility workarounds for Apache Arrow with libc++
 */

#ifndef STREAMING_ARROW_COMPATIBILITY_H_
#define STREAMING_ARROW_COMPATIBILITY_H_

#include <string>
#include <string_view>

// Workaround for Arrow + libc++ char_traits issue with unsigned char
// Arrow's buffer.h uses util::bytes_view which requires std::basic_string_view<unsigned char>
// but libc++ doesn't provide char_traits<unsigned char> by default
namespace std {
template <>
struct char_traits<unsigned char> : public char_traits<char>
{
    using char_type = unsigned char;
    using int_type = int;
    using off_type = streamoff;
    using pos_type = streampos;
    using state_type = mbstate_t;
};
} // namespace std

#endif // STREAMING_ARROW_COMPATIBILITY_H_
