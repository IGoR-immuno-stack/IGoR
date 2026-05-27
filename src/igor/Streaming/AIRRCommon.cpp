/*
 * AIRRCommon.cpp
 *
 *  Created on: Feb 4, 2026
 *      Author: IGoR Development Team
 *
 *  This source code is distributed as part of the IGoR software.
 */

#include <igor/Streaming/AIRRCommon.h>

#include <fstream>
#include <algorithm>
#include <stdexcept>

namespace igor::airr {

char delimiter_char(Delimiter delimiter)
{
    switch (delimiter) {
        case Delimiter::TAB:
            return '\t';
        case Delimiter::COMMA:
            return ',';
        case Delimiter::AUTO:
            return '\t'; // Shouldn't be called with AUTO
    }
    return '\t';
}

Delimiter detect_delimiter(const std::string& filepath)
{
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("airr::detect_delimiter: Cannot open file: " + filepath);
    }

    std::string first_line;
    if (!std::getline(file, first_line)) {
        throw std::runtime_error("airr::detect_delimiter: File is empty: " + filepath);
    }

    // Count tabs and commas in header line
    size_t tab_count = std::count(first_line.begin(), first_line.end(), '\t');
    size_t comma_count = std::count(first_line.begin(), first_line.end(), ',');

    // AIRR standard is TSV, so prefer tabs if both are present
    if (tab_count >= comma_count && tab_count > 0) {
        return Delimiter::TAB;
    } else if (comma_count > 0) {
        return Delimiter::COMMA;
    }

    // Default to TSV (AIRR standard)
    return Delimiter::TAB;
}

} // namespace igor::airr
