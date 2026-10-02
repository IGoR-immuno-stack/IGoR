#pragma once

/**
 * \file JsonDetail.h
 * \brief Validation helpers shared by the events' JSON constructors.
 *
 * Internal header: it includes the full json.hpp, so it is included by .cpp files only and is
 * deliberately absent from Core's installed FILE_SET. An installed header that included this
 * one would put 25k lines of json.hpp in front of every consumer.
 *
 * The point of these helpers is that a JSON document moves errors from compile time to run
 * time, so every one of them has to be caught and named. A missing key, an unknown key and a
 * gap in the realization indices all raise, rather than producing an event that looks built
 * and is not.
 */

#include <nlohmann/json.hpp>

#include <algorithm>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <vector>

namespace igor::core::legacy {}
namespace igor::alignment::legacy {}
namespace igor::model::legacy::json_detail {
using namespace igor::core::legacy;
using namespace igor::alignment::legacy;

/// Keys every event node may carry. "parents" is model-level, written and read by
/// ModelJson, but it is listed here so that an event constructor does not reject it.
inline constexpr std::initializer_list<const char *> kEventKeys = {
        "type", "gene_class", "seq_type", "side", "priority", "nickname", "realizations", "parents"
};

inline const nlohmann::json &require(const nlohmann::json &node, const char *key)
{
    const auto it = node.find(key);
    if (it == node.end())
        throw std::runtime_error(std::string("event json: missing key \"") + key + "\"");
    return *it;
}

inline void reject_unknown_keys(const nlohmann::json &node,
                                std::initializer_list<const char *> allowed)
{
    for (const auto &entry : node.items()) {
        const bool known = std::any_of(allowed.begin(), allowed.end(),
                                       [&](const char *a) { return entry.key() == a; });
        if (!known)
            throw std::runtime_error("event json: unknown key \"" + entry.key() + "\"");
    }
}

inline void expect_type(const nlohmann::json &node, const char *expected)
{
    const auto type = require(node, "type").get<std::string>();
    if (type != expected)
        throw std::runtime_error("event json: expected type \"" + std::string(expected)
                                 + "\", got \"" + type + "\"");
}

/**
 * \brief The node's realizations, ordered by index.
 *
 * Also enforces the schema's one numbering rule: indices are exactly 0 to n-1. The rule is
 * not decoration. Every concrete add_realization() assigns the index itself, from the current
 * realization count, so a document with a gap or a duplicate would be silently renumbered and
 * would no longer match the marginal array it was written for. All six shipped models satisfy
 * the rule, so enforcing it costs nothing and closes that hole.
 */
inline std::vector<const nlohmann::json *> realizations_in_index_order(const nlohmann::json &node)
{
    const nlohmann::json &array = require(node, "realizations");
    std::vector<const nlohmann::json *> out;
    out.reserve(array.size());
    for (const auto &realization : array)
        out.push_back(&realization);

    std::sort(out.begin(), out.end(),
              [](const nlohmann::json *a, const nlohmann::json *b) {
                  return require(*a, "index").get<int>() < require(*b, "index").get<int>();
              });

    for (std::size_t position = 0; position != out.size(); ++position) {
        const int index = require(*out[position], "index").get<int>();
        if (index != static_cast<int>(position))
            throw std::runtime_error(
                    "event json: realization indices must run from 0 without gaps, found "
                    + std::to_string(index) + " where " + std::to_string(position)
                    + " was expected in event \""
                    + require(node, "nickname").get<std::string>() + "\"");
    }
    return out;
}

}  // namespace igor::model::legacy::json_detail
