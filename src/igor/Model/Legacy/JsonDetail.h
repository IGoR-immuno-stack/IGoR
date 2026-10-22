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

#include <igor/Core/Types.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <initializer_list>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace igor::model::legacy::json_detail {

/// Keys every event node may carry. The nickname is not one: it is the node's key in the
/// document's "events" object, and is handed to the constructor beside the node.
inline constexpr std::initializer_list<const char *> kEventKeys = { "type", "priority", "targets",
                                                                    "realizations" };

/// A gene choice's keys: the common ones plus "library", the template library it draws from.
/// No other event has a library of its own, so a node of another type carrying one is
/// refused like any other unknown key.
inline constexpr std::initializer_list<const char *> kGeneChoiceKeys = { "type", "library", "priority",
                                                                         "targets", "realizations" };

/// A realization's keys: its index, the tensor coordinate, and its label, its identity.
inline constexpr std::initializer_list<const char *> kRealizationKeys = { "index", "label" };

/// A gene choice's realization also carries its template. In a model document the template sits
/// in the "library" section instead, and the reader puts it back before the node reaches here.
inline constexpr std::initializer_list<const char *> kTemplatedRealizationKeys = { "index", "label", "sequence" };

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

/// "5'" and "3'", as the model document spells an end. An event that acts on no particular end
/// leaves the side out.
inline std::string side_to_json(core::SeqSide side)
{
    return (side == core::Five_prime) ? "5'" : "3'";
}

/**
 * \brief The one segment, and end, an event acts on.
 *
 * "targets" is a list because an event of the new model may act on several (REVIEW §7.1). The
 * legacy events act on exactly one, so any other count is refused.
 */
inline std::pair<core::SeqTypeString, core::SeqSide> single_target(const nlohmann::json &node)
{
    const nlohmann::json &targets = require(node, "targets");
    if (!targets.is_array() || targets.size() != 1)
        throw std::runtime_error("event json: \"targets\" must list exactly one target, as a legacy event "
                                 "acts on one segment");
    const nlohmann::json &target = targets.front();
    reject_unknown_keys(target, { "segment", "side" });
    const auto segment = require(target, "segment").get<core::SeqTypeString>();
    if (segment.empty())
        throw std::runtime_error("event json: a target's \"segment\" is empty");

    const auto side = target.find("side");
    if (side == target.end())
        return { segment, core::Undefined_side };
    if (*side == "5'")
        return { segment, core::Five_prime };
    if (*side == "3'")
        return { segment, core::Three_prime };
    throw std::runtime_error("event json: side " + side->dump() + " is neither \"5'\" nor \"3'\"");
}

/// The label of a realization whose identity is a number: a deletion or an insertion.
inline int integer_label(const nlohmann::json &realization)
{
    const nlohmann::json &label = require(realization, "label");
    if (!label.is_number_integer())
        throw std::runtime_error("event json: label " + label.dump() + " is not an integer");
    return label.get<int>();
}

/// The label of a realization whose identity is a name: a gene or a nucleotide.
inline std::string text_label(const nlohmann::json &realization)
{
    const nlohmann::json &label = require(realization, "label");
    if (!label.is_string() || label.get_ref<const std::string &>().empty())
        throw std::runtime_error("event json: label " + label.dump() + " is not a non-empty string");
    return label.get<std::string>();
}

/**
 * \brief The node's realizations, ordered by index.
 *
 * Enforces the schema's two rules: labels are unique, since a label is a realization's identity;
 * and indices are exactly 0 to n-1. The second is not decoration. Every concrete
 * add_realization() assigns the index itself, from the current realization count, so a document
 * with a gap or a duplicate would be silently renumbered and would no longer match the marginal
 * array it was written for.
 */
inline std::vector<const nlohmann::json *> realizations_in_index_order(
        const nlohmann::json &node, std::initializer_list<const char *> allowed = kRealizationKeys)
{
    const nlohmann::json &array = require(node, "realizations");
    std::vector<const nlohmann::json *> out;
    out.reserve(array.size());
    std::set<std::string> labels;
    for (const auto &realization : array) {
        reject_unknown_keys(realization, allowed);
        if (!labels.insert(require(realization, "label").dump()).second)
            throw std::runtime_error("event json: two realizations are labelled " + realization.at("label").dump());
        out.push_back(&realization);
    }

    std::sort(out.begin(), out.end(),
              [](const nlohmann::json *a, const nlohmann::json *b) {
                  return require(*a, "index").get<int>() < require(*b, "index").get<int>();
              });

    for (std::size_t position = 0; position != out.size(); ++position) {
        const int index = require(*out[position], "index").get<int>();
        if (index != static_cast<int>(position))
            throw std::runtime_error("event json: realization indices must run from 0 without gaps, found "
                                     + std::to_string(index) + " where " + std::to_string(position)
                                     + " was expected");
    }
    return out;
}

}  // namespace igor::model::legacy::json_detail
