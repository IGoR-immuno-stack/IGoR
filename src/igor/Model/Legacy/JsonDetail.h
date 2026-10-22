#pragma once

/**
 * \file JsonDetail.h
 * \brief Helpers shared by the events' JSON constructors.
 *
 * Internal header: it includes the full json.hpp, so it is included by .cpp files only and is
 * deliberately absent from Core's installed FILE_SET. An installed header that included this
 * one would put 25k lines of json.hpp in front of every consumer.
 *
 * A constructor first checks its node against its type's schema (ModelSchema.h), which settles
 * the node's shape: its keys, their types, the side's spelling, one target. What is left here is
 * what a schema cannot say: that indices run from 0 without gaps, and that labels are unique.
 */

#include <igor/Core/Types.h>
#include <igor/Model/ModelSchema.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace igor::model::legacy::json_detail {

/// \a node, once checked against the schema of the event type \a type: what a constructor reads.
inline const nlohmann::json &validated(const nlohmann::json &node, const std::string &type)
{
    igor::model::schema::validate(node, "event-" + type + ".schema.json", "the " + type + " node");
    return node;
}

/// "5'" and "3'", as the model document spells an end. An event that acts on no particular end
/// leaves the side out.
inline std::string side_to_json(core::SeqSide side)
{
    return (side == core::Five_prime) ? "5'" : "3'";
}

/// The one segment, and end, a validated node's event acts on: the legacy events act on one.
inline std::pair<core::SeqTypeString, core::SeqSide> single_target(const nlohmann::json &node)
{
    const nlohmann::json &target = node.at("targets").front();
    const auto side = target.find("side");
    const core::SeqSide parsed = (side == target.end()) ? core::Undefined_side
                                 : (*side == "5'")      ? core::Five_prime
                                                        : core::Three_prime;
    return { target.at("segment").get<core::SeqTypeString>(), parsed };
}

/// single_target() of \a node once validated: for a constructor's base-class arguments, whose order
/// of evaluation is unspecified, so that each validates before it reads.
inline std::pair<core::SeqTypeString, core::SeqSide> validated_target(const nlohmann::json &node,
                                                                       const std::string &type)
{
    return single_target(validated(node, type));
}

/**
 * \brief A validated node's realizations, ordered by index.
 *
 * Enforces the two rules a schema cannot: labels are unique, since a label is a realization's
 * identity; and indices are exactly 0 to n-1. The second is not decoration. Every concrete
 * add_realization() assigns the index itself, from the current realization count, so a document
 * with a gap or a duplicate would be silently renumbered and would no longer match the marginal
 * array it was written for.
 */
inline std::vector<const nlohmann::json *> realizations_in_index_order(const nlohmann::json &node)
{
    std::vector<const nlohmann::json *> out;
    std::set<std::string> labels;
    for (const auto &realization : node.at("realizations")) {
        if (!labels.insert(realization.at("label").dump()).second)
            throw std::runtime_error("event json: two realizations are labelled " + realization.at("label").dump());
        out.push_back(&realization);
    }

    std::sort(out.begin(), out.end(), [](const nlohmann::json *a, const nlohmann::json *b) {
        return a->at("index").get<int>() < b->at("index").get<int>();
    });

    for (std::size_t position = 0; position != out.size(); ++position) {
        const int index = out[position]->at("index").get<int>();
        if (index != static_cast<int>(position))
            throw std::runtime_error("event json: realization indices must run from 0 without gaps, found "
                                     + std::to_string(index) + " where " + std::to_string(position)
                                     + " was expected");
    }
    return out;
}

}  // namespace igor::model::legacy::json_detail
