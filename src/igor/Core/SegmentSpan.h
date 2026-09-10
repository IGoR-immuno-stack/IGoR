/*
 * SegmentSpan.h
 *
 *  This source code is distributed as part of the IGoR software.
 *  IGoR (Inference and Generation of Repertoires) is a versatile software to analyze and model immune receptors
 *  generation, selection, mutation and all other processes.
 *   Copyright (C) 2017  Quentin Marcou
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <igor/Core/CoreEnums.h>
#include <igor/Core/SeqTypeRegistry.h>

#include <stdexcept>
#include <string>

/**
 * \brief An ordered, anchor-exclusive range of the registry ordering.
 *
 * The addressing unit for junction-length bounds (G5), for the effect queries
 * Rec_Event::affects_length_of() / affects_proba_of(), and -- once it lands -- for Phase D's
 * cluster boundaries. One object, named once, in all three.
 *
 * \a left and \a right are the *anchor* segments bounding the span; neither belongs to it.
 * A span's length is the read-space distance between \a left's and \a right's **as-created**
 * facing boundaries, so an event that trims one of those boundaries (a Deletion) widens the
 * span rather than shortening the anchor. See section 2.5 of
 * docs/ITERATE_GENERIC_REWRITE_PLAN.md for the frame and the composition algebra.
 */
struct SegmentSpan {
    SeqTypeId left = kNoSeqType;
    SeqTypeId right = kNoSeqType;

    friend bool operator==(SegmentSpan, SegmentSpan) = default;
};

/**
 * \brief The span a legacy junction seq_type names.
 *
 * The Len_proba machinery addresses junctions by the Seq_type enum, where one value carries
 * two meanings: VJ_ins_seq is a *segment* in a VJ model but the whole V->J *span* in a VDJ
 * one. This function resolves that ambiguity in the only direction the legacy code needs --
 * the argument is always the span reading -- and is the single place the legacy VDJ topology
 * is hardcoded. It goes away with S4c, when the length maps become span-keyed.
 *
 * \throws std::invalid_argument for a seq_type that names no junction.
 */
inline SegmentSpan legacy_span_of(Seq_type junction)
{
    switch (junction) {
    case VD_ins_seq:
        return {static_cast<SeqTypeId>(V_gene_seq), static_cast<SeqTypeId>(D_gene_seq)};
    case DJ_ins_seq:
        return {static_cast<SeqTypeId>(D_gene_seq), static_cast<SeqTypeId>(J_gene_seq)};
    case VJ_ins_seq:
        return {static_cast<SeqTypeId>(V_gene_seq), static_cast<SeqTypeId>(J_gene_seq)};
    default:
        throw std::invalid_argument("legacy_span_of: seq_type " + std::to_string(junction)
                                    + " names no junction span");
    }
}

/**
 * \brief Inverse of legacy_span_of(), for predicates whose tables are still enum-keyed.
 *
 * Every current affects_length_of() override reduces to its pre-S4a has_effect_on() table,
 * which is written over the enum. Rather than generalise those tables -- a semantic change,
 * and S4b's job once the traversal carries the registry ordering -- S4a maps back here so the
 * answers are unchanged by construction.
 *
 * \throws std::invalid_argument for a span no legacy junction names.
 */
inline Seq_type legacy_junction_of(SegmentSpan span)
{
    if (span == legacy_span_of(VD_ins_seq)) { return VD_ins_seq; }
    if (span == legacy_span_of(DJ_ins_seq)) { return DJ_ins_seq; }
    if (span == legacy_span_of(VJ_ins_seq)) { return VJ_ins_seq; }
    throw std::invalid_argument("legacy_junction_of: span (" + std::to_string(span.left) + ","
                                + std::to_string(span.right) + ") is not a legacy junction");
}
