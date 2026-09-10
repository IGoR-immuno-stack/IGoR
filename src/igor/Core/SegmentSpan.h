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
#include <igor/Core/StdTypedefs.h>

#include <stdexcept>
#include <string>

/**
 * \brief One end of one segment: the addressing unit for everything side-taking.
 *
 * Already the de facto key in about twenty signatures -- Seq_offsets_map's eight accessors,
 * Rec_Event::get_offset_role() and get_offset_delta_bounds() on all four subclasses, and
 * PendingModifierBounds, which flattens the pair by hand into `end_index(id, side)`. Naming it
 * is what lets a span address a segment's own extent rather than only the gap between two
 * segments (decision O11).
 *
 * ### The coordinate convention, stated once
 *
 * A boundary is a **half-open range bound**, not a nucleotide index: cut_position() turns it
 * into the `begin` or `end` you would slice the read with. So a span's two boundaries are
 * exactly a `[begin, end)` pair, its length is `end - begin`, and two adjacent spans compose
 * because `[a, b)` and `[b, c)` make `[a, c)` -- the ordinary property, not a new algebra.
 *
 * This is why boundaries are not simply offsets. `Seq_Offset` is a *closed*-interval index,
 * and closed intervals do not compose by subtraction: taking a segment as `[5', 3']` and the
 * gap after it as `(3', 5')` mixes two conventions, so every measurement needs its own +/-1
 * and the two cannot be added. Sections 7.1 and 7.8 are both off-by-one bugs in this area,
 * which is why the conversion lives in one function rather than at each call site. Section
 * 2.5 of docs/ITERATE_GENERIC_REWRITE_PLAN.md has the derivation.
 */
struct SegmentBoundary {
    SeqTypeId id = kNoSeqType;
    Seq_side side = Undefined_side;

    friend bool operator==(SegmentBoundary, SegmentBoundary) = default;
};

/**
 * \brief Turn a segment end's offset into a half-open range bound -- the index you slice with.
 *
 * `Seq_Offset` names a nucleotide: `Five_prime` is the segment's first, `Three_prime` its last.
 * Slicing wants `[begin, end)`, so the conversion is side-dependent:
 *
 *     5' -> the offset itself   (the first nucleotide is included: it is `begin`)
 *     3' -> the offset plus one (`end` is one past the last nucleotide)
 *
 * Use it in pairs. For a span, `cut_position()` of the left boundary is `begin` and of the
 * right boundary is `end`, so the region is `str.substr(begin, end - begin)` for a
 * std::string, or `{seq.begin() + begin, seq.begin() + end}` for an Int_Str, and the span's
 * length is `end - begin`.
 *
 * With V at read positions 0..9 and D starting at 15:
 *
 *     the V segment : [cut(5',0), cut(3',9)) == [0, 10)  -> 10 nt
 *     the VD gap    : [cut(3',9), cut(5',15)) == [10, 15) -> 5 nt, positions 10..14
 *     V and the gap : [0, 15)                             -> 15 nt, with no correction term
 *
 * \throws std::invalid_argument for Undefined_side, which names no end of anything.
 */
inline Seq_Offset cut_position(Seq_side side, Seq_Offset segment_end_offset)
{
    switch (side) {
    case Five_prime:
        return segment_end_offset;
    case Three_prime:
        return segment_end_offset + 1;
    default:
        throw std::invalid_argument("cut_position: Undefined_side names no cut");
    }
}

/**
 * \brief An ordered range of the registry ordering, delimited by two boundaries.
 *
 * The addressing unit for junction-length bounds (G5), for the effect queries
 * Rec_Event::affects_length_of() / affects_proba_of(), and -- once it lands -- for Phase D's
 * cluster boundaries. One object, named once, in all three.
 *
 * ### What a span's length is, and is not
 *
 * A span does **not** have a length. Until every event contributing to it has chosen a
 * realization it has a *profile*: the set of achievable lengths with a best probability each,
 * which is what the `map<int,double>` bound objects are and why their composition operator is
 * a convolution. The single number is the **observed distance** -- what the scenario has
 * already committed to at the consumer's position -- and the bound is one used as a key into
 * the other.
 *
 * That also means a span's profile is relative to *where it is read*: a boundary is a slot
 * that a `Creates` event writes and a `Modifies` event moves, so the same span has different
 * profiles at different points in the priority ordering. See section 2.5 of
 * docs/ITERATE_GENERIC_REWRITE_PLAN.md.
 *
 * ### Constructing one
 *
 * Use gap(), which is the only shape any current caller needs. The general two-boundary form
 * is representable so that Phase D and 5b need no type migration, but it has no consumer yet
 * and the query semantics for it are deliberately not implemented (O11, step 3).
 */
struct SegmentSpan {
    SegmentBoundary left{};
    SegmentBoundary right{};

    friend bool operator==(SegmentSpan, SegmentSpan) = default;

    /**
     * The span *between* two segments, excluding both: from \a l's 3' end to \a r's 5' end.
     * Today's only meaning, and what every caller wants -- an event that trims one of those
     * two boundaries widens the span rather than shortening the anchor.
     */
    static constexpr SegmentSpan gap(SeqTypeId l, SeqTypeId r)
    {
        return SegmentSpan{SegmentBoundary{l, Three_prime}, SegmentBoundary{r, Five_prime}};
    }
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
        return SegmentSpan::gap(static_cast<SeqTypeId>(V_gene_seq), static_cast<SeqTypeId>(D_gene_seq));
    case DJ_ins_seq:
        return SegmentSpan::gap(static_cast<SeqTypeId>(D_gene_seq), static_cast<SeqTypeId>(J_gene_seq));
    case VJ_ins_seq:
        return SegmentSpan::gap(static_cast<SeqTypeId>(V_gene_seq), static_cast<SeqTypeId>(J_gene_seq));
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
    throw std::invalid_argument("legacy_junction_of: span (" + std::to_string(span.left.id) + ","
                                + std::to_string(span.right.id) + ") is not a legacy junction");
}
