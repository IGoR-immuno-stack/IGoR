/*
 * GenerationState.h
 *
 *  The segments of one generated sequence, by SeqTypeId.
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

#include <igor/Core/Legacy/DynamicSequenceMap.h>
#include <igor/Core/Legacy/SeqTypeRegistry.h>

#include <string>
#include <utility>
#include <vector>

namespace igor::model::legacy {
using namespace igor::core::legacy;

/**
 * \class GenerationState GenerationState.h
 * \brief What the legacy generator builds one sequence in: a segment per seq_type.
 *
 * The "generation state" an event's `apply` writes into (ARCHITECTURE_SYNTHESIS §5; here
 * Rec_Event::construct_realization()). See docs/GENERATION_REWRITE_PLAN.md §3.
 *
 * The segments live in the container inference keeps its constructed sequences in, a
 * DynamicSequenceMap, so that the two sides already agree on what "not built" means when they
 * come to share a ScenarioContext layer. Generation needs one layer, since it never backtracks,
 * and it stores pointers into buffers it owns, as inference stores pointers into the events'
 * buffers: the values stay small and trivially copyable, and a modifier edits in place.
 *
 * Three states, as in inference: **not built** (no event has created the segment: it does not
 * exist, and reading it throws std::out_of_range), **built and empty** (an insertion of length
 * zero, a template deleted away), **built**.
 *
 * Neither copyable nor movable: the map points into the buffers.
 */
class GenerationState {
public:
    /// \param registry frozen; must outlive this state
    explicit GenerationState(const SeqTypeRegistry &registry)
        : segments_(registry), buffers_(registry.total_count())
    {
        //One layer per segment, claimed up front and written by whichever event creates it.
        for (std::size_t id = 0; id != segments_.count(); ++id) {
            segments_.request_layer(id);
        }
    }

    GenerationState(const GenerationState &) = delete;
    GenerationState &operator=(const GenerationState &) = delete;

    /// Create the segment `id` with `content`, as an event that constructs it does. Replaces
    /// whatever was there.
    std::string &create(SeqTypeId id, std::string content)
    {
        std::string &buffer = buffers_.at(id);
        buffer = std::move(content);
        segments_.set(id, &buffer, 0);
        return buffer;
    }

    /// A segment some event has created, to be edited in place.
    /// \throws std::out_of_range if no event has created it
    std::string &modify(SeqTypeId id) { return *segments_.get(id); }

    /// \throws std::out_of_range if no event has created it
    const std::string &read(SeqTypeId id) const { return *segments_.get(id); }

    bool built(SeqTypeId id) const { return segments_.exists(id); }

    /// The segments in the registry's 5'->3' order, skipping those no event built.
    std::string assemble() const
    {
        std::string sequence;
        for (const SeqTypeId id : segments_.registry().ordering()) {
            if (built(id)) {
                sequence += read(id);
            }
        }
        return sequence;
    }

    const SeqTypeRegistry &registry() const { return segments_.registry(); }

private:
    DynamicSequenceMap<std::string *> segments_;
    std::vector<std::string> buffers_;
};

} // namespace igor::model::legacy
