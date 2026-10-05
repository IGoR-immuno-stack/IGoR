/*
 * LayerClaim.h
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

/**
 * \brief Claim layer 0 for a key a fixture is about to write.
 *
 * Since R3b a write lands only on a layer its key has claimed, and layer 0 is no exception:
 * it is not a free base layer, it is simply the layer the first `request_layer()` grants
 * (see `LayeredArray`, decision O10). Production code claims in `initialize_event()`; a
 * fixture that writes a map directly is standing in for those events and has to stand in
 * for their claims too.
 *
 * Idempotent on purpose. Presetting the same key twice keeps the claim it already has
 * rather than walking the ownership mark up, which would move the baseline the layer and
 * ownership harnesses snapshot before `iterate()` and quietly widen what they accept.
 *
 * Variadic over the key so it works for every map shape the tests touch: one key for
 * `DynamicSequenceMap`, a `(SeqTypeId, Seq_side)` pair for `Seq_offsets_map`, a `SafetyCell`
 * for `SafetyMatrix`.
 */
template <typename Map, typename... Key>
void claim_layer_zero(Map &map, Key... key)
{
    if (map.claimed_layer(key...) < 0) {
        map.request_layer(key...);
    }
}
