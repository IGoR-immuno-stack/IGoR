/*
 * DynamicSequenceMap.h
 *
 *  Created on: Nov 26, 2025
 *      Author: IGoR Team
 *
 *  This source code is distributed as part of the IGoR software.
 *  IGoR (Inference and Generation of Repertoires) is a versatile software to
 * analyze and model immune receptors generation, selection, mutation and all
 * other processes. Copyright (C) 2017  Quentin Marcou
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
 *
 */

#pragma once

#include <igor/Core/FastMemoryMap.h>
#include <igor/Core/SequenceTypes.h>

/**
 * \class DynamicSequenceMap
 * \brief A wrapper around Enum_fast_memory_map that automatically sizes itself
 * based on the SequenceTypeRegistry.
 *
 * This class is designed to replace Enum_fast_memory_map<Seq_type, V> to
 * support dynamic sequence types (like tandem D genes) while maintaining the
 * performance characteristics of the original implementation. It uses 'int' as
 * the key type to accommodate extended TypeIds.
 */
template <typename V>
class DynamicSequenceMap : public Enum_fast_memory_map<int, V>
{
public:
    /**
   * \brief Default constructor.
   * Initializes the map with the size from the global SequenceTypeRegistry.
   */
    DynamicSequenceMap() : Enum_fast_memory_map<int, V>(get_sequence_type_registry().size()) { }

    /**
   * \brief Constructor with explicit size.
   * \param size The size of the map (number of keys).
   */
    DynamicSequenceMap(int size) : Enum_fast_memory_map<int, V>(size) { }

    virtual ~DynamicSequenceMap() = default;
};
