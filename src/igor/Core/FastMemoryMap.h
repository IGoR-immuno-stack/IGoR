/*
 * FastMemoryMap.h
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

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

/*
 * This class provides a fast alternative to unordered_map<Seq_type,string*> for
 * the constructed_sequences objects Change this and give some kind of matrix
 * with memory levels Create a 0 size at first? Get rid of it in the deletions
 */
template <typename K, typename V>
class Enum_fast_memory_map
{
public:
    Enum_fast_memory_map(int defined_range) : max_layer(0), range(defined_range)
    {
        value_ptr_arr = new V[range];
        memory_layer_ptr = new int[range];
        for (size_t i = 0; i != range; ++i) {
            this->memory_layer_ptr[i] = -1;
        }
    }
    virtual ~Enum_fast_memory_map()
    {
        delete[] value_ptr_arr;
        delete[] memory_layer_ptr;
    }

    // Accessors
    V &operator[](const K &key)
    {
        if (key > range - 1) {
            throw std::out_of_range("Unknown key in Enum_fast_memory_map::operator(Seq_type)");
        }
        if (memory_layer_ptr[key] > -1) {
            return (*(value_ptr_arr + key + memory_layer_ptr[key] * range));
        } else {
            memory_layer_ptr[key] = 0;
            return (*(value_ptr_arr + key));
        }
    }

    V &at(const K &key)
    {
        if (key > range - 1) {
            throw std::out_of_range("Unknown key in Enum_fast_memory_map::operator(K key)");
        } else {
            if (memory_layer_ptr[key] > -1) {
                return (*(value_ptr_arr + key + memory_layer_ptr[key] * range));
            } else {
                throw std::out_of_range("Trying to access uninitialized position in "
                                        "Enum_fast_memory_map::at(const K& key)");
            }
        }
    }

    const V &at(const K &key) const
    {
        if (key > range - 1) {
            throw std::out_of_range("Unknown key in Enum_fast_memory_map::operator(K key)");
        } else {
            if (memory_layer_ptr[key] > -1) {
                return (*(value_ptr_arr + key + memory_layer_ptr[key] * range));
            } else {
                throw std::out_of_range("Trying to access uninitialized position in "
                                        "Enum_fast_memory_map::at(const K& key)");
            }
        }
    }

    V &at(const K &key, int memory_layer)
    {
        if (key > range - 1) {
            throw std::out_of_range("Unknown seq type in Enum_fast_memory_map::operator(Seq_type)");
        } else {
            if (memory_layer <= (memory_layer_ptr[key] + 1)) {
                memory_layer_ptr[key] = memory_layer;
                return (*(value_ptr_arr + key + memory_layer * range));
            } else {
                throw std::out_of_range(
                        "Trying to access uninitialized position in "
                        "Enum_fast_memory_map::at( const K& key, int memory_layer)");
            }
        }
    }

    const V &at(const K &key, int memory_layer) const
    {
        if (key > range - 1) {
            throw std::out_of_range("Unknown seq type in Enum_fast_memory_map::operator(Seq_type)");
        } else {
            if (memory_layer <= (memory_layer_ptr[key] + 1)) {
                memory_layer_ptr[key] = memory_layer;
                return (*(value_ptr_arr + key + memory_layer * range));
            } else {
                throw std::out_of_range(
                        "Trying to access uninitialized position in "
                        "Enum_fast_memory_map::at( const K& key, int memory_layer)");
            }
        }
    }

    int get_current_memory_layer(const K &key) { return memory_layer_ptr[key]; }

    void get_all_current_memory_layer(int *memory_layers_recipient)
    {
        for (size_t i = 0; i != range; ++i) {
            memory_layers_recipient[i] = memory_layer_ptr[i];
        }
    }

    bool exist(const K &key) { return memory_layer_ptr[key] > -1; }

    void request_memory_layer(const K &key)
    {
        if (key > range - 1) {
            throw std::out_of_range("Unknown key in Enum_fast_memory_map::request_memory_layer()");
        }
        // Get current memory layer at this position
        if (memory_layer_ptr[key] < max_layer) {
            ++memory_layer_ptr[key];
        } else {
            ++max_layer;
            V *new_value_ptr = new V[range * (max_layer + 1)];
            for (size_t i = 0; i != range; ++i) {
                for (size_t j = 0; j != (max_layer); ++j) {
                    (*(new_value_ptr + i + j * range)) = (*(value_ptr_arr + i + j * range));
                }
            }
            delete[] value_ptr_arr;
            value_ptr_arr = new_value_ptr;
            ++memory_layer_ptr[key];
        }
    }

    // Setters
    void set_value(const K &key, const V &value, int memory_layer)
    {
        if (key > range - 1) {
            throw std::out_of_range("Unknown seq type in Seq_type_str_p_map::operator(Seq_type)");
        }
        // Cannot fill memory layer without filling the ones downstream
        if (memory_layer <= (memory_layer_ptr[key] + 1)) {
            (*(value_ptr_arr + key + memory_layer * range)) = value;
            // Setting a value at a given layer invalidate upper layers
            memory_layer_ptr[key] = memory_layer;
        } else {
            throw std::out_of_range("Trying to access incorrect memory layer in "
                                    "Enum_fast_memory_map::set_value()");
        }
    }

    void multiply_all(double &prod_operand, int *memory_adresses)
    {
        for (size_t i = 0; i != range; ++i) {
            prod_operand *= value_ptr_arr[i + memory_adresses[i] * range];
        }
    }

    void reset()
    {
        for (size_t i = 0; i != range; ++i) {
            if (memory_layer_ptr[i] > -1) {
                memory_layer_ptr[i] = 0;
            }
        }
    }

    void init_first_layer(V value)
    {
        for (size_t i = 0; i != range; ++i) {
            if (memory_layer_ptr[i] > -1) {
                throw std::runtime_error("First memory layer already initialized for key "
                                         + std::to_string(i)
                                         + " in Enum_fast_memory_map::init_first_layer");
            } else {
                value_ptr_arr[i] = value;
                memory_layer_ptr[i] = 0;
            }
        }
    }

protected:
    V *value_ptr_arr;
    int *memory_layer_ptr;
    int max_layer;
    size_t range;
};

/*
 * This class provides a fast alternative to unordered_map<Seq_type,string*> for
 * the constructed_sequences objects Change this and give some kind of matrix
 * with memory levels Create a 0 size at first? Get rid of it in the deletions
 */
template <typename K1, typename K2, typename V>
class Enum_fast_memory_dual_key_map
{
public:
    Enum_fast_memory_dual_key_map(size_t Key1_range, size_t Key2_range)
        : max_layer(0), range_key1(Key1_range), range_key2(Key2_range)
    {
        total_range = range_key1 * range_key2;
        value_ptr_arr = new V[total_range];
        memory_layer_ptr = new int[total_range];
        for (size_t i = 0; i != total_range; ++i) {
            this->memory_layer_ptr[i] = -1;
        }
    }
    virtual ~Enum_fast_memory_dual_key_map()
    {
        delete[] value_ptr_arr;
        delete[] memory_layer_ptr;
    }

    V &at(const K1 &key1, const K2 &key2)
    {
        if (key1 > range_key1 - 1) {
            throw std::out_of_range("Unknown key1 in Enum_fast_memory_dual_key_map::at()");
        } else if (key2 > range_key2 - 1) {
            throw std::out_of_range("Unknown key2 in Enum_fast_memory_dual_key__map::at()");
        } else {
            if (memory_layer_ptr[key1 + range_key1 * key2] > -1) {
                return (*(value_ptr_arr + key1 + range_key1 * key2
                          + memory_layer_ptr[key1 + range_key1 * key2] * total_range));
            } else {
                throw std::out_of_range("Trying to access uninitialized position in "
                                        "Enum_fast_memory_dual_key__map::at()");
            }
        }
    }

    const V &at(const K1 &key1, const K2 &key2) const
    {
        if (key1 > range_key1 - 1) {
            throw std::out_of_range("Unknown key1 in Enum_fast_memory_dual_key_map::at()");
        } else if (key2 > range_key2 - 1) {
            throw std::out_of_range("Unknown key2 in Enum_fast_memory_dual_key__map::at()");
        } else {
            if (memory_layer_ptr[key1 + range_key1 * key2] > -1) {
                return (*(value_ptr_arr + key1 + range_key1 * key2
                          + memory_layer_ptr[key1 + range_key1 * key2] * total_range));
            } else {
                throw std::out_of_range("Trying to access uninitialized position in "
                                        "Enum_fast_memory_dual_key__map::at()");
            }
        }
    }

    V &at(const K1 &key1, const K2 &key2, int memory_layer)
    {
        if (key1 > range_key1 - 1) {
            throw std::out_of_range("Unknown key1 in Enum_fast_memory_dual_key__map::at()");
        } else if (key2 > range_key2 - 1) {
            throw std::out_of_range("Unknown key2 in Enum_fast_memory_dual_key__map::at()");
        } else {
            if (memory_layer <= (memory_layer_ptr[key1 + range_key1 * key2] + 1)) {
                memory_layer_ptr[key1 + range_key1 * key2] = memory_layer;
                return (*(value_ptr_arr + key1 + range_key1 * key2 + memory_layer * total_range));
            } else {
                throw std::out_of_range("Trying to access uninitialized position in "
                                        "Enum_fast_memory_dual_key__map::at()");
            }
        }
    }

    const V &at(const K1 &key1, const K2 &key2, int memory_layer) const
    {
        if (key1 > range_key1 - 1) {
            throw std::out_of_range("Unknown key1 in Enum_fast_memory_dual_key__map::at()");
        } else if (key2 > range_key2 - 1) {
            throw std::out_of_range("Unknown key2 in Enum_fast_memory_dual_key__map::at()");
        } else {
            if (memory_layer <= (memory_layer_ptr[key1 + range_key1 * key2] + 1)) {
                memory_layer_ptr[key1 + range_key1 * key2] = memory_layer;
                return (*(value_ptr_arr + key1 + range_key1 * key2 + memory_layer * total_range));
            } else {
                throw std::out_of_range("Trying to access uninitialized position in "
                                        "Enum_fast_memory_dual_key__map::at()");
            }
        }
    }

    int get_current_memory_layer(const K1 &key1, const K2 &key2)
    {
        return memory_layer_ptr[key1 + range_key1 * key2];
    }

    void request_memory_layer(const K1 &key1, const K2 &key2)
    {
        if (key1 > range_key1 - 1) {
            throw std::out_of_range("Unknown key1 in "
                                    "Enum_fast_memory_dual_key__map::request_memory_layer()");
        }
        if (key2 > range_key2 - 1) {
            throw std::out_of_range("Unknown key2 in "
                                    "Enum_fast_memory_dual_key__map::request_memory_layer()");
        }
        // Get current memory layer at this position
        if (memory_layer_ptr[key1 + range_key1 * key2] < max_layer) {
            ++memory_layer_ptr[key1 + range_key1 * key2];
        } else {
            ++max_layer;
            V *new_value_ptr = new V[total_range * (max_layer + 1)];
            for (size_t i = 0; i != total_range; ++i) {
                for (size_t j = 0; j != (max_layer); ++j) {
                    (*(new_value_ptr + i + j * total_range)) =
                            (*(value_ptr_arr + i + j * total_range));
                }
            }
            delete[] value_ptr_arr;
            value_ptr_arr = new_value_ptr;
            ++memory_layer_ptr[key1 + range_key1 * key2];
        }
    }

    // Setters
    void set_value(const K1 &key1, const K2 &key2, V value, int memory_layer)
    {
        if (key1 > range_key1 - 1) {
            throw std::out_of_range("Unknown key1 in Seq_type_str_p_map::set_value()");
        }
        if (key2 > range_key2 - 1) {
            throw std::out_of_range("Unknown key2 in Seq_type_str_p_map::set_value()");
        }
        // Cannot fill memory layer without filling the ones downstream
        if (memory_layer <= (memory_layer_ptr[key1 + range_key1 * key2] + 1)) {
            (*(value_ptr_arr + key1 + range_key1 * key2 + memory_layer * total_range)) = value;
            // Setting a value at a given layer invalidate upper layers
            memory_layer_ptr[key1 + range_key1 * key2] = memory_layer;
        } else {
            throw std::out_of_range("Trying to access incorrect memory layer in "
                                    "Enum_fast_memory_dual_key__map::set_value()");
        }
    }

protected:
    V *value_ptr_arr;
    int *memory_layer_ptr;
    int max_layer;
    size_t range_key1;
    size_t range_key2;
    size_t total_range;
};
