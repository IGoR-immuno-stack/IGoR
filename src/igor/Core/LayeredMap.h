/*
 * LayeredMap.h
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
 * \file LayeredMap.h
 * \brief Layered, integer-keyed containers backing the scenario state maps.
 *
 * A "memory layer" is one level of the scenario traversal: writing a value at layer n
 * leaves the value at layer n-1 intact, so backtracking is a decrement rather than a copy.
 * Layers are tracked per key independently.
 *
 * The templates here are key-type agnostic; the concrete instantiations used by the
 * scenario contexts are declared as typedefs in Utils.h.
 *
 * These are the legacy containers. They are being replaced by LayeredArray/DynamicSequenceMap
 * (see docs/REC_EVENT_CAPABILITY_REFACTORING_PLAN.md, decision D5) and will be removed once
 * their last consumer is ported.
 */

#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <string>

/*
 * This class provides a fast alternative to unordered_map<Seq_type,string*> for the constructed_sequences objects
 * Change this and give some kind of matrix with memory levels
 * Create a 0 size at first?
 * Get rid of it in the deletions
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
        /*for(size_t i = 0 ; i!=range ; i++){
				str_ptr_arr[i] = nullptr;
			}*/
    }
    virtual ~Enum_fast_memory_map()
    {
        delete[] value_ptr_arr;
        delete[] memory_layer_ptr;
    }

    //Accessors
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
                throw std::out_of_range(
                        "Trying to access uninitialized position in Enum_fast_memory_map::at(const K& key)");
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
                throw std::out_of_range("Trying to access uninitialized position in Enum_fast_memory_map::at( const K& "
                                        "key, int memory_layer)");
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
                throw std::out_of_range("Trying to access uninitialized position in Enum_fast_memory_map::at( const K& "
                                        "key, int memory_layer)");
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
        /*std::cout<<key<<std::endl;
			std::cout<<memory_layer_ptr[key]<<std::endl;*/
        if (key > range - 1) {
            throw std::out_of_range("Unknown key in Enum_fast_memory_map::request_memory_layer()");
        }
        //Get current memory layer at this position
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

    //Setters
    void set_value(const K &key, const V &value, int memory_layer)
    {
        assert(key <= range - 1);
        //Cannot fill memory layer without filling the ones downstream. A key whose
        //layers have never been written (memory_layer_ptr[key] == -1) may therefore
        //only be written at layer 0. Writing it at a higher layer means the caller is
        //using a layer index that belongs to some other map's sequence, which this
        //assertion exists to catch -- do not relax it.
        assert(memory_layer <= (memory_layer_ptr[key] + 1));
        //Grow the buffer if this layer lies beyond what is currently allocated.
        //request_memory_layer() grows max_layer as it goes, so a map driven through it
        //never gets here. A map driven purely by set_value() walking 0,1,2,... does:
        //the constructor allocates only layer 0, and set_value() has no other growth
        //path (the ones in at()/operator[] are not on this route). Without this, such
        //a write lands past the end of value_ptr_arr.
        if (memory_layer >= max_layer) {
            int new_max_layer = memory_layer + 1;
            V *new_value_ptr = new V[range * (new_max_layer + 1)];
            for (size_t i = 0; i != range; ++i) {
                for (size_t j = 0; j != static_cast<size_t>(max_layer); ++j) {
                    (*(new_value_ptr + i + j * range)) = (*(value_ptr_arr + i + j * range));
                }
            }
            delete[] value_ptr_arr;
            value_ptr_arr = new_value_ptr;
            max_layer = new_max_layer;
        }
        (*(value_ptr_arr + key + memory_layer * range)) = value;
        //Setting a value at a given layer invalidate upper layers
        memory_layer_ptr[key] = memory_layer;
    }

    void multiply_all(double &prod_operand, int *memory_adresses)
    {
        for (size_t i = 0; i != range; ++i) {
            /*				std::cout<<i<<std::endl;
				std::cout<<(*(value_ptr_arr + i))<<std::endl;
				std::cout<<(*(memory_adresses + i))*range<<std::endl;*/
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
                throw std::runtime_error("First memory layer already initialized for key " + std::to_string(i)
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
    size_t range; //= 6; //number of outcomes in Seq_type
};

/*
 * This class provides a fast alternative to unordered_map<Seq_type,string*> for the constructed_sequences objects
 * Change this and give some kind of matrix with memory levels
 * Create a 0 size at first?
 * Get rid of it in the deletions
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
        /*for(size_t i = 0 ; i!=range ; i++){
				str_ptr_arr[i] = nullptr;
			}*/
    }
    virtual ~Enum_fast_memory_dual_key_map()
    {
        delete[] value_ptr_arr;
        delete[] memory_layer_ptr;
    }

    //Accessors
    //Cannot use [] with more than one argument
    /*V& operator[](const K1& key1 , const K2& key2){
			if(key1>range_key1-1){throw std::out_of_range("Unknown key1 in Enum_fast_memory_map::operator(Seq_type)");}
			if(key2>range_key2-1){throw std::out_of_range("Unknown key2 in Enum_fast_memory_map::operator(Seq_type)");}
			if(memory_layer_ptr[key1+range_key1*key2]>-1){
				return (*(value_ptr_arr+ key1+range_key1*key2 + memory_layer_ptr[key1+range_key1*key2]*total_range));
			}
			else{
				memory_layer_ptr[key1+range_key1*key2] = 0;
				return (*(value_ptr_arr+key1+range_key1*key2));
			}

		}*/

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
                throw std::out_of_range(
                        "Trying to access uninitialized position in Enum_fast_memory_dual_key__map::at()");
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
                throw std::out_of_range(
                        "Trying to access uninitialized position in Enum_fast_memory_dual_key__map::at()");
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
                throw std::out_of_range(
                        "Trying to access uninitialized position in Enum_fast_memory_dual_key__map::at()");
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
                throw std::out_of_range(
                        "Trying to access uninitialized position in Enum_fast_memory_dual_key__map::at()");
            }
        }
    }

    int get_current_memory_layer(const K1 &key1, const K2 &key2) { return memory_layer_ptr[key1 + range_key1 * key2]; }

    bool exist(const K1 &key1, const K2 &key2) { return memory_layer_ptr[key1 + range_key1 * key2] > -1; }

    void request_memory_layer(const K1 &key1, const K2 &key2)
    {
        /*std::cout<<key<<std::endl;
			std::cout<<memory_layer_ptr[key]<<std::endl;*/
        if (key1 > range_key1 - 1) {
            throw std::out_of_range("Unknown key1 in Enum_fast_memory_dual_key__map::request_memory_layer()");
        }
        if (key2 > range_key2 - 1) {
            throw std::out_of_range("Unknown key2 in Enum_fast_memory_dual_key__map::request_memory_layer()");
        }
        //Get current memory layer at this position
        if (memory_layer_ptr[key1 + range_key1 * key2] < max_layer) {
            ++memory_layer_ptr[key1 + range_key1 * key2];
        } else {
            ++max_layer;
            V *new_value_ptr = new V[total_range * (max_layer + 1)];
            for (size_t i = 0; i != total_range; ++i) {
                for (size_t j = 0; j != (max_layer); ++j) {
                    (*(new_value_ptr + i + j * total_range)) = (*(value_ptr_arr + i + j * total_range));
                }
            }
            delete[] value_ptr_arr;
            value_ptr_arr = new_value_ptr;
            ++memory_layer_ptr[key1 + range_key1 * key2];
        }
    }

    //Setters
    void set_value(const K1 &key1, const K2 &key2, V value, int memory_layer)
    {
        if (key1 > range_key1 - 1) {
            throw std::out_of_range("Unknown key1 in Seq_type_str_p_map::set_value()");
        }
        if (key2 > range_key2 - 1) {
            throw std::out_of_range("Unknown key2 in Seq_type_str_p_map::set_value()");
        }
        //Cannot fill memory layer without filling the ones downstream
        if (memory_layer <= (memory_layer_ptr[key1 + range_key1 * key2] + 1)) {
            (*(value_ptr_arr + key1 + range_key1 * key2 + memory_layer * total_range)) = value;
            //Setting a value at a given layer invalidate upper layers
            memory_layer_ptr[key1 + range_key1 * key2] = memory_layer;
        } else {
            throw std::out_of_range(
                    "Trying to access incorrect memory layer in Enum_fast_memory_dual_key__map::set_value()");
        }
    }

protected:
    V *value_ptr_arr;
    int *memory_layer_ptr;
    int max_layer;
    size_t range_key1; //= 6; //number of outcomes in Seq_type
    size_t range_key2;
    size_t total_range;
};
