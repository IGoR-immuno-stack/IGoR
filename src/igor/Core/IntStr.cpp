/*
 * IntStr.cpp
 *
 *  Created on: Jul 21, 2016
 *      Author: Quentin Marcou
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

 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include <igor/Core/IntStr.h>


namespace igor::core {

using namespace std;

IntStr &IntStr::operator+=(const IntStr &other)
{
    this->insert(this->end(), other.begin(), other.end());
    return *this;
}

IntStr &IntStr::operator+=(const int &a)
{
    this->push_back(a);
    return *this;
}

IntStr &IntStr::operator+=(int &&a)
{
    this->push_back(a);
    return *this;
}

IntStr &IntStr::append(const IntStr &other)
{
    return (*this) += other;
}

IntStr &IntStr::append(const int &a)
{
    return (*this) += a;
}

IntStr IntStr::operator+(const IntStr &other) const
{
    IntStr new_int_str;
    new_int_str.reserve(this->size() + other.size());
    new_int_str += (*this);
    new_int_str += other;
    return new_int_str;
}

void IntStr::substr(IntStr &new_int_str, size_t pos /*= 0*/, size_t len /*= 99999999999999*/) const
{
    IntStr::const_iterator begin = this->begin() + pos;
    IntStr::const_iterator last;
    if (len == npos) {
        last = this->end();
    } else if (pos + len >= this->size()) {
        last = this->end();
    } else {
        last = begin + len;
    }
    new_int_str.assign(begin, last);
}

IntStr IntStr::substr(size_t pos /*= 0*/, size_t len /*= 99999999999999*/) const
{
    IntStr::const_iterator begin = this->begin() + pos;
    IntStr::const_iterator last;
    if (len == npos) {
        last = this->end();
    } else if (pos + len >= this->size()) {
        last = this->end();
    } else {
        last = begin + len;
    }
    IntStr new_int_str;
    new_int_str.assign(begin, last);
    return new_int_str;
}

IntStr &IntStr::erase(size_t pos, size_t len)
{
    IntStr::iterator begin = this->begin() + pos;
    IntStr::iterator last;
    if (len == npos) {
        last = this->end();
    } else if (pos + len >= this->size()) {
        last = this->end();
    } else {
        last = begin + len;
    }
    this->erase(begin, last);
    return *this;
}

ostream &operator<<(ostream &out, const IntStr &int_str)
{
    for (IntStr::const_iterator iter = int_str.begin(); iter != int_str.end(); ++iter) {
        out << (*iter);
    }
    return out;
}

/*
IntStr::IntStr(): int_vector() {
}

IntStr::IntStr(const IntStr& other){
	this->int_vector = other.int_vector;
}

IntStr::~IntStr() {
	// TODO Auto-generated destructor stub
}

IntStr& IntStr::operator=(const IntStr& other){
	this->int_vector = other.int_vector;
	return this;
}

size_t IntStr::max_size() const noexcept{
	return this->int_vector.max_size();
}

size_t IntStr::capacity() const noexcept{
	return this->int_vector.capacity();
}

void IntStr::clear(){

}*/

} // namespace igor::core
