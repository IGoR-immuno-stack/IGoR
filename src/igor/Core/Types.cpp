/*
 * Types.cpp
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

#include <igor/Core/Types.h>

#include <ostream>
#include <stdexcept>

// Moved verbatim from Legacy/Utils.cpp (step 1c): the textual forms of the vocabulary enums live
// in igor::core next to the enums, where argument-dependent lookup finds them.
using std::invalid_argument;
using std::ostream;
using std::string;

namespace igor::core {

ostream &operator<<(ostream &os, SeqSide ss)
{
    switch (ss) {
    case Five_prime:
        os << "Five_prime";
        break;
    case Three_prime:
        os << "Three_prime";
        break;
    case Undefined_side:
        os << "Undefined_side";
        break;

    default:
        throw invalid_argument("Unknown Seq_side in operator << ");
    }
    return os;
}

string operator+(const string &str, SeqSide ss)
{
    string next_str;
    switch (ss) {
    case Five_prime:
        next_str = "Five_prime";
        break;
    case Three_prime:
        next_str = "Three_prime";
        break;
    case Undefined_side:
        next_str = "Undefined_side";
        break;
    }
    return str + next_str;
}

string operator+(const string &str, EventType et)
{
    string next_str;
    switch (et) {
    case GeneChoice_t:
        next_str = "GeneChoice";
        break;
    case Deletion_t:
        next_str = "Deletion";
        break;
    case Insertion_t:
        next_str = "Insertion";
        break;
    case Dinuclmarkov_t:
        next_str = "DinucMarkov";
        break;
    }
    return str + next_str;
}

string to_string(const SeqType st)
{
    switch (st) {
    case V_gene_seq: return "V_gene_seq";
    case VD_ins_seq: return "VD_ins_seq";
    case D_gene_seq: return "D_gene_seq";
    case DJ_ins_seq: return "DJ_ins_seq";
    case J_gene_seq: return "J_gene_seq";
    case VJ_ins_seq: return "VJ_ins_seq";
    default: throw invalid_argument("Unknown Seq_type in to_string(Seq_type)");
    }
}

string to_string(const SeqSide ss)
{
    switch (ss) {
    case Five_prime:
        return "Five_prime";
    case Three_prime:
        return "Three_prime";
    case Undefined_side:
        return "Undefined_side";

    default:
        throw invalid_argument("Unknown Seq_side in to_string(const Seq_side).");
    }
}

} // namespace igor::core
