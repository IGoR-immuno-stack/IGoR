/*
 * CoreEnums.h
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
 * \file CoreEnums.h
 * \brief The small dependency-free enums describing events and sequence segments.
 *
 * Split out of Utils.h so the container headers can name them without including Utils.h,
 * which now includes those containers in turn. Nothing here depends on any other IGoR
 * header, which is what makes the split possible; the Gene_class enums stay in Utils.h for
 * now because they carry export annotations and a cluster of conversion functions.
 */

enum Event_type { GeneChoice_t, Deletion_t, Insertion_t, Dinuclmarkov_t, Undefined_t };

enum Event_safety { VD_safe = 0, DJ_safe = 1, VJ_safe = 2 };

/// Which end of a constructed sequence segment an offset or a deletion refers to.
enum Seq_side { Five_prime = 0, Three_prime = 1, Undefined_side = 2 };

/**
 * \brief The six sequence segment kinds of a standard VDJ/VJ model.
 *
 * Superseded at runtime by SeqTypeId, which a SeqTypeRegistry allocates and which can
 * describe any topology. These values are pinned to the matching SeqTypeIds by
 * Model_Parms::read_model_parms(), so enum-keyed code still addresses the right slots
 * while it is migrated.
 */
enum Seq_type { V_gene_seq = 0, VD_ins_seq = 1, D_gene_seq = 2, DJ_ins_seq = 3, J_gene_seq = 4, VJ_ins_seq = 5 };
