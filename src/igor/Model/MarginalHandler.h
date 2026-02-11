#pragma once

#include <igor/Math/Tensor.h>

#include <string>
#include <ostream>
#include <istream>
#include <vector>

namespace igor::model {

// ─── Event Descriptor ──────────────────────────────────────────────────

/// Describes a recombination event for handler construction.
/// Built from Model_Parms + Rec_Event at runtime.
struct EventDescriptor {
    std::string name;               // Rec_Event_name (e.g. "v_choice")
    int type;                       // Event_type enum value (0=GeneChoice_t, 1=Insertion_t, 2=Deletion_t, 3=Dinuclmarkov_t)
    int gene_class;                 // Gene_class enum value (e.g., V_gene, D_gene, J_gene, VD_genes, DJ_genes)
    int side;                       // Seq_side enum value (Five_prime, Three_prime, or Undefined_side)
    std::vector<std::size_t> shape; // Tensor dimensions from Rec_Event::size() [n_realizations, parent1, parent2, ...]
};

// ─── Handler Interface ─────────────────────────────────────────────────

template <typename T = double>
class MarginalHandler {
public:
    using scalar_type = T;

    virtual ~MarginalHandler() = default;

    // Identity
    const std::string& name() const { return name_; }

    // Tensor access
    virtual const math::Tensor<T>& parameters() const = 0;
    virtual math::Tensor<T>& parameters() = 0;
    virtual const math::Tensor<T>& accumulator() const = 0;
    virtual math::Tensor<T>& accumulator() = 0;

    // EM operations
    virtual void reset_accumulator() = 0;
    virtual void maximize_likelihood() = 0;

    // I/O (backward compatible with legacy txt format)
    virtual void write_parameters(std::ostream& out) const = 0;
    virtual void read_parameters(std::istream& in) = 0;

protected:
    explicit MarginalHandler(std::string name) : name_(std::move(name)) {}
    std::string name_;
};

} // namespace igor::model

