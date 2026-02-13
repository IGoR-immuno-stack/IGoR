/*
 * Errorscounter.cpp
 */

#include <igor/Core/Errorscounter.h>
#include <igor/Core/SequenceTypes.h>

using namespace std;

Errors_counter::Errors_counter() : Counter()
{
    counter_size = 100;
    updated_counter_pointer = new long double[counter_size]();
}

Errors_counter::Errors_counter(size_t size) : Counter()
{
    counter_size = size;
    updated_counter_pointer = new long double[counter_size]();
}

Errors_counter::Errors_counter(size_t size, std::string) : Counter()
{
    counter_size = size;
    updated_counter_pointer = new long double[counter_size]();
}

Errors_counter::~Errors_counter()
{
    delete[] updated_counter_pointer;
}

void Errors_counter::count_scenario(
        long double scenario_error_w_proba, double scenario_proba, const std::string &sequence,
        Seq_type_str_p_map &constructed_sequences, const Seq_offsets_map &seq_offsets,
        const std::unordered_map<std::tuple<Event_type, int, Seq_side>, std::shared_ptr<Rec_Event>> &events_map,
        Mismatch_vectors_map &mismatches_lists)
{
    int n_errors = 0;
    auto &registry = get_sequence_type_registry();
    for (size_t i = 0; i < registry.size(); ++i) {
        auto type_info = registry.get_type_info((SequenceTypeRegistry::TypeId)i);
        if (type_info.is_gene) {
            // We need to know which memory layer to use.
            // For now, let's use the current one if it exists.
            int layer = mismatches_lists.get_current_memory_layer(i);
            if (layer != -1) {
                const std::vector<int> *mismatches = mismatches_lists.at(i, layer);
                if (mismatches) {
                    n_errors += mismatches->size();
                }
            }
        }
    }

    if (n_errors < (int)this->counter_size) {
        this->updated_counter_pointer[n_errors] += scenario_error_w_proba;
    }
}

void Errors_counter::initialize_counter(const Model_Parms &model_parms, const Model_marginals &model_marginals) { }

void Errors_counter::add_checked(std::shared_ptr<Counter> other)
{
    auto other_ptr = dynamic_pointer_cast<Errors_counter>(other);
    if (other_ptr) {
        for (size_t i = 0; i < counter_size && i < other_ptr->counter_size; ++i) {
            this->updated_counter_pointer[i] += other_ptr->updated_counter_pointer[i];
        }
    }
}

shared_ptr<Counter> Errors_counter::copy() const
{
    auto new_ptr = make_shared<Errors_counter>(this->counter_size);
    return new_ptr;
}
