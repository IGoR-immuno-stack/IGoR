/*
 * Errorscounter.h
 */

#pragma once

#include <igor/Core/Counter.h>
#include <map>
#include <igor/Core/Export.h>

class CORE_EXPORT Errors_counter : public Counter
{
public:
    Errors_counter();
    Errors_counter(size_t);
    Errors_counter(size_t, std::string); // for demo output path
    virtual ~Errors_counter();

    std::string type() const override { return "ErrorsCounter"; };

    void count_scenario(long double, double, const std::string &, Seq_type_str_p_map &, const Seq_offsets_map &,
                        const std::unordered_map<std::tuple<Event_type, int, Seq_side>, std::shared_ptr<Rec_Event>> &,
                        Mismatch_vectors_map &) override;

    void initialize_counter(const Model_Parms &, const Model_marginals &) override;

    void add_checked(std::shared_ptr<Counter>) override;

    std::shared_ptr<Counter> copy() const override;

private:
    std::map<int, int> memory_layer_mismatches_map;
    size_t counter_size;
    long double *updated_counter_pointer;
};
