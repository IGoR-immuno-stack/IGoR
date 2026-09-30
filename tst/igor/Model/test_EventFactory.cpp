/*
 * test_EventFactory.cpp
 *
 * Tests for EventFactory registration and two-phase initialization.
 *
 * Two-phase init pattern:
 *   Phase 1 – EventFactory::create(type)  → shared_ptr<Rec_Event>  (blank object)
 *   Phase 2 – configure via Rec_Event base-class setters:
 *              set_event_class(), set_event_side(), set_nickname(), set_priority()
 *              (+ derived-class extras where absolutely needed, e.g. Gene_choice templates)
 *
 * The Rec_Event base-class API is preferred throughout.  Downcasts to derived
 * types appear only where a derived-specific setter/getter is strictly required,
 * and in the direct-construction comparison targets.
 */

#include <catch2/catch_test_macros.hpp>
#include <igor/Core/Dinuclmarkov.h>
#include <igor/Core/Deletion.h>
#include <igor/Core/Genechoice.h>
#include <igor/Core/Insertion.h>
#include <igor/Model/EventFactory.h>

// ─── Registration & basic creation ──────────────────────────────────────────

TEST_CASE("EventFactory - Registration and Creation", "[factory]")
{
    using namespace igor::model::event_factory;

    SECTION("all concrete types are registered")
    {
        REQUIRE(is_registered(GeneChoice_t));
        REQUIRE(is_registered(Deletion_t));
        REQUIRE(is_registered(Insertion_t));
        REQUIRE(is_registered(Dinuclmarkov_t));
        REQUIRE_FALSE(is_registered(Undefined_t));
    }

    SECTION("create returns non-null with correct type - GeneChoice")
    {
        auto ev = create(GeneChoice_t);
        REQUIRE(ev != nullptr);
        REQUIRE(ev->get_type() == GeneChoice_t);
        REQUIRE(std::dynamic_pointer_cast<Gene_choice>(ev) != nullptr);
    }

    SECTION("create returns non-null with correct type - Deletion")
    {
        auto ev = create(Deletion_t);
        REQUIRE(ev != nullptr);
        REQUIRE(ev->get_type() == Deletion_t);
        REQUIRE(std::dynamic_pointer_cast<Deletion>(ev) != nullptr);
    }

    SECTION("create returns non-null with correct type - Insertion")
    {
        auto ev = create(Insertion_t);
        REQUIRE(ev != nullptr);
        REQUIRE(ev->get_type() == Insertion_t);
        REQUIRE(std::dynamic_pointer_cast<Insertion>(ev) != nullptr);
    }

    SECTION("create returns non-null with correct type - Dinuclmarkov")
    {
        auto ev = create(Dinuclmarkov_t);
        REQUIRE(ev != nullptr);
        REQUIRE(ev->get_type() == Dinuclmarkov_t);
        REQUIRE(std::dynamic_pointer_cast<Dinucl_markov>(ev) != nullptr);
    }

    SECTION("create throws for unregistered type")
    {
        REQUIRE_THROWS_AS(create(Undefined_t), std::runtime_error);
    }
}

// ─── Phase-1 blank state ─────────────────────────────────────────────────────
//
// A freshly-created event (Phase 1 only) must match the default-constructed
// instance of the same class in every base-class-visible property.

TEST_CASE("EventFactory - Phase 1 blank state matches default constructor", "[factory][2phase]")
{
    using namespace igor::model::event_factory;

    SECTION("Gene_choice blank state")
    {
        auto ev = create(GeneChoice_t);
        Gene_choice direct;

        CHECK(ev->get_type()  == direct.get_type());
        CHECK(ev->get_class() == direct.get_class());
        CHECK(ev->get_side()  == direct.get_side());
        CHECK(ev->size()      == direct.size());
    }

    SECTION("Deletion blank state")
    {
        auto ev = create(Deletion_t);
        Deletion direct;

        CHECK(ev->get_type()  == direct.get_type());
        CHECK(ev->get_class() == direct.get_class());
        CHECK(ev->get_side()  == direct.get_side());
        CHECK(ev->size()      == direct.size());
    }

    SECTION("Insertion blank state")
    {
        auto ev = create(Insertion_t);
        Insertion direct;

        CHECK(ev->get_type()  == direct.get_type());
        CHECK(ev->get_class() == direct.get_class());
        CHECK(ev->get_side()  == direct.get_side());
        CHECK(ev->size()      == direct.size());
    }

    SECTION("Dinucl_markov blank state")
    {
        auto ev = create(Dinuclmarkov_t);
        Dinucl_markov direct;

        CHECK(ev->get_type()  == direct.get_type());
        CHECK(ev->get_class() == direct.get_class());
        CHECK(ev->get_side()  == direct.get_side());
    }
}

// ─── Two-phase init via Rec_Event API ────────────────────────────────────────
//
// For each realistic VDJ event configuration:
//   1. Create via factory (Phase 1)
//   2. Configure via Rec_Event base-class API (Phase 2)
//   3. Verify base-class getters
//   4. Compare against a directly-constructed equivalent (sanity check)

TEST_CASE("EventFactory - Two-Phase Init via Rec_Event API", "[factory][2phase][api]")
{
    using namespace igor::model::event_factory;

    // Shared helper: configure common properties through the base class only.
    auto configure = [](std::shared_ptr<Rec_Event> ev,
                        Gene_class gene_class, Seq_side side,
                        const std::string& nick, int prio)
    {
        ev->set_event_class(gene_class);
        ev->set_event_side(side);
        ev->set_nickname(nick);
        ev->set_priority(prio);
    };

    // Shared helper: assert all base-class getters.
    auto check = [](const std::shared_ptr<Rec_Event>& ev,
                    Event_type type, Gene_class gene_class, Seq_side side,
                    const std::string& nick, int prio)
    {
        CHECK(ev->get_type()     == type);
        CHECK(ev->get_class()    == gene_class);
        CHECK(ev->get_side()     == side);
        CHECK(ev->get_nickname() == nick);
        CHECK(ev->get_priority() == prio);
    };

    // Shared helper: base-class comparison against a directly-constructed reference.
    auto matches = [](const std::shared_ptr<Rec_Event>& ev, const Rec_Event& ref)
    {
        CHECK(ev->get_type()     == ref.get_type());
        CHECK(ev->get_class()    == ref.get_class());
        CHECK(ev->get_side()     == ref.get_side());
        CHECK(ev->get_nickname() == ref.get_nickname());
        CHECK(ev->get_priority() == ref.get_priority());
    };

    // ── Gene_choice ─────────────────────────────────────────────────────────

    SECTION("Gene_choice - V gene / Undefined_side / prio 7")
    {
        auto ev = create(GeneChoice_t);
        configure(ev, V_gene, Undefined_side, "v_choice", 7);
        check(ev, GeneChoice_t, V_gene, Undefined_side, "v_choice", 7);

        // Derived-class extra: load genomic templates (no base equivalent)
        auto gc = std::static_pointer_cast<Gene_choice>(ev);
        gc->set_genomic_templates({{"TRBV1", "ATGC"}, {"TRBV2", "GCTA"}});
        REQUIRE(gc->get_realizations_map().size() == 2);
        REQUIRE(gc->get_realizations_map().at("TRBV1").value_str == "ATGC");

        // Reference comparison (base API only)
        Gene_choice ref(V_gene);
        ref.set_nickname("v_choice");
        ref.set_priority(7);
        matches(ev, ref);
    }

    SECTION("Gene_choice - D gene / Undefined_side / prio 6")
    {
        auto ev = create(GeneChoice_t);
        configure(ev, D_gene, Undefined_side, "d_gene", 6);
        check(ev, GeneChoice_t, D_gene, Undefined_side, "d_gene", 6);

        auto gc = std::static_pointer_cast<Gene_choice>(ev);
        gc->set_genomic_templates({{"TRBD1", "AGGT"}, {"TRBD2", "TGGG"}});
        REQUIRE(gc->get_realizations_map().size() == 2);

        Gene_choice ref(D_gene);
        ref.set_nickname("d_gene");
        ref.set_priority(6);
        matches(ev, ref);
    }

    SECTION("Gene_choice - J gene / Undefined_side / prio 7")
    {
        auto ev = create(GeneChoice_t);
        configure(ev, J_gene, Undefined_side, "j_choice", 7);
        check(ev, GeneChoice_t, J_gene, Undefined_side, "j_choice", 7);

        auto gc = std::static_pointer_cast<Gene_choice>(ev);
        gc->set_genomic_templates({{"TRBJ1-1", "TGG"}, {"TRBJ1-2", "GGG"}, {"TRBJ2-1", "AGG"}});
        REQUIRE(gc->get_realizations_map().size() == 3);
        REQUIRE(gc->get_realizations_map().at("TRBJ1-1").value_str == "TGG");

        Gene_choice ref(J_gene);
        ref.set_nickname("j_choice");
        ref.set_priority(7);
        matches(ev, ref);
    }

    // ── Deletion – four standard configurations ──────────────────────────────

    SECTION("Deletion - V gene / Three_prime / prio 5")
    {
        auto ev = create(Deletion_t);
        configure(ev, V_gene, Three_prime, "v_3_del", 5);
        check(ev, Deletion_t, V_gene, Three_prime, "v_3_del", 5);

        Deletion ref(V_gene, Three_prime);
        ref.set_nickname("v_3_del");
        ref.set_priority(5);
        matches(ev, ref);
    }

    SECTION("Deletion - D gene / Three_prime / prio 5")
    {
        auto ev = create(Deletion_t);
        configure(ev, D_gene, Three_prime, "d_3_del", 5);
        check(ev, Deletion_t, D_gene, Three_prime, "d_3_del", 5);

        Deletion ref(D_gene, Three_prime);
        ref.set_nickname("d_3_del");
        ref.set_priority(5);
        matches(ev, ref);
    }

    SECTION("Deletion - D gene / Five_prime / prio 5")
    {
        auto ev = create(Deletion_t);
        configure(ev, D_gene, Five_prime, "d_5_del", 5);
        check(ev, Deletion_t, D_gene, Five_prime, "d_5_del", 5);

        Deletion ref(D_gene, Five_prime);
        ref.set_nickname("d_5_del");
        ref.set_priority(5);
        matches(ev, ref);
    }

    SECTION("Deletion - J gene / Five_prime / prio 5")
    {
        auto ev = create(Deletion_t);
        configure(ev, J_gene, Five_prime, "j_5_del", 5);
        check(ev, Deletion_t, J_gene, Five_prime, "j_5_del", 5);

        Deletion ref(J_gene, Five_prime);
        ref.set_nickname("j_5_del");
        ref.set_priority(5);
        matches(ev, ref);
    }

    // ── Insertion – VD and DJ junctions ─────────────────────────────────────

    SECTION("Insertion - VD genes / Undefined_side / prio 4")
    {
        auto ev = create(Insertion_t);
        configure(ev, VD_genes, Undefined_side, "vd_ins", 4);
        check(ev, Insertion_t, VD_genes, Undefined_side, "vd_ins", 4);

        Insertion ref(VD_genes);
        ref.set_nickname("vd_ins");
        ref.set_priority(4);
        matches(ev, ref);
    }

    SECTION("Insertion - DJ genes / Undefined_side / prio 2")
    {
        auto ev = create(Insertion_t);
        configure(ev, DJ_genes, Undefined_side, "dj_ins", 2);
        check(ev, Insertion_t, DJ_genes, Undefined_side, "dj_ins", 2);

        Insertion ref(DJ_genes);
        ref.set_nickname("dj_ins");
        ref.set_priority(2);
        matches(ev, ref);
    }

    // ── Dinucl_markov – VD and DJ junctions ─────────────────────────────────

    SECTION("Dinucl_markov - VD genes / Undefined_side / prio 3")
    {
        auto ev = create(Dinuclmarkov_t);
        configure(ev, VD_genes, Undefined_side, "vd_dinucl", 3);
        check(ev, Dinuclmarkov_t, VD_genes, Undefined_side, "vd_dinucl", 3);

        Dinucl_markov ref(VD_genes);
        ref.set_nickname("vd_dinucl");
        ref.set_priority(3);
        matches(ev, ref);
    }

    SECTION("Dinucl_markov - DJ genes / Undefined_side / prio 1")
    {
        auto ev = create(Dinuclmarkov_t);
        configure(ev, DJ_genes, Undefined_side, "dj_dinucl", 1);
        check(ev, Dinuclmarkov_t, DJ_genes, Undefined_side, "dj_dinucl", 1);

        Dinucl_markov ref(DJ_genes);
        ref.set_nickname("dj_dinucl");
        ref.set_priority(1);
        matches(ev, ref);
    }
}
