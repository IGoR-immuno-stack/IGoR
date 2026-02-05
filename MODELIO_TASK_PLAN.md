# ModelIO Task Plan

## GitHub Issue
- **Issue:** [#32 - Update model parms and marginals file format](https://github.com/IGoR-immuno-stack/IGoR/issues/32)
- **Labels:** enhancement, refactoring
- **Branch:** `feature/model-io`

## Overview

Modernize IGoR's model storage from the legacy banner format to JSON, separating immutable model structure from inferred weights and improving editability, comment support, and interoperability.

## Goals

1. Use JSON as the modern format (with JSON Schema validation)
2. Split into `model_struct.json` (immutable structure) and `model_weights.json` (learned probabilities)
3. Use `event_nickname` as the unique, stable event identifier
4. Provide legacy ↔ JSON converters for OLGA/ALICE compatibility
5. Follow namespace + free function pattern (like AIRR modules)
6. **Event-centric parsing:** Introduce `Rec_Event` factory methods that accept structured objects
7. Delegate event parsing to each event type for maintainability

## Non-Goals

- Full round-trip compatibility for advanced features (tandem-D) back to legacy format
- Breaking changes to existing Model_Parms/Model_Marginals classes (we wrap, not replace)
- Immediate replacement of existing legacy I/O (gradual migration path)

---

## Architecture

### Two-Layer Design

The ModelIO module provides a **bridge architecture** between modern JSON format and existing Core classes:

1. **Layer 1: Standalone I/O (Current Phase)**
   - `igor::modelio` namespace with free functions
   - Independent JSON ↔ data structure conversion
   - No dependency on Core event classes
   - Useful for: CLI tools, format conversion, validation

2. **Layer 2: Event Factory Integration (Future Phase)**
   - Factory methods in `Rec_Event` base class: `from_json()`, `to_json()`
   - Delegate to concrete event types: `Gene_choice`, `Deletion`, `Insertion`, `Dinucl_markov`
   - Bidirectional conversion: `Rec_Event* ↔ EventData`
   - Integration with `Model_Parms` and `Model_Marginals`

### Current Implementation (Layer 1)

```cpp
namespace igor::modelio {
    // Common types
    struct EventData { ... };
    struct EdgeData { ... };
    struct ModelData { ... };
    struct WeightsData { ... };

    namespace detail {
        // Private helper functions (following sparrow pattern)
        std::string generate_nickname(const EventData& event);
        std::string generate_legacy_name(const EventData& event);
    }

    namespace structure {
        ModelData read_json(const std::string& filepath);
        void write_json(const std::string& filepath, const ModelData& model, bool pretty = false);
        ModelData read_legacy(const std::string& filepath);
        void write_legacy(const std::string& filepath, const ModelData& model);
    }

    namespace weights {
        WeightsData read_json(const std::string& filepath);
        void write_json(const std::string& filepath, const WeightsData& weights, bool pretty = false);
        WeightsData read_legacy(const std::string& filepath);
        void write_legacy(const std::string& filepath, const WeightsData& weights);
    }
}
```

### Future Integration (Layer 2)

```cpp
// In Rec_Event.h (Core module)
class Rec_Event {
public:
    // Factory method: construct event from JSON/structured data
    static std::shared_ptr<Rec_Event> from_json(const igor::modelio::EventData& data);

    // Serialize event to structured data
    virtual igor::modelio::EventData to_json() const = 0;

    // Existing methods...
    virtual void write2txt(std::ofstream&, std::queue<int>&) = 0;
};

// In Gene_choice.h
class Gene_choice : public Rec_Event {
public:
    static std::shared_ptr<Gene_choice> from_json(const igor::modelio::EventData& data);
    igor::modelio::EventData to_json() const override;
    // ...
};
```

### File Layout

```
src/igor/ModelIO/
├── CMakeLists.txt
├── ModelIOCommon.h      # Common types (EventData, ModelData, WeightsData)
├── ModelReader.h/.cpp   # structure::read_*, weights::read_* (JSON + legacy)
├── ModelWriter.h/.cpp   # structure::write_*, weights::write_* (JSON + legacy)
└── README.md

app/igor-model/
├── CMakeLists.txt
└── igor-model.cpp       # CLI tool for model conversion

tst/igor/ModelIO/
├── CMakeLists.txt
├── test_model_structure.cpp  # Tests for structure I/O
└── test_data/
    ├── TRB_model_struct.json
    ├── TRB_model_weights.json
    ├── TRB_model_parms.txt     # Legacy for testing
    └── TRB_model_marginals.txt # Legacy for testing
```
---

## JSON Schema Design

### model_struct.json

```json
{
  "format_version": "1.0.0",
  "metadata": {
    "species": "Homo sapiens",
    "chain": "TRB",
    "model_type": "VDJ",
    "created": "2026-02-05T12:00:00Z",
    "description": "Human TRB recombination model"
  },
  "events": [
    {
      "id": "v_choice",
      "legacy_name": "GeneChoice_V_gene_Undefined_side_prio7_size89",
      "type": "GeneChoice",
      "gene_class": "V",
      "seq_side": "Undefined_side",
      "priority": 7,
      "realizations": [
        {
          "index": 0,
          "name": "TRBV1*01",
          "sequence": "GATACTGGA..."
        }
      ]
    }
  ],
  "edges": [
    { "parent": "v_choice", "child": "v_3_del" },
    { "parent": "j_choice", "child": "d_gene" }
  ],
  "sequence_types": {
    "order": ["V_gene", "VD_ins", "D_gene", "DJ_ins", "J_gene"],
    "definitions": {
      "V_gene": { "id": 0 },
      "VD_ins": { "id": 1, "parents": ["V_gene"] },
      "D_gene": { "id": 2, "parents": ["VD_ins"] },
      "DJ_ins": { "id": 3, "parents": ["D_gene"] },
      "J_gene": { "id": 4, "parents": ["DJ_ins"] }
    }
  }
}
```

### model_weights.json

```json
{
  "format_version": "1.0.0",
  "struct_version": "1.0.0",
  "events": [
    {
      "event_id": "v_choice",
      "dims": [89],
      "conditioning_events": [],
      "values": [0.011236, 0.011236, ...],
      "normalized": true
    },
    {
      "event_id": "d_gene",
      "dims": [15, 3],
      "conditioning_events": ["j_choice"],
      "values": [0.333, 0.333, 0.333, ...],
      "normalized": true
    }
  ],
  "error_rate": {
    "type": "SingleErrorRate",
    "rate": 0.001,
    "learn_on": "seq",
    "apply_on": "seq"
  }
}
```

---

## Task Breakdown

### Phase 0: Dependencies and Setup
**Estimated:** 0.5 days | **Status:** ✅ Complete

- [x] **Task 0.1:** Add nlohmann/json to pixi.toml
- [x] **Task 0.2:** Create ModelIO module directory structure
- [x] **Task 0.3:** Create CMakeLists.txt for ModelIO library
- [x] **Task 0.4:** Create empty header/source files with license headers

**Acceptance Criteria:**
- ✅ nlohmann/json available in build
- ✅ ModelIO module builds
- ✅ Tests can link to ModelIO

**Implementation Notes:**
- Namespace: `igor::modelio` (consistency with module name)
- Using `detail` namespace for private helpers (following sparrow pattern)
- Helper functions defined before use (no forward declarations needed)

---

### Phase 1: Common Types and Factory | **Status:** 🟡 In Progress
**Estimated:** 1.5 days

- [x] **Task 1.1:** Define `EventData` struct
- [x] **Task 1.2:** Define `EdgeData` struct
- [x] **Task 1.3:** Define `ModelData` struct
- [x] **Task 1.4:** Define `WeightsData` struct
- [x] **Task 1.5:** Implement JSON Schema validation helper
- [x] **Task 1.6:** Write unit tests for common types
- [ ] **Task 1.7:** Create `EventFactory` to turn `EventData` into `Rec_Event` objects

**Acceptance Criteria:**
- ✅ All structs compile and are copyable/movable
- ✅ JSON serialization works for simple types
- ✅ 8 test cases pass (26 assertions)
- [ ] Factory correctly instantiates `Gene_choice`, `Deletion`, `Insertion`, and `Dinucl_markov`

**Implementation Notes:**
- All types defined in `ModelIOCommon.h`
- Factory will live in `EventFactory.h/.cpp`
- Test file: `test_model_structure.cpp` covers basic types


---

### Phase 2: Model Structure I/O
**Estimated:** 3 days

#### Task 2.1: JSON Reading
- [ ] Implement `structure::read_json()`
- [ ] Parse events array with type-specific realization handling
- [ ] Parse edges array
- [ ] Parse sequence_types
- [ ] Validate format_version
- [ ] Handle missing optional fields gracefully

#### Task 2.2: JSON Writing
- [ ] Implement `structure::write_json()`
- [ ] Pretty-print with 2-space indentation
- [ ] Write metadata section
- [ ] Write events with realizations
- [ ] Write edges
- [ ] Write sequence_types

#### Task 2.3: Legacy Reading
- [ ] Implement `structure::read_legacy()`
- [ ] Parse @Event_list section
- [ ] Parse #GeneChoice, #Deletion, #Insertion, #DinucMarkov
- [ ] Parse @Edges section
- [ ] Parse @ErrorRate section (extract type only, not weights)
- [ ] Map legacy names to nicknames

#### Task 2.4: Legacy Writing
- [ ] Implement `structure::write_legacy()`
- [ ] Generate @Event_list section
- [ ] Generate event realizations
- [ ] Generate @Edges section
- [ ] Generate @ErrorRate section (type only)

#### Task 2.5: Integration with Model_Parms
- [ ] Add `to_model_data()` method or free function
- [ ] Add `from_model_data()` method or free function
- [ ] Preserve backward compatibility

**Acceptance Criteria:**
- Round-trip: JSON → ModelData → JSON preserves all data
- Round-trip: Legacy → ModelData → Legacy preserves all data
- Cross-format: Legacy → JSON → Legacy round-trips
- 15+ test cases pass

---

### Phase 3: Model We | **Status:** ✅ Complete (Core functionality)

#### Task 2.1: JSON Reading ✅ Complete
- [x] Implement `structure::read_json()`
- [x] Parse events array with type-specific realization handling
- [x] Parse edges array
- [x] Parse sequence_types
- [x] Validate format_version
- [x] Handle missing optional fields gracefully

#### Task 2.2: JSON Writing ✅ Complete
- [x] Implement `structure::write_json()`
- [x] Pretty-print with 2-space indentation
- [x] Write metadata section
- [x] Write events with realizations
- [x] Write edges
- [x] Write sequence_types

#### Task 2.3: Legacy Reading ✅ Complete
- [x] Implement `structure::read_legacy()`
- [x] Parse @Event_list section
- [x] Parse #GeneChoice, #Deletion, #Insertion, #DinucMarkov
- [x] Parse @Edges section
- [x] Parse @ErrorRate section (extract type only, not weights)
- [x] Map legacy names to nicknames

#### Task 2.4: Legacy Writing ✅ Complete
- [x] Implement `structure::write_legacy()`
- [x] Generate @Event_list section
- [x] Generate event realizations
- [x] Generate @Edges section
- [x] Generate @ErrorRate section (type only)

#### Task 2.5: Integration with Model_Parms ⬜ Not Started
- [ ] Add `to_model_data()` method or free function
- [ ] Add `from_model_data()` method or free function
- [ ] Preserve backward compatibility

**Acceptance Criteria:**
- ✅ Round-trip: JSON → ModelData → JSON preserves all data
- ✅ Round-trip: Legacy → ModelData → Legacy preserves all data (basic)
- ⬜ Cross-format: Legacy → JSON → Legacy round-trips (needs real data testing)
- ✅ 8+ test cases pass

**Implementation Notes:**
- `detail::generate_nickname()` creates short IDs from event properties
- `detail::generate_legacy_name()` reconstructs original format names
- Error handling with | **Status:** ✅ Complete (Core functionality)

#### Task 3.1: JSON Reading ✅ Complete
- [x] Implement `weights::read_json()`
- [x] Parse event weights with dims and values
- [x] Parse error_rate section
- [x] Validate tensor sizes match dims

#### Task 3.2: JSON Writing ✅ Complete
- [x] Implement `weights::write_json()`
- [x] Write event weights in order
- [x] Write error_rate with all parameters
- [x] Handle different error rate types

#### Task 3.3: Legacy Reading ✅ Complete
- [x] Implement `weights::read_legacy()`
- [x] Parse @nickname sections
- [x] Parse $Dim[...] dimensions
- [x] Parse #[conditioning] headers
- [x] Parse %probability values
- [x] Handle multi-dimensional tensors

#### Task 3.4: Legacy Writing ✅ Complete
- [x] Implement `weights::write_legacy()`
- [x] Generate @nickname headers
- [x] Generate $Dim[...] lines
- [x] Generate conditioning headers
- [x] Generate probability rows

#### Task 3.5: Integration with Model_Marginals ⬜ Not Started
- [ ] Add conversion functions
- [ ] Preserve normalization

**Acceptance Criteria:**
- ✅ Round-trip: JSON → WeightsData → JSON preserves all data (basic)
- ✅ Round-trip: Legacy → WeightsData → Legacy preserves all data (basic)
- ⬜ Tensor dimensions validated (needs comprehensive testing)
- ⬜ 15+ test cases pass (currently 8, need weights-specific tests)

**Implementation Notes:**
- Legacy format parsing handles comma-separated values
- Multi-dimensional tensor support via flattened arrays
- Error rate types: SingleErrorRate, HypermutationGlobal, HypermutationFullNmer
- Conditioning events tracked per event
- Need to add dedicated test file: `test_model_weights.cpp`

---

### Phase 4: Core Integration (Rec_Event ↔ EventData)
**Estimated:** 3-5 days | **Status:** ⬜ Not Started | **Priority:** HIGH

**Goal:** Integrate ModelIO with existing Core event classes for bidirectional conversion

**✅ Design Pattern Chosen:** **Option 7 (Core Factory)** - Domain-driven design

**Architecture:** Core provides `igor::event::Factory`; ModelIO readers call it to populate Core structures directly.

#### Task 4.1: Implement Core Factory
- [ ] Create `src/igor/Core/Event/Factory.h` and `.cpp`
- [ ] Define `igor::event::Factory` class
- [ ] Implement `create(const EventData& data)` method (type dispatch based on data.type)
- [ ] Implement type-specific creation methods:
  - `create_gene_choice(const EventData&)`
  - `create_deletion(const EventData&)`
  - `create_insertion(const EventData&)`
  - `create_dinucl_markov(const EventData&)`
- [ ] Implement `extract(const Rec_Event& event)` for serialization (Rec_Event → EventData)

#### Task 4.2: Update Core Event Classes (if needed)
- [ ] Add public accessors to `Gene_choice` (if needed for `extract()`)
- [ ] Add public accessors to `Deletion` (if needed for `extract()`)
- [ ] Add public accessors to `Insertion` (if needed for `extract()`)
- [ ] Add public accessors to `Dinucl_markov` (if needed for `extract()`)
- [ ] **Goal:** Minimal changes; only add getters, not serialization logic

#### Task 4.3: ModelIO Integration Functions
- [ ] Create `src/igor/ModelIO/CoreIntegration.h` and `.cpp`
- [ ] Implement `Model_Parms model_parms_from_json(const std::string& filepath)`
  - Read JSON → EventData
  - Use `igor::event::Factory` to create Rec_Event objects
  - Populate Model_Parms with events
  - Reconstruct edges from JSON
- [ ] Implement `void model_parms_to_json(const Model_Parms&, const std::string&, bool pretty)`
  - Extract EventData from each Rec_Event using `Factory::extract()`
  - Write to JSON
- [ ] Handle sequence type mapping (coordinate with #15)

#### Task 4.4: Model_Marginals Integration
- [ ] Implement `Model_Marginals model_marginals_from_json(const std::string& filepath)`
  - Read JSON → WeightsData (already implemented)
  - Populate Model_Marginals from WeightsData
- [ ] Implement `void model_marginals_to_json(const Model_Marginals&, const std::string&, bool pretty)`
  - Extract WeightsData from Model_Marginals
  - Write to JSON (already implemented)
- [ ] Handle tensor dimension validation
- [ ] Preserve normalization

#### Task 4.5: Optional Convenience Methods in Core
- [ ] Add `Model_Parms::save_json(const std::string& path, bool pretty = false)`
  - Delegates to `modelio::model_parms_to_json(*this, path, pretty)`
- [ ] Add static `Model_Parms::load_json(const std::string& path)`
  - Delegates to `modelio::model_parms_from_json(path)`
- [ ] Add `Model_Marginals::save_json(const std::string& path, bool pretty = false)`
- [ ] Add static `Model_Marginals::load_json(const std::string& path)`

#### Task 4.6: Comprehensive Testing
- [ ] Unit tests for `igor::event::Factory` creation methods
- [ ] Unit tests for `Factory::extract()` method
- [ ] Test Core → ModelIO → Core round-trips
- [ ] Test with real TRB/TRA models from `demo/`
- [ ] Verify inference results are identical (JSON vs legacy)
- [ ] Test with tandem-D models (if supported)
- [ ] Performance benchmarking vs legacy I/O
- [ ] Test edge cases (empty models, missing fields, etc.)

**Acceptance Criteria:**
- All Core event types can convert to/from EventData via Factory
- Model_Parms can load from/save to JSON
- Model_Marginals can load from/save to JSON
- Core controls object creation (domain-centric design)
- Existing Core tests still pass (backward compatibility)
- 20+ new integration tests pass
- No performance regression vs legacy format

**Implementation Notes:**
- **Factory lives in Core:** `src/igor/Core/Event/Factory.{h,cpp}`
- **Dependency direction:** ModelIO → Core (ModelIO depends on Core)
- **Core depends on EventData:** Must `#include <igor/ModelIO/ModelIOCommon.h>`
- **Simple and direct:** No registry, no abstract factory complexity
- **Type dispatch:** Factory uses `if/else` or `std::map<std::string, creator_fn>` on `event.type`
- **Extract method:** Uses RTTI (`dynamic_cast`) or visitor pattern to determine event type
- **Coordinate with #15:** Sequence type handling
- **Coordinate with #16:** Tensor/mdspan handling in marginals

**Trade-offs Accepted:**
- ❌ Core now depends on ModelIO headers (EventData)
- ❌ ModelIO cannot be used standalone without Core
- ✅ Simple, straightforward implementation
- ✅ Core controls its object creation (encapsulation)
- ✅ Natural dependency direction for domain-driven design

**Dependencies:**
- Requires Phase 2 and Phase 3 complete ✅
- Requires Core module changes (Factory class)
- May need minimal Core event accessors
- Coordinate with #15 (sequence types) and #16 (marginals modernization)

---

### Phase 5: CLI Tool (igor-model)
**Estimated:** 1 day | **Status:** ⬜ Not Started

- [ ] **Task 5.1:** Create `app/igor-model/igor-model.cpp`
- [ ] **Task 5.2:** Implement command structure:
  ```
  igor-model convert <input> <output>
  igor-model validate <file>
  igor-model info <file>
  ```
- [ ] **Task 5.3:** Auto-detect format from extension (.json, .txt)
- [ ] **Task 5.4:** Support conversion between formats:
  - Legacy parms → JSON struct
  - Legacy marginals → JSON weights
  - JSON struct → Legacy parms
  - JSON weights → Legacy marginals
- [ ] **Task 5.5:** Add --pretty flag for formatted output
- [ ] **Task 5.6:** Add --validate flag for schema validation

**Acceptance Criteria:**
- Convert legacy TRB model to JSON
- Convert JSON back to legacy
- Help message is informative
- Error messages are clear

---

### Phase 6: Documentation and Polish
**Estimated:** 1 day | **Status:** ⬜ Not Started

- [ ] **Task 6.1:** Write README.md for ModelIO module
- [ ] **Task 6.2:** Document JSON schema with examples
- [ ] **Task 6.3:** Add Doxygen comments to all public functions
- [ ] **Task 6.4:** Create migration guide from legacy format
- [ ] **Task 6.5:** Update main IGoR README with ModelIO info
- [ ] **Task 6.6:** Create example model files

**Acceptance Criteria:**
- README explains usage with examples
- All public API documented
- Migration guide complete

---

## Design Pattern Options for Core Integration

### Overview

Phase 4 requires bidirectional conversion between Core classes (`Rec_Event`, `Model_Parms`, `Model_Marginals`) and ModelIO data structures (`EventData`, `ModelData`, `WeightsData`). Several design patterns are viable; each has trade-offs.

---

### Option 1: Factory Method Pattern (Member Functions)

**Concept:** Each `Rec_Event` subclass implements serialization/deserialization as member functions.

#### Implementation

```cpp
// In Rec_Event.h (Core module)
class Rec_Event {
public:
    // Virtual method: serialize to EventData
    virtual igor::modelio::EventData to_event_data() const = 0;

    // Static factory: create from EventData
    static std::shared_ptr<Rec_Event> from_event_data(const igor::modelio::EventData& data);

    // Existing legacy I/O (preserved)
    virtual void write2txt(std::ofstream&, std::queue<int>&) = 0;
};

// In Gene_choice.h
class Gene_choice : public Rec_Event {
public:
    igor::modelio::EventData to_event_data() const override;

    // Helper for base class factory
    static std::shared_ptr<Gene_choice> from_event_data(const igor::modelio::EventData& data);
};
```

#### Pros
- ✅ **Encapsulation:** Each event type knows its own serialization
- ✅ **Single Responsibility:** Event logic + serialization in one place
- ✅ **Simple to understand:** Straightforward virtual function pattern
- ✅ **Easy to extend:** New event types just implement the interface
- ✅ **Type safety:** Compile-time checking of conversions

#### Cons
- ❌ **Tight coupling:** Core depends on ModelIO headers
- ❌ **Recompilation:** Changes to ModelIO require Core rebuild
- ❌ **Mixing concerns:** Domain logic mixed with serialization
- ❌ **Hard to test:** Need full Core objects for serialization tests
- ❌ **Dependency direction:** Core → ModelIO (undesirable if ModelIO is "utility")

#### Complexity: ⭐⭐ (Low-Medium)

---

### Option 2: Abstract Factory Pattern

**Concept:** Separate factory class hierarchy handles all conversions, keeping Core and ModelIO independent.

#### Implementation

```cpp
// In ModelIO module: EventFactory.h
namespace igor::modelio {

class EventFactory {
public:
    virtual ~EventFactory() = default;

    // Create Core event from EventData
    virtual std::shared_ptr<Rec_Event> create_event(const EventData& data) const = 0;

    // Extract EventData from Core event
    virtual EventData to_event_data(const Rec_Event& event) const = 0;

    // Check if this factory handles the event type
    virtual bool handles(const std::string& event_type) const = 0;
};

class GeneChoiceFactory : public EventFactory {
public:
    std::shared_ptr<Rec_Event> create_event(const EventData& data) const override;
    EventData to_event_data(const Rec_Event& event) const override;
    bool handles(const std::string& event_type) const override { return event_type == "GeneChoice"; }
};

class DeletionFactory : public EventFactory { /* ... */ };
class InsertionFactory : public EventFactory { /* ... */ };
class DinuclMarkovFactory : public EventFactory { /* ... */ };

// Registry for all factories
class EventFactoryRegistry {
public:
    void register_factory(std::unique_ptr<EventFactory> factory);

    std::shared_ptr<Rec_Event> create_event(const EventData& data) const;
    EventData to_event_data(const Rec_Event& event) const;

private:
    std::vector<std::unique_ptr<EventFactory>> factories_;
};

} // namespace igor::modelio
```

#### Usage

```cpp
// Setup (once, in main or module init)
EventFactoryRegistry registry;
registry.register_factory(std::make_unique<GeneChoiceFactory>());
registry.register_factory(std::make_unique<DeletionFactory>());
registry.register_factory(std::make_unique<InsertionFactory>());
registry.register_factory(std::make_unique<DinuclMarkovFactory>());

// Convert ModelData → Model_Parms
ModelData json_model = structure::read_json("model.json");
std::vector<std::shared_ptr<Rec_Event>> events;
for (const auto& event_data : json_model.events) {
    events.push_back(registry.create_event(event_data));
}

// Convert Model_Parms → ModelData
ModelData output_model;
for (const auto& event : model_parms.get_events()) {
    output_model.events.push_back(registry.to_event_data(*event));
}
```

#### Pros
- ✅ **Loose coupling:** Core and ModelIO remain independent
- ✅ **Testability:** Factory logic tested separately from Core
- ✅ **Dependency inversion:** Both depend on abstractions
- ✅ **Extensibility:** Register new factories at runtime
- ✅ **Clear separation:** Serialization logic isolated from domain logic
- ✅ **Flexible deployment:** Can have different factory implementations

#### Cons
- ❌ **More complex:** Extra abstraction layer to maintain
- ❌ **Boilerplate:** Need factory class for each event type
- ❌ **Runtime registration:** Must remember to register all factories
- ❌ **Indirection:** More objects and function calls
- ❌ **Type casting:** Factories may need dynamic_cast for event-specific operations

#### Complexity: ⭐⭐⭐⭐ (Medium-High)

---

### Option 3: Visitor Pattern

**Concept:** Separate visitor classes traverse Core objects and produce EventData without Core knowing about serialization.

#### Implementation

```cpp
// In ModelIO module: EventDataSerializer.h
namespace igor::modelio {

class EventDataSerializer {
public:
    EventData serialize(const Rec_Event& event) const;

    // Visit each event type
    EventData visit(const Gene_choice& event) const;
    EventData visit(const Deletion& event) const;
    EventData visit(const Insertion& event) const;
    EventData visit(const Dinucl_markov& event) const;
};

class EventDataDeserializer {
public:
    std::shared_ptr<Rec_Event> deserialize(const EventData& data) const;
};

} // namespace igor::modelio
```

#### Core changes (minimal)

```cpp
// Add accept method to Rec_Event (only change to Core)
class Rec_Event {
public:
    // Accept visitor for double dispatch
    virtual void accept(class EventDataSerializer& visitor) const = 0;
};

class Gene_choice : public Rec_Event {
public:
    void accept(EventDataSerializer& visitor) const override {
        visitor.visit(*this);
    }
};
```

#### Pros
- ✅ **Open/Closed Principle:** Add new operations without changing Core
- ✅ **Separation of concerns:** Serialization logic in ModelIO
- ✅ **Single place for format logic:** All serialization in one visitor class
- ✅ **Easier to add new formats:** Create new visitor (e.g., YAMLSerializer)
- ✅ **Core stays clean:** Minimal intrusion (just accept() method)

#### Cons
- ❌ **Double dispatch:** Requires accept() method in Core
- ❌ **Still coupled:** Core needs forward declaration of visitor
- ❌ **Hard to extend types:** Adding new event types requires visitor changes
- ❌ **Circular dependency:** Core includes visitor, visitor includes Core
- ❌ **Not idiomatic for serialization:** Better suited for operations on object graphs

#### Complexity: ⭐⭐⭐ (Medium)

---

### Option 4: Adapter/Bridge Pattern

**Concept:** ModelIO provides adapter functions that operate on Core objects via their public interface only.

#### Implementation

```cpp
// In ModelIO module: CoreAdapters.h
namespace igor::modelio {

// No changes to Core classes needed!

// Adapters use only public interfaces
EventData adapt_to_event_data(const Rec_Event& event);
std::shared_ptr<Rec_Event> adapt_from_event_data(const EventData& data);

// Specialize for each type using type traits or tag dispatch
template<typename EventType>
EventData adapt_event(const EventType& event) {
    // Use only public methods: get_type(), get_realizations(), etc.
}

} // namespace igor::modelio
```

#### Pros
- ✅ **Zero Core changes:** Uses existing public interface
- ✅ **Complete decoupling:** Core unaware of ModelIO
- ✅ **Testable:** Adapter logic independent
- ✅ **Non-invasive:** Can be added as external module

#### Cons
- ❌ **Limited access:** Can only use public interfaces
- ❌ **May be incomplete:** If Core lacks necessary accessors
- ❌ **Duplication:** Adapter logic duplicates domain knowledge
- ❌ **Maintenance burden:** Must track Core API changes
- ❌ **Type dispatch:** May need RTTI (dynamic_cast/typeid) for event types

#### Complexity: ⭐⭐⭐ (Medium)

---

### Option 5: Serialization Traits (Policy-Based Design)

**Concept:** Use template specialization to define serialization behavior per type, no runtime polymorphism.

#### Implementation

```cpp
// In ModelIO module: SerializationTraits.h
namespace igor::modelio {

// Default trait (not implemented)
template<typename T>
struct SerializationTraits;

// Specialize for each event type
template<>
struct SerializationTraits<Gene_choice> {
    static EventData to_event_data(const Gene_choice& event) {
        EventData data;
        data.type = "GeneChoice";
        data.gene_class = event.get_gene_class();
        // ... use public accessors
        return data;
    }

    static std::unique_ptr<Gene_choice> from_event_data(const EventData& data) {
        // Construct Gene_choice from data
        return std::make_unique<Gene_choice>(/* ... */);
    }
};

// Generic serialize/deserialize using traits
template<typename EventType>
EventData serialize(const EventType& event) {
    return SerializationTraits<EventType>::to_event_data(event);
}

template<typename EventType>
std::unique_ptr<EventType> deserialize(const EventData& data) {
    return SerializationTraits<EventType>::from_event_data(data);
}

} // namespace igor::modelio
```

#### Pros
- ✅ **No virtual calls:** Compile-time dispatch
- ✅ **Type-safe:** Compile-time type checking
- ✅ **Zero Core changes:** Pure template specialization
- ✅ **Performance:** No runtime overhead
- ✅ **Flexible:** Can specialize behavior per type easily

#### Cons
- ❌ **Header-only or explicit instantiation:** Template code in headers
- ❌ **No runtime polymorphism:** Can't serialize `Rec_Event*` directly
- ❌ **Type erasure needed:** Need wrapper for runtime dispatch
- ❌ **Compile-time bloat:** Template instantiation for each type
- ❌ **Limited to public API:** Same issue as Adapter pattern

#### Complexity: ⭐⭐⭐⭐ (Medium-High, due to template metaprogramming)

---

### Option 6: Current Hybrid Approach (Standalone + Manual Bridge)

**Concept:** Keep ModelIO standalone; users manually bridge Core ↔ ModelIO when needed.

#### Implementation

```cpp
// User code manually converts
void save_model_as_json(const Model_Parms& parms, const std::string& path) {
    igor::modelio::ModelData data;

    // Manual conversion for each event
    for (const auto& event : parms.get_events()) {
        igor::modelio::EventData ed;

        if (auto* gc = dynamic_cast<Gene_choice*>(event.get())) {
            ed.type = "GeneChoice";
            ed.gene_class = gc->get_gene_class();
            // ... manual field copying
        } else if (auto* del = dynamic_cast<Deletion*>(event.get())) {
            ed.type = "Deletion";
            // ...
        }
        // ... handle other types

        data.events.push_back(ed);
    }

    igor::modelio::structure::write_json(path, data);
}
```

#### Pros
- ✅ **Simple:** No pattern, just procedural code
- ✅ **Flexible:** User controls conversion logic
- ✅ **No framework:** Minimal abstraction
- ✅ **Complete independence:** ModelIO and Core fully decoupled

#### Cons
- ❌ **Manual work:** Users write boilerplate for every use case
- ❌ **Error-prone:** Easy to forget fields or types
- ❌ **Not reusable:** Conversion logic scattered
- ❌ **Maintenance:** Must update all conversion sites when format changes
- ❌ **No standardization:** Different users may convert differently

#### Complexity: ⭐ (Very Low, but pushes complexity to users)

---

### Option 7: Core-Owned Factory (Domain-Driven Design)

**Concept:** Core module provides factory interface; ModelIO readers call it to create Rec_Event objects directly.

#### Implementation

```cpp
// In Core module: Event/Factory.h
namespace igor::event {

class Factory {
public:
    // Create Rec_Event from structured data
    std::shared_ptr<Rec_Event> create(const igor::modelio::EventData& data) const;

    // Or with type dispatch
    std::shared_ptr<Rec_Event> create_gene_choice(const igor::modelio::EventData& data) const;
    std::shared_ptr<Rec_Event> create_deletion(const igor::modelio::EventData& data) const;
    std::shared_ptr<Rec_Event> create_insertion(const igor::modelio::EventData& data) const;
    std::shared_ptr<Rec_Event> create_dinucl_markov(const igor::modelio::EventData& data) const;
};

} // namespace igor::event
```

#### Usage in ModelIO

```cpp
// In ModelIO: ModelReader.cpp
namespace igor::modelio::structure {

ModelData read_json(const std::string& filepath) {
    // Parse JSON file to EventData
    nlohmann::json j = /* ... */;
    ModelData model;

    // Core's factory creates the actual events
    igor::event::Factory factory;

    for (const auto& event_json : j["events"]) {
        EventData event_data = /* parse */;

        // Factory creates Rec_Event based on type
        auto rec_event = factory.create(event_data);

        // Store in model or convert back to EventData for ModelData
        model.events.push_back(event_data);
    }

    return model;
}

} // namespace igor::modelio::structure
```

#### Alternative: Direct Population

```cpp
// ModelIO readers directly populate Core structures
namespace igor::modelio::structure {

Model_Parms read_json_to_parms(const std::string& filepath) {
    nlohmann::json j = /* ... */;

    igor::event::Factory factory;
    Model_Parms parms;

    for (const auto& event_json : j["events"]) {
        EventData event_data = /* parse */;
        auto rec_event = factory.create(event_data);
        parms.add_event(std::move(rec_event));
    }

    return parms;
}

} // namespace igor::modelio::structure
```

#### Pros
- ✅ **Natural dependency direction:** I/O layer depends on domain (standard architecture)
- ✅ **Core controls creation:** Encapsulation - Core knows how to build its objects
- ✅ **Single source of truth:** One place for object construction logic
- ✅ **Domain-driven design:** Domain layer provides services to infrastructure layer
- ✅ **Testable:** Factory logic tested in Core; I/O tested in ModelIO
- ✅ **Extensible:** Core can evolve creation logic without touching ModelIO
- ✅ **Clean separation:** ModelIO does I/O; Core does domain logic

#### Cons
- ❌ **Core depends on ModelIO headers:** Must include EventData definition
- ❌ **Bidirectional dependency:** Core imports ModelIO types (coupling)
- ❌ **Recompilation:** Changes to EventData require Core rebuild
- ❌ **Module boundary blur:** Core must understand external data format
- ⚠️ **Serialization direction:** Creating events from JSON is clear, but extracting EventData from Rec_Event needs separate mechanism

#### Complexity: ⭐⭐ (Low-Medium)

#### Comparison with Option 2

| Aspect | Option 2 (ModelIO Factory) | Option 7 (Core Factory) |
|--------|---------------------------|------------------------|
| **Dependency** | Core ← ModelIO (independent) | Core → ModelIO (coupled) |
| **Who owns factory?** | ModelIO | Core |
| **Who knows formats?** | ModelIO | Core (via EventData) |
| **Best for...** | Pluggable serialization | Domain-centric design |
| **Serialization/Deserialization** | Symmetric (both in ModelIO) | Asymmetric (create in Core, extract in ModelIO?) |

#### When to Use

**Prefer Option 7 if:**
- Core is the primary module, ModelIO is infrastructure
- You want domain layer to control object creation
- Core team owns both Core and ModelIO
- Single application (no external tools using ModelIO standalone)

**Prefer Option 2 if:**
- Need ModelIO standalone (CLI tools, external apps)
- Want to keep Core independent of file formats
- Multiple serialization targets (JSON, YAML, XML)
- Different teams maintain Core vs ModelIO

---

## Recommendation Matrix

| Criterion | Option 1<br/>Factory Method | Option 2<br/>Abstract Factory | Option 3<br/>Visitor | Option 4<br/>Adapter | Option 5<br/>Traits | Option 6<br/>Manual | Option 7<br/>Core Factory |
|-----------|---------|------------|---------|---------|--------|--------|--------|
| **Coupling** | ❌ Tight | ✅ Loose | ⚠️ Medium | ✅ Loose | ✅ Loose | ✅ None | ❌ Bidirectional |
| **Dependency Direction** | Core → ModelIO | ModelIO → Core | Core ↔ ModelIO | ModelIO → Core | ModelIO → Core | None | Core → ModelIO |
| **Core Changes** | ❌ Significant | ✅ None | ⚠️ Minimal | ✅ None | ✅ None | ✅ None | ⚠️ Add Factory |
| **Complexity** | ✅ Low | ❌ High | ⚠️ Medium | ⚠️ Medium | ❌ High | ✅ Very Low | ✅ Low |
| **Extensibility** | ✅ Good | ✅ Excellent | ⚠️ Fair | ⚠️ Fair | ✅ Good | ❌ Poor | ✅ Good |
| **Testability** | ⚠️ Medium | ✅ Excellent | ✅ Good | ✅ Good | ✅ Good | ❌ Poor | ✅ Good |
| **Performance** | ✅ Excellent | ⚠️ Good | ⚠️ Good | ⚠️ Good | ✅ Excellent | ✅ Excellent | ✅ Excellent |
| **Type Safety** | ✅ Excellent | ⚠️ Good | ✅ Good | ⚠️ Fair | ✅ Excellent | ❌ Poor | ✅ Excellent |
| **Maintainability** | ✅ Good | ⚠️ Medium | ✅ Good | ⚠️ Medium | ⚠️ Medium | ❌ Poor | ✅ Good |
| **Idiomatic C++** | ✅ Yes | ✅ Yes | ⚠️ Unusual | ⚠️ Unusual | ⚠️ Advanced | ✅ Yes | ✅ Yes |
| **Standalone ModelIO** | ❌ No | ✅ Yes | ❌ No | ✅ Yes | ✅ Yes | ✅ Yes | ❌ No |
| **Domain-Centric** | ⚠️ Medium | ❌ No | ❌ No | ❌ No | ❌ No | N/A | ✅ Yes |

---

## Recommended Approach

### User's Proposal: **Option 7 (Core Factory)** ⭐

**Rationale:**
> "To me, the readers directly populate the Core structure."

This is a **domain-driven design** approach where:
1. **Core owns object creation:** `igor::event::Factory` in Core module
2. **ModelIO readers call Core factory:** `structure::read_json()` → `Factory::create(EventData)`
3. **Natural dependency flow:** ModelIO (infrastructure) → Core (domain)
4. **Single responsibility:** ModelIO does I/O; Core does business logic

#### Implementation Strategy

```cpp
// Core module: src/igor/Core/Event/Factory.h
namespace igor::event {

class Factory {
public:
    // Create Rec_Event from EventData
    std::shared_ptr<Rec_Event> create(const igor::modelio::EventData& data) const;

    // Type-specific creation methods
    std::shared_ptr<Gene_choice> create_gene_choice(const igor::modelio::EventData& data) const;
    std::shared_ptr<Deletion> create_deletion(const igor::modelio::EventData& data) const;
    std::shared_ptr<Insertion> create_insertion(const igor::modelio::EventData& data) const;
    std::shared_ptr<Dinucl_markov> create_dinucl_markov(const igor::modelio::EventData& data) const;

    // Extract EventData from Rec_Event (for serialization)
    igor::modelio::EventData extract(const Rec_Event& event) const;
};

} // namespace igor::event
```

```cpp
// ModelIO module: src/igor/ModelIO/CoreIntegration.h
namespace igor::modelio {

// Conversion functions using Core's factory
Model_Parms model_parms_from_json(const std::string& filepath);
void model_parms_to_json(const Model_Parms& parms, const std::string& filepath, bool pretty = false);

} // namespace igor::modelio
```

**Trade-offs:**
- ✅ Core controls its object creation (encapsulation)
- ✅ Simple and straightforward
- ✅ Domain-centric design
- ❌ Core depends on EventData (ModelIO headers)
- ❌ Cannot use ModelIO standalone without Core

**When to use:**
- Core and ModelIO are part of the same project
- Domain layer should control object lifecycle
- No need for standalone ModelIO tools

---

### Alternative: **Option 2 (Abstract Factory in ModelIO)**

**If you need ModelIO to be standalone** (e.g., for external CLI tools, Python bindings):

**Rationale:**
1. **ModelIO owns serialization:** Factories live in ModelIO module
2. **No Core dependency in ModelIO data structures:** Core depends on ModelIO only at factory boundary
3. **Pluggable:** Can swap factory implementations

#### Implementation Strategy (from earlier)

```cpp
// ModelIO provides the abstract interface
namespace igor::modelio {

class EventFactory {
public:
    virtual ~EventFactory() = default;
    virtual std::shared_ptr<Rec_Event> create(const EventData&) const = 0;
    virtual EventData extract(const Rec_Event&) const = 0;
    virtual bool handles(const std::string& type) const = 0;
};

// Concrete factories
class GeneChoiceFactory : public EventFactory { /* ... */ };
class DeletionFactory : public EventFactory { /* ... */ };

class EventFactoryRegistry {
public:
    void register_factory(std::unique_ptr<EventFactory> factory);
    std::shared_ptr<Rec_Event> create_event(const EventData& data) const;
    EventData extract_data(const Rec_Event& event) const;
};

} // namespace igor::modelio
```

**Trade-offs:**
- ✅ ModelIO can be used without Core (standalone tools)
- ✅ Highly extensible and testable
- ✅ Clean separation of concerns
- ❌ More complex (registry, abstract factories)
- ❌ Runtime dispatch overhead

**When to use:**
- Need standalone ModelIO library
- Multiple serialization targets
- Different teams maintain Core vs ModelIO

---

## Decision Summary

| Requirement | Recommended Option |
|-------------|-------------------|
| **Core and ModelIO tightly integrated** | Option 7 (Core Factory) ⭐ |
| **Domain-driven design** | Option 7 (Core Factory) ⭐ |
| **Standalone ModelIO tools needed** | Option 2 (Abstract Factory) |
| **Maximum flexibility/extensibility** | Option 2 (Abstract Factory) |
| **Simplest implementation** | Option 7 (Core Factory) ⭐ |
| **Need to avoid Core ↔ ModelIO coupling** | Option 2 (Abstract Factory) |

**Based on user's statement, proceed with Option 7 (Core Factory).**

---

### Rejected Options (with rationale)

- **Option 3 (Visitor):** Overly complex for serialization use case; better suited for AST traversal
- **Option 5 (Traits):** Too template-heavy; loses runtime polymorphism benefits
- **Option 6 (Manual):** Pushes burden to users; not suitable for library

---

## Performance Targets

| Operation | Target | Notes |
|-----------|--------|-------|
| JSON read (TRB model) | < 50ms | ~100 events, ~90 realizations |
| JSON write (TRB model) | < 50ms | Pretty-printed |
| Legacy read | < 100ms | Existing performance baseline |
| Legacy write | < 100ms | Existing performance baseline |
| Schema validation | < 10ms | Optional, disabled by default |

---

## Test Strategy

### Unit Tests
- Common type serialization/deserialization
- Edge cases: empty events, missing fields, malformed JSON
- Tensor dimension validation

### Integration Tests
- Full round-trip tests with real TRB/TRA models
- Cross-format conversion tests
- CLI tool tests

### Compatibility Tests
- Load converted models in existing IGoR code
- Verify inference produces same results
- Test with OLGA/ALICE (if available)

---

## Dependencies

| Dependency | Version | Purpose |
|------------|---------|---------|
| nlohmann/json | >= 3.11 | JSON parsing/generation |
| Catch2 | >= 3.0 | Unit testing |

---

## Open Questions

1. **Weights storage:** Per-event vs central container?
   - **Decision:** Per-event (matches JSON structure, simpler)

2. **Error rate in struct or weights?**
   - **Decision:** Type in struct, parameters in weights

3. **Sequence types format?**
   - **Decision:** Need to align with #15 refactoring branch

4. **Format versioning strategy?**
   - **Decision:** Semantic versioning (1.0.0), check format_version on read

5. **Should ModelIO directly integrate with Rec_Event classes or be standalone?**
   - **Decision:** Two-layer architecture
     - **Layer 1 (Current):** Standalone module for format conversion, CLI tools, validation
     - **Layer 2 (Phase 6):** Factory integration with Core classes for in-memory conversion
   - **Rationale:**
     - Separation of concerns: format I/O vs domain logic
     - Enables standalone tools without Core dependencies
     - Allows gradual migration path
     - Factory pattern provides clean extension point

6. **How to handle event-centric parsing mentioned in the issue?**
   - **Decision:** Deferred to Phase 6
   - Each `Rec_Event` subclass will implement:
     - `static std::shared_ptr<Rec_Event> from_json(const EventData&)`
     - `EventData to_json() const`
   - Maintains single responsibility: events know their own serialization
   - Current ModelIO provides the data structures for these methods to consume/produce

---

## Timeline

| Phase | Duration | Dependencies |
|-------|----------|--------------|
| Phase 0: Setup | 0.5 days | None |
| Phase 1: Common Types | 1 day | Phase 0 |
| Phase 2: Structure I/O | 3 days | Phase 1 |
| Phase 3: Weights I/O | 2 days | Phase 1 |
| **Phase 4: Core Integration** | **3-5 days** | **Phases 2, 3; Optional coordination with #15, #16** |
| Phase 5: CLI Tool | 1 day | Phases 2, 3, 4 |
| Phase 6: Documentation | 1 day | Phase 5 |
| **Total** | **11.5-13.5 days** | **Complete ModelIO with Core integration** |

**Critical Path:** Phases 0 → 1 → 2/3 (parallel) → 4 → 5 → 6

---

## Status Tracking

### Current Phase: Phase 2/3 Complete, Phase 4 Next (Design Decision Required)

| Phase | Status | Notes |
|-------|--------|-------|
| Phase 0: Setup | ✅ Complete | All infrastructure in place |
| Phase 1: Common Types | ✅ Complete | All structs defined, 10 test cases (51 assertions) |
| Phase 2: Structure I/O | ✅ Complete | JSON + Legacy I/O working |
| Phase 3: Weights I/O | ✅ Complete | JSON + Legacy I/O working |
| Phase 4: Core Integration | ⬜ **NEXT** | **Design pattern decision required (see options above)** |
| Phase 5: CLI Tool | ⬜ Not Started | Depends on Phase 4 |
| Phase 6: Documentation | ⬜ Not Started | |

### Progress Summary

**Completed:**
- ✅ Module setup and build configuration
- ✅ All data structures (ModelData, WeightsData, etc.)
- ✅ JSON reading/writing for structure and weights
- ✅ Legacy reading/writing for structure and weights
- ✅ Comprehensive test suite (10 test cases, 51 assertions)
- ✅ Modern C++ idioms (detail namespace, RAII test fixtures, Catch2 matchers)
- ✅ Test utilities following Streaming module pattern

**Current Decision Point:**
- ⚠️ **Phase 4 design pattern selection:** Choose between 6 options for Core integration
  - **Recommended:** Abstract Factory + Adapter (loose coupling, testable, extensible)
  - **Alternative:** Factory Method (simpler but tighter coupling)
  - See "Design Pattern Options" section for detailed comparison

**Next Steps:**
1. **Choose design pattern** for Phase 4 (Core integration)
2. Implement factories for event conversion
3. Add Model_Parms/Model_Marginals conversion functions
4. Build CLI tool (Phase 5)
5. Documentation (Phase 6)

**Future Work (after Phase 4):**
- ⬜ CLI tool (igor-model) for format conversion
- ⬜ Documentation and examples
- ⬜ Performance benchmarking
- ⬜ Integration with #15 (sequence types) and #16 (marginals modernization)

---

## References

- [Issue #32](https://github.com/IGoR-immuno-stack/IGoR/issues/32) - Original issue (Model format modernization plan)
- [Issue #15](https://github.com/IGoR-immuno-stack/IGoR/pull/15) - Sequence type refactoring
- [Issue #16](https://github.com/IGoR-immuno-stack/IGoR/issues/16) - Marginals modernization
- [Issue #9](https://github.com/IGoR-immuno-stack/IGoR/issues/9) - Inferred genomic info
- [nlohmann/json](https://github.com/nlohmann/json) - JSON library
- [JSON Schema](https://json-schema.org/) - Schema specification

## Architecture Notes

### Two-Layer Design Rationale

The ModelIO module adopts a **two-layer architecture** to balance immediate utility with long-term integration:

**Layer 1 (Current Implementation - Phases 0-5):**
- **Purpose:** Standalone format conversion and validation
- **Benefits:**
  - No dependencies on Core event classes
  - Can be used by external tools (CLI, Python bindings, etc.)
  - Easier to test and maintain
  - Provides immediate value for model migration
- **Use cases:**
  - `igor-model convert` CLI tool
  - Format validation and schema checking
  - Model inspection and debugging
  - Integration with external tools (OLGA, ALICE)

**Layer 2 (Future - Phase 6):**
- **Purpose:** Deep integration with Core module via factory pattern
- **Benefits:**
  - Event-centric serialization (each event type handles its own format)
  - In-memory conversion between Core objects and JSON
  - Enables gradual migration of existing codebase
  - Maintains backward compatibility with legacy I/O
- **Use cases:**
  - `Model_Parms::read_model_parms()` auto-detects JSON vs legacy
  - Inference pipeline can load models from either format
  - Model building tools can output both formats

### Event-Centric Parsing (Phase 6)

Per the updated issue #32 proposal, event parsing should be delegated to each event type:

```cpp
// Each Rec_Event subclass knows how to serialize itself
class Gene_choice : public Rec_Event {
public:
    // Factory: construct from JSON
    static std::shared_ptr<Gene_choice> from_json(const igor::modelio::EventData& data);

    // Serialize to JSON
    igor::modelio::EventData to_json() const override;

    // Existing legacy format (preserved for compatibility)
    void write2txt(std::ofstream&, std::queue<int>&) override;
};
```

This approach:
- Maintains single responsibility (events know their serialization)
- Allows extension for custom event types
- Keeps format knowledge close to the domain logic
- Enables progressive enhancement (add JSON support per event type)

### Coordination with Other Refactorings

- **Issue #15 (Sequence Types):** Phase 6 must align sequence type definitions between ModelIO and Core
- **Issue #16 (Marginals/mdspan):** Tensor dimension handling should coordinate with mdspan adoption
- **Issue #9 (Inferred Genomic Info):** Decide whether inferred data belongs in structure or weights

---

## Next Action

**Immediate (Phase 4):** Create CLI tool for model conversion and validation
- Provides immediate utility for users migrating to JSON format
- Tests the standalone ModelIO API in a real-world context
- No Core dependencies needed

**Near-term (Phase 5):** Documentation and examples
- Enable users to understand and adopt the new format
- Provide migration guide from legacy format

**Long-term (Phase 6):** Factory integration with Core module
- Coordinate with #15 and #16 refactorings
- Requires design review for Core module changes
- Enables in-memory format conversion
