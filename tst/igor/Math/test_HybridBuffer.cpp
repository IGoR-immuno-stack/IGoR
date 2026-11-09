#include <catch2/catch_test_macros.hpp>
#include <igor/Math/HybridBuffer.h>

using namespace igor::math;

TEST_CASE("HybridBuffer Ownership Semantics", "[Math][HybridBuffer]") {
    
    SECTION("Owning Constructor allocates memory") {
        HybridBuffer<double> buf(10);
        REQUIRE(buf.size() == 10);
        REQUIRE(buf.data() != nullptr);
        REQUIRE(buf.is_owning() == true);
        
        // Write access
        buf[0] = 3.14;
        REQUIRE(buf[0] == 3.14);
    }

    SECTION("Borrowing Constructor wraps memory") {
        double raw[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
        HybridBuffer<double> buf(raw, 5);
        
        REQUIRE(buf.size() == 5);
        REQUIRE(buf.data() == raw);
        REQUIRE(buf.is_owning() == false);
        
        // Check content
        REQUIRE(buf[2] == 3.0);
        
        // Modify via buffer
        buf[0] = 99.0;
        REQUIRE(raw[0] == 99.0);
    }

    SECTION("Move Semantics transfer ownership") {
        HybridBuffer<int> source(100);
        int* original_ptr = source.data();
        
        HybridBuffer<int> dest = std::move(source);
        
        REQUIRE(dest.size() == 100);
        REQUIRE(dest.data() == original_ptr);
        REQUIRE(dest.is_owning() == true);
        
        REQUIRE(source.size() == 0);
        REQUIRE(source.data() == nullptr);
        REQUIRE(source.is_owning() == false);
    }
}
