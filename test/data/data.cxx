#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "llrte/grids/regular.h"
#include "llrte/data.h"
#include "llrte/maths.h"


TEMPLATE_TEST_CASE("Data, Ownership",
                   "Allocation and operators",
                   float, double) {
    using llrte::maths::small;
    auto array = llrte::Data<TestType>(10);
    REQUIRE(small(array[0] - 0.0));
    array.fill(1.0);
    REQUIRE(small(array[0] - 1.0));
    array += 1.0;
    REQUIRE(small(array[0] - 2.0));
    array *= 2.0;
    REQUIRE(small(array[0] - 4.0));
    array /= 4.0;
    REQUIRE(small(array[0] - 1.0));
}
