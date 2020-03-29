#define CATCH_CONFIG_MAIN
#include "llrte/grids/regular.h"

#include "catch.hpp"
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

TEMPLATE_TEST_CASE("Data, Ownership",
                   "Ownership",
                   float, double) {
    using llrte::maths::small;
    auto array1 = llrte::Data<TestType>(10);
    auto array2 = new llrte::Data<TestType>(array1.get_data_pointer(), 10);
    auto array3 = llrte::Data<TestType>(10);
    array1.fill(1.0);
    array3.fill(3.0);
    REQUIRE(array1[0] == (*array2)[0]);
    (*array2) = array3;
    delete array2;
    REQUIRE(array1[0] == array3[0]);
}
