#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "llrte/grids/regular.h"
#include "llrte/data.h"
#include "llrte/maths.h"
#include "llrte/common.h"


template<typename T>
__global__ void scale_array(T t) {
    for (size_t i = 0; i < t.size(); ++i) {
        t[i] *= 2.0;
    }
}

TEMPLATE_TEST_CASE("Data, Cuda",
                   "Basics",
                   float, double) {
    using llrte::maths::small;
    auto array = llrte::Data<TestType>(10);
    array.fill(1.0);
    array.device();
    scale_array<<<1,1>>>(array);
    array.host();
    REQUIRE(small(array[0] - 2.0));
}
