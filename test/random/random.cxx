#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "llrte/random.h"
#include "llrte/io/netcdf.h"
#include "llrte/data.h"
#include "llrte/maths.h"

TEMPLATE_TEST_CASE("Generator",
                   "Threaded numbers",
                   float, double) {
    using llrte::maths::small;
    size_t n = 2048;
    auto data = llrte::Tensor<TestType, 2>({n, 8});
    auto generator = llrte::Generator<float>();

    #pragma omp parallel for
    for (size_t i = 0; i < 8; ++i) {
        for (size_t j = 0; j < n; ++j) {
            data(j, i) = generator.sample_tau();
        }
    }

    REQUIRE(data(0u, 0u) != 0u);
    REQUIRE(data(0u, 0u) != data(0u, 1u));

    llrte::io::NetCDFFile file("data_exp.nc", true);
    file.add_dimension("x", n);
    file.add_dimension("y", 8);
    file.store_variable(data, "data", {"x", "y"});

    #pragma omp parallel for
    for (size_t i = 0; i < 8; ++i) {
        for (size_t j = 0; j < n; ++j) {
            data(j, i) = generator.sample_angle_uniform();
        }
    }

    REQUIRE(data(0u,0u) != 0);
    REQUIRE(data(0u, 0u) != data(0u, 1u));

    file = llrte::io::NetCDFFile("data_phi.nc", true);
    file.add_dimension("x", n);
    file.add_dimension("y", 8);
    file.store_variable(data, "data", {"x", "y"});

#pragma omp parallel for
    for (size_t i = 0; i < 8; ++i) {
        for (size_t j = 0; j < n; ++j) {
            data(j, i) = generator.sample_zenith_angle();
        }
    }

    REQUIRE(data(0u, 0u) != 0);
    REQUIRE(data(0u, 0u) != data(0u, 1u));

    file = llrte::io::NetCDFFile("data_theta.nc", true);
    file.add_dimension("x", n);
    file.add_dimension("y", 8);
    file.store_variable(data, "data", {"x", "y"});
}
