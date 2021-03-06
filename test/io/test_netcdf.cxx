#include "llrte/io/netcdf.h"
#include "catch.hpp"

template <typename T>
bool read_write() {

    auto t = llrte::eigen::Tensor<T, 3>(3, 2, 1);
    t.setConstant(3.0);

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            t.coeffRef(i, j, 0) = i + j;
        }
    }

    llrte::io::NetCDFFile file("file.nc", true);
    file.add_dimension("x", 3);
    file.add_dimension("y", 2);
    file.store_variable(t, "data", {"x", "y"});

    file = llrte::io::NetCDFFile("file.nc", false);
    auto t2 = file.load_variable<T, 3>("data");

    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            std::cout << t.coeffRef(i, j, 0) << std::endl;
        }
    }

    return true;
}

int main(int /*argc*/, char **/*argv*/) {

    auto t = llrte::eigen::Tensor<double, 2>(3, 3);
    t.setConstant(1.0);

    read_write<float>();
}
