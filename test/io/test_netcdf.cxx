#include "llrte/io/netcdf.h"
#include "llrte/data.h"
#include "catch2/catch.hpp"

int main(int /*argc*/, char **/*argv*/) {

    auto t = llrte::Tensor<double, 2>({3, 3});
    t.fill(1.0);

    llrte::io::NetCDFFile file("file.nc", true);
    file.add_dimension("x", 3);
    file.add_dimension("y", 3);
    file.store_variable(t, "data", {"x", "y"});
}
