#include "llrte/data.h"
#include "llrte/io/netcdf.h"
#include "llrte/random.h"


template<typename FloatType>
void sample_angle() {
    auto data = llrte::Array<FloatType>(1000000);
    auto generator = llrte::Generator<FloatType>();
    auto sample = [&generator]() {return generator.sample_zenith_angle();};
    data.template map(sample);

    llrte::io::NetCDFFile file("file.nc", true);
    file.add_dimension("x", data.size());
    file.store_variable(data, "data_full", {"x"});

    auto sample2 = [&generator]() {return generator.sample_zenith_angle(0.0,
                                                                        llrte::Constants<FloatType>::pi / 2.0);};
    data.template map(sample2);
    file.store_variable(data, "data_half", {"x"});
}

int main(int /*argc*/, const char **/*argv*/) {
    sample_angle<float>();
}
