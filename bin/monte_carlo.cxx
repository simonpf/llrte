#include <llrte/solvers/forward.h>
#include <llrte/grids.h>
#include <llrte/absorption.h>
#include <llrte/atmosphere.h>
#include <llrte/tracers.h>
#include <llrte/sources.h>
#include <llrte/photons.h>
#include <llrte/types/vector.h>
#include <memory>

template<typename F>
std::shared_ptr<F[]> make_linear_vector(F start,
                                        F stop,
                                        size_t steps) {


    std::shared_ptr<F[]> v{new F[steps]};

    F d = (stop - start) / (steps - 1);
    F x = start;
    for (size_t i = 0; i < steps;  ++i) {
        v[i] = x;
        x = x + d;
    }
    return v;
}

int main(int /*argc*/, const char **/***argv*/) {


    using V3 = llrte::Vector<3, float>;
    using Float = float;
    using Grid = llrte::RegularGrid<Float>;
    using AbsorptionModel = llrte::ConstantAbsorption<Float>;
    using ScatteringModel = llrte::NoScattering<Float>
    using Atmosphere = llrte::Atmosphere<Grid, AbsorptionModel, ScatteringModel>;
    using Tracer = llrte::Histogram<Grid>;
    using Photon = llrte::Photon<V3, Tracer>;
    using Source = llrte::BeamSource<Photon>;
    using Solver = llrte::ForwardSolver<Atmosphere, Source>;

    auto source_position = V3{};
    source_position[0] = 0.0;
    source_position[1] = 0.0;
    source_position[2] = 0.0;

    auto source_direction = V3{};
    source_direction[0] = 1.0;
    source_direction[1] = 0.0;
    source_direction[2] = 0.0;

    auto source = Source(source_position, source_direction);

    float start = 0.0e3;
    float stop = 10.0e3;
    auto x = make_linear_vector<Float>(start, stop, 101);
    auto y = make_linear_vector<Float>(-0.5, 0.5, 2);
    auto z = make_linear_vector<Float>(-0.5, 0.5, 2);
    size_t shape[3] = {101, 2, 2};

    auto grid = Grid{shape, x, y, z};
    auto absorption_model = llrte::ConstantAbsorption<float>(1e-3);
    auto scattering_model = NoScattering();
    auto atmosphere = Atmosphere{grid, absorption_model, scattering_model};

    auto solver = Solver(atmosphere, source);

    Tracer::initialize(grid);
    for (size_t i = 0; i < 10000000; i++) {
        solver.sample_photon();
    }

    Tracer::dump("results.bin");
}
