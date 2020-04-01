#include <llrte/data.h>
#include <llrte/grids.h>
#include <llrte/absorption.h>
#include <llrte/atmosphere.h>
#include <llrte/tracers.h>
#include <llrte/sources.h>
#include <llrte/photons.h>
#include <llrte/scattering.h>
#include <llrte/monte_carlo.h>
#include <llrte/types/vector.h>
#include <memory>

template<typename Float = float>
void run_experiment(size_t n_grid_cells,
                    size_t n_photons,
                    float optical_depth,
                    float fb_ratio,
                    float ssa,
                    std::string filename) {

    using V3 = llrte::Vector3<Float>;
    using Grid = llrte::RegularGrid<Float>;
    using AbsorptionModel = llrte::ConstantAbsorption<Float>;
    using ScatteringModel = llrte::BidirectionalScattering<Float>;
    using Atmosphere = llrte::Atmosphere<Grid, AbsorptionModel, ScatteringModel>;
    using Tracer = llrte::tracers::AbsorptionTracer<Grid>;
    using Photon = llrte::Photon<V3, llrte::GridPosition>;
    using Source = llrte::BeamSource<Photon>;
    using Generator = llrte::Generator<Float>;
    using Solver = llrte::MonteCarlo<Atmosphere, Generator, Tracer>;

    auto source_position = V3{0.0, 0.0, 0.0};
    auto source_direction = V3{1.0, 0.0, 0.0};
    auto source = Source(source_position, source_direction);

    Float start = 0.0e3;
    Float stop = 10.0e3;
    auto x = llrte::Array<Float>::fill_linear(start, stop, n_grid_cells + 1);
    auto y = llrte::Array<Float>::fill_linear(-0.5, 0.5, 2);
    auto z = llrte::Array<Float>::fill_linear(-0.5, 0.5, 2);

    auto grid = Grid{x, y, z};
    auto absorption_model = llrte::ConstantAbsorption<Float>((1.0 - ssa) * optical_depth / 1e4);
    auto scattering_model = llrte::BidirectionalScattering<Float>(ssa * optical_depth / 1e4, fb_ratio);
    auto atmosphere = Atmosphere{grid, absorption_model, scattering_model};

    Tracer tracer{grid};
    Generator generator{};
    auto solver = Solver(atmosphere, generator, tracer);

    solver.initialize();
    for (size_t i = 0; i < n_photons; i++) {
        solver.forward(source);
    }
    tracer.save(filename);
}

int main(int /*argc*/, const char **/***argv*/) {
    run_experiment(10,   10000,   1.0, 0.8, 0.8, "results_2_a_1.nc");
    run_experiment(100,  10000,   1.0, 0.8, 0.8, "results_2_a_2.nc");
    run_experiment(1000, 10000,   1.0, 0.8, 0.8, "results_2_a_3.nc");
    run_experiment(100,  100,     1.0, 0.8, 0.8, "results_2_a_4.nc");
    run_experiment(100,  10000,   1.0, 0.8, 0.8, "results_2_a_5.nc");
    run_experiment(100,  1000000, 1.0, 0.8, 0.8, "results_2_a_6.nc");
    run_experiment(100,  10000,    1.0, 0.8, 0.8, "results_2_a_7.nc");
    run_experiment(100,  10000,    2.0, 0.8, 0.8, "results_2_a_8.nc");
    run_experiment(100,  10000,    3.0, 0.8, 0.8, "results_2_a_9.nc");
    run_experiment(100,  10000,    1.0, 0.2, 0.8, "results_2_a_10.nc");
    run_experiment(100,  10000,    1.0, 0.5, 0.8, "results_2_a_11.nc");
    run_experiment(100,  10000,    1.0, 0.8, 0.8, "results_2_a_12.nc");
    run_experiment(100,  10000,    1.0, 0.8, 0.2, "results_2_a_13.nc");
    run_experiment(100,  10000,    1.0, 0.8, 0.5, "results_2_a_14.nc");
    run_experiment(100,  10000,    1.0, 0.8, 0.8, "results_2_a_15.nc");
}
