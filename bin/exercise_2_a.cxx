#include <llrte/grids.h>
#include <llrte/absorption.h>
#include <llrte/atmosphere.h>
#include <llrte/tracers.h>
#include <llrte/sources.h>
#include <llrte/photons.h>
#include <llrte/scattering.h>
#include <llrte/solvers/monte_carlo.h>
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

void run_experiment(size_t n_grid_cells,
                    size_t n_photons,
                    std::string filename) {

    using Float = float;
    using V3 = llrte::Vector<3, Float>;
    using Grid = llrte::RegularGrid<Float>;
    using AbsorptionModel = llrte::ConstantAbsorption<Float>;
    using ScatteringModel = llrte::BidirectionalScattering<Float>;
    using Atmosphere = llrte::Atmosphere<Grid, AbsorptionModel, ScatteringModel>;
    using Tracer = llrte::AbsorptionTracer<Grid>;
    using Photon = llrte::FixedEnergyPhoton<V3, Tracer>;
    using Source = llrte::BeamSource<Photon>;
    using Solver = llrte::MonteCarloSolver<Atmosphere, Source>;

    auto source_position = V3{};
    source_position[0] = 0.0;
    source_position[1] = 0.0;
    source_position[2] = 0.0;

    auto source_direction = V3{};
    source_direction[0] = 1.0;
    source_direction[1] = 0.0;
    source_direction[2] = 0.0;

    auto source = Source(source_position, source_direction);

    Float start = 0.0e3;
    Float stop = 10.0e3;
    auto x = make_linear_vector<Float>(start, stop, n_grid_cells + 1);
    auto y = make_linear_vector<Float>(-0.5, 0.5, 2);
    auto z = make_linear_vector<Float>(-0.5, 0.5, 2);
    size_t shape[3] = {n_grid_cells + 1, 2, 2};

    auto grid = Grid{shape, x, y, z};
    auto absorption_model = llrte::ConstantAbsorption<Float>(2e-5);
    auto scattering_model = llrte::BidirectionalScattering<Float>(8e-5, 0.8);
    auto atmosphere = Atmosphere{grid, absorption_model, scattering_model};

    auto solver = Solver(atmosphere, source);

    Tracer::initialize(grid);
    for (size_t i = 0; i < n_photons; i++) {
        solver.sample_photon();
    }
    Tracer::dump(filename);
}

int main(int /*argc*/, const char **/***argv*/) {
    run_experiment(10, 10000, "results_2_a_1.bin");
    run_experiment(100, 10000, "results_2_a_2.bin");
    run_experiment(1000, 10000, "results_2_a_3.bin");
    run_experiment(100, 100, "results_2_a_4.bin");
    run_experiment(100, 10000, "results_2_a_5.bin");
    run_experiment(100, 1000000, "results_2_a_6.bin");
}
