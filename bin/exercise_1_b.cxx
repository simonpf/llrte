#include <llrte/grids.h>
#include <llrte/absorption.h>
#include <llrte/atmosphere.h>
#include <llrte/tracers.h>
#include <llrte/data.h>
#include <llrte/sources.h>
#include <llrte/photons.h>
#include <llrte/scattering.h>
#include <llrte/monte_carlo.h>
#include <llrte/types/vector.h>
#include <memory>

template <typename Float>
class HeterogeneousAbsorption {

public:

    HeterogeneousAbsorption(Float abs_1,
                            Float abs_2,
                            Float boundary) :
        abs_1_(abs_1),
        abs_2_(abs_2),
        boundary_(boundary) {
        // Nothing to do.
    }

    template <typename Grid, typename Position>
    Float get_absorption_coefficient(Grid /*grid*/, Position position) {
        if (position.x < boundary_) {
            return abs_1_;
        } else {
            return abs_2_;
        }
    }

private:

    Float abs_1_, abs_2_, boundary_;

};

template<typename Float=float>
void run_experiment(size_t n_grid_cells,
                    size_t n_photons,
                    std::string filename) {

    using V3 = llrte::Vector3<Float>;
    using Grid = llrte::RegularGrid<Float>;
    using Grid = llrte::RegularGrid<Float>;
    using AbsorptionModel = HeterogeneousAbsorption<Float>;
    using ScatteringModel = llrte::NoScattering<Float>;
    using Atmosphere = llrte::Atmosphere<Grid, AbsorptionModel, ScatteringModel>;
    using Tracer = llrte::tracers::Histogram<Grid>;
    using Photon = llrte::Photon<V3, llrte::GridPosition>;
    using Source = llrte::BeamSource<Photon>;
    using Generator = llrte::Generator<Float>;
    using Solver = llrte::MonteCarlo<Atmosphere, Generator, Tracer>;

    auto source_position = V3{0.0, 0.0, 0.0};
    auto source_direction = V3{1.0, 0.0, 0.0};
    auto source = Source(source_position, source_direction);

    float start = 0.0e3;
    float stop = 10.0e3;
    auto x = llrte::Array<Float>::fill_linear(start, stop, n_grid_cells + 1);
    auto y = llrte::Array<Float>::fill_linear(-0.5, 0.5, 2);
    auto z = llrte::Array<Float>::fill_linear(-0.5, 0.5, 2);

    auto grid = Grid{x, y, z};
    auto absorption_model = HeterogeneousAbsorption<Float>(0.5e-4, 1.5e-4, 5e3);
    auto scattering_model = llrte::NoScattering<Float>();
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
    run_experiment(10, 10000, "results_1_b_1.nc");
    run_experiment(100, 10000, "results_1_b_2.nc");
    run_experiment(1000, 10000, "results_1_b_3.nc");
    run_experiment(100, 100, "results_1_b_4.nc");
    run_experiment(100, 10000, "results_1_b_5.nc");
    run_experiment(100, 1000000, "results_1_b_6.nc");
}
