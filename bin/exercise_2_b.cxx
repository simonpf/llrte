#include <llrte/grids.h>
#include <llrte/absorption.h>
#include <llrte/atmosphere.h>
#include <llrte/tracers.h>
#include <llrte/data.h>
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
        if (position.x() < boundary_) {
            return abs_1_;
        } else {
            return abs_2_;
        }
    }

private:

    Float abs_1_, abs_2_, boundary_;

};

template <typename F>
class HeterogeneousScattering : public llrte::BidirectionalScattering<F> {

public:

    using Float = F;
    using PhaseFunction = typename llrte::BidirectionalScattering<Float>::PhaseFunction;

    HeterogeneousScattering(Float sigma_s_1,
                            Float sigma_2_2,
                            Float fb_ratio_1,
                            Float fb_ratio_2,
                            Float boundary) :
        llrte::BidirectionalScattering<Float>(0.0, 0.0),
        sigma_s_1_(sigma_s_1),
        sigma_s_2_(sigma_2_2),
        fb_ratio_1_(fb_ratio_1),
        fb_ratio_2_(fb_ratio_2),
        boundary_(boundary)
    {
        // Nothing to do.
    }

    template <typename Grid, typename Position>
    Float get_scattering_coefficient(Grid /*grid*/, Position position) {
        if (position.x() < boundary_) {
            return sigma_s_1_;
        } else {
            return sigma_s_2_;
        }
    }

    template <typename Grid, typename Position>
    PhaseFunction get_phase_function(Grid /*grid*/, Position position) {
        if (position.x() < boundary_) {
            return PhaseFunction{fb_ratio_1_};
        } else {
            return PhaseFunction{fb_ratio_2_};
        }
    }

private:

    F sigma_s_1_, sigma_s_2_;
    F fb_ratio_1_, fb_ratio_2_;
    F boundary_;

};

int main(int /*argc*/, const char **/***argv*/) {

    size_t n_grid_cells = 100;
    size_t n_photons = 10000;
    std::string filename{"results_2_b.bin"};

    using Float = float;
    using V3 = llrte::Vector<3, Float>;
    using Grid = llrte::RegularGrid<Float>;
    using AbsorptionModel = HeterogeneousAbsorption<Float>;
    using ScatteringModel = HeterogeneousScattering<Float>;
    using Atmosphere = llrte::Atmosphere<Grid, AbsorptionModel, ScatteringModel>;
    using Tracer = llrte::AbsorptionTracer<Grid>;
    using Photon = llrte::Photon<V3>;
    using Source = llrte::BeamSource<Photon>;
    using Solver = llrte::ForwardSolver<Atmosphere, Source, Tracer>;

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
    auto x = llrte::Array<Float>::fill_linear(start, stop, n_grid_cells + 1);
    auto y = llrte::Array<Float>::fill_linear(-0.5, 0.5, 2);
    auto z = llrte::Array<Float>::fill_linear(-0.5, 0.5, 2);

    auto grid = Grid{x, y, z};
    auto absorption_model = HeterogeneousAbsorption<Float>(0.2e-4, 0.2e-4, 5e3);
    auto scattering_model = HeterogeneousScattering<Float>(0.8e-4, 0.8e-4, 0.8, 0.2, 5e3);
    auto atmosphere = Atmosphere{grid, absorption_model, scattering_model};

    auto solver = Solver(atmosphere, source);

    Tracer::initialize(grid);
    for (size_t i = 0; i < n_photons; i++) {
        solver.sample_photon();
    }
    Tracer::dump(filename);
}
