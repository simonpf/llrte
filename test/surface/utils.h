#ifndef _LLRTE_TEST_SURFACE_UTILS_H_
#define _LLRTE_TEST_SURFACE_UTILS_H_

#include <llrte/grids.h>
#include <llrte/absorption.h>
#include <llrte/atmosphere.h>
#include <llrte/tracers.h>
#include <llrte/sources.h>
#include <llrte/photons.h>
#include <llrte/scattering.h>
#include <llrte/solvers/monte_carlo.h>
#include <llrte/types/vector.h>

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

template<typename Boundaries>
auto make_test_atmosphere(Boundaries &boundaries) {

    using V3 = llrte::Vector<3, float>;
    using Float = float;
    using Grid = llrte::RegularGrid<Float>;
    using AbsorptionModel = llrte::NoAbsorption<Float>;
    using ScatteringModel = llrte::NoScattering<Float>;
    using Atmosphere = llrte::Atmosphere<Grid, AbsorptionModel, ScatteringModel, Boundaries>;
    using Photon = llrte::FixedEnergyPhoton<V3>;
    using Source = llrte::BeamSource<Photon>;

    auto source_position = V3{};
    source_position[0] = 0.0;
    source_position[1] = 0.0;
    source_position[2] = 0.0;

    auto source_direction = V3{};
    source_direction[0] = 1.0;
    source_direction[1] = 0.0;
    source_direction[2] = 0.0;

    auto source = Source(source_position, source_direction);

    auto x = make_linear_vector<Float>(0.0, 1.0, 2);
    auto y = make_linear_vector<Float>(-0.5, 0.5, 2);
    auto z = make_linear_vector<Float>(-0.5, 0.5, 2);
    size_t shape[3] = {2, 2, 2};

    auto grid = Grid{shape, x, y, z};
    auto absorption_model = llrte::NoAbsorption<Float>();
    auto scattering_model = llrte::NoScattering<Float>();
    auto atmosphere = Atmosphere{grid,
                                 absorption_model,
                                 scattering_model,
                                 boundaries};

    return std::make_pair(atmosphere, source);
}


#endif
