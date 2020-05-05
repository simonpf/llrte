#ifndef _LLRTE_TEST_SURFACE_UTILS_H_
#define _LLRTE_TEST_SURFACE_UTILS_H_

#include <llrte/grids.h>
#include <llrte/absorption.h>
#include <llrte/atmosphere.h>
#include <llrte/tracers.h>
#include <llrte/sources.h>
#include <llrte/photons.h>
#include <llrte/tracers.h>
#include <llrte/scattering.h>
#include <llrte/monte_carlo.h>
#include <llrte/types/vector.h>

template<typename Boundaries, typename Tracer>
auto make_test_atmosphere(Boundaries &boundaries,
                          Tracer tracer) {

    using V3 = llrte::Vector3<float>;
    using Float = float;
    using Grid = llrte::RegularGrid<Float>;
    using AbsorptionModel = llrte::NoAbsorption<Float>;
    using ScatteringModel = llrte::NoScattering<Float>;
    using Atmosphere = llrte::Atmosphere<Grid, AbsorptionModel, ScatteringModel, Boundaries>;
    using Photon = llrte::Photon<V3, llrte::GridPosition>;
    using Source = llrte::BeamSource<Photon>;
    using Generator = llrte::Generator<Float>;
    using Solver = llrte::MonteCarlo<Atmosphere, Generator, Tracer>;

    auto source_position = V3{0.0, 0.0, 0.0};
    auto source_direction = V3{1.0, 0.0, 0.0};
    auto source = Source(source_position, source_direction);

    auto x = llrte::Array<Float>::fill_linear(0.0, 1.0, 2);
    auto y = llrte::Array<Float>::fill_linear(-0.5, 0.5, 2);
    auto z = llrte::Array<Float>::fill_linear(-0.5, 0.5, 2);

    auto grid = Grid{std::move(x),
                     std::move(y),
                     std::move(z)};
    auto absorption_model = llrte::NoAbsorption<Float>();
    auto scattering_model = llrte::NoScattering<Float>();
    auto atmosphere = Atmosphere{std::move(grid),
                                 absorption_model,
                                 scattering_model,
                                 boundaries};
    Generator generator{};
    auto solver = Solver(std::move(atmosphere), generator, std::move(tracer));

    return std::make_tuple(std::move(solver), source, grid);
}


#endif
