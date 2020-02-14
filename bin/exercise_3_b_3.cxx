#include <llrte/absorption.h>
#include <llrte/atmosphere.h>
#include <llrte/grids.h>
#include <llrte/photons.h>
#include <llrte/scattering.h>
#include <llrte/solvers/monte_carlo.h>
#include <llrte/sources.h>
#include <llrte/surfaces.h>
#include <llrte/tracers.h>
#include <llrte/types/vector.h>

#include <memory>

template <typename F>
std::shared_ptr<F[]> make_linear_vector(F start, F stop, size_t steps) {
  std::shared_ptr<F[]> v{new F[steps]};

  F d = (stop - start) / (steps - 1);
  F x = start;
  for (size_t i = 0; i < steps; ++i) {
    v[i] = x;
    x = x + d;
  }
  return v;
}

std::tuple<float, float, float> run_experiment(size_t n_grid_cells,
                                               size_t n_photons,
                                               float optical_depth) {
  using V3 = llrte::Vector<3, float>;
  using Float = float;
  using Grid = llrte::RegularGrid<Float>;
  using AbsorptionModel = llrte::ConstantAbsorption<Float>;
  using ScatteringPlane =
      llrte::maths::geometry::FixedScatteringPlane<2>;
  using ScatteringModel =
      llrte::RayleighScattering<Float, ScatteringPlane>;
  using Tracer = llrte::AbsorptionTracer<Grid>;
  using Photon = llrte::Photon<V3>;
  using Source = llrte::RandomOffset<llrte::BeamSource<Photon>>;

  //////////////////////////////////////////////////////////////////////
  // Surface
  //////////////////////////////////////////////////////////////////////

  // Setup black boundary.
  auto base_b = V3{};
  base_b[0] = 0.0;
  base_b[1] = 0.0;
  base_b[2] = 0.0;

  auto normal_b = V3{};
  normal_b[0] = 1.0;
  normal_b[1] = 0.0;
  normal_b[2] = 0.0;

  // Setup periodic boundary.
  auto base_1_p = V3{};
  base_1_p[0] = 0.0;
  base_1_p[1] = 0.0;
  base_1_p[2] = 0.0;

  auto base_2_p = V3{};
  base_2_p[0] = 0.0;
  base_2_p[1] = 10e3;
  base_2_p[2] = 0.0;

  auto normal_p = V3{};
  normal_p[0] = 0.0;
  normal_p[1] = -1.0;
  normal_p[2] = 0.0;

  using ReflectingSurface =
      llrte::surfaces::ReflectingPlane<V3,
                                       llrte::surfaces::Lambertian<ScatteringPlane>>;

  auto surfaces = std::make_tuple(
      ReflectingSurface(base_b, normal_b, 0.8),
      llrte::surfaces::PeriodicBoundary<V3>(base_1_p, base_2_p, normal_p));
  using Surfaces = decltype(surfaces);
  using Atmosphere =
      llrte::Atmosphere<Grid, AbsorptionModel, ScatteringModel, Surfaces>;
  using Solver = llrte::MonteCarloSolver<Atmosphere &, Source &, Tracer>;

  //////////////////////////////////////////////////////////////////////
  // Source
  //////////////////////////////////////////////////////////////////////

  auto source_position = V3{};
  source_position[0] = 100e3;
  source_position[1] = 0.0;
  source_position[2] = 0.0;

  auto source_direction = V3{};
  source_direction[0] = -1.0;
  source_direction[1] = 0.0;
  source_direction[2] = 0.0;

  auto source_offset = V3{};
  source_offset[0] = 0.0;
  source_offset[1] = 1.0;
  source_offset[2] = 0.0;

  auto source =
      Source(source_offset, 0, 10e3, source_position, source_direction);

  //////////////////////////////////////////////////////////////////////
  // Domain
  //////////////////////////////////////////////////////////////////////

  float start = 0.0e3;
  float stop = 100.0e3;
  auto x = make_linear_vector<Float>(start, stop, n_grid_cells + 1);
  auto y = make_linear_vector<Float>(start, stop, n_grid_cells + 1);
  auto z = make_linear_vector<Float>(-0.5, 0.5, 2);
  size_t shape[3] = {n_grid_cells + 1, n_grid_cells + 1, 2};

  auto grid = Grid{shape, x, y, z};
  auto absorption_model = llrte::ConstantAbsorption<Float>(0.5 * optical_depth * 1e-5);
  auto scattering_model = ScatteringModel(0.5 * optical_depth * 1e-5, 180);
  auto atmosphere =
      Atmosphere{grid, absorption_model, scattering_model, surfaces};
  auto solver = Solver(atmosphere, source);

  Tracer::initialize(grid);
  for (size_t i = 0; i < n_photons; i++) {
    solver.sample_photon();
  }

  float upwelling = Tracer::get_total_leaving_counts(1) / n_photons;
  float absorbed = Tracer::get_total_absorption_counts() / n_photons;
  float surface =
      atmosphere.get_boundary<0>().get_absorbed_energy() / n_photons;

  return std::make_tuple(upwelling, absorbed, surface);
}

int main(int /*argc*/, const char ** /***argv*/) {
  std::vector<float> ods = {0.2, 0.4, 0.6, 0.8, 1.0,
                                1.2, 1.4, 1.6, 1.8, 2.0};
  for (size_t i = 0; i < ods.size(); ++i) {
    auto results = run_experiment(100, 10000, ods[i]);
    std::cout << std::get<0>(results) << " " << std::get<1>(results);
    std::cout << " " << std::get<2>(results) << std::endl;
  }
}
