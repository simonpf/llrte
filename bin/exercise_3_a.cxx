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

void run_experiment(size_t n_grid_cells, size_t n_photons,
                    std::string filename) {
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
  base_b[0] = 1.0;
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
  source_position[0] = 10e3;
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
  float stop = 10.0e3;
  auto x = make_linear_vector<Float>(start, stop, n_grid_cells + 1);
  auto y = make_linear_vector<Float>(start, stop, n_grid_cells + 1);
  auto z = make_linear_vector<Float>(-0.5, 0.5, 2);
  size_t shape[3] = {n_grid_cells + 1, n_grid_cells + 1, 2};

  auto grid = Grid{shape, x, y, z};
  auto absorption_model = llrte::ConstantAbsorption<Float>(0.5e-4);
  auto scattering_model = ScatteringModel(0.5e-4, 100);
  auto atmosphere =
      Atmosphere{grid, absorption_model, scattering_model, surfaces};
  auto solver = Solver(atmosphere, source);

  Tracer::initialize(grid);
  for (size_t i = 0; i < n_photons; i++) {
    solver.sample_photon();
  }
  Tracer::dump(filename);

  std::cout << "surface energy: " << atmosphere.get_boundary<0>().get_absorbed_energy() << std::endl;
}

int main(int /*argc*/, const char ** /***argv*/) {
  run_experiment(100, 100000, "./results_3_a.bin");
}
