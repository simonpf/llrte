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

void run_experiment(size_t n_grid_cells, size_t n_photons, float optical_depth,
                    float fb_ratio, float ssa, float surface_albedo,
                    std::string filename) {
  using Float = float;
  using V3 = llrte::Vector<3, Float>;
  using Grid = llrte::RegularGrid<Float>;
  using AbsorptionModel = llrte::ConstantAbsorption<Float>;
  using ScatteringModel = llrte::BidirectionalScattering<Float>;
  using Tracer = llrte::AbsorptionTracer<Grid>;
  using Photon = llrte::Photon<V3>;
  using Source = llrte::BeamSource<Photon>;

  //////////////////////////////////////////////////////////////////////
  // Surface
  //////////////////////////////////////////////////////////////////////

  // Setup black boundary.
  auto base_b = V3{};
  base_b[0] = 10e3;
  base_b[1] = 0.0;
  base_b[2] = 0.0;

  auto normal_b = V3{};
  normal_b[0] = -1.0;
  normal_b[1] = 0.0;
  normal_b[2] = 0.0;

  using ReflectingSurface =
      llrte::surfaces::ReflectingPlane<V3, llrte::surfaces::BackwardsDirection>;

  auto surfaces =
      std::make_tuple(ReflectingSurface(base_b, normal_b, surface_albedo));
  using Surfaces = decltype(surfaces);
  using Atmosphere =
      llrte::Atmosphere<Grid, AbsorptionModel, ScatteringModel, Surfaces>;
  using Solver = llrte::MonteCarloSolver<Atmosphere &, Source &, Tracer>;

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
  auto absorption_model =
      llrte::ConstantAbsorption<Float>((1.0 - ssa) * optical_depth / 1e4);
  auto scattering_model = llrte::BidirectionalScattering<Float>(
      ssa * optical_depth / 1e4, fb_ratio);
  auto atmosphere =
      Atmosphere{grid, absorption_model, scattering_model, surfaces};

  auto solver = Solver(atmosphere, source);

  Tracer::initialize(grid);
  for (size_t i = 0; i < n_photons; i++) {
    solver.sample_photon();
  }
  Tracer::dump(filename);

  std::cout << "Upwelling intensity:           ";
  std::cout << Tracer::get_total_leaving_counts(0) / n_photons << std::endl;
  std::cout << "Total absorbed intensity:      ";
  std::cout << Tracer::get_total_absorption_counts() / n_photons << std::endl;
  std::cout << "Intensity absorbed by surface: ";
  std::cout << atmosphere.get_boundary<0>().get_absorbed_energy() / n_photons
            << std::endl;
}

int main(int /*argc*/, const char ** /***argv*/) {
  std::cout << "Surface albedo A = 0.3:" << std::endl;
  run_experiment(100, 100000, 1.0, 0.7, 0.8, 0.3, "results_2_c.bin");
  std::cout << "Surface albedo A = 0.0:" << std::endl;
  run_experiment(100, 100000, 1.0, 0.7, 0.8, 0.0, "results_2_c_r.bin");
}
