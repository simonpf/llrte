#include <llrte/absorption.h>
#include <llrte/atmosphere.h>
#include <llrte/grids.h>
#include <llrte/photons.h>
#include <llrte/scattering.h>
#include <llrte/solvers/monte_carlo.h>
#include <llrte/data.h>
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

void run_experiment(size_t n_photons,
                    std::string filename) {
  using Float = double;
  using V3 = llrte::Vector<3, Float>;
  using Grid = llrte::RegularGrid<Float>;
  using AbsorptionModel = llrte::ConstantAbsorption<Float>;
  using ScatteringPlane =
      llrte::maths::geometry::FixedScatteringPlane<2>;
  using ScatteringModel =
      llrte::HenyeyGreenstein<Float, ScatteringPlane>;
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
  base_2_p[1] = 4000;
  base_2_p[2] = 0.0;

  auto normal_p = V3{};
  normal_p[0] = 0.0;
  normal_p[1] = -1.0;
  normal_p[2] = 0.0;

  using ReflectingSurface =
      llrte::surfaces::ReflectingPlane<V3,
                                       llrte::surfaces::Lambertian<ScatteringPlane>>;

  Float sa = 0.7;
  auto surfaces = std::make_tuple(
      ReflectingSurface(base_b, normal_b, sa),
      llrte::surfaces::PeriodicBoundary<V3>(base_1_p, base_2_p, normal_p));
  using Surfaces = decltype(surfaces);
  using Atmosphere =
      llrte::Atmosphere<Grid, AbsorptionModel, ScatteringModel, Surfaces>;
  using Solver = llrte::ForwardSolver<Atmosphere &, Source &, Tracer>;

  //////////////////////////////////////////////////////////////////////
  // Source
  //////////////////////////////////////////////////////////////////////

  auto source_position = V3{};
  source_position[0] = 100;
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
      Source(source_offset, 0, 4000.0, source_position, source_direction);

  //////////////////////////////////////////////////////////////////////
  // Domain
  //////////////////////////////////////////////////////////////////////

  auto x = llrte::Array<Float>::fill_linear(0.0, 100.0, 20);
  auto y = llrte::Array<Float>::fill_linear(0.0, 4000.0, 1000);
  auto z = llrte::Array<Float>::fill_linear(-0.5, 0.5, 2);

  Float od = 5.0;
  Float ssa = 0.8;
  auto grid = Grid{x, y, z};
  auto absorption_model = llrte::ConstantAbsorption<Float>((1.0 - ssa) * od / 100.0);
  auto scattering_model = ScatteringModel(ssa * od / 100.0, 0.9, 180);
  auto atmosphere =
      Atmosphere{grid, absorption_model, scattering_model, surfaces};
  auto solver = Solver(atmosphere, source);

  Tracer::initialize(grid);
  for (size_t i = 0; i < n_photons; i++) {
    solver.sample_photon();
  }
  Tracer::dump(filename);


  std::cout << "Upwelling intensity:           ";
  std::cout << Tracer::get_total_leaving_counts(1) / n_photons << std::endl;
  std::cout << "Total absorbed intensity:      ";
  std::cout << Tracer::get_total_absorption_counts() / n_photons << std::endl;
  std::cout << "Intensity absorbed by surface: ";
  std::cout << atmosphere.get_boundary<0>().get_absorbed_energy() / n_photons << std::endl;

}

int main(int /*argc*/, const char ** /***argv*/) {
  run_experiment(10000, "./results_4_a.bin");
}
