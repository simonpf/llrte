#include <llrte/absorption.h>
#include <llrte/atmosphere.h>
#include <llrte/grids.h>
#include <llrte/photons.h>
#include <llrte/scattering.h>
#include <llrte/solvers/forward.h>
#include <llrte/sources.h>
#include <llrte/surfaces.h>
#include <llrte/tracers.h>
#include <llrte/types/vector.h>
#include <llrte/data.h>
#include <llrte/utils/data.h>

#include <memory>

using Float = double;

std::tuple<Float, Float, Float> run_experiment(size_t n_photons, Float g) {
  using V3 = llrte::Vector<3, Float>;
  using Grid = llrte::RegularGrid<Float>;
  using AbsorptionModel = llrte::ConstantAbsorption<Float>;
  using ScatteringPlane = llrte::maths::geometry::FixedScatteringPlane<2>;
  using ScatteringModel = llrte::HenyeyGreenstein<Float, ScatteringPlane>;
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

  using ReflectingSurface = llrte::surfaces::ReflectingPlane<
      V3, llrte::surfaces::Lambertian<ScatteringPlane>>;

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
  auto absorption_model =
      llrte::ConstantAbsorption<Float>((1.0 - ssa) * od / 100.0);
  auto scattering_model = ScatteringModel(ssa * od / 100.0, g, 180);
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
  std::vector<float> g = {-0.9, -0.8, -0.6, -0.4, -0.2, 0.0,
                          0.2,  0.4,  0.6,  0.8,  0.9};
  for (size_t i = 0; i < g.size(); ++i) {
    auto results = run_experiment(100000, g[i]);
    std::cout << g[i] << " " << std::get<0>(results) << " " << std::get<1>(results);
    std::cout << " " << std::get<2>(results) << std::endl;
  }
}
