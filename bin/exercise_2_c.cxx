#include <llrte/absorption.h>
#include <llrte/atmosphere.h>
#include <llrte/grids.h>
#include <llrte/photons.h>
#include <llrte/scattering.h>
#include <llrte/monte_carlo.h>
#include <llrte/sources.h>
#include <llrte/data.h>
#include <llrte/random.h>
#include <llrte/surfaces.h>
#include <llrte/tracers.h>
#include <llrte/types/vector.h>
#include <memory>

template<typename Float>
void run_experiment(size_t n_grid_cells,
                    size_t n_photons,
                    float optical_depth,
                    float fb_ratio,
                    float ssa,
                    float surface_albedo,
                    std::string filename) {

  using V3 = llrte::Vector3<Float>;
  using Grid = llrte::RegularGrid<Float>;
  using AbsorptionModel = llrte::ConstantAbsorption<Float>;
  using ScatteringModel = llrte::BidirectionalScattering<Float>;
  using Tracer = llrte::tracers::AbsorptionTracer<Grid>;
  using Photon = llrte::Photon<V3, llrte::GridPosition>;
  using Source = llrte::BeamSource<Photon>;

  //////////////////////////////////////////////////////////////////////
  // Surface
  //////////////////////////////////////////////////////////////////////

  // Setup black boundary.
  auto base_b = V3{10e3, 0.0, 0.0};
  auto normal_b = V3{-1.0, 0.0, 0.0};

  using ReflectingSurface =
      llrte::surfaces::ReflectingPlane<V3, llrte::surfaces::BackwardsDirection>;

  auto surfaces =
      std::make_tuple(ReflectingSurface(base_b, normal_b, surface_albedo));
  using Surfaces = decltype(surfaces);
  using Atmosphere =
      llrte::Atmosphere<Grid, AbsorptionModel, ScatteringModel, Surfaces>;
  using Generator = llrte::Generator<Float>;
  using Solver = llrte::MonteCarlo<Atmosphere, Generator, Tracer>;

  auto source_position = V3{0.0, 0.0, 0.0};
  auto source_direction = V3{1.0, 0.0, 0.0};
  auto source = Source(source_position, source_direction);

  Float start = 0.0e3;
  Float stop = 10.0e3;
  auto x = llrte::Array<Float>::fill_linear(start, stop, n_grid_cells + 1);
  auto y = llrte::Array<Float>::fill_linear(-0.5, 0.5, 2);
  auto z = llrte::Array<Float>::fill_linear(-0.5, 0.5, 2);

  auto grid = Grid{x, y, z};
  auto absorption_model =
      llrte::ConstantAbsorption<Float>((1.0 - ssa) * optical_depth / 1e4);
  auto scattering_model = llrte::BidirectionalScattering<Float>(
      ssa * optical_depth / 1e4, fb_ratio);
  auto atmosphere =
      Atmosphere{grid, absorption_model, scattering_model, surfaces};

  Tracer tracer{grid};
  Generator generator{};
  auto solver = Solver(atmosphere, generator, tracer);

  solver.initialize();
  for (size_t i = 0; i < n_photons; i++) {
      solver.forward(source);
  }
  tracer.save(filename);

  std::cout << "Upwelling intensity:           ";
  std::cout << tracer.get_total_leaving_counts(0) / n_photons << std::endl;
  std::cout << "Total absorbed intensity:      ";
  std::cout << tracer.get_total_absorption_counts() / n_photons << std::endl;
  std::cout << "Intensity absorbed by surface: ";
  std::cout << solver.atmosphere().template get_boundary<0>().get_absorbed_energy() / n_photons
            << std::endl;
}

int main(int /*argc*/, const char ** /***argv*/) {
  std::cout << "Surface albedo A = 0.3:" << std::endl;
  run_experiment<float>(100, 100000, 1.0, 0.7, 0.8, 0.3, "results_2_c.bin");
  std::cout << "Surface albedo A = 0.0:" << std::endl;
  run_experiment<float>(100, 100000, 1.0, 0.7, 0.8, 0.0, "results_2_c_r.bin");
}
