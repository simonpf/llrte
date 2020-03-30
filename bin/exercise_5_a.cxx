#include <llrte/absorption.h>
#include <llrte/atmosphere.h>
#include <llrte/data.h>
#include <llrte/grids.h>
#include <llrte/photons.h>
#include <llrte/scattering.h>
#include <llrte/solvers/forward.h>
#include <llrte/sources.h>
#include <llrte/surfaces.h>
#include <llrte/tracers.h>
#include <llrte/types/vector.h>

#include <memory>

////////////////////////////////////////////////////////////////////////////////
// Heterogeneous scattering
////////////////////////////////////////////////////////////////////////////////

template <typename F>
class HeterogeneousScattering {
public:
    using Float = F;
    using ScatteringPlane =
        llrte::maths::geometry::FixedScatteringPlane<2>;
    using PhaseFunction = llrte::NumericPhaseFunction<F, ScatteringPlane>;

    HeterogeneousScattering(Float sigma_1,
                            Float sigma_2,
                            Float x_0,
                            Float x_1,
                            Float g)
        : hg_(sigma_2, g, 180), rayleigh_(sigma_1, 180), x_0_(x_0), x_1_(x_1)
        {}

    template <typename Grid, typename GridPos>
    Float get_scattering_coefficient(const Grid &g, const GridPos &gp) {
        if ((gp.x() > x_0_) && (gp.x() <= x_1_)) {
            return hg_.get_scattering_coefficient(g, gp);
        } else {
            return rayleigh_.get_scattering_coefficient(g, gp);
        }
    }

    template <typename Grid, typename GridPos>
    PhaseFunction get_phase_function(const Grid &g, const GridPos &gp) {
        if ((gp.x() > x_0_) && (gp.x() <= x_1_)) {
            return hg_.get_phase_function(g, gp);
        } else {
            return rayleigh_.get_phase_function(g, gp);
        }
    }

private:
    llrte::HenyeyGreenstein<Float, ScatteringPlane> hg_;
    llrte::RayleighScattering<Float, ScatteringPlane> rayleigh_;
    Float x_0_;
    Float x_1_;
};

////////////////////////////////////////////////////////////////////////////////
// Heterogeneous absorption
////////////////////////////////////////////////////////////////////////////////

template <typename F>
class HeterogeneousAbsorption {
public:
    using Float = F;
    HeterogeneousAbsorption(Float sigma_1,
                            Float sigma_2,
                            Float x_0,
                            Float x_1)
        :  sigma_1_(sigma_1), sigma_2_(sigma_2), x_0_(x_0), x_1_(x_1)
        {}

    template <typename Grid, typename GridPos>
    constexpr F get_absorption_coefficient(const Grid &/*g*/, const GridPos &gp) {
        if ((gp.x() > x_0_) && (gp.x() <= x_1_)) {
            return sigma_2_;
        } else {
            return sigma_1_;
        }
    }

private:
    Float sigma_1_;
    Float sigma_2_;
    Float x_0_;
    Float x_1_;
};

////////////////////////////////////////////////////////////////////////////////
// The simulation
////////////////////////////////////////////////////////////////////////////////

void run_experiment(size_t n_photons,
                    std::string filename) {
  using Float = float;
  using V3 = llrte::Vector<3, Float>;
  using Grid = llrte::RegularGrid<Float>;
  using AbsorptionModel = HeterogeneousAbsorption<Float>;
  using ScatteringModel = HeterogeneousScattering<Float>;

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
  base_2_p[1] = 60e3;
  base_2_p[2] = 0.0;

  auto normal_p = V3{};
  normal_p[0] = 0.0;
  normal_p[1] = -1.0;
  normal_p[2] = 0.0;

  using ScatteringPlane =
      llrte::maths::geometry::FixedScatteringPlane<2>;
  using ReflectingSurface =
      llrte::surfaces::ReflectingPlane<V3,
                                       llrte::surfaces::Lambertian<ScatteringPlane>>;

  Float sa = 0.07;
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
  source_position[0] = 20e3;
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
      Source(source_offset, 0, 60e3, source_position, source_direction);

  //////////////////////////////////////////////////////////////////////
  // Domain
  //////////////////////////////////////////////////////////////////////

  auto x = llrte::Array<Float>::fill_linear(0.0, 20e3, 200);
  auto y = llrte::Array<Float>::fill_linear(0.0, 60e3, 600);
  auto z = llrte::Array<Float>::fill_linear(-0.5, 0.5, 2);

  auto grid = Grid{x, y, z};
  auto absorption_model = AbsorptionModel(0.0, 0.1 * 1e-3, 2e3, 7e3);
  auto scattering_model = ScatteringModel(0.0025 * 1e-3, 0.9 * 1e-3, 2e3, 7e3, 0.9);
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
  run_experiment(1000000, "./results_5_a.bin");
}
