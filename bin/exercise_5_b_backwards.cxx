#include <llrte/io/netcdf.h>
#include <llrte/absorption.h>
#include <llrte/atmosphere.h>
#include <llrte/grids.h>
#include <llrte/photons.h>
#include <llrte/scattering.h>
#include <llrte/solvers/backward.h>
#include <llrte/sources.h>
#include <llrte/data.h>
#include <llrte/surfaces.h>
#include <llrte/random.h>
#include <llrte/types/vector.h>
#include <llrte/sensor.h>

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
                            Float y_0,
                            Float y_1,
                            Float g)
        : hg_(sigma_2, g, 180), rayleigh_(sigma_1, 180), x_0_(x_0), x_1_(x_1),
          y_0_(y_0), y_1_(y_1)
        {}

    template <typename Grid, typename GridPos>
    Float get_scattering_coefficient(const Grid &g, const GridPos &gp) {
        if ((gp.x() > x_0_) && (gp.x() <= x_1_) && (gp.y() > y_0_) && (gp.y() <= y_1_)) {
            return hg_.get_scattering_coefficient(g, gp);
        } else {
            return rayleigh_.get_scattering_coefficient(g, gp);
        }
    }

    template <typename Grid, typename GridPos>
    PhaseFunction get_phase_function(const Grid &g, const GridPos &gp) {
        if ((gp.x() > x_0_) && (gp.x() <= x_1_) && (gp.y() > y_0_) && (gp.y() <= y_1_)) {
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
    Float y_0_;
    Float y_1_;
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
                            Float x_1,
                            Float y_0,
                            Float y_1)
        :  sigma_1_(sigma_1), sigma_2_(sigma_2), x_0_(x_0), x_1_(x_1),
           y_0_(y_0), y_1_(y_1)
        {}

    template <typename Grid, typename GridPos>
    constexpr F get_absorption_coefficient(const Grid &/*g*/, const GridPos &gp) {
        if ((gp.x() > x_0_) && (gp.x() <= x_1_) && (gp.y() > y_0_) && (gp.y() <= y_1_)) {
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
    Float y_0_;
    Float y_1_;
};

////////////////////////////////////////////////////////////////////////////////
// The simulation
////////////////////////////////////////////////////////////////////////////////

void run_experiment() {
  using Float = float;
  using V3 = llrte::Vector<3, Float>;
  using Grid = llrte::RegularGrid<Float>;
  using AbsorptionModel = HeterogeneousAbsorption<Float>;
  using ScatteringModel = HeterogeneousScattering<Float>;
  using C = llrte::Constants<Float>;

  using Photon = llrte::Photon<V3>;
  using Source = llrte::PlanarSource<Photon>;

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
  using Solver = llrte::BackwardSolver<Atmosphere &, Source &>;

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

  auto source_normal = V3{{-1.0, 0.0, 0.0}};
  auto source = Source(1.0, source_normal);

  //////////////////////////////////////////////////////////////////////
  // Sensor
  //////////////////////////////////////////////////////////////////////

  auto sensor_position = V3{{1.0, 30e3, 0.0}};
  auto sensor_x = V3{{0.0, 1.0, 0.0}};
  auto sensor_y = V3{{0.0, 0.0, 1.0}};
  auto sensor_dx = llrte::Array<Float>::fill_linear(-30e3, 30e3, 240);
  auto sensor_dy = llrte::Array<Float>::fill_linear(0.0, 0.0, 1);
  auto sensor_dz = llrte::Array<Float>::fill_linear(-0.5 * C::pi, 0.5 * C::pi, 2);
  auto sensor_da = llrte::Array<Float>::fill_linear(0.0, 0.0, 1);

  llrte::SensorArray<V3> sensor{sensor_position,
                                sensor_x,
                                sensor_y,
                                sensor_dx,
                                sensor_dy,
                                sensor_dz,
                                sensor_da};


  //////////////////////////////////////////////////////////////////////
  // Domain
  //////////////////////////////////////////////////////////////////////

  auto x = llrte::Array<Float>::fill_linear(0.0, 20e3, 100);
  auto y = llrte::Array<Float>::fill_linear(0.0, 60e3, 1000);
  auto z = llrte::Array<Float>::fill_linear(-0.5, 0.5, 2);

  auto grid = Grid{x, y, z};
  auto absorption_model = AbsorptionModel(0.0, 0.1 * 1e-3, 2e3, 7e3, 20e3, 40e3);
  auto scattering_model = ScatteringModel(0.0025 * 1e-3, 0.9 * 1e-3, 2e3, 7e3, 20e3, 40e3, 0.9);
  auto atmosphere =
      Atmosphere{grid, absorption_model, scattering_model, surfaces};

  llrte::Generator<Float> generator{};
  auto solver = Solver(atmosphere, source);
  sensor.sample(generator, solver, 1000);

  sensor.dump("exercise_5_b_backwards.nc");

}

int main(int /*argc*/, const char ** /***argv*/) {
  run_experiment();
}
