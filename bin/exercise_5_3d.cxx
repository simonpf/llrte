#include <llrte/io/netcdf.h>
#include <llrte/absorption.h>
#include <llrte/atmosphere.h>
#include <llrte/grids.h>
#include <llrte/photons.h>
#include <llrte/scattering.h>
#include <llrte/solvers/monte_carlo.h>
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
    using ScatteringPlane = llrte::maths::geometry::RandomPlane;
    using PhaseFunction = llrte::NumericPhaseFunction<F, ScatteringPlane>;

    HeterogeneousScattering(Float sigma_1,
                            Float sigma_2,
                            Float x_0,
                            Float x_1,
                            Float y_0,
                            Float y_1,
                            Float z_0,
                            Float z_1,
                            Float g)
        : hg_(sigma_2, g, 180), rayleigh_(sigma_1, 180), x_0_(x_0), x_1_(x_1),
          y_0_(y_0), y_1_(y_1), z_0_(z_0), z_1_(z_1)
        {}

    template <typename Grid, typename GridPos>
    Float get_scattering_coefficient(const Grid &g, const GridPos &gp) {
        if ((gp.x() > x_0_) && (gp.x() <= x_1_)
            && (gp.y() > y_0_) && (gp.y() <= y_1_)
            && (gp.z() > z_0_) && (gp.z() <= z_1_)) {
            return hg_.get_scattering_coefficient(g, gp);
        } else {
            return rayleigh_.get_scattering_coefficient(g, gp);
        }
    }

    template <typename Grid, typename GridPos>
    PhaseFunction get_phase_function(const Grid &g, const GridPos &gp) {
        if ((gp.x() > x_0_) && (gp.x() <= x_1_)
            && (gp.y() > y_0_) && (gp.y() <= y_1_)
            && (gp.z() > z_0_) && (gp.z() <= z_1_)) {
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
    Float z_0_;
    Float z_1_;
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
                            Float y_1,
                            Float z_0,
                            Float z_1)
        :  sigma_1_(sigma_1), sigma_2_(sigma_2), x_0_(x_0), x_1_(x_1),
           y_0_(y_0), y_1_(y_1), z_0_(z_0), z_1_(z_1)
        {}

    template <typename Grid, typename GridPos>
    constexpr F get_absorption_coefficient(const Grid &/*g*/, const GridPos &gp) {
        if ((gp.x() > x_0_) && (gp.x() <= x_1_)
            && (gp.y() > y_0_) && (gp.y() <= y_1_)
            && (gp.z() > z_0_) && (gp.z() <= z_1_)) {
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
    Float z_0_;
    Float z_1_;
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
  auto base_1_p_y = V3{{0.0, 0.0, 0.0}};
  auto base_2_p_y = V3{{0.0, 60e3, 0.0}};
  auto normal_p_y = V3{{0.0, -1.0, 0.0}};
  auto base_1_p_z = V3{{0.0, 0.0, 0.0}};
  auto base_2_p_z = V3{{0.0,  0.0, 30e3}};
  auto normal_p_z = V3{{0.0,  0.0, -1.0}};

  using ScatteringPlane =
      llrte::maths::geometry::FixedScatteringPlane<2>;
  using ReflectingSurface =
      llrte::surfaces::ReflectingPlane<V3,
                                       llrte::surfaces::Lambertian<ScatteringPlane>>;

  Float sa = 0.07;
  auto surfaces = std::make_tuple(
      ReflectingSurface(base_b, normal_b, sa),
      llrte::surfaces::PeriodicBoundary<V3>(base_1_p_y, base_2_p_y, normal_p_y),
      llrte::surfaces::PeriodicBoundary<V3>(base_1_p_z, base_2_p_z, normal_p_z)
      );
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
  auto source = Source(1.0, source_normal, C::pi / 180.0);

  //////////////////////////////////////////////////////////////////////
  // Sensor
  //////////////////////////////////////////////////////////////////////

  auto sensor_position = V3{{20e3, 10e3, 15e3}};
  auto sensor_x = V3{{-1.0, -1.0, 0.0}}.normed();
  auto sensor_y = V3{{0.0, 0.0, 1.0}};
  auto sensor_dx = llrte::Array<Float>::fill_linear(-10e3, 10e3, 100);
  auto sensor_dy = llrte::Array<Float>::fill_linear(-10e3, 10e3, 100);
  auto sensor_dz = llrte::Array<Float>::fill_linear(0.0, 0.0, 1);
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
  auto y = llrte::Array<Float>::fill_linear(0.0, 60e3, 200);
  auto z = llrte::Array<Float>::fill_linear(0.0, 30e3, 100);

  auto grid = Grid{x, y, z};
  auto absorption_model = AbsorptionModel(0.0, 0.1 * 1e-3, 2e3, 7e3, 20e3, 40e3, 10e3, 20e3);
  auto scattering_model = ScatteringModel(0.0025 * 1e-3, 0.9 * 1e-3, 2e3, 7e3, 20e3, 40e3, 10e3, 20e3, 0.9);
  auto atmosphere =
      Atmosphere{grid, absorption_model, scattering_model, surfaces};

  llrte::Generator<Float> generator{};
  auto solver = Solver(atmosphere, source);
  sensor.sample(generator, solver, 200);

  sensor.dump("exercise_5_3d.nc");

}

int main(int /*argc*/, const char ** /***argv*/) {
  run_experiment();
}
