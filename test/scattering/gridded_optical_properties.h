#ifndef _LLRTE_TEST_GRIDDED_H_
#define _LLRTE_TEST_GRIDDED_H_

#include <llrte/grids.h>
#include <llrte/absorption.h>
#include <llrte/atmosphere.h>
#include <llrte/constants.h>
#include <llrte/tracers.h>
#include <llrte/sources.h>
#include <llrte/photons.h>
#include <llrte/tracers.h>
#include <llrte/scattering.h>
#include <llrte/monte_carlo.h>
#include <llrte/types/vector.h>

// pxx :: export
class PointScattering {
 public:
  using V3 = llrte::Vector3<float>;
  using Float = float;
  using Grid = llrte::RegularGrid<Float>;
  using OpticalProperties = llrte::GriddedOpticalProperties<Grid>;
  using AbsorptionModel = OpticalProperties::AbsorptionModel;
  using ScatteringModel = OpticalProperties::ScatteringModel;
  using Atmosphere =
      llrte::Atmosphere<Grid&, AbsorptionModel, ScatteringModel, std::tuple<>>;
  using Photon = llrte::Photon<V3, llrte::GridPosition>;
  using Source = llrte::BeamSource<Photon>;
  using Generator = llrte::Generator<Float>;
  using Tracer = llrte::tracers::PhotonTracer<Float>;
  using Solver = llrte::MonteCarlo<Atmosphere, Generator, Tracer&>;

  using Vector = llrte::eigen::Vector<Float>;
  template<size_t rank>
      using Tensor = llrte::eigen::Tensor<Float, rank>;
  using PhaseFunctionTensor = llrte::eigen::Tensor<Float, 4>;
  using CoefficientTensor = llrte::eigen::Tensor<Float, 3>;
  using Results = std::tuple<Tensor<2>, Tensor<2>, Tensor<2>, Tensor<2>, llrte::eigen::Tensor<int, 1>>;

  PointScattering(Vector x, Vector y, Vector z, Vector scattering_angles,
                  CoefficientTensor absorption_coefficient,
                  CoefficientTensor scattering_coefficient,
                  PhaseFunctionTensor phase_function)
      :

        grid_(x, y, z),
        optical_properties_(grid_, scattering_angles, absorption_coefficient,
                            scattering_coefficient, phase_function),
        atmosphere_(grid_, optical_properties_.get_absorption_model(),
                    optical_properties_.get_scattering_model()),
        solver_(atmosphere_, Generator{}, tracer_),
        source_(V3{0.0, 0.0, 0.0}, V3{1.0, 0.0, 0.0}) {}

  Results run(size_t n) {
    tracer_.reset(n);
    for (size_t i = 0; i < n; ++i) {
      solver_.forward(source_);
    }
    return std::make_tuple(tracer_.get_incoming_positions(),
                           tracer_.get_incoming_directions(),
                           tracer_.get_outgoing_positions(),
                           tracer_.get_outgoing_directions(),
                           tracer_.get_scattering_frequencies());
  }

 private:
  Grid grid_;
  llrte::GriddedOpticalProperties<Grid> optical_properties_;
  Atmosphere atmosphere_;
  Solver solver_;
  Source source_;
  Tracer tracer_{0};
};

#endif
