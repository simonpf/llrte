#ifndef _LLRTE_SOLVERS_BACKWARD_H_
#define _LLRTE_SOLVERS_BACKWARD_H_

#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <tuple>

#include <llrte/constants.h>
#include <llrte/random.h>
#include <math.h>

template <
    typename Atmosphere,
    typename Source
    >
class BackwardSolver {
 public:
  using Vector = typename std::remove_reference_t<Source>::Vector;
  using Float = typename std::remove_reference_t<Source>::Float;
  using Photon = typename std::remove_reference_t<Source>::Photon;

  BackwardSolver(Atmosphere atmosphere,
                 Source source)
      : atmosphere_(atmosphere), source_(source) {
    // Nothing to do here.
  }

  template <typename Generator>
  Photon propagate_photon(Generator &generator, Photon &photon) {

    auto position = atmosphere_.get_grid_position(photon.get_position());
    auto tau = generator.sample_tau();

    while (true) {
      auto absorption_xc = atmosphere_.get_absorption_coefficient(position);
      auto scattering_xc = atmosphere_.get_scattering_coefficient(position);

      auto l = tau / scattering_xc;
      auto intersection =
          atmosphere_.get_intersection(position, photon.get_direction(), l);
      auto d = std::get<0>(intersection);
      tau -= d * scattering_xc;

      // Check if left atmosphere.
      if (d <= -1.0) {
        photon.scale_energy(source_.get_intensity(photon));
        break;
      }

      // Handle absorption.
      auto f_abs = exp(-absorption_xc * d);
      photon.scale_energy(f_abs);
      if (photon.get_energy() < minimum_energy_) {
        break;
      }

      position = std::get<1>(intersection);

      // Check if scattering event.
      if (l <= d) {
        auto phase_function = atmosphere_.get_phase_function(position);
        photon.scatter(generator, phase_function);
        tau = generator.sample_tau();
      } else {
        bool hit = atmosphere_.apply_boundaries(generator_, photon);
        if (hit) {
          position = atmosphere_.get_grid_position(photon.get_position());
        }
      }
    }
    return photon;
  }

 private:
  size_t n_photons = 0;

  Atmosphere atmosphere_;
  Source source_;
  Generator<Float> generator_;
  Float minimum_energy_ = 1e-6;
};

}  // namespace llrte
#endif
