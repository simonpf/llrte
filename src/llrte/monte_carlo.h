#ifndef _LLRTE_SOLVERS_MONTE_CARLO_H_
#define _LLRTE_SOLVERS_MONTE_CARLO_H_

#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <tuple>

#include <llrte/constants.h>
#include <llrte/random.h>
#include <llrte/tracers.h>
#include <math.h>

namespace llrte {

/** MonteCarlo
 *
 * The MonteCarlo implements the basic Monte Carlo functionality
 * for MC simulations of radiative transfer. This class provides an abstract
 * implementation of the Monte Carlo algorithm, which consists of propagating
 * a photon through an atmosphere.
 *
 * @tparam Atmosphere The atmosphere type providing basic atmospheric
 * @tparam Photon The photon type to use.
 * properties such as absorption and scattering properties.
 * @tparam Tracer A tracer type that is used to store results of the simulation.
 */
template <
    typename Atmosphere,
    typename Generator,
    typename Tracer = tracers::NoTrace
>
class MonteCarlo {
 public:

 using Float = typename Atmosphere::Float;

 MonteCarlo(Atmosphere atmosphere,
            Generator generator = Generator{},
            Tracer tracer = Tracer{})
 : atmosphere_(atmosphere), tracer_(tracer), generator_(generator) {
     // Nothing to do here.
  }

 void initialize() {
     generator_.initialize();
 }

 /**
  * Propagate photon through atmosphere.
  *
  * Given a photon in the atmosphere propagate the photon until it leaves the
  * atmosphere or its energy falls below a threshold.
  *
  * @param photon The photon to propagate through the atmosphere.
  * @tparam Photon The photon type used.
  */
 template <typename Photon>
  Photon propagate_photon(Photon photon) {


    auto tau = generator_.sample_tau();

    while (true) {
      auto absorption_xc = atmosphere_.get_absorption_coefficient(photon.position);
      auto scattering_xc = atmosphere_.get_scattering_coefficient(photon.position);

      auto l = tau / scattering_xc;
      auto d = atmosphere_.step(photon, l);
      tau -= d * scattering_xc;

      // Check if left atmosphere.
      if (d <= -1.0) {
        tracer_.left_atmosphere(photon);
        break;
      }

      // Handle absorption.
      auto f_abs = exp(-absorption_xc * d);
      tracer_.absorption(photon, photon.get_energy() * (1.0 - f_abs));
      photon.scale_energy(f_abs);
      if (photon.get_energy() < minimum_energy_) {
        tracer_.out_of_energy(photon);
        break;
      }

      // Check if scattering event.
      if (l <= d) {
        auto phase_function = atmosphere_.get_phase_function(photon.position);
        photon.scatter(generator_, phase_function);
        tau = generator_.sample_tau();
        tracer_.scattering(photon);
      } else {
        bool hit = atmosphere_.apply_boundaries(generator_, photon);
        if (hit) {
            photon = atmosphere_.place_on_grid(photon.position,
                                               photon.direction);
        }
      }
    }
     return photon;
  }

 /**
  * Run forward solver for single photon.
  *
  * This method calls the sample_photon function of the provided source and then
  * propagates the resulting photon through the atmosphere. One call to this function
  * corresponds to computing a single sample for the MC method.
  *
  * @param source Source object to provide the photon
  * @tparam Source The Source type.
  */
 template<typename Source>
 void forward(Source &source) {
    auto photon = source.sample_photon(generator_, atmosphere_);
    tracer_.created(photon);
    propagate_photon(photon);
  }
 /**
  * Run backward solver for single photon.
  *
  * This method calls the sample_photon function of the provided source and then
  * propagates the resulting photon through the atmosphere. One call to this function
  * corresponds to computing a single sample for the MC method.
  *
  * @param source Source object to provide the photon
  * @tparam Source The Source type.
  */

 template<typename Photon,
          typename Vector,
          typename Source>
 void backward(Vector position,
               Vector direction,
               Source &source) {
     Photon photon(atmosphere_.place_on_grid(position, direction));
     tracer_.created(photon);
     photon = propagate_photon(generator_, photon);
     photon.scale_energy(source.get_energy(photon));
     return photon;
 }

 private:

  Atmosphere atmosphere_;
  Tracer tracer_;
  Generator generator_;
  Float minimum_energy_ = 1e-6;

};

}  // namespace llrte
#endif
