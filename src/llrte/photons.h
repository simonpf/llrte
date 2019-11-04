#ifndef _LLRTE_PHOTONS_
#define _LLRTE_PHOTONS_

#include "llrte/tracers.h"

namespace llrte {}

/**
 * Simple photon.
 *
 * Simple photon that propagates through space until it is
 * absorbed.
 */
template <typename V, typename Tracer = NoTrace>
class Photon {
 public:
  using Vector = V;
  using Float = typename V::Float;

  Photon(Vector position, Vector direction)
 : position_(position), direction_(direction), n_scattered_(0) {
    // Nothing to do here.
  }

  size_t get_scattering_events() const {
     return n_scattered_;
  }


  template <typename Atmosphere, typename Random>
  void propagate(Atmosphere atmosphere, Random& generator) {
    auto position = atmosphere.get_grid_position(position_);
    auto tau = generator.sample_tau();
    Tracer::trace(*this, position, Event::step);

    while (true) {
      auto absorption_xc = atmosphere.get_absorption_coefficient(position);
      auto scattering_xc = atmosphere.get_scattering_coefficient(position);

      auto l = tau / (absorption_xc + scattering_xc);
      auto intersection = atmosphere.get_intersection(position, direction_, l);
      auto d = std::get<0>(intersection);
      tau -= d * (absorption_xc + scattering_xc);

      position = std::get<1>(intersection);

      // Check if left atmosphere
      if (!atmosphere.is_inside(position)) {
        Tracer::trace(*this, position, Event::left_domain);
        break;
      }

      // Scattering or absorption event.
      if (l <= d) {
        auto uniform = generator.sample_uniform();
        if (uniform < scattering_xc / (scattering_xc + absorption_xc)) {
          auto phase_function = atmosphere.get_phase_function(position);
          direction_ = phase_function.get_direction(generator, direction_);
          n_scattered_++;
          tau = generator.sample_tau();
          Tracer::trace(*this, position, Event::scattering);
        } else {
          Tracer::trace(*this, position, Event::absorption);
          break;
        }
      }
      Tracer::trace(*this, position, Event::step);
    }
  }

 private:
  Vector position_;
  Vector direction_;
  size_t n_scattered_ = 0;

};

/**
 * Simple photon.
 *
 * Simple photon that propagates through space until it is
 * absorbed.
 */
template <typename V,
          typename Tracer = NoTrace>
class FixedEnergyPhoton {
 public:
  using Vector = V;
  using Float = typename V::Float;

  FixedEnergyPhoton(Vector position, Vector direction)
 : position_(position), direction_(direction), n_scattered_(0), energy_(1.0) {
    // Nothing to do here.
  }

  size_t get_scattering_events() const {
     return n_scattered_;
  }


  template <typename Atmosphere, typename Random>
  void propagate(Atmosphere atmosphere, Random& generator) {
    auto position = atmosphere.get_grid_position(position_);
    auto tau = generator.sample_tau();

    while (true) {
      auto absorption_xc = atmosphere.get_absorption_coefficient(position);
      auto scattering_xc = atmosphere.get_scattering_coefficient(position);

      auto l = tau / scattering_xc;
      auto intersection = atmosphere.get_intersection(position, direction_, l);
      auto d = std::get<0>(intersection);
      tau -= d * scattering_xc;

      // Handle absorption.
      auto f_abs = exp(-absorption_xc * d);
      Tracer::trace(*this, position, energy_ * (1.0 - f_abs), Event::absorption);
      energy_ *= f_abs;
      if (energy_ < minimum_energy_) {
        break;
      }

      position = std::get<1>(intersection);

      // Check if left atmosphere.
      if (!atmosphere.is_inside(position)) {
        Tracer::trace(*this, position, energy_, Event::left_domain);
        break;
      }

      // Check if scattering event.
      if (l <= d) {
        auto phase_function = atmosphere.get_phase_function(position);
        direction_ = phase_function.get_direction(generator, direction_);
        n_scattered_++;
        tau = generator.sample_tau();
        Tracer::trace(*this, position, energy_, Event::scattering);
      }

      Tracer::trace(*this, position, Event::step);
    }
  }

 private:
  Vector position_;
  Vector direction_;
  size_t n_scattered_ = 0;
  Float minimum_energy_ = 1e-6;
  Float energy_ = 1.0;
};

// template <typename V>
// class BackwardsRay {
// public:
//
//    using Vector = V;
//    using Float = typename Vector::Float;
//
//    Ray(Vector position,
//        Vector direction)
//    : position_(position), direction_(direction) {
//        // Nothing to do here.
//    }
//
//    Vector get_position() {
//        return position_;
//    }
//
//    Vector get_direction() {
//        return direction_;
//    }
//
//    template<typename Atmosphere, typename Random>
//        void propagate(Float maximum_optical_depth,
//                       Atmosphere atmosphere,
//                       Random & generator) {
//
//        auto position = atmosphere.get_grid_position(position_);
//
//        auto absorption_xs_old = atmosphere.get_absorption(position);
//        auto scattering_xs_old = atmosphere.get_scattering(position);
//
//        while (true) {
//
//
//            auto p = generator.sample_path_length(1.0 / scattering_xs);
//            auto intersection = atmosphere.get_intersection(position,
//            direction_, p);
//
//            auto d = std::get<0>(intersection);
//            auto posistion_new = std::get<1>(intersection);
//            auto absorption_xs_new = atmosphere.get_absorption(position_new);
//            auto scattering_xs_new = atmosphere.get_scattering(position_new);
//
//            optical_depth += 0.5 * (absorption_xs_old + absorption_xs_new) *
//            d; absorption_xs_old = absorptions_xs_new; scattering_xs_old =
//            scattering_xs_new;
//
//            if (p > d) {
//                // No scattering event has occured
//            } else {
//                // Scattering event
//                auto phase_matrix = atmosphere.get_phase_matrix(position_new);
//                direction_ = phase_matrix.sample_direction(generator);
//            }
//            if (d > p) {
//                position = position +
//                Tracer::trace(position, Event::scattering);
//                break;
//            } else {
//                Tracer::trace(position, Event::step);
//                if (!atmosphere.is_inside(std::get<1>(intersection))) {
//                    break;
//                }
//                position = std::get<1>(intersection);
//            }
//        }
//    }
//};
//
// template <typename V>
// class OpticalDepthIntegral {
// public:
//
//    using Vector = V;
//    using Float = typename Vector::Float;
//
//    OpticalDepthIntegral(Vector start,
//                         Vector end)
//    : position_(position), direction_(direction) {
//        // Nothing to do here.
//    }
//
//    Vector get_position() {
//        return position_;
//    }
//
//    Vector get_direction() {
//        return direction_;
//    }
//
//    template<typename Atmosphere, typename Random>
//        void propagate(Float maximum_optical_depth,
//                       Atmosphere atmosphere,
//                       Random & generator) {
//
//        auto position = atmosphere.get_grid_position(position_);
//
//        auto absorption_xs_old = atmosphere.get_absorption(position);
//        auto scattering_xs_old = atmosphere.get_scattering(position);
//
//        while (true) {
//
//
//            auto p = generator.sample_path_length(1.0 / scattering_xs);
//            auto intersection = atmosphere.get_intersection(position,
//            direction_, p);
//
//            auto d = std::get<0>(intersection);
//            auto posistion_new = std::get<1>(intersection);
//            auto absorption_xs_new = atmosphere.get_absorption(position_new);
//            auto scattering_xs_new = atmosphere.get_scattering(position_new);
//
//            optical_depth += 0.5 * (absorption_xs_old + absorption_xs_new) *
//            d; absorption_xs_old = absorptions_xs_new; scattering_xs_old =
//            scattering_xs_new;
//
//            if (p > d) {
//                // No scattering event has occured
//            } else {
//                // Scattering event
//                auto phase_matrix = atmosphere.get_phase_matrix(position_new);
//                direction_ = phase_matrix.sample_direction(generator);
//            }
//            if (d > p) {
//                position = position +
//                Tracer::trace(position, Event::scattering);
//                break;
//            } else {
//                Tracer::trace(position, Event::step);
//                if (!atmosphere.is_inside(std::get<1>(intersection))) {
//                    break;
//                }
//                position = std::get<1>(intersection);
//            }
//        }
//    }
//};
//
//
// private:
//
//    Float optical_depth_;
//    Vector position_;
//    Vector direction_;
//
//};
}
#endif
