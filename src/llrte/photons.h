#ifndef _LLRTE_PHOTONS_
#define _LLRTE_PHOTONS_

#include "llrte/tracers.h"

namespace llrte {

}

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

    Photon(Vector position,
           Vector direction)
    : position_(position), direction_(direction) {
        // Nothing to do here.
    }

template<typename Atmosphere, typename Random>
void propagate(Atmosphere atmosphere, Random & generator) {

        auto position = atmosphere.get_grid_position(position_);
        auto tau = generator.sample_tau();

        while (true) {
            auto absorption_xc = atmosphere.get_absorption_coefficient(position);
            auto scattering_xc = atmosphere.get_scattering_coefficient(position);

            auto l = tau / (absorption_xc + scattering_xc);
            auto intersection = atmosphere.get_intersection(position, direction_, l);
            auto d = std::get<0>(intersection);

            // Scattering event.
            if (l <= d) {
                auto uniform = generator.sample_uniform();
                if (uniform < scattering_xc / (scattering_xc + absorption_xc)) {
                    auto phase_function = atmosphere.get_phase_function(position);
                    direction_ = phase_function.get_direction(generator, direction_);
                    Tracer::trace(position, Event::scattering);
                } else {
                    Tracer::trace(position, Event::absorption);
                    break;
                }
            // Stepping on.
            } else {
                if (!atmosphere.is_inside(std::get<1>(intersection))) {
                    Tracer::trace(position, Event::left_domain);
                    break;
                }
                Tracer::trace(position, Event::step);
                position = std::get<1>(intersection);
                tau -= d * (absorption_xc + scattering_xc);
            }
        }
    }


private:
    Vector position_;
    Vector direction_;
};

//template <typename V>
//class BackwardsRay {
//public:
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
//            auto intersection = atmosphere.get_intersection(position, direction_, p);
//
//            auto d = std::get<0>(intersection);
//            auto posistion_new = std::get<1>(intersection);
//            auto absorption_xs_new = atmosphere.get_absorption(position_new);
//            auto scattering_xs_new = atmosphere.get_scattering(position_new);
//
//            optical_depth += 0.5 * (absorption_xs_old + absorption_xs_new) * d;
//            absorption_xs_old = absorptions_xs_new;
//            scattering_xs_old = scattering_xs_new;
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
//template <typename V>
//class OpticalDepthIntegral {
//public:
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
//            auto intersection = atmosphere.get_intersection(position, direction_, p);
//
//            auto d = std::get<0>(intersection);
//            auto posistion_new = std::get<1>(intersection);
//            auto absorption_xs_new = atmosphere.get_absorption(position_new);
//            auto scattering_xs_new = atmosphere.get_scattering(position_new);
//
//            optical_depth += 0.5 * (absorption_xs_old + absorption_xs_new) * d;
//            absorption_xs_old = absorptions_xs_new;
//            scattering_xs_old = scattering_xs_new;
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
//private:
//
//    Float optical_depth_;
//    Vector position_;
//    Vector direction_;
//
//};

}
#endif
