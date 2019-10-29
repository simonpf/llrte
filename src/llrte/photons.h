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

        while (true) {
            auto absorption = atmosphere.get_absorption(position);
            auto intersection = atmosphere.get_intersection(position, direction_);
            auto d = std::get<0>(intersection);

            auto p = generator.sample_path_length(1.0 / absorption);
            if (d > p) {
                Tracer::trace(position, Event::scattering);
                break;
            } else {
                Tracer::trace(position, Event::step);
                if (!atmosphere.is_inside(std::get<1>(intersection))) {
                    break;
                }
                position = std::get<1>(intersection);
            }
        }
    }


private:
    Vector position_;
    Vector direction_;
};

template <typename V>
class Ray {
public:

    using Vector = V;
    using Float = typename Vector::Float;

    Ray(Vector position,
        Vector direction)
    : position_(position), direction_(direction) {
        // Nothing to do here.
    }

    Vector get_position() {
        return position_;
    }

    Vector get_direction() {
        return direction_;
    }


private:

    Float optical_depth_;
    Vector position_;
    Vector direction_;
};

}
#endif
