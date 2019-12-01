#ifndef _LLRTE_SOURCES_H_
#define _LLRTE_SOURCES_H_

#include <random>

#include "llrte/constants.h"

namespace llrte {

/**
 * Isotropic point source
 *
 * As the name suggests this point source emits photons in all direction
 * with equal probability.
 */
template <typename P>
class PointSource {
 public:
  using Photon = P;
  using Vector = typename Photon::Vector;
  using Float = typename Vector::Float;
  using C = Constants<Float>;

  PointSource(Vector position) : position_(position) {}

  Photon sample_photon() {
    auto theta = theta_d_(generator_);
    auto phi = asin(2.0 * phi_d_(generator_) - 1.0);
    Vector v = Vector{};
    v[0] = cos(phi) * cos(theta);
    v[1] = cos(phi) * sin(theta);
    v[2] = sin(phi);
    return Photon(position_, v);
  }

 private:
  Vector position_;

  std::default_random_engine generator_;
  std::uniform_real_distribution<Float> theta_d_{-C::pi, C::pi};
  std::uniform_real_distribution<Float> phi_d_{0, 1};
};

/**
 * A beam of photons.
 *
 * Beam source that emits a perfectly collimated beam of photons.
 */
template <typename P>
class BeamSource {
 public:
  using Photon = P;
  using Vector = typename Photon::Vector;
  using Float = typename Vector::Float;
  using C = Constants<Float>;

  BeamSource(Vector position, Vector direction)
      : position_(position), direction_(direction) {
    // Nothing to do here.
  }

  void set_direction(const Vector &v) { direction_ = v; }
  const Vector & get_direction(const Vector &v) const { direction_ = v; }

  template <typename Generator>
  Photon sample_photon(Generator &) { return Photon(position_, direction_); }

 private:
  Vector position_, direction_;
};

template <typename Source>
class RandomOffset : public Source {
public:

  using Float   = typename Source::Float;
  using Vector = typename Source::Vector;
  using Photon = typename Source::Photon;

  template <typename ... Ts>
  RandomOffset(Vector direction,
               Float limit_low,
               Float limit_high,
               Ts ... ts)
  : Source(ts ...), direction_(direction), limit_low_(limit_low), limit_high_(limit_high)
  {}

  template <typename Generator>
  Photon sample_photon(Generator &g) {
    auto p = Source::sample_photon(g);
    auto s = limit_low_ + (limit_high_ - limit_low_) * g.sample_uniform();
    auto pp = p.get_position();
    p.set_position(pp + s * direction_);
    return p;
    }

 private:
  Vector direction_;
  Float limit_low_ = 0.0;
  Float limit_high_ = 0.0;
};

}  // namespace llrte
#endif
