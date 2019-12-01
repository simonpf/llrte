#ifndef _LLRTE_PHOTONS_
#define _LLRTE_PHOTONS_

#include "llrte/tracers.h"

namespace llrte {

/**
 * Simple photon.
 *
 * Simple photon that propagates through space until it is
 * absorbed.
 */
template <typename V>
class Photon {
 public:
  using Vector = V;
  using Float = typename V::Float;

  Photon(Vector position, Vector direction)
      : position_(position),
        direction_(direction),
        n_scattered_(0),
        energy_(1.0) {
    // Nothing to do here.
  }

  template <typename Generator, typename PhaseFunction>
  void scatter(Generator& generator, const PhaseFunction& p) {
    p.scatter(generator, *this);
    n_scattered_++;
  }

  size_t get_scattering_events() const { return n_scattered_; }

  Float get_energy() const { return energy_; }
  void set_energy(Float e) { energy_ = e; }
  void scale_energy(Float f) { energy_ *= f; }

  const Vector& get_position() const { return position_; }
  Vector& get_position() { return position_; }
  void set_position(Vector p) { position_ = p; }

  const Vector& get_direction() const { return direction_; }
  void set_direction(Vector d) { direction_ = d; }

 private:
  Vector position_;
  Vector direction_;
  size_t n_scattered_ = 0;
  Float energy_ = 1.0;
};

template <typename Vector>
std::ostream& operator<<(std::ostream& os, const Photon<Vector>& p) {
    os << "Photon ::" << std::endl;
    os << "pos: " << p.get_position() << std::endl;
    os << "dir: " << p.get_direction() << std::endl;
    os << "e:   " << p.get_energy() << std::endl;
    return os;
}

}  // namespace llrte
#endif
