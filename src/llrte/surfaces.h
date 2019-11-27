#ifndef _LLRTE_SURFACES_
#define _LLRTE_SURFACES_

#include <iostream>
#include "llrte/tracers.h"

namespace llrte::surfaces {

template <typename Vector>
class Plane {
 protected:
  Plane(const Vector &base, const Vector &normal)
      : base_(base.normed()), normal_(normal.normed()) {
    // Nothing to do here.
  }

 public:
  bool has_crossed(const Vector &v) {
    Vector dv = v - base_;
    return dot(dv, normal_) <= 0.0;
  }

 private:
  Vector base_;
  Vector normal_;
};

template <typename Vector>
class BlackPlane : public Plane<Vector> {
 public:
  using Float = typename Vector::Float;
  BlackPlane(const Vector &base, const Vector &normal)
      : Plane<Vector>(base, normal) {
    // Nothing to do here.
  }

  template <typename Photon>
  void apply(Photon &p) {
    absorbed_energy_ += p.get_energy();
    p.set_energy(0.0);
  }

  Float get_absorbed_energy() const { return absorbed_energy_; }

 private:
  Float absorbed_energy_ = 0.0;
};

template <typename Vector>
class AbsorbingPlane : public Plane<Vector> {
 public:
  using Float = typename Vector::Float;
  AbsorbingPlane(const Vector &base, const Vector &normal, Float albedo)
      : albedo_(albedo), Plane<Vector>(base, normal) {
    // Nothing to do here.
  }

  template <typename Photon>
  void apply(Photon &p) {
    auto e = p.get_energy();
    absorbed_energy_ += (1.0 - albedo_) * e;
    p.set_energy(albedo_ * e);

    auto d = p.get_direction();
    auto dn = normal_ * dot(normal_, d) - d;
    p.set_direction(normal_ + dn);
  }

  Float get_absorbed_energy() const { return absorbed_energy_; }

 private:
  Float absorbed_energy_ = 0.0;
  Float albedo_ = 0.0;
};

template <typename Vector>
class PeriodicBoundary {
 public:
  using Float = typename Vector::Float;
  PeriodicBoundary(const Vector &base_1, const Vector &base_2,
                   const Vector &normal)
      : base_1_(base_1),
        base_2_(base_2),
        normal_(dot(normal, base_2 - base_1) > 0.0
                    ? normal
                    : static_cast<Float>(-1.0) * normal) {
    // Nothing to do here.
  }

  bool has_crossed(const Vector &v) {
    if (dot(v - base_1_, normal_) <= 0.0) {
      return true;
    }
    if (dot(v - base_2_, normal_) >= 0.0) {
      return true;
    }
    return false;
  }

  template <typename Photon>
  void apply(Photon &p) {
    using Float = typename Photon::Float;
    auto position = p.get_position();
    decltype(position) db{};
    if (dot(position - base_1_, normal_) <= 0.0) {
      db = base_2_ - base_1_;
    } else if (dot(position - base_2_, normal_) >= 0.0) {
      db = base_1_ - base_2_;
    }
    p.set_position(position +
                   static_cast<Float>((1.0 - 1e-6) / normal_.length()) *
                       dot(normal_, db) * normal_);
  }

 private:
  Vector base_1_, base_2_, normal_;
};

}  // namespace llrte::surfaces
#endif
