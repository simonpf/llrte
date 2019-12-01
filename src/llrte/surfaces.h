#ifndef _LLRTE_SURFACES_
#define _LLRTE_SURFACES_

#include <iostream>

#include "llrte/configurations.h"
#include "llrte/tracers.h"
#include "llrte/constants.h"
#include "llrte/rotations.h"
#include "llrte/maths/geometry.h"

namespace llrte::surfaces {

template <typename Vector>
class Plane {
 protected:
  Plane(const Vector &base, const Vector &normal)
      : base_(base.normed()), normal_(normal.normed()) {
    // Nothing to do here.
  }

 public:
  template <typename Photon>
  bool has_crossed(const Photon &p) {
    using Float = typename Photon::Float;
    auto dir = p.get_direction();
    if (dot(dir, normal_) > 0.0) {
      return false;
    }
    Vector dv = p.get_position() - base_;
    return dot(dv, normal_) <= Limits<Float>::eps;
  }

 protected:
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

  template <typename Generator, typename Photon>
  void apply(Generator &/*g*/, Photon &p) {
    absorbed_energy_ += p.get_energy();
    p.set_energy(0.0);
  }

  Float get_absorbed_energy() const { return absorbed_energy_; }

 private:
  Float absorbed_energy_ = 0.0;
};

struct Specular {
  template <typename Generator, typename Vector>
  static Vector get_outgoing_direction(Generator &/*g*/,
                                       const Vector &d,
                                       const Vector &n) {
    using Float = typename Vector::Float;
    auto ns = n * dot(n, d);
    auto dn = n * dot(n, d) - d;
    return static_cast<Float>(-1.0) * (ns + dn);
  }
};

template <typename ReflectionPlane = maths::geometry::RandomPlane>
struct Lambertian {
  template <typename Generator, typename Vector>
  static Vector get_outgoing_direction(Generator &g,
                                       const Vector &/*d*/,
                                       const Vector &n) {
    using Float = typename Vector::Float;
    auto ns = ReflectionPlane::get_normal(g, n);
    auto phi = g.sample_angle_uniform() - Constants<Float>::pi / 2.0 ;
    auto dn = rotations::rotate(n, ns, phi);
    return dn;
  }
};

template <typename Vector, typename Type>
class ReflectingPlane : public Plane<Vector> {
 public:
  using Float = typename Vector::Float;
  using Plane<Vector>::normal_;
  using Plane<Vector>::base_;

  ReflectingPlane(const Vector &base, const Vector &normal, Float albedo)
      : Plane<Vector>(base, normal), absorbed_energy_(0.0), albedo_(albedo)  {
    // Nothing to do here.
  }

  template <typename Generator, typename Photon>
      void apply(Generator &g, Photon &p) {
    auto e = p.get_energy();
    absorbed_energy_ += (1.0 - albedo_) * e;
    p.set_energy(albedo_ * e);

    auto d = Type::get_outgoing_direction(g, p.get_direction(), normal_);
    p.set_direction(d);
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

  template <typename Photon>
  bool has_crossed(const Photon &p) {
    auto pos = p.get_position();
    if (dot(pos - base_1_, normal_) <= 0.0) {
      return true;
    }
    if (dot(pos - base_2_, normal_) >= 0.0) {
      return true;
    }
    return false;
  }

  template <typename Generator, typename Photon>
  void apply(Generator &/*g*/, Photon &p) {
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
