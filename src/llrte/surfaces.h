#ifndef _LLRTE_SURFACES_
#define _LLRTE_SURFACES_

#include <iostream>

#include "llrte/common.h"
#include "llrte/configurations.h"
#include "llrte/constants.h"
#include "llrte/geometry.h"
#include "llrte/rotations.h"
#include "llrte/tracers.h"

namespace llrte::surfaces {

/**ReflectingPlane
 *
 * Provides a base class for reflecting planes. A reflecting plane reflects a
 * photon in a given direction and may absorb parts of its energy.
 *
 * The Reflecting surface keeps track of the total energy it absorbs.
 *
 * @tparam Vector The vector type used to represent the position of the plane.
 * @tparam Type Type trait defining the reflection type.
 */
template <typename Vector, typename Type>
    class ReflectingPlane : public geometry::Plane<Vector> {
public:
    using Float = typename Vector::Float;
    using geometry::Plane<Vector>::has_crossed;
    using geometry::Plane<Vector>::normal_;
    using geometry::Plane<Vector>::base_;

    /**
     * Determine whether a photon has crossed the surface.
     */
    template <typename Photon>
        bool has_crossed(const Photon &photon) {
        std::cout << "has crossed " << photon << std::endl;
        return has_crossed(photon.position, photon.direction);
    }

    /**
    * Create reflecting plane with a given albedo.
    * @param base Base vector describing the position of any one point in the plane
    * @param normal The normal vector describing the direction of the plane
    * @param albedo The fraction of incoming radiative energy that is reflected.
    */
    ReflectingPlane(const Vector &base,
                    const Vector &normal,
                    Float albedo)
        : geometry::Plane<Vector>(base, normal), absorbed_energy_(0.0), albedo_(albedo) {
            // Nothing to do here.
        }

    /**
     * Apply reflection to photon.
     * @param generator Random number generator to use.
     * @param photon The photon being reflected.
     */
    template <typename Generator, typename Photon>
        void apply(Generator &generator, Photon &photon) {
        auto e = photon.get_energy();
        absorbed_energy_ += (1.0 - albedo_) * e;
        photon.set_energy(albedo_ * e);

        auto d = Type::get_outgoing_direction(generator, photon.direction, normal_);
        photon.change_direction(d);
    }

    /** Get total energy absorbed by surface.
     * @return The absorbed energy.
     */
    Float get_absorbed_energy() const {return absorbed_energy_; }

private:
    Float absorbed_energy_ = 0.0;
    Float albedo_ = 0.0;
};

/**
 * Type trait representing specular reflection.
 */
struct Specular {
  template <typename Generator, typename Vector>
  static Vector get_outgoing_direction(Generator & /*g*/,
                                       const Vector &d,
                                       const Vector &n) {
      std::cout << "dumb " << std::endl;
    using Float = typename Vector::Float;
    auto ns = n * dot(n, d);
    auto dn = n * dot(n, d) - d;
    return static_cast<Float>(-1.0) * (ns + dn);
  }
};

/**
 * Reflection in the backwards direction.
 */
struct BackwardsDirection {
  template <typename Generator, typename Vector>
  static Vector get_outgoing_direction(Generator & /*g*/,
                                       const Vector &d,
                                       const Vector & /*n*/) {
      std::cout << "wtf" << std::endl;
    using Float = typename Vector::Float;
    return static_cast<Float>(-1.0) * d;
  }
};

/**
 * Type trait representing lambertian reflection.
 * @tparam ReflectionPlane Type trait defining the plane in which reflection will
 * take place.
 */
template <typename ReflectionPlane = geometry::RandomPlane>
struct Lambertian {
  template <typename Generator, typename Vector>
  static Vector get_outgoing_direction(Generator &g, const Vector & /*d*/,
                                       const Vector &n) {
      std::cout << "lamb" << std::endl;
    using Float = typename Vector::Float;
    auto ns = ReflectionPlane::get_normal(g, n);
    auto phi = g.sample_angle_uniform() - Constants<Float>::pi / 2.0;
    auto dn = rotations::rotate(n, ns, phi);
    return dn;
  }
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
    auto pos = p.position;
    if (dot(pos - base_1_, normal_) <= 0.0) {
      return true;
    }
    if (dot(pos - base_2_, normal_) >= 0.0) {
      return true;
    }
    return false;
  }

  template <typename Generator, typename Photon>
  void apply(Generator & /*g*/, Photon &p) {
    using Float = typename Photon::Float;
    auto position = p.position;
    decltype(position) db{};
    if (dot(position - base_1_, normal_) <= 0.0) {
      db = base_2_ - base_1_;
    } else if (dot(position - base_2_, normal_) >= 0.0) {
      db = base_1_ - base_2_;
    }
    p.position = position + static_cast<Float>((1.0 - 1e-6) / normal_.length()) * dot(normal_, db) * normal_;
  }

 private:
  Vector base_1_, base_2_, normal_;
};

}  // namespace llrte::surfaces
#endif
