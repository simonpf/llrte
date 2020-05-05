#ifndef _LLRTE_MATHS_GEOMETRY_H_
#define _LLRTE_MATHS_GEOMETRY_H_

#include "llrte/configurations.h"
#include "llrte/maths.h"
#include "llrte/rotations.h"

namespace llrte::geometry {

////////////////////////////////////////////////////////////////////////////////
// Vectors and matrices
////////////////////////////////////////////////////////////////////////////////

/**
 * Return unit vector.
 * @tparam Vector The vector type
 * @tparam The index of the coordinate axis
 */
template <typename Vector, size_t d>
Vector unit_vector() {
  if (d == 0) {
    return Vector{1.0, 0.0, 0.0};
  }
  if (d == 1) {
    return Vector{0.0, 1.0, 0.0};
  }
  if (d == 2) {
    return Vector{0.0, 0.0, 1.0};
  }
}

/**
 * Determine whether two vectors are close. Vectors
 * are considered close when the Euclidean norm of
 * the difference is small.
 *
 * @param u The first vector
 * @param y The second vector
 */
template <typename U, typename V>
bool is_close(U u, V v) {
  auto d = u - v;
  auto eps = Limits<typename decltype(d)::Float>::eps;
  return maths::small(d.length());
}

/**
 * Compute perpendicular vector.
 * @param v The vector to which the returned vector will be perpendicular.
 */
template <typename Vector>
Vector perpendicular(const Vector &v) {
  Vector w;
  auto l = v.length();
  decltype(l) vx, vy, vz;
  vx = v.x / l;
  vy = v.y / l;
  vz = v.z / l;

  if (abs(vx) > 0.5) {
    l = sqrt(vz * vz + vx * vx);
    w.x = -vz / l;
    w.z = vx / l;
  } else {
    l = sqrt(vy * vy + vz * vz);
    w.y = vz / l;
    w.z = -vy / l;
  }
  return w;
}

/**
 * Compute angle between two vectors in radians.
 * @param v1 The first vector
 * @param v2 The other vector
 */
template <typename Vector>
typename Vector::Float angle(const Vector &v1, const Vector &v2) {
  auto theta = dot(v1, v2) / v1.length() / v2.length();
  return acos(theta);
}

/**
 * Return matrix containing an orthonormal basis that includes
 * the vector.
 * @param v The vector around which to form the orthonormal basis.
 */
template <typename Matrix, typename Vector>
Matrix orthonormal_basis(const Vector &v) {
  Matrix b{};

  auto l = v.length();
  decltype(l) vx, vy, vz, vmax;
  vx = v[0] / l;
  vy = v[1] / l;
  vz = v[2] / l;

  vmax = vx;
  if (abs(vy) > abs(vmax)) vmax = vy;
  if (abs(vz) > abs(vmax)) vmax = vz;

  if (vmax < 0.0) {
    vx *= -1.0;
    vy *= -1.0;
    vz *= -1.0;
  }

  b(0, 0) = vx;
  b(1, 0) = vy;
  b(2, 0) = vz;

  if (abs(vx) > 0.5) {
    l = sqrt(vz * vz + vx * vx);
    b(0, 1) = -vz / l;
    b(2, 1) = vx / l;
  } else {
    l = sqrt(vy * vy + vz * vz);
    b(1, 1) = vz / l;
    b(2, 1) = -vy / l;
  }

  b(0, 2) = b(1, 0) * b(2, 1) - b(2, 0) * b(1, 1);
  b(1, 2) = b(2, 0) * b(0, 1) - b(0, 0) * b(2, 1);
  b(2, 2) = b(0, 0) * b(1, 1) - b(1, 0) * b(0, 1);

  return b;
}

////////////////////////////////////////////////////////////////////////////////
// Planes
////////////////////////////////////////////////////////////////////////////////

/** Generic plane
 *
 * Generic representation of planes using normal and base vector. The plane
 * is assumed to be directed, with the front side pointing in the direction
 * of the normal vector.
 * @tparam Vector The vector type to use to represent the plane.
 */
template <typename Vector>
class Plane {
 public:
  using Float = Vector::Float;
  /** Create a plane.
   * @base The base vector of the plane
   * @base The normal vector of the plane.
   */
  Plane(const Vector &base, const Vector &normal)
      : base_(base), normal_(normal.normed()) {
    // Nothing to do here.
  }

  /**
   * Determine which side of the plane a position lies.
   * @return The returned float will be 1.0 if the position lies on the
   * front side of the plane and -1.0 otherwise.
   */
  Float get_side(const Vector &position) {
    Vector d = position - base_;
    if (dot(d, normal_) > 0.0) {
      return 1.0;
    }
    return -1.0;
  }

  /**
   * Determine whether an object has crossed a boundary. An object
   * is considered to have crossed the plane when it is on or behind
   * the plane and moving in the opposite direction of the normal vector.
   *
   * @position position Current-position vector of the object.
   * @position direction Direction of the objectk
   */
  bool has_crossed(const Vector &position, const Vector &direction) {
    if (dot(direction, normal_) >= 0.0) {
      return false;
    }
    return get_side(position) < 0.0;
  }

 protected:
  Vector base_;
  Vector normal_;
};

/** Random scattering plane
 * Type trait representing a random scattering plane. This can be in scattering
 * models to extend scattering into a third dimension.
 */
struct RandomPlane {
  /**
   * Returns normal vector of a random scattering plane which includes
   * the vector d.
   * @param g Random number generator to use.
   * @param d Orientation of the reference frame, which the scattering plane
   * should include.
   */
  template <typename Generator, typename Vector>
  static Vector get_normal(Generator &g, const Vector &d) {
    auto n = geometry::perpendicular(d);
    auto phi = g.sample_angle_uniform();
    n = rotations::rotate(n, d, phi);
    return n;
  }
};

/** Fixed scattering plane
 *
 * Restrict scattering plane to a given coordinate plane.
 * @tparam size_t Index of the coordinate axis taken as normal to the scattering
 * plane
 */
template <size_t i>
struct FixedScatteringPlane {
  template <typename Generator, typename Vector>
  static Vector get_normal(Generator &g, const Vector & /*v*/) {
    using Float = typename Vector::Float;
    auto r = g.sample_uniform();
    Vector d = geometry::unit_vector<Vector, i>();
    if (r > 0.5) {
      d = static_cast<Float>(-1.0) * d;
    }
    return d;
  }
};

}  // namespace llrte::geometry

#endif
