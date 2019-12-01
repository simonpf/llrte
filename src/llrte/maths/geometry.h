#ifndef _LLRTE_MATHS_GEOMETRY_H_
#define _LLRTE_MATHS_GEOMETRY_H_

#include "llrte/configurations.h"
#include "llrte/rotations.h"

namespace llrte::maths::geometry {

template <typename Vector, size_t d>
Vector unit_vector() {
  Vector u{};
  u[d] = 1.0;
  return u;
}

template <typename U, typename V>
bool is_close(U u, V v) {
  auto d = u - v;
  auto eps = Limits<typename decltype(d)::Float>::eps;
  if (d.length() > eps) {
    return false;
  }
  return true;
}

template<typename Vector>
Vector perpendicular(const Vector &v) {
  Vector w;
  auto l = v.length();
  decltype(l) vx, vy, vz;
  vx = v[0] / l;
  vy = v[1] / l;
  vz = v[2] / l;

  if (abs(vx) > 0.5) {
    l = sqrt(vz * vz + vx * vx);
    w[0] = -vz / l;
    w[2] = vx / l;
  } else {
    l = sqrt(vy * vy + vz * vz);
    w[1] = vz / l;
    w[2] = - vy / l;
  }
  return w;
}

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
    b(2, 1) = - vy / l;
  }

  b(0, 2) = b(1, 0) * b(2, 1) - b(2, 0) * b(1, 1);
  b(1, 2) = b(2, 0) * b(0, 1) - b(0, 0) * b(2, 1);
  b(2, 2) = b(0, 0) * b(1, 1) - b(1, 0) * b(0, 1);

  return b;
}

template <size_t i>
struct FixedScatteringPlane {
  template<typename Generator, typename Vector>
  static Vector get_normal(Generator &g, const Vector &/*v*/) {
    using Float = typename Vector::Float;
    auto r = g.sample_uniform();
    Vector d = maths::geometry::unit_vector<Vector, i>();
    if (r > 0.5) {
      d = static_cast<Float>(-1.0) * d;
    }
    return d;
  }
};

struct RandomPlane {
  template<typename Generator, typename Vector>
  static Vector get_normal(Generator &g, const Vector &d) {
    auto n = maths::geometry::perpendicular(d);
    auto phi = g.sample_angle_uniform();
    n = rotations::rotate(n, d, phi);
    return n;
  }
};


}  // namespace llrte::maths::geometry

#endif
