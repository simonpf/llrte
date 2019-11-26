#ifndef _LLRTE_MATHS_GEOMETRY_H_
#define _LLRTE_MATHS_GEOMETRY_H_

#include "llrte/configurations.h"

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
  std::cout << d.length() << std::endl;
  if (d.length() > eps) {
    return false;
  }
  return true;
}

template <typename Matrix, typename Vector>
Matrix orthonormal_basis(const Vector &v) {

  Matrix b{};

  auto l = v.length();
  decltype(l) vx, vy, vz;
  vx = v[0] / l;
  vy = v[1] / l;
  vz = v[2] / l;

  b(0, 0) = vx;
  b(1, 0) = vy;
  b(2, 0) = vz;

  if (vx > 0.5) {
    l = sqrt(vz * vz + vx * vx);
    b(0, 1) = -vz / l;
    b(2, 1) = vx / l;
    l = sqrt(vx * vx + vy * vy);
    b(0, 2) = vy / l;
    b(1, 2) = vx / l;
  } else if (vy > 0.5) {
    l = sqrt(vx * vx + vy * vy);
    b(0, 1) = vy / l;
    b(1, 1) = vx / l;
    l = sqrt(vy * vy + vz * vz);
    b(1, 2) = vz / l;
    b(2, 2) = -vy / l;
  } else if (vz > 0.5) {
    l = sqrt(vy * vy + vz * vz);
    b(1, 1) = vz / l;
    b(2, 1) = - vy / l;
    l = sqrt(vx * vx + vz * vz);
    b(0, 2) = - vz / l;
    b(2, 2) = vx / l;
  }
  return b;
}

}  // namespace llrte::maths::geometry

#endif
