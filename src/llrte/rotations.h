#ifndef _LLRTE_ROTATIONS_H_
#define _LLRTE_ROTATIONS_H_

#include "llrte/types/matrix.h"

namespace llrte::rotations {

  template <template <size_t, size_t, typename F> typename Matrix, typename F>
  Matrix<3, 3, F> make_rotation_x(F theta) {
    Matrix<3, 3, F> r{};
    r(0, 0) = 1.0;
    r(1, 1) = cos(theta);
    r(1, 2) = -sin(theta);
    r(2, 1) = sin(theta);
    r(2, 2) = cos(theta);
    return r;
  }

  template <template <size_t, size_t, typename F> typename Matrix, typename F>
  Matrix<3, 3, F> make_rotation_y(F theta) {
    Matrix<3, 3, F> r{};
    r(0, 0) = cos(theta);
    r(0, 2) = sin(theta);
    r(1, 1) = 1.0;
    r(2, 0) = -sin(theta);
    r(2, 2) = cos(theta);
    return r;
  }

  template <typename Vector, typename Float>
  auto rotate_y(const Vector &u, Float theta) -> Vector {
    auto v = Vector{};
    auto ct = cos(theta);
    auto st = sin(theta);
    v.x = ct * u.x + st * u.z;
    v.y = u.y;
    v.z = -st * u.x + ct * u.z;
    return v;
  }

  template <template <size_t, size_t, typename F> typename Matrix, typename F>
  Matrix<3, 3, F> make_rotation_z(F theta) {
    Matrix<3, 3, F> r{};
    r(0, 0) = cos(theta);
    r(0, 1) = -sin(theta);
    r(1, 0) = sin(theta);
    r(1, 1) = cos(theta);
    r(2, 2) = 1.0;
    return r;
  }

  template <typename Matrix, typename Vector>
  Matrix make_rotation(const Vector &v,
                       typename Vector::Float  theta) {
    auto ct = cos(theta);
    auto st = sin(theta);
    Matrix r{};

    r(0, 0) = ct + v.x * v.x * (1.0 - ct);
    r(0, 1) = v.x * v.y * (1.0 - ct) - v.z * st;
    r(0, 2) = v.x * v.z * (1.0 - ct) + v.y * st;

    r(1, 0) = v.y * v.x * (1.0 - ct) + v.z * st;
    r(1, 1) = ct + v.y * v.y * (1.0 - ct);
    r(1, 2) = v.y * v.z * (1.0 - ct) - v.x * st;

    r(2, 0) = v.z * v.x * (1 - ct) - v.y * st;
    r(2, 1) = v.z * v.y * (1 - ct) + v.x * st;
    r(2, 2) = ct + v.z * v.z * (1 - ct);

    return r;
  }

  template <typename Vector>
  Vector rotate(const Vector &v,
                const Vector &n,
                typename Vector::Float  theta) {
    using Float = typename Vector::Float;
    return make_rotation<Matrix<3, 3, Float>>(n, theta) * v;
  }

  template <typename Vector, typename Float>
  auto rotate_z(const Vector &u, Float theta) -> Vector {
    auto v = Vector{};
    auto ct = cos(theta);
    auto st = sin(theta);
    v.x = ct * u.x - st * u.y;
    v.y = st * u.x + ct * u.y;
    v.z = u.z;
    return v;
  }

  template <typename Vector, typename Float>
  Vector rotate(const Vector &v, Float phi, Float theta) {
    return rotate_z(rotate_y(v, theta), phi);
  }

  template <typename Matrix, typename Vector, typename Float>
  Vector rotate(const Matrix reference_frame, const Vector &v, Float phi, Float theta) {
    return reference_frame & rotate_z(rotate_y(reference_frame * v, theta), phi);
  }

} // namespace llrte::rotations
#endif
