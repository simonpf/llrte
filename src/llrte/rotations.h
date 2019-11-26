#ifndef _LLRTE_ROTATIONS_H_
#define _LLRTE_ROTATIONS_H_

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
    v[0] = ct * u[0] + st * u[2];
    v[1] = u[1];
    v[2] = -st * u[0] + ct * u[2];
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

  template <typename Vector, typename Float>
  auto rotate_z(const Vector &u, Float theta) -> Vector {
    auto v = Vector{};
    auto ct = cos(theta);
    auto st = sin(theta);
    v[0] = ct * u[0] - st * u[1];
    v[1] = st * u[0] + ct * u[1];
    v[2] = u[2];
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
