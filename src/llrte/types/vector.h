#ifndef _LLRTE_TYPES_VECTOR_H_
#define _LLRTE_TYPES_VECTOR_H_

#include <math.h>
#include <iostream>

namespace llrte {

/**
 * 3D Vector
 */
template <size_t N, typename F>
class Vector {
 public:
  using Float = F;

  Vector() {
    for (size_t i = 0; i < N; ++i) {
      elements_[i] = 0.0;
    }
  }

  Vector(const std::array<F, 3> &data) {
      for (size_t i = 0; i < 3; ++i) {
          elements_[i] = data[i];
      }
  }

  Vector (const Vector& other) {
      for (size_t i = 0; i < 3; ++i) {
          elements_[i] = other.elements_[i];
      }
  }

  Vector& operator=(const Vector &other) {
      for (size_t i = 0; i < 3; ++i) {
          elements_[i] = other.elements_[i];
      }
      return *this;
  }

  Float operator[](size_t i) const { return elements_[i]; }
  Float& operator[](size_t i) { return elements_[i]; }

  Vector operator+(const Vector& v) const {
    Vector w;
    for (size_t i = 0; i < N; ++i) {
      w[i] = v[i] + elements_[i];
    }
    return w;
  }

  Vector operator*(const Float& v) const {
    Vector w;
    for (size_t i = 0; i < N; ++i) {
      w[i] = v * elements_[i];
    }
    return w;
  }

  Vector operator-(const Vector& v) const {
    Vector w;
    for (size_t i = 0; i < N; ++i) {
      w[i] = elements_[i] - v[i];
    }
    return w;
  }

  Vector operator==(const Vector& v) const {
    bool same = false;
    for (size_t i = 0; i < 3; ++i) {
      if (v[i] != elements_[i]) {
        return false;
      }
    }
    return true;
  }

  Float length() const {
    Float s = 0.0;
    for (size_t i = 0; i < N; ++i) {
      s += elements_[i] * elements_[i];
    }
    return sqrt(s);
  }

  Vector normed() const { return (*this) * (1.0 / length()); }

 public:
  Float elements_[N];
};

template <size_t N, typename Float>
Float dot(Vector<N, Float> v, Vector<N, Float> w) {
  Float d = 0.0;
  for (size_t i = 0; i < N; ++i) {
    d += v[i] * w[i];
  }
  return d;
}

template <typename Float>
Vector<3, Float> cross(Vector<3, Float> u, Vector<3, Float> v) {
  Vector<3, Float> w{};
  w[0] = u[1] * v[2] - u[2] * v[1];
  w[1] = u[2] * v[0] - u[0] * v[2];
  w[2] = u[0] * v[1] - u[1] * v[0];
  return w;
}

template <typename Vector>
    typename Vector::Float angle(const Vector &v1, const Vector &v2) {
    auto theta = dot(v1, v2) / v1.length() / v2.length();
    return acos(theta);
}

template <size_t N, typename Float>
Vector<N, Float> operator*(Float c, Vector<N, Float> v) {
  Vector<N, Float> w{};
  for (size_t i = 0; i < N; ++i) {
    w[i] = c * v[i];
  }
  return w;
}

template <size_t N, typename Real>
std::ostream& operator<<(std::ostream& os, const Vector<N, Real>& v) {
  os << "[";
  for (size_t i = 0; i < N - 1; ++i) {
    os << v[i] << ",";
  }
  os << v[N - 1] << "]" << std::endl;
  return os;
}

}  // namespace llrte
#endif
