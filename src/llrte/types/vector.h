#ifndef _LLRTE_TYPES_VECTOR_H_
#define _LLRTE_TYPES_VECTOR_H_

#include <iostream>
#include <math.h>

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

  Float length() {
    Float s = 0.0;
    for (size_t i = 0; i < N; ++i) {
      s += elements_[i] * elements_[i];
    }
    return sqrt(s);
  }

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
