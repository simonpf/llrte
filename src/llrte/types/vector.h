#ifndef _LLRTE_TYPES_VECTOR_H_
#define _LLRTE_TYPES_VECTOR_H_

#include <math.h>
#include <iostream>

namespace llrte {

/**
 * 3D Vector
 */
template <typename F>
struct Vector3 {
    using Float = F;

    Float x;
    Float y;
    Float z;

  Vector3() {
      x = 0.0;
      y = 0.0;
      z = 0.0;
  }

  Vector3(Float x_, Float y_, Float z_) : x(x_), y(y_), z(z_) {}

  Vector3 operator+(const Vector3 &other) const {
      return Vector3{x + other.x, y + other.y, z + other.z};
  }

    Vector3 operator*(const Vector3 &other) const {
        return Vector3{x * other.x, y * other.y, z * other.z};
    }

    Vector3 operator-(const Vector3 &other) const {
        return Vector3{x - other.x, y - other.y, z - other.z};
    }

    Vector3 operator==(const Vector3 &other)  const {
        return (x == other.x) && (y == other.y) && (z == other.z);
    }


  Float length() const {
      return sqrt(x * x + y * y + z * z);
  }

  Vector3 normed() const { return (*this) * (1.0 / length()); }

};

template <typename Float>
inline Float dot(Vector3<Float> v, Vector3<Float> w) {
    return v.x * w.x + v.y * w.y + v.z * w.z;
}

template <typename Float>
inline Vector3<Float> cross(Vector3<Float> u, Vector3<Float> v) {
  return Vector3<Float>{
          u.y * v.z - u.z * v.y,
          u.z * v.x - u.x * v.z,
          u.x * v.y - u.y * v.x,
  };
}

template <typename Float>
Vector3<Float> operator*(Float c, Vector3<Float> v) {
    return Vector3<Float>{c * v.x, c * v.y, c * v.z};
}

template <typename Real>
std::ostream& operator<<(std::ostream& os, const Vector3<Real>& v) {
  os << "[";
  os << v.x << ", ";
  os << v.y << ", ";
  os << v.z << "]";
  return os;
}

}  // namespace llrte
#endif
