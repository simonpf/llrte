#ifndef _LLRTE_TYPES_VECTOR_H_
#define _LLRTE_TYPES_VECTOR_H_

#include <math.h>
#include <iostream>

#include "llrte/common.h"

namespace llrte {

/** 3D Vector
 *
 * The Vector3 class provides a simple implementation of 3-dimensional
 * vectors.
 */
template <typename F>
struct Vector3 {
    using Float = F;

    /* x-component of vector */
    Float x;
    /* y-component of vector */
    Float y;
    /* z-component of vector */
    Float z;

  /** Create vector with all components 0.0 */
  Vector3() {
      x = 0.0;
      y = 0.0;
      z = 0.0;
  }

/** Create a vector from given components
 * @param x_ The x component
 * @param y_ The y component
 * @param z_ The z component
 */
  Vector3(Float x_, Float y_, Float z_) : x(x_), y(y_), z(z_) {}

    /**
     * Add two vectors
     * @param other The vector to add to this
     * @return Vector representing the sum of the two.
     */
  __DEV__ Vector3 operator+(const Vector3 &other) const {
      return Vector3{x + other.x, y + other.y, z + other.z};
  }

    /**
     * Multiply two vectors
     * @param other The vector to multiply with this
     * @return Vector representing the product of the two.
     */
    __DEV__ Vector3 operator*(const Vector3 &other) const {
        return Vector3{x * other.x, y * other.y, z * other.z};
    }

    /**
     * Multiply two vectors
     * @param other The vector to multiply with this
     * @return Vector representing the product of the two.
     */
    __DEV__ Vector3 operator*(Float c) const {
        return Vector3{x * c, y * c, z * c};
    }

    /**
     * Subtract vector
     * @param other The vector to subtract
     * @return Vector representing the difference.
     */
    __DEV__ Vector3 operator-(const Vector3 &other) const {
        return Vector3{x - other.x, y - other.y, z - other.z};
    }

    /**
     * Vector equality
     * @param other The vector to compare with
     * @return true if all components are equal
     */
    __DEV__ Vector3 operator==(const Vector3 &other)  const {
        return (x == other.x) && (y == other.y) && (z == other.z);
    }

    /**
     * The length of the vector
     * @return The Euclidean norm of the vector
     */
  __DEV__ Float length() const {
      return sqrt(x * x + y * y + z * z);
  }

    /**
     * Norm vector
     * @return Vector pointing in same direction but with unit length.
     */
  __DEV__ Vector3 normed() const { return (*this) * (static_cast<Float>(1.0) / length()); }

};

/**
 * Dot product of two vectors
 * @param One vector
 * @param Other vector
 */
template <typename Float>
__DEV__ inline Float dot(Vector3<Float> v, Vector3<Float> w) {
    return v.x * w.x + v.y * w.y + v.z * w.z;
}

/**
 * Cross product of two vectors
 * @param One vector
 * @param Other vector
 */
template <typename Float>
__DEV__ inline Vector3<Float> cross(Vector3<Float> u, Vector3<Float> v) {
  return Vector3<Float>{
          u.y * v.z - u.z * v.y,
          u.z * v.x - u.x * v.z,
          u.x * v.y - u.y * v.x,
  };
}

/**
 * Scale vector
 * @param Scaling factor
 * @param Other vector
 */
template <typename Float>
__DEV__ Vector3<Float> operator*(Float c, Vector3<Float> v) {
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
