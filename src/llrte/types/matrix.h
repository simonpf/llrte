#ifndef _LLRTE_TYPES_MATRIX_H_
#define _LLRTE_TYPES_MATRIX_H_

#include <math.h>

#include <iostream>

namespace llrte {

/**
 * 3D Vector
 */
template <size_t M, size_t N, typename F>
class Matrix {
 public:
  using Float = F;

  Matrix() {
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        (*this)(i, j) = 0.0;
      }
    }
  }

  Matrix(const Matrix& other) {
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        (*this)(i, j) = other(i, j);
      }
    }
  }

  Matrix& operator=(const Matrix &other) {
      for (size_t i = 0; i < M; ++i) {
          for (size_t j = 0; j < N; ++j) {
              (*this)(i, j) = other(i, j);
          }
      }
      return *this;
  }

  Float operator()(size_t i, size_t j) const { return elements_[i * N + j]; }
  Float& operator()(size_t i, size_t j) { return elements_[i * N + j]; }

  Matrix operator+(const Matrix& b) const {
    Matrix c;
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        c(i, j) = (*this)(i, j) + b(i, j);
      }
    }
    return c;
  }

  Matrix operator-(const Matrix& b) const {
    Matrix c;
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < N; ++j) {
        c(i, j) = (*this)(i, j) - b(i, j);
      }
    }
    return c;
  }

  template <size_t O>
  Matrix operator*(const Matrix<N, O, F>& b) const {
    Matrix<M, O, F> c{};
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < O; ++j) {
        for (size_t k = 0; k < N; k++) {
          c(i, j) += (*this)(i, k) * b(k, j);
        }
      }
    }
    return c;
  }

  template <size_t O>
  Matrix operator&(const Matrix<O, N, F>& b) const {
    Matrix<M, O, F> c{};
    for (size_t i = 0; i < M; ++i) {
      for (size_t j = 0; j < O; ++j) {
        for (size_t k = 0; k < N; k++) {
          c(i, j) += (*this)(i, k) * b(j, k);
        }
      }
    }
    return c;
  }


 private:
  Float elements_[M * N];
};

template <typename F,
    template <typename> typename Vector>
Vector<F> operator*(const Matrix<3, 3, F> mat, const Vector<F>& b) {
    return Vector<F>{mat(0, 0) * b.x + mat(0, 1) * b.y + mat(0, 2) * b.z,
                     mat(1, 0) * b.x + mat(1, 1) * b.y + mat(1, 2) * b.z,
                     mat(2, 0) * b.x + mat(2, 1) * b.y + mat(2, 2) * b.z};
}

template <size_t M, size_t N, typename Float>
Matrix<M, N, Float> operator*(Float c, const Matrix<M, N, Float>& a) {
  Matrix<M, N, Float> b{};
  for (size_t i = 0; i < M; ++i) {
    for (size_t j = 0; j < N; ++j) {
      b(i, j) = c * a(i, j);
    }
  }
  return b;
}

template <typename Vector, size_t M, size_t N, typename Float>
Vector column(const Matrix<M, N, Float>& m, size_t j) {
  Vector v{};
  for (size_t i = 0; i < M; ++i) {
    v[i] = m(i, j);
  }
  return v;
}

template <size_t M, size_t N, typename Real>
std::ostream& operator<<(std::ostream& os, const Matrix<M, N, Real>& a) {
  os << "[";
  for (size_t i = 0; i < M - 1; ++i) {
    os << "[";
    for (size_t j = 0; j < N - 1; ++j) {
      os << a(i, j) << ",";
    }
    os << a(i, N - 1) << "]," << std::endl << " ";
  }
  os << "[";
  for (size_t j = 0; j < N - 1; ++j) {
    os << a(M - 1, j) << ",";
  }
  os << a(M - 1, N - 1) << "]]" << std::endl << " ";
  return os;
}

}  // namespace llrte
#endif
