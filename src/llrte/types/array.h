#ifndef _LLRTE_TYPES_ARRAY_H_
#define _LLRTE_TYPES_ARRAY_H_

#include <math.h>

#include <iostream>

namespace llrte {

template <typename F>
class Array {
 public:
  using Float = F;

  Array() : n_(0), data_(nullptr) {};
  Array(size_t n) : n_(n) { data_ = new Float[n_]; }

  Array(const Array& other) : n_(other.n_) {
    data_ = new Float[n_];
    std::copy(other.data_, other.data_ + n_, data_);
  }

  Array(Array&& other) : n_(other.n_), data_(other.data_) {
    other.data_ = nullptr;
  }

  Array& operator=(const Array& other) {
    n_ = other.n_;
    data_ = new Float[n_];
    std::copy(other.data_, other.data_ + n_, data_);
  }

  Array& operator=(Array&& other) {
    n_ = other.n_;
    data_ = other.data_;
    other.data_ = nullptr;
    return *this;
  }

  Float& operator[](size_t i) { return data_[i]; }
  Float operator[](size_t i) const { return data_[i]; }

  size_t size() const { return n_; }

  ~Array() {
    if (data_) {
      delete data_;
    }
  }

 private:
  size_t n_ = 0;
  F* data_ = nullptr;
};

template <typename Real>
std::ostream& operator<<(std::ostream& os, const Array<Real>& array) {
  os << "[";
  size_t n = array.size();
  for (size_t i = 0; i < n - 1; ++i) {
    os << array[i] << ",";
  }
  os << array[n - 1] << "]" << std::endl;
  return os;
}

}  // namespace llrte
#endif
