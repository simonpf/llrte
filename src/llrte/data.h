#ifndef _LLRTE_DATA_H_
#define _LLRTE_DATA_H_

#include <initializer_list>
#include "llrte/utils/array.h"
#include <iostream>


namespace llrte {

using utils::array::zip, utils::array::reduce, utils::array::Multiply, utils::array::Add;

template <typename T>
class Data {
public:
Data(size_t n) : data_(new T[n])
    {}

const T& operator()(size_t i) const { return data_[i]; }

T& operator()(size_t i) { return data_[i]; }

  ~Data() {
    if (data_) {
        delete[] data_;
    }
  }

 private:
  T* data_;
};

template <typename T, size_t rank>
class Tensor {
public:
    template <typename ... Ts>
    Tensor(Ts ... sizes)
        : shape_{sizes ...}, data_{reduce<Multiply>(shape_)} {
    // Nothing to do here.
  }

  template <typename ... Ts>
  T& operator()(Ts ... sizes) {
    std::array<size_t, rank> sizes_(sizes...);
    return data_(reduce<Add>(zip<Multiply>(shape_, sizes_)));
  }

  template <typename ... Ts>
  const T& operator()(Ts ... sizes) const {
    std::array<size_t, rank> sizes_(sizes...);
    return data_(reduce<Add>(zip<Multiply>(shape_, sizes_)));
  }

  void fill(const T &t) {
      for (size_t i = 0; i <  reduce<Multiply>(shape_); ++i) {
          data_[i] = t;
      }
  }

  private :
  std::array<size_t, rank> shape_;
  size_t bytes_;
  Data<T> data_;
};

template <typename T, size_t rank>
    std::ostream& operator<<(std::ostream& os, const Tensor<T, rank> &tensor) {
    std::cout << "Tensor(" << rank << ")" << std::endl;
}


}  // namespace llrte
#endif
