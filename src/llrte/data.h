#ifndef _LLRTE_DATA_H_
#define _LLRTE_DATA_H_

#include <iostream>

#include "llrte/utils/array.h"

namespace llrte {

using utils::array::zip, utils::array::reduce, utils::array::index,
    utils::array::Multiply, utils::array::Add;

//******************************************************************************
// Array data
//******************************************************************************
/**
 * \brief Array data
 *
 * The Data class provides linear array-type storage for a given
 * data type. It should be used to store "heavy-weight" data.
 * A Data object owns the data it contains and will free all its
 * memory when destroyed.
 *
 * \tparam T The type of the data to store.
 */
template <typename T>
class Data {
 public:
  /** Create new array holding n elements. */
  Data(size_t n) : data_(new T[n]) {}
    /** Destroys array and frees memory. */
    ~Data() {if (data_) {
            delete[] data_;
        }
    }
  /** Access element i. */
  const T& operator[](size_t i) const { return data_[i]; }
  /** Access element i. */
  T& operator[](size_t i) { return data_[i]; }
  /** Get pointer to raw data.*/
  const T *get_data_pointer() const {
      return data_;
  }
  T *get_data_pointer() {
      return data_;
  }

 private:
  T* data_;
};

//******************************************************************************
// Tensors
//******************************************************************************
/**
 * \brief Multi-dimensional data
 *
 * The Tensor class provides functionality for storing and manipulating
 * multidimensional, gridded data.
 *
 * \tparam T The type of the data to store.
 * \tparam rank The rank k of the tensor
 */
template <typename T, size_t rank>
class Tensor {
 public:
  template <typename... Ts>
  Tensor(Ts... sizes) : shape_{sizes...}, data_{reduce<Multiply>(shape_)} {
    // Nothing to do here.
  }

  const std::array<size_t, rank> & shape() const {
      return shape_;
  }

  template <typename... Ts>
  T& operator()(Ts... indices) {
    return data_[index({indices...}, shape_)];
  }

  template <typename... Ts>
  const T& operator()(Ts... indices) const {
    return data_[index({indices...}, shape_)];
  }

  void fill(const T& t) {
    for (size_t i = 0; i < reduce<Multiply>(shape_); ++i) {
      data_[i] = t;
    }
  }

  const T *get_data_pointer() const {
      return data_.get_data_pointer();
  }

   T *get_data_pointer() {
      return data_.get_data_pointer();
  }

  friend std::ostream& operator<<(std::ostream& os, const Tensor& tensor) {
    if (rank > 2) {
      std::cout << "Tensor(" << tensor.shape_ << ")" << std::endl;
    } else if (rank == 2) {
      std::cout << "[";
      for (size_t i = 0; i < tensor.shape_[0] - 1; ++i) {
        std::cout << "[";
        for (size_t j = 0; j < tensor.shape_[1] - 1; ++j) {
          std::cout << tensor(i, j) << ", ";
        }
        std::cout << tensor(i, tensor.shape_[1] - 1) << "]," << std::endl
                  << " ";
      }
      std::cout << "[";
      for (size_t j = 0; j < tensor.shape_[1] - 1; ++j) {
        std::cout << tensor(tensor.shape_[0] - 1, j) << ", ";
      }
      std::cout << tensor(tensor.shape_[0] - 1, tensor.shape_[1] - 1) << "]]"
                << std::endl;
    } else if (rank == 1) {
      std::cout << "[";
      for (size_t i = 0; i < tensor.shape_[0] - 1; ++i) {
        std::cout << tensor(i) << ", ";
      }
      std::cout << tensor(tensor.shape_[0] - 1) << "]" << std::endl;
    }
    return os;
  }

 private:
  std::array<size_t, rank> shape_;
  size_t bytes_;
  Data<T> data_;
};
}  // namespace llrte
#endif
