#ifndef _LLRTE_DATA_H_
#define _LLRTE_DATA_H_

#include <iostream>
#include "assert.h"

#include "llrte/utils/array.h"

namespace llrte {

using utils::array::zip, utils::array::reduce, utils::array::index,
    utils::array::Multiply, utils::array::Add, utils::array::tail;

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
  Data(size_t n) : data_(new T[n]), size_(n), owner_(true) {
    }
  Data(T* ptr, size_t n) : data_(ptr), size_(n), owner_(false) {
    }

  Data(const Data& other)
      : data_(new T[other.size_]), size_(other.size_), owner_(true) {
    std::copy(other.data_, other.data_ + size_, data_);
  }

  Data(Data&& other)
      : data_(other.data_), size_(other.size_), owner_(other.owner_) {
    other.data_ = nullptr;
    other.owner_ = false;
    other.size_ = 0;
  }

  Data& operator=(const Data& other) {
    if (data_ && owner_) delete[] data_;
    data_ = new T[size_];
    size_ = other.size_;
    owner_ = true;
    std::copy(other.data_, other.data_ + size_, data_);
  }

  Data& operator=(Data&& other) {
    if (data_ && owner_) delete[] data_;
    data_ = other.data_;
    size_ = other.size_;
    owner_ = other.owner_;
    other.data_ = nullptr;
    other.owner_ = false;
    other.size_ = 0;
  }

  void copy(const Data &other) {
      std::copy(other.get_data_pointer(),
                other.get_data_pointer(),
                data_);
  }

  /** Destroys array and frees memory. */
  ~Data() {
      if (owner_ && data_) {
          delete[] data_;
          data_=nullptr;
      }
  }
  /** Access element i. */
  const T& operator[](size_t i) const {return data_[i]; }
  /** Access element i. */
  T& operator[](size_t i) {return data_[i]; }

  /** Get pointer to raw data.*/
  const T* get_data_pointer() const { return data_; }
  T* get_data_pointer() { return data_; }


////////////////////////////////////////////////////////////////////////////////
// Basic operators
////////////////////////////////////////////////////////////////////////////////

  Data& operator*=(const T &t) {
      for (size_t i = 0; i < size_; ++i) {
          data_[i] *= t;
      }
  }

  Data& operator+=(const T &t) {
      for (size_t i = 0; i < size_; ++i) {
          data_[i] += t;
      }
  }

  Data& operator-=(const T &t) {
      for (size_t i = 0; i < size_; ++i) {
          data_[i] -= t;
      }
  }

  template <typename F>
  void map(F f){
      for (size_t i = 0; i < size_; ++i) {
          data_[i] = f();
      }
  }

 private:
  T* data_;
  size_t size_;
  bool owner_;
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
  Tensor(std::array<size_t, rank> shape)
      : shape_{shape}, data_{reduce<Multiply>(shape)} {
    // Nothing to do here.
  }
  Tensor(T* ptr, std::array<size_t, rank> shape)
      : shape_{shape}, data_{ptr, reduce<Multiply>(shape)} {}

  const std::array<size_t, rank>& shape() const { return shape_; }

  template <typename... Ts>
  auto operator()(Ts... indices)
      -> std::conditional<sizeof...(indices) < rank,
                          Tensor<T, rank - sizeof...(indices)>, T&>::type {
    constexpr size_t n = sizeof...(indices);
    std::array<size_t, n> index_array{indices...};
    if constexpr (n < rank) {
      return Tensor<T, rank - n>(
          get_data_pointer() + index(index_array, shape_), tail<n>(shape_));
    } else {
      return data_[index(index_array, shape_)];
    }
  }

  template <typename... Ts>
  const T& operator()(Ts... indices) const {
    std::array<size_t, sizeof...(indices)> index_array{indices...};
    return data_[index(index_array, shape_)];
  }

  void fill(const T& t) {
    for (size_t i = 0; i < reduce<Multiply>(shape_); ++i) {
      data_[i] = t;
    }
  }

  void copy(const Tensor &other) {
      data_.copy(other.get_data());
  }

  const T* get_data_pointer() const { return data_.get_data_pointer(); }
  T* get_data_pointer() { return data_.get_data_pointer(); }

  const Data<T> & get_data() const {return data_;}
  Data<T> & get_data() {return data_;}

  Tensor& operator*=(const T &t) {
      data_ *= t;
  }

  Tensor& operator+=(const T &t) {
      data_ += t;
  }

  Tensor& operator-=(const T &t) {
      data_ -= t;
  }

  template <typename F>
  void map(F f) {
      data_.template map(f);
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

 protected:
  std::array<size_t, rank> shape_;
  Data<T> data_;
};

template <typename FloatType>
class Array : public Tensor<FloatType, 1> {
 public:
  using Tensor<FloatType, 1>::shape_;
  using Tensor<FloatType, 1>::operator();

  Array(size_t size) : Tensor<FloatType, 1>({size}) {}

  Array(const Array& other) = default;
  Array(Array&& other) = default;
  Array& operator=(const Array& other) = default;
  Array& operator=(Array&& other) = default;

  FloatType& operator[](size_t i) { return this->operator()(i); }
  FloatType operator[](size_t i) const { return this->operator()(i); }

  Array cumulative_integral(const Array dx,
                            FloatType c = static_cast<FloatType>(0.0)) const {
    assert(dx.size() == size());
    Array result(size() + 1);
    FloatType integral = c;
    result[0] = c;
    for (size_t i = 0; i < size(); ++i) {
        result[i + 1] = result[i] + this->operator[](i)* dx[i];
    }
    return result;
  }

  Array diff() const {
    Array output{size() - 1};
    for (size_t i = 0; i < size() - 1; ++i) {
        output[i] = this->operator[](i + 1) - this->operator[](i);
    }
    return output;
  }

  Array centers() const {
      Array output{size() - 1};
      for (size_t i = 0; i < size() - 1; ++i) {
          output[i] = 0.5(this->operator[](i + 1) + this->operator[](i - 1));
      }
      return output;
  }

  size_t size() const { return shape_[0]; }

  static Array fill_linear(FloatType start, FloatType stop, size_t n_elements) {
    Array<FloatType> result{n_elements};
    FloatType d = (stop - start) / (n_elements - 1);
    FloatType x = start;
    for (size_t i = 0; i < n_elements; ++i) {
      result[i] = x;
      x += d;
    }
    return result;
  }

  static Array range(size_t n_elements) {
      Array<FloatType> result{n_elements};
      for (size_t i = 0; i < n_elements; ++i) {
          result[i] = static_cast<FloatType>(i);
      }
      return result;
  }
};

}  // namespace llrte
#endif
