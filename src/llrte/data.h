#ifndef _LLRTE_DATA_H_
#define _LLRTE_DATA_H_

#include <iostream>
#include "assert.h"
#include "common.h"

#include "llrte/utils/array.h"

namespace llrte {

using utils::array::zip, utils::array::reduce, utils::array::index,
    utils::array::Multiply, utils::array::Add, utils::array::tail;

////////////////////////////////////////////////////////////////////////////////
// Array data
////////////////////////////////////////////////////////////////////////////////
/** Data
 *
 * The Data class provides linear array-type storage for a given
 * data type. It should be used to store "heavy-weight" data.
 * A Data object in own the data that it allocated itself during construction
 * but can also be used as a view on existing data.
 *
 * \tparam T The type of the data to store.
 */
template <typename T>
class Data {
 public:
  /** Create new array holding n elements.
   *
   * Create a new array and initializes data to 0. The created array
   * will own the data it created and destroy it when it is destroyed.
   *
   * @param size Size of the array to create.
   */
  Data(size_t size) : data_(new T[size]), size_(size), owner_(true) {
    fill(0.0);
  }

  /** Create array from existing data.
   *
   * When called with this constructor the array will not release memory
   * when it goes out of scope.
   *
   * @param ptr Pointer to the existing array
   * @param size Size of the array.
   */
  Data(T* ptr, size_t size) : data_(ptr), size_(size), owner_(false) {}

  /** Copy constructor
   *
   * Performs a deep copy of the data of the other array.
   * @param other The array to copy.
   */

    Data(const Data& other)
        : data_(other.data_), device_data_(other.device_data_), size_(other.size_), owner_(false) {}

  /**
   * Take over memory from other array. No data is copied.
   *
   * @param other The array to take over.
   */
  Data(Data&& other)
      : data_(other.data_), size_(other.size_), owner_(other.owner_) {
    other.data_ = nullptr;
    other.owner_ = false;
    other.size_ = 0;
  }

  /**
   * Performs a deep copy of the data.
   * @param other The array to copy from.
   */
  Data& operator=(const Data& other) {
    assert(other.size_ == size_);
    std::copy(other.data_, other.data_ + size_, data_);
    return *this;
  }


  /**
   * Move over data from other array. No data is copied.
   * @param other The array to move.
   */
  Data& operator=(Data&& other) {
    if (data_ && owner_) delete[] data_;
    data_ = other.data_;
    size_ = other.size_;
    owner_ = other.owner_;
    other.data_ = nullptr;
    other.owner_ = false;
    other.size_ = 0;
  }

  /** Destroys array and frees memory. */
  ~Data() {
    if (owner_ && data_) {
      delete[] data_;
      data_ = nullptr;
    }
    #ifdef __CUDA__
    if (owner_ && device_data_){
        cudaFree(reinterpret_cast<void*>(device_data_));
        device_data_=nullptr;
    }
    #endif
  }

  /** Access element i.*/
  const T& operator[](size_t i) const { return data_[i]; }
  /** Access element i. */
  T& operator[](size_t i) { return data_[i]; }
  /** Get pointer to raw data.*/

  const T* get_data_pointer() const { return data_; }
  /** Get pointer to raw data.*/
  T* get_data_pointer() { return data_; }

  size_t size() { return size_; }

  /** Scale numeric data in array. */
  Data& operator*=(const T& t) {
    for (size_t i = 0; i < size_; ++i) {
      data_[i] *= t;
    }
    return *this;
  }

  /** Scale numeric data in array. */
  Data& operator/=(const T& t) {
      auto ti = 1.0 / t;
      for (size_t i = 0; i < size_; ++i) {
          data_[i] *= ti;
      }
      return *this;
  }

  /** In-place add to numeric data in array. */
  Data& operator+=(const T& t) {
    for (size_t i = 0; i < size_; ++i) {
      data_[i] += t;
    }
    return *this;
  }

  /** In-place subtract to numeric data in array. */
  Data& operator-=(const T& t) {
    for (size_t i = 0; i < size_; ++i) {
      data_[i] -= t;
    }
    return *this;
  }

  /** Map function over data. */
  template <typename F>
  void map(F f) {
    for (size_t i = 0; i < size_; ++i) {
      data_[i] = f();
    }
  }

  /**
   *Fill array with value
   * @param t Value to fill array with.
   */
  void fill(T t) {
    for (size_t i = 0; i < size_; ++i) {
        data_[i] = t;
    }
  }

  #ifdef __CUDACC__
  __device__ const T& operator[](size_t i) const { return device_data_[i]; }
  __device__ T& operator[](size_t i) { return device_data_[i]; }
  __device__ size_t size() { return size_; }
  __device__ T* get_ptr() { return data_; }
  #endif 
  #ifdef __CUDA__
__device__ Data(const Data& other)
    : data_(other.data_), device_data_(other.device_data_), size_(other.size_), owner_(false) {
  }
  void device() {
      if (device_data_) {
          cudaFree(reinterpret_cast<void*>(device_data_));
      }
      int count;
      cudaGetDeviceCount(&count);
      CUDAERROR(cudaMalloc(&device_data_, size_ * sizeof(T)));
      CUDAERROR(cudaMemcpy(device_data_,
                           data_,
                           (int) size_ * sizeof(T),
                           cudaMemcpyHostToDevice));
  }
  void host() {
      if (device_data_) {
          cudaMemcpy(data_,
                     device_data_,
                     size_ * sizeof(T),
                     cudaMemcpyDeviceToHost);
      }
      cudaDeviceSynchronize();
  }
  #endif


  T* data_;
  T* device_data_;
  size_t size_;
  bool owner_;
};

////////////////////////////////////////////////////////////////////////////////
// Tensors
////////////////////////////////////////////////////////////////////////////////
/**
 * Rank-k tensor.
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
    /**
     * Create a tensor with a given shape. When called with this constructor the
     * tensor will allocate and take care of all required memory.
     * @param shape Array with rank elements specifying the number of elements along
     * each rank.
     */
  Tensor(std::array<size_t, rank> shape)
      : shape_{shape}, data_{reduce<Multiply>(shape)} {}

    /**
     * Create a tensor with a given shape from existing data.
     * When created with this constructor the tensor does not take ownership
     * of the data.
     * @param ptr Pointer to the data holding the elements of the tensor.
     * @param shape Array with rank elements specifying the number of elements along
     * each rank.
     */
  Tensor(T* ptr, std::array<size_t, rank> shape)
      : shape_{shape}, data_{ptr, reduce<Multiply>(shape)} {}

    /** @return Array containing the shape of the tensor. */
  const std::array<size_t, rank>& shape() const { return shape_; }

  /**
   * Return element or sub-tensor. If the length of the provided index sequence
   * is less than the rank of the tensor a sub-tensor view is returned. Other-
   * wise a double is returned.
   * @param Ts Sequence of indices of the elements to extract along each tank.
   * @return Depending on length of input sequence, reference to element of
   * view on sub-tensor.
   */
  template <typename... Ts>
  auto operator()(Ts... indices)
      -> typename std::conditional<sizeof...(indices) < rank,
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

  /**
   * Return element from tensor
   * @param Ts Sequence of indices of the elements to extract along each tank.
   */
  template <typename... Ts>
  const T& operator()(Ts... indices) const {
    std::array<size_t, sizeof...(indices)> index_array{indices...};
    return data_[index(index_array, shape_)];
  }

  /**
   * Fill tensor.
   * @param t Element to fill tensor with.
   */
  void fill(const T& t) {
      data_.fill(t);
  }

  /**
   * Return pointers to tensor data. Data is stored following C conventions,
   * i.e. with the last rank being continuous in memory.
   * @return Pointer to first element of tensor.
   */
  const T* get_data_pointer() const { return data_.get_data_pointer(); }
  T* get_data_pointer() { return data_.get_data_pointer(); }

  /**
   * Return data as array.
   * @return Data object holding the tensor data.
   */
  const Data<T> & get_data() const {return data_;}
  Data<T> & get_data() {return data_;}

  /** In-place scaling of tensor elements. */
  Tensor& operator*=(const T &t) {
      data_ *= t;
  }

  /** In-place addition to tensor elements. */
  Tensor& operator+=(const T &t) {
      data_ += t;
  }

  /** In-place subtraction from elements. */
  Tensor& operator-=(const T &t) {
      data_ -= t;
  }

  /** Map function over elements. */
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

////////////////////////////////////////////////////////////////////////////////
// Tensors
////////////////////////////////////////////////////////////////////////////////

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
    result[0] = c;
    for (size_t i = 0; i < size(); ++i) {
        result[i + 1] = result[i] + this->operator[](i)* dx[i];
    }
    return result;
  }

  FloatType first() const {return this->operator[](0);}
  FloatType & first() {return this->operator[](0);}
  FloatType last() const {return this->operator[](size() - 1);}
  FloatType & last() {return this->operator[](size() - 1);}

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
