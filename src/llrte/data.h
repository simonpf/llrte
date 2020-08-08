#ifndef _LLRTE_DATA_H_
#define _LLRTE_DATA_H_

#include <iostream>
#include <memory>
#include "assert.h"
#include "common.h"

#include "llrte/common.h"
#include "llrte/utils/array.h"

namespace llrte {

using utils::array::zip, utils::array::reduce, utils::array::index,
    utils::array::Multiply, utils::array::Add, utils::array::tail;

template <typename Float>
    std::shared_ptr<Float[]> make_shared(size_t n) {
    return std::shared_ptr<Float[]>{new Float[n], [](Float *p)
        {
            if (p) delete [] p;
        }
    };
}

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
  Data() : size_(0), offset_(0) {};
  /** Create new array holding n elements.
   *
   * Create a new array and initializes data to 0. The created array
   * will own the data it created and destroy it when it is destroyed.
   *
   * @param size Size of the array to create.
   */
    Data(size_t size) : data_(make_shared<T>(size)), offset_(0), size_(size) {
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
  Data(std::shared_ptr<T[]> ptr, size_t offset, size_t size) : data_(ptr), offset_(offset), size_(size) {}

  /** Copy constructor
   *
   * Performs a deep copy of the data of the other array.
   * @param other The array to copy.
   */
  Data(const Data& other) :
    data_(make_shared<T>(other.get_size())), offset_(0), size_(other.get_size()) {
        copy(other);
    }

  /**
   * Take over memory from other array. No data is copied.
   *
   * @param other The array to take over.
   */
    Data(Data&& other) = default;

  /**
   * Performs a deep copy of the data.
   * @param other The array to copy from.
   */
    Data& operator=(const Data& other) { copy(other); }

    /**
     * Move over data from other array. No data is copied.
     * @param other The array to move.
     */
    Data& operator=(Data&& other) = default;

  /** Destroys array and frees memory. */
  ~Data() {
    #ifdef __CUDA__
    if (device_data_){
        printf("Freeing device mem: %p" , device_data_);
        cudaFree(reinterpret_cast<void*>(device_data_));
        device_data_=nullptr;
    }
    #endif
  }

  /** Access element i.*/
  const T& operator[](size_t i) const { return data_[i + offset_]; }
  /** Access element i. */
  T& operator[](size_t i) { return data_[i + offset_]; }
  /** Get pointer to raw data.*/
  const T* get_data_pointer() const { return data_.get() + offset_; }
  /** Get pointer to raw data.*/
  T* get_data_pointer() { return data_.get() + offset_; }
  /** Get shared pointer to holding the data.*/
  std::shared_ptr<T[]> get_shared_pointer() const { return data_; }

  /** Get size of data array. */
  size_t get_size() const {return size_;}

  size_t size() { return size_; }

  /** Scale numeric data in array. */
  Data& operator*=(const T& t) {
    for (size_t i = 0; i < size_; ++i) {
      data_[i + offset_] *= t;
    }
    return *this;
  }

  /** Element-wise sum. */
  Data& operator+=(const Data &d) {
      if (size_ != d.get_size()) {
          throw std::runtime_error("Can't add arrays of different sizes.");
      }
      for (size_t i = 0; i < size_; ++i) {
          data_[i + offset_] += d[i];
      }
      return *this;
  }

  /** Scale numeric data in array. */
  Data& operator/=(const T& t) {
      auto ti = 1.0 / t;
      for (size_t i = 0; i < size_; ++i) {
          data_[i + offset_] *= ti;
      }
      return *this;
  }

  /** In-place add to numeric data in array. */
  Data& operator+=(const T& t) {
    for (size_t i = 0; i < size_; ++i) {
      data_[i + offset_] += t;
    }
    return *this;
  }

  /** In-place subtract to numeric data in array. */
  Data& operator-=(const T& t) {
    for (size_t i = 0; i < size_; ++i) {
      data_[i + offset_] -= t;
    }
    return *this;
  }

  /** Map function over data. */
  template <typename F>
  void map(F f) {
    for (size_t i = 0; i < size_; ++i) {
      data_[i + offset_] = f();
    }
  }

  /**
   *Fill array with value
   * @param t Value to fill array with.
   */
  void fill(T t) {
    for (size_t i = 0; i < size_; ++i) {
        data_[i + offset_] = t;
    }
  }

  /**
   * Copy data from other array.
   * @param other Array to copy data from.
   */
  void copy(const Data &other) {
    if (size_ != other.get_size()) {
      throw std::runtime_error("Sizes of array to copy don't match.");
    }
    std::copy(other.get_data_pointer(), other.get_data_pointer() + size_, get_data_pointer());
  }

  /**
   *Fill array with value
   * @param t Value to fill array with.
   */
  T sum() {
      T t = static_cast<T>(0.0);
      for (size_t i = 0; i < size_; ++i) {
          t += data_[i + offset_];
      }
      return t;
  }

  void atomic_add(size_t i, T t) {
      OPENMP_ATOMIC(data_[i + offset_] = data_[i + offset_] + t;)
  }

  #ifdef __CUDACC__
  __device__ void atomic_add(size_t i, T t) {atomicAdd(device_data_ + offset_ + i, t);}
  __device__ const T& operator[](size_t i) const { return device_data_[i + offset_]; }
  __device__ T& operator[](size_t i) { return device_data_[i + offest_]; }
  __device__ size_t size() { return size_; }
  __device__ T* get_ptr() { return device_data_ + offset_; }
  #endif 
  #ifdef __CUDA__
  __device__ Data(const Data& other) = default;
  void device() {
      if (offset_ > 0) {
          throw std::runtime_error("Sub-tensors cannot be copied to any device.");
      }
      if (device_data_) {
          CUDA_CALL(cudaFree(reinterpret_cast<void*>(device_data_)));
          CUDA_CALL(cudaDeviceSynchronize());
      }
      CUDA_CALL(cudaMalloc(&device_data_, size_ * sizeof(T)));
      CUDA_CALL(cudaMemcpy(device_data_,
                           data_,
                           (int) size_ * sizeof(T),
                           cudaMemcpyHostToDevice));
  }
  void host() {
      if (device_data_) {
          CUDA_CALL(cudaMemcpy(data_,
                               device_data_ + offset_,
                               size_ * sizeof(T),
                               cudaMemcpyDeviceToHost));
      }
      CUDA_CALL(cudaDeviceSynchronize());
  }
  #endif

  friend Data operator*(T t, Data u) {
      Data results(u.size_);
      for (size_t i = 0; i < u.size_; ++i) {
          results[i] = t * u[i];
      }
      return results;
  }

  friend Data operator/(Data u, T t) {
      Data results(u.size_);
      for (size_t i = 0; i < u.size_; ++i) {
          results[i] = u[i] / t;
      }
      return results;
  }

  friend Data operator+(T t, Data u) {
      Data results(u.size_);
      for (size_t i = 0; i < u.size_; ++i) {
          results[i] = u[i] + t;
      }
  }

  friend Data operator+(Data u, T t) {
      Data results(u.size_);
      for (size_t i = 0; i < u.size_; ++i) {
          results[i] = u[i] + t;
      }
  }

  friend Data operator-(T t, Data u) {
      Data results(u.size_);
      for (size_t i = 0; i < u.size_; ++i) {
          results[i] = t - u[i];
      }
  }

  friend Data operator-(Data u, T t) {
      Data results(u.size_);
      for (size_t i = 0; i < u.size_; ++i) {
          results[i] = u[i] - t;
      }
  }


protected:

  std::shared_ptr<T[]> data_ = nullptr;
  T* device_data_ = nullptr;
  size_t offset_, size_;
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
    Tensor() {}
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
    Tensor(std::shared_ptr<T[]> ptr, size_t offset, std::array<size_t, rank> shape)
      : shape_{shape}, data_{ptr, offset, reduce<Multiply>(shape)} {}

    /** @return Array containing the shape of the tensor. */
  const std::array<size_t, rank>& shape() const { return shape_; }

  Tensor(const Tensor &) = default;
  Tensor(Tensor &&) = default;
  Tensor& operator =(const Tensor &) = default;
  Tensor& operator =(Tensor &&) = default;

  /**
   * Return element or sub-tensor. If the length of the provided index sequence
   * is less than the rank of the tensor a sub-tensor view is returned. Other-
   * wise a double is returned.
   * @param Ts Sequence of indices of the elements to extract along each tank.
   * @return Depending on length of input sequence, reference to element of
   * view on sub-tensor.
   */
  template <typename... Ts>
  __DEV__ auto operator()(Ts... indices)
      -> typename std::conditional<sizeof...(indices) < rank,
                          Tensor<T, rank - sizeof...(indices)>, T&>::type {
    constexpr size_t n = sizeof...(indices);
    std::array<size_t, n> index_array{indices...};
    if constexpr (n < rank) {
      return Tensor<T, rank - n>(
          get_shared_pointer(), index(index_array, shape_), tail<n>(shape_));
    } else {
      return data_[index(index_array, shape_)];
    }
  }

  /**
   * Return element from tensor
   * @param Ts Sequence of indices of the elements to extract along each tank.
   */
  template <typename... Ts>
      __DEV__ auto operator()(Ts... indices) const
      -> typename std::conditional<sizeof...(indices) < rank,
      const Tensor<T, rank - sizeof...(indices)>, T>::type {
      constexpr size_t n = sizeof...(indices);
      std::array<size_t, n> index_array{indices...};
      if constexpr (n < rank) {
              return Tensor<T, rank - n>(
                  get_shared_pointer(), index(index_array, shape_), tail<n>(shape_));
          } else {
          return data_[index(index_array, shape_)];
      }
  }

  /**
   * Fill tensor.
   * @param t Element to fill tensor with.
   */
  void fill(const T& t) {
      data_.fill(t);
  }

  /**
   * Copy data from other tensor..
   * @param t Tensor to copy data from.
   */
  void copy(const Tensor& t) {
      data_.copy(t.get_data());
  }

  Tensor copy() const {
      Tensor tensor(shape_);
      tensor.copy(*this);
      return tensor;
  }

  /**
   * Return pointers to tensor data. Data is stored following C conventions,
   * i.e. with the last rank being continuous in memory.
   * @return Pointer to first element of tensor.
   */
  const T* get_data_pointer() const { return data_.get_data_pointer(); }
  T* get_data_pointer() { return data_.get_data_pointer(); }

  std::shared_ptr<T[]> get_shared_pointer() const { return data_.get_shared_pointer(); }
  /**
   * Return data as array.
   * @return Data object holding the tensor data.
   */
  const Data<T> & get_data() const {return data_;}
  Data<T> & get_data() {return data_;}

  /** In-place scaling of tensor elements. */
  Tensor& operator*=(const T &t) {
      data_ *= t;
      return *this;
  }

  /** In-place addition to tensor elements. */
  Tensor& operator+=(const T &t) {
      data_ += t;
      return *this;
  }

  /** Element-wise addition to tensor elements. */
  Tensor& operator+=(const Tensor &t) {
      data_ += t.data_;
      return *this;
  }

  /** In-place subtraction from elements. */
  Tensor& operator-=(const T &t) {
      data_ -= t;
      return *this;
  }

  /** Map function over elements. */
  template <typename F>
  void map(F f) {
      data_.template map(f);
  }

  T sum() {
      return data_.sum();
  }

  template <typename... Ts>
      __DEV__ void atomic_add(T t, Ts... indices) {
      std::array<size_t, sizeof...(indices)> index_array{indices...};
      data_.atomic_add(index(index_array, shape_), t);
  }

    #ifdef __CUDA__
    void device() {data_.device();}
    void host() {data_.host();}
    #endif

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

  friend Tensor operator*(T t, Tensor u) {
    Tensor results(u.shape_);
    results.data_ = t * u.data_;
    return results;
  }

  friend Tensor operator/(Tensor t, T u) {
    Tensor results(t.shape_);
    results.data_ = u.data_ / u.data_;
    return results;
  }

  friend Tensor operator+(T t, Tensor u) {
    Tensor results(u.shape_);
    results.data_ = t + u.data_;
    return results;
  }

  friend Tensor operator+(Tensor u, T t) {
      Tensor results(u.shape_);
      results.data_ = t + u.data_;
      return results;
  }

  friend Tensor operator-(T t, Tensor u) {
      Tensor results(u.shape_);
      results.data_ = t - u.data_;
      return results;
  }

  friend Tensor operator-(Tensor u, T t) {
      Tensor results(u.shape_);
      results.data_ = u.data_ - t;
      return results;
  }

 protected:
  std::array<size_t, rank> shape_;
  Data<T> data_;
};

////////////////////////////////////////////////////////////////////////////////
// Arrays
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

  __DEV__ FloatType& operator[](size_t i) { return this->operator()(i); }
  __DEV__ FloatType operator[](size_t i) const { return this->operator()(i); }

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

  __DEV__ size_t size() const { return shape_[0]; }

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
