#ifndef _LLRTE_UTILS_ARRAY_H_
#define _LLRTE_UTILS_ARRAY_H_

#include <array>
#include <iostream>

namespace llrte::utils::array {

    struct Multiply {

        template <typename T>
        static constexpr T identity() {
            return static_cast<T>(1.0);
        }

        template<typename T>
        static T apply(const T &a, const T &b) {
            return a * b;
        }
    };

    struct Add {

        template <typename T>
        static constexpr T identity() {
            return static_cast<T>(0.0);
        }

        template<typename T>
        static T apply(const T &a, const T &b) {
            return a + b;
        }
    };

    template <typename Operation, typename T, size_t N>
    T reduce(const std::array<T, N> &array) {
        T acc = Operation::template identity<T>();
      for (size_t i = 0; i < N; ++i) {
          acc = Operation::apply(acc, array[i]);
      }
      return acc;
    }

    template <typename T, size_t N>
    T index(const std::array<T, N> &indices, const std::array<T, N> &sizes) {
      auto i = static_cast<T>(0);
      for (size_t i = 0; i < N; ++i) {
        i += indices[i];
        if (i < N - 1) {
            i *= sizes[i];
        }
      }
      return i;
    }

    template <typename Operation, typename T, size_t N>
    std::array<T, N> zip(const std::array<T, N> &as,
          const std::array<T, N> &bs) {
      std::array<T, N> cs;
      for (size_t i = 0; i < N; ++i) {
          cs[i] = Operation::apply(as[i], bs[i]);
      }
      return cs;
    }

}  // namespace llrte::utils::array

template<typename T, size_t rank>
    std::ostream& operator<<(std::ostream& os, const std::array<T, rank> &array) {
    os << "[";
    for (size_t i = 0; i < rank - 1; ++i) {
        os << array[i] << ", ";
    }
    os << array[rank - 1] << "]";
    return os;
}

#endif
