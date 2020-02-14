#ifndef _LLRTE_UTILS_ARRAY_H_
#define _LLRTE_UTILS_ARRAY_H_

#include <array>

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
    T reduce(std::array<T, N> &array) {
        T acc = Operation::template identity<T>();
      for (size_t i = 0; i < N; ++i) {
          acc = Operation::apply(acc, array[i]);
      }
      return acc;
    }

    template <typename Operation, typename T, size_t N>
    T zip(const std::array<T, N> &as,
          const std::array<T, N> &bs) {
      std::array<T, N> cs;
      for (size_t i = 0; i < N; ++i) {
          cs[i] = Operation::apply(as[i], bs[i]);
      }
      return cs;
    }
}

#endif
