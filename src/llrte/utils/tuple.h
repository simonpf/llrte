#ifndef _LLRTE_UTILS_TUPLE_H_
#define _LLRTE_UTILS_TUPLE_H_

#include <iostream>
#include <tuple>

namespace llrte::tuple {
namespace detail {

template <typename Tuple>
struct tuple_size {
  static constexpr size_t value = std::tuple_size<Tuple>::value;
};

template <>
struct tuple_size<std::tuple<>> {
  static constexpr size_t value = 0;
};

template <typename F, class Tuple, size_t i>
struct MapImpl {
  template <typename... Args>
__DEV__  static void apply(F &f, Tuple &t, Args & ...args) {
    constexpr size_t ind =
        std::tuple_size<std::remove_reference_t<Tuple>>::value - i;
    auto& tt = std::get<ind>(t);
    f.apply(tt, args...);
    MapImpl<F, Tuple, i - 1>::apply(f, t, args...);
  }
};
template <typename F, class Tuple>
struct MapImpl<F, Tuple, 0> {
  template <typename... Args>
  __DEV__ static void apply(F &, Tuple &, Args & ...) {}
};
template <typename F, class Tuple, size_t i>
struct LoopImpl {
  template <typename... Args>
  __DEV__ static void apply(Tuple t, F f, Args... args) {
    constexpr size_t ind = std::tuple_size<Tuple>::value - i;
    auto tt = std::get<ind>(t);
    F::apply(tt, f, i, args...);
    LoopImpl<F, Tuple, ind>::apply(t, f, args...);
  }
};
template <typename F, class Tuple>
struct LoopImpl<F, Tuple, 0> {
  template <typename... Args>
  __DEV__ static void apply(Tuple t, F f, Args...) {}
};

}  // namespace detail

template <typename F, typename Tuple, typename... Args>
__DEV__ void map(F &f, Tuple &t, Args & ...args) {
  constexpr size_t tuple_size =
      detail::tuple_size<std::remove_reference_t<Tuple>>::value;
  detail::MapImpl<F, Tuple, tuple_size>::apply(f, t, args...);
}

template <typename F, typename Tuple, typename... Args>
void loop(Tuple&& t, F f, Args... args) {
  detail::LoopImpl<F, Tuple, std::tuple_size<Tuple>::value>::apply(f, t,
                                                                   args...);
}

}  // namespace llrte::tuple
#endif
