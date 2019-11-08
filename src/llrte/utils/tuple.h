#ifndef _LLRTE_UTILS_TUPLE_H_
#define _LLRTE_UTILS_TUPLE_H_

#include <tuple>
#include <iostream>

namespace llrte::tuple {

namespace detail {
    template <
        template<typename T> typename F,
        class Tuple,
        size_t i
        >
        struct MapImpl
    {
        static void apply(Tuple t) {
            constexpr size_t ind = std::tuple_size<Tuple>::value - i;
            auto tt = std::get<ind>(t);
            F<decltype(tt)>::apply(tt);
            MapImpl<F, Tuple, i - 1>::apply(t);
        }
    };
    template <
        template<typename T> typename F,
        class Tuple
        >
        struct MapImpl<F, Tuple, 0>
        {
            static void apply(Tuple t) {}
        };
}  // namespace detail

template <
    template <typename T> typename F,
    class Tuple
    >
    void map(Tuple&& t)
{
    return detail::MapImpl<F, Tuple, std::tuple_size<Tuple>::value>::apply(t);
}

}
#endif
