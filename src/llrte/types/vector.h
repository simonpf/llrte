#ifndef _LLRTE_TYPES_VECTOR_H_

#include <iostream>

namespace llrte {

/**
 * 3D Vector
 */
template <size_t N, typename F>
class Vector {
public:

    using Float = F;

    Vector () {
        for (size_t i = 0; i < N; ++i) {
            elements_[i] = 0.0;
        }
    }

    template <typename ... Ts>
    Vector(Ts ... ts) {
        auto t = std::make_tuple(ts ...);
        for (size_t i = 0; i < N; ++i) {
            //elements_[i] = std::get<i>(t);
        }
    }

    Float operator[](size_t i) const {
        return elements_[i];
    }

    Float& operator[](size_t i) {
        return elements_[i];
    }

    Vector operator+(const Vector &v) {
        Vector w;
        for (size_t i = 0; i < N; ++i) {
            w[i] = v[i] + elements_[i];
        }
        return w;
    }

    Vector operator*(const Float &v) {
        Vector w;
        for (size_t i = 0; i < N; ++i) {
            w[i] = v * elements_[i];
        }
        return w;
    }

    Vector operator-(const Vector &v) {
        Vector w;
        for (size_t i = 0; i < N; ++i) {
            w[i] = v[i] - elements_[i];
        }
        return w;
    }

    Float length() {
        Float s = 0.0;
        for (size_t  i = 0; i < N; ++i) {
            s += elements_[i] * elements_[i];
        }
        return sqrt(s);
    }

public:
    Float elements_[N];
};

template <size_t N, typename Float>
Vector<N, Float> operator*(Float c, Vector<N, Float> v) {
    Vector <N, Float> w{};
    for (size_t i = 0; i < N; ++i) {
    w[i] = c * v[i];
    }
}

template <size_t N, typename Real>
    std::ostream& operator<<(std::ostream& os, const Vector<N, Real>& v)
{
    os << "[";
    for (size_t i = 0; i < N - 1; ++i) {
        os << v[i] << ",";
    }
    os << v[N - 1] << "]" << std::endl;
    return os;
}

}
#endif
