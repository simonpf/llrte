#include "llrte/data.h"


template<typename FloatType>
bool test_fill() {
    auto t = llrte::Tensor<FloatType, 2>({3, 3});

    auto t2 = t(1);
    t2.fill(2.0);

    for (size_t i = 0; i < 3; ++i) {
        if (t(1, i) != 2.0) {
            return false;
        }
    }

    std::cout << t2 << std::endl;
    return true;
}

template<typename FloatType>
bool test_diff() {
    auto t = llrte::Array<FloatType>::range(5);
    auto d = t.diff();
    std::cout << t << std::endl;
    std::cout << d << std::endl;
    return true;
}

template<typename FloatType>
bool test_integral() {
    auto t = llrte::Array<FloatType>::range(5);
    auto d = t.diff();
    auto i = d.cumulative_integral(d);
    
    std::cout << t << std::endl;
    std::cout << d << std::endl;
    std::cout << i << std::endl;
    return true;
}

int main(int /*argc*/, char **/*argv*/) {
    test_diff<double>();

    test_integral<float>();
    test_integral<double>();

}
