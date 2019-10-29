#ifndef _LLRTE_RANDOM_
#define _LLRTE_RANDOM_

template<typename V>
class Generator {

    public:

    using Vector = V;
    using Float = typename Vector::Float;

    Generator() {
        // Nothing to do here.
    }

    Float sample_path_length(Float m) {
        auto y = distribution_(generator_);
        return -m * log(y);
    }


public:
    std::default_random_engine generator_{};
    std::uniform_real_distribution<Float> distribution_{0.0,1.0};
};

#endif