#ifndef _LLRTE_RANDOM_
#define _LLRTE_RANDOM_

template <typename V>
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

  Float sample_uniform() { return distribution_(generator_); }

  Float sample_tau() {
    auto y = distribution_(generator_);
    return -log(y);
  }

  template<typename Vector>
  Vector random_direction() {
      using Float = typename Vector::Float;
      Vector v{};
      Float phi, theta;

  }

 public:
  std::default_random_engine generator_{};
  std::uniform_real_distribution<Float> distribution_{0.0, 1.0};
};

#endif
