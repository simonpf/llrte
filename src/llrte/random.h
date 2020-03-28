#ifndef _LLRTE_RANDOM_
#define _LLRTE_RANDOM_

#include <random>

#include "llrte/constants.h"
#include "llrte/rotations.h"

namespace llrte {

template <typename F>
class Generator {
 public:
  using Float = F;
  using C = Constants<Float>;

  Generator() {
    // Nothing to do here.
  }

  Float sample_uniform() {return distribution_(generator_);}
  Float sample_uniform(Float a, Float b) {
      return a  + (b - a) * distribution_(generator_);
  }
  Float sample_angle_uniform() { return sample_uniform(-C::pi, C::pi);}

  Float sample_zenith_angle(float a = 0.0,
                            float b = Constants<Float>::pi){
      Float a_i = -cos(a);
      Float b_i = -cos(b);
      Float r = sample_uniform(a_i, b_i);
      return acos(-r);
  }

  Float sample_path_length(Float m) {
    auto y = distribution_(generator_);
    return -m * log(y);
  }

  Float sample_tau() {
    auto y = distribution_(generator_);
    return -log(y);
  }

 public:
  std::default_random_engine generator_{};
  std::uniform_real_distribution<Float> distribution_{0.0, 1.0};
};

template <typename Vector, typename Generator>
Vector random_direction(Generator &generator) {
  using Float = typename Vector::Float;
  Float phi, theta;

  Vector v{};
  v[0] = 0;
  v[1] = 0;
  v[2] = 1.0;

  phi = 2.0 * Constants<Float>::pi * generator.sample_uniform();
  theta = Constants<Float>::pi * generator.sample_uniform();

  return rotations::rotate(v, phi, theta);
}

}  // namespace llrte
#endif
