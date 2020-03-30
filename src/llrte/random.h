#ifndef _LLRTE_RANDOM_
#define _LLRTE_RANDOM_

#include <random>
#include <chrono>

#include "llrte/common.h"
#include "llrte/constants.h"
#include "llrte/rotations.h"

namespace llrte {

////////////////////////////////////////////////////////////////////////////////
// Random number generator
////////////////////////////////////////////////////////////////////////////////

/** Generator
 *
 * The Generator class providees n architecture agnostic, thread-safe number
 * generator.
 *
 * @tparam F The floating point type to use to represent floating point numbers.
 */
template <typename F>
class Generator {

 public:

  using Float = F;
  using C = Constants<Float>;

  /** Create a random number generator. */
  Generator() {
      std::chrono::seconds sec(1);
      seed_ = std::chrono::microseconds(sec).count();
  }

  ~Generator() {
    #ifdef CUDA
    if (curand_state_) {
        delete curand_state_;
        curand_state_=nullptr;
    }
    #endif
  }

  /** Initialize random number generator.
   *
   * This needs to be called before any random number is generated.
   */
  void initialize() {}

  /**
   * Generate a sample from a uniform distribution.
   * @return A uniform sample from the range [0, 1] of the given floating point precision
   */
  Float sample_uniform() {return distribution_(generator_);}

  /**
   * Generate uniform sample in given range.
   * @param a Left end of the range
   * @param b Right end of the range
   * @return A uniform sample from the range [0, 1] of the given floating point precision
   */
  __DEV__ Float sample_uniform(Float a, Float b) {
      return a  + (b - a) * sample_uniform();
  }

  /**
   * Sample angle uniformly from range [-PI, +PI]
   * @return A uniform sample from the range [-PI, +PI]
   */
  __DEV__ Float sample_angle_uniform() { return sample_uniform(-C::pi, C::pi);}
  __DEV__ Float sample_uniform_angle() { return sample_uniform(-C::pi, C::pi);}

  /**
   * Sample zenith angle so that points are uniformly distributed on sphere.
   * @param a Lower bound (default 0)
   * @param b Upper bound (default PI)
   * @return A uniform sample from the range [a, b]
   */
  __DEV__ Float sample_zenith_angle(float a = 0.0,
                                    float b = Constants<Float>::pi){
      Float a_i = -cos(a);
      Float b_i = -cos(b);
      Float r = sample_uniform(a_i, b_i);
      return acos(-r);
  }

  /**
   * Sample exponentially distributed path length.
   * @param m The mean path length
   * @return Sample from exponential distribution with given path length.
   */
  __DEV__ Float sample_path_length(Float m) {
      auto y = distribution_(generator_);
      return -m * log(y);
  }

  /**
   * Sample exponentially distributed optical depth.
   * @return Sample from exponential distribution with path length 1.
   */
  __DEV__ Float sample_tau() {
      auto y = sample_uniform();
      return -log(y);
  }

  //
  // Device code
  //

  #ifdef __CUDACC__
  __device__ void initialize() {
      curand_state_ = new curandState;
      int idx = threadIdx.x + blockIdx.x * blockDim.x;
      curand_init(seed_, idx, 0, curand_state_);
  }
  __device__ Float sample_uniform() {return curand_uniform(curand_state_);}
  #endif

 private:
  unsigned seed_;
  std::default_random_engine generator_{};
  std::uniform_real_distribution<Float> distribution_{0.0, 1.0};
  #ifdef CUDA
  curandState_t *curand_state_ = nullptr;
  #endif

};

////////////////////////////////////////////////////////////////////////////////
// Utility functions.
////////////////////////////////////////////////////////////////////////////////

/**
 * Sample direction uniformly from sphere.
 *
 * @tparam Vector The vector type to represent 3D vector.
 * @tparam Generator The random number generator to use.
 * @param generator Random number generate to use to
 * generate random numbers.
 */
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
