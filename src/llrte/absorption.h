#ifndef _LLRTE_ABSORPTION_H_
#define _LLRTE_ABSORPTION_H_

namespace llrte {

/**
 * The simplest absorption model: No absorption anywhere in the atmosphere.
 */
template <typename Float>
class NoAbsorption {
 public:
  /** Create no absorption model. */
  NoAbsorption() {}

  /**
   * This is the interface function to the Monte Carlo solver.
   */
  template <typename... Ts>
  __DEV__ constexpr Float get_absorption_coefficient(Ts &&...) {
    return 0.0;
  }
};

/**
 * The next simplest absorption model, constant absorption everywhere in the
 * atmosphere.
 */
template <typename Float>
class ConstantAbsorption {
 public:
    /**
     * Create absorption model with given constant absorption.
     *
     * @param absorption The absorption cross section. Note that units
     * should match that of the atmosphere.
     */
  ConstantAbsorption(Float absorption) : absorption_coefficient_(absorption) {}


  /**
  * Returns the constant absorption coefficient.
  */
  template <typename... Ts>
  __DEV__ Float get_absorption_coefficient(Ts && ...) {
    return absorption_coefficient_;
  }

 private:
  Float absorption_coefficient_;
};

}  // namespace llrte
#endif
