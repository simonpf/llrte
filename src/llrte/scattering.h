#ifndef _LLRTE_SCATTERING_H_
#define _LLRTE_SCATTERING_H_

namespace llrte {

template <typename F>
class NoScattering {
 public:
  struct PhaseFunction {
    PhaseFunction() {}
    template <typename G, typename T>
    constexpr const T& get_direction(G& /*g*/, const T& t) {
      return t;
    }
  };

  NoScattering() {}

  template <typename... Ts>
  constexpr F get_scattering_coefficient(Ts...) {
    return 0.0;
  }

  template <typename... Ts>
  constexpr PhaseFunction get_phase_function(Ts...) {
    return PhaseFunction();
  }
};

template <typename F>
class BidirectionalScattering {
 public:
  using Float = F;

  struct PhaseFunction {
    PhaseFunction() {}
    template <typename G, typename T>
    T get_direction(G& g, const T& t) {
      auto x = g.sample_uniform();
      if (x < 0.5) {
        return t;
      } else {
        return static_cast<Float>(-1.0) * t;
      }
    }
  };

  BidirectionalScattering(Float scattering_coefficient)
      : scattering_coefficient_(scattering_coefficient) {
    // Nothing to do here.
  }

  template <typename... Ts>
  constexpr F get_scattering_coefficient(Ts...) {
    return scattering_coefficient_;
  }

  template <typename... Ts>
  constexpr PhaseFunction get_phase_function(Ts...) {
    return PhaseFunction();
  }

 private:
  F scattering_coefficient_;
};

}  // namespace llrte
#endif
