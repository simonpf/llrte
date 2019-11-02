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
    template <typename G, typename T>
    T get_direction(G& g, const T& t) {
      auto x = g.sample_uniform();
      if (x > fb_ratio) {
        return static_cast<Float>(-1.0) * t;
      } else {
          return t;
      }
    }
    Float fb_ratio;
  };

  BidirectionalScattering(Float scattering_coefficient, Float fb_ratio)
      : scattering_coefficient_(scattering_coefficient), fb_ratio_(fb_ratio) {
    // Nothing to do here.
  }

  template <typename... Ts>
  constexpr F get_scattering_coefficient(Ts...) {
    return scattering_coefficient_;
  }

  template <typename... Ts>
  constexpr PhaseFunction get_phase_function(Ts...) {
    return PhaseFunction{fb_ratio_};
  }

 private:
  F scattering_coefficient_;
  F fb_ratio_;
};

}  // namespace llrte
#endif
