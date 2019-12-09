#ifndef _LLRTE_SCATTERING_H_
#define _LLRTE_SCATTERING_H_

#include <math.h>

#include <llrte/types/array.h>
#include <llrte/maths/geometry.h>
#include <llrte/constants.h>
#include <llrte/rotations.h>

namespace llrte {

template <typename F>
class NoScattering {
 public:
  struct PhaseFunction {
    PhaseFunction() {}
    template <typename G, typename P>
    constexpr void scatter(G& /*g*/, const P& /*t*/) const {}
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
    template <typename G, typename P>
    void scatter(G& g, P& p) const {
      auto x = g.sample_uniform();
      auto d = p.get_direction();
      if (x > fb_ratio) {
        p.set_direction(static_cast<Float>(-1.0) * d);
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

template <typename F, typename ScatteringPlane = maths::geometry::RandomPlane>
class NumericPhaseFunction {
 public:
  NumericPhaseFunction(const Array<F> &sa,
                       const Array<F> &p_int)
 : p_int_(p_int), sa_(sa) {}

  template <typename Generator, typename Photon>
  void scatter(Generator &g, Photon &p) const {
    auto r = g.sample_uniform();

    size_t i = 0;
    while (r > p_int_[i]) {
      i++;
    }

    // Determine scattering angle.
    if (i > 0) {
      if ((r - p_int_[i - 1]) < (p_int_[i] - r)) {
        i--;
      }
    }
    auto sa = sa_[i];

    // Determine normal to scattering plane.
    auto d = p.get_direction();
    auto n = ScatteringPlane::get_normal(g, d);
    auto nd = rotations::rotate(d, n, sa);

    p.set_direction(nd);
  }

  private:
  const Array<F> &p_int_;
  const Array<F> &sa_;
};

template <typename F, typename ScatteringPlane = maths::geometry::RandomPlane>
class RayleighScattering {
public:
  using Float = F;
  using PhaseFunction = NumericPhaseFunction<Float, ScatteringPlane>;
  RayleighScattering(Float scattering_coefficient,
                     size_t n_steps)
  : scattering_coefficient_(scattering_coefficient) {
    sa_ = Array<Float>(n_steps + 1);
    p_int_ = Array<Float>(n_steps + 1);
    Float da = Constants<Float>::pi / n_steps;
    Float x = 0.0;
    Float s = 0.0;

    for (size_t i = 0; i < n_steps; ++i) {

      auto xl = x - 0.5 * da;
      auto xr = x + 0.5 * da;
      auto yl = 0.75 * (1.0 + cos(xl) * cos(xl)) * sin(xl);
      auto yr = 0.75 * (1.0 + cos(xr) * cos(xr)) * sin(xr);

      sa_[i] = x;
      p_int_[i] = s;
      x += da;
      s += 0.5 * (yl + yr) * da;
    }
    sa_[n_steps] = x;
    p_int_[n_steps] = s;

    for (size_t i = 0; i < n_steps + 1; ++i) {
      p_int_[i] /= s;
    }
  }

  template <typename... Ts>
  constexpr F get_scattering_coefficient(Ts...) {
    return scattering_coefficient_;
  }

  template <typename... Ts>
  PhaseFunction get_phase_function(Ts...) {
    return PhaseFunction(sa_, p_int_);
  }

  private:
  Float scattering_coefficient_ = 0.0;
  Array<Float> sa_ = Array<Float>();
  Array<Float> p_int_ = Array<Float>();
};

template<typename Float>
Float henyey_greenstein(Float g, Float theta) {
    return (static_cast<Float>(1.0) - g * g) / pow(1.0 +  g * g - g * cos(theta) * 2.0, 1.5);
}

template <typename F, typename ScatteringPlane = maths::geometry::RandomPlane>
class HenyeyGreenstein {
public:
  using Float = F;
  using PhaseFunction = NumericPhaseFunction<Float, ScatteringPlane>;
  HenyeyGreenstein(Float scattering_coefficient,
                   Float g,
                   size_t n_steps)
  : scattering_coefficient_(scattering_coefficient) {
    sa_ = Array<Float>(n_steps + 1);
    p_int_ = Array<Float>(n_steps + 1);
    Float da = Constants<Float>::pi / n_steps;
    Float x = 0.0;
    Float s = 0.0;

    for (size_t i = 0; i < n_steps; ++i) {

      Float xl = x - 0.5 * da;
      Float xr = x + 0.5 * da;
      Float yl = henyey_greenstein(g, xl);
      Float yr = henyey_greenstein(g, xr);

      sa_[i] = x;
      p_int_[i] = s;
      x += da;
      s += 0.5 * (yl + yr) * da;
    }
    sa_[n_steps] = x;
    p_int_[n_steps] = s;

    for (size_t i = 0; i < n_steps + 1; ++i) {
      p_int_[i] /= s;
    }
  }

  template <typename... Ts>
  constexpr F get_scattering_coefficient(Ts...) {
    return scattering_coefficient_;
  }

  template <typename... Ts>
  PhaseFunction get_phase_function(Ts...) {
    return PhaseFunction(sa_, p_int_);
  }

  private:
  Float scattering_coefficient_ = 0.0;
  Array<Float> sa_ = Array<Float>();
  Array<Float> p_int_ = Array<Float>();
};
}  // namespace llrte
#endif
