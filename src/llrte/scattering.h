#ifndef _LLRTE_SCATTERING_H_
#define _LLRTE_SCATTERING_H_

#include <math.h>

#include <llrte/eigen.h>
#include <llrte/geometry.h>
#include <llrte/constants.h>
#include <llrte/rotations.h>

namespace llrte {

using eigen::Vector;
using eigen::Index;
using eigen::ConstVectorMap;
using eigen::Tensor;
using eigen::array_cast;

/**
 * The simplest scattering class: No scattering at all. This class basically
 * does nothing, but implementing the interface to the Monte Carlo solver with
 * NOOPs.
 */
template <typename F>
class NoScattering {
 public:

  struct PhaseFunction {
    __DEV__ PhaseFunction() {}
    template <typename G, typename P>
    constexpr void scatter(G& /*g*/, const P& /*t*/) const {}
  };

  /** I got 99 problems but scattering ain't one.*/
  NoScattering() {}

  template <typename... Ts>
  __DEV__ constexpr F get_scattering_coefficient(Ts && ...) {
    return 0.0;
  }

  template <typename... Ts>
  constexpr PhaseFunction get_phase_function(Ts && ...) {
    return PhaseFunction();
  }
};

////////////////////////////////////////////////////////////////////////////////
// Bi-directional scattering
////////////////////////////////////////////////////////////////////////////////

/** Bidirectional scattering
 *
 * This class provides scattering in only the forward or backward
 * direcitons.
 */
template <typename F>
class BidirectionalScattering {
 public:
  using Float = F;

  struct PhaseFunction {
    template <typename G, typename P>
    void scatter(G& g, P& p) const {
      auto x = g.sample_uniform();
      auto d = p.direction;
      if (x > fb_ratio) {
          p.change_direction(static_cast<Float>(-1.0) * d);
      }
    }
    Float fb_ratio;
  };

  /**
   * Create bi-directional scattering model with given scattering cross-section
   * and foward-to-backward scattering ratio.
   * @param scattering_coefficient The scattering cross-section
   * @param fb_ratio The forward-to-backscattering ratio
   */
  BidirectionalScattering(Float scattering_coefficient,
                          Float fb_ratio)
      : scattering_coefficient_(scattering_coefficient), fb_ratio_(fb_ratio) {
    // Nothing to do here.
  }

  /**
   * Monte Calor interface.
   */
  template <typename... Ts>
  constexpr F get_scattering_coefficient(Ts...) {
    return scattering_coefficient_;
  }

  /**
   * Monte Calor interface.
   */
  template <typename... Ts>
  constexpr PhaseFunction get_phase_function(Ts...) {
    return PhaseFunction{fb_ratio_};
  }

 private:
  F scattering_coefficient_;
  F fb_ratio_;
};

template <typename T>
using ConstVectorRef = const Vector<T> &;

template <typename F, typename ScatteringPlane = geometry::RandomPlane,
        template <typename> typename VectorType = ConstVectorRef>
class NumericPhaseFunction {
 public:
  NumericPhaseFunction(VectorType<F> sa, VectorType<F> p_int)
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
    auto d = p.direction;
    auto n = ScatteringPlane::get_normal(g, d);
    auto nd = rotations::rotate(d, n, sa);

    p.change_direction(nd);
  }

 private:
  VectorType<F> p_int_;
  VectorType<F> sa_;
};

template <typename F, typename ScatteringPlane = geometry::RandomPlane>
class RayleighScattering {
public:
  using Float = F;
  using PhaseFunction = NumericPhaseFunction<Float, ScatteringPlane>;
  RayleighScattering(Float scattering_coefficient,
                     size_t n_steps)
: scattering_coefficient_(scattering_coefficient),
  sa_(Vector<Float>(n_steps + 1)),
  p_int_(Vector<Float>(n_steps + 1)) {
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
  Vector<Float> sa_;
  Vector<Float> p_int_;
};

template<typename Float>
Float henyey_greenstein(Float g, Float theta) {
    return (static_cast<Float>(1.0) - g * g) / pow(1.0 +  g * g - g * cos(theta) * 2.0, 1.5);
}

template <typename F, typename ScatteringPlane = geometry::RandomPlane>
class HenyeyGreenstein {
public:
  using Float = F;
  using PhaseFunction = NumericPhaseFunction<Float, ScatteringPlane>;
  HenyeyGreenstein(Float scattering_coefficient,
                   Float g,
                   size_t n_steps)
: scattering_coefficient_(scattering_coefficient),
  sa_(Vector<Float>(n_steps + 1)),
  p_int_(Vector<Float>(n_steps + 1))
{
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
  Vector<Float> sa_;
  Vector<Float> p_int_;
};

template <typename Grid>
class GriddedOpticalProperties {
 using Float = typename std::remove_reference<typename Grid::Float>::type;
 public:

    using PhaseFun =
        NumericPhaseFunction<Float, geometry::RandomPlane, eigen::ConstVectorMap>;

  class AbsorptionModel {
   public:
    AbsorptionModel(GriddedOpticalProperties *opt_props)
        : opt_props_(opt_props) {}

    template <typename GridPosition>
    Float get_absorption_coefficient(Grid &/*grid*/, GridPosition gp) {
      return opt_props_->get_absorption_coefficient(gp);
    }

   private:
    GriddedOpticalProperties *opt_props_;
  };

  AbsorptionModel get_absorption_model() {
      return AbsorptionModel(this);
  }

  class ScatteringModel {
   public:
      using PhaseFunction = PhaseFun;
    ScatteringModel(GriddedOpticalProperties *opt_props)
        : opt_props_(opt_props) {}

    template <typename GridPosition>
    Float get_scattering_coefficient(Grid & /*grid*/, GridPosition gp) {
      return opt_props_->get_scattering_coefficient(gp);
    }

    template <typename GridPosition>
    PhaseFunction get_phase_function(Grid & /*grid*/, GridPosition gp) {
      return opt_props_->get_phase_function(gp);
    }

   private:
    GriddedOpticalProperties *opt_props_;
  };

  ScatteringModel get_scattering_model() {
      return ScatteringModel(this);
  }

  GriddedOpticalProperties(Grid &grid, Vector<Float> scattering_angles,
                           Tensor<Float, 3> absorption_coefficient,
                           Tensor<Float, 3> scattering_coefficient,
                           Tensor<Float, 4> phase_function)
      : grid_(grid),
        scattering_angles_(scattering_angles),
        scattering_angle_map_(scattering_angles_.data(),
                              scattering_angles_.size()),
        absorption_coefficient_(absorption_coefficient),
      scattering_coefficient_(scattering_coefficient),
        phase_function_(phase_function) {}

  template <typename GridPosition>
  Float get_scattering_coefficient(GridPosition gp) {
      return scattering_coefficient_(gp.template get_index_array<Index>());
  }

  template <typename GridPosition>
  Float get_absorption_coefficient(GridPosition gp) {
      return absorption_coefficient_(gp.template get_index_array<Index>());
  }

  template <typename GridPosition>
  PhaseFun get_phase_function(GridPosition gp) {
      return PhaseFun(scattering_angle_map_, phase_function_(gp.template get_index_array<Index>()));
  }

 private:
  Grid &grid_;
  Vector<Float> scattering_angles_;
  ConstVectorMap<Float> scattering_angle_map_;
  Tensor<Float, 3> absorption_coefficient_;
  Tensor<Float, 3> scattering_coefficient_;
  Tensor<Float, 4> phase_function_;
};

}  // namespace llrte
#endif
