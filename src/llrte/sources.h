#ifndef _LLRTE_SOURCES_H_
#define _LLRTE_SOURCES_H_

#include <random>

#include "llrte/constants.h"
#include "llrte/rotations.h"
#include "llrte/geometry.h"

namespace llrte {

////////////////////////////////////////////////////////////////////////////////
// Isotropic point source
////////////////////////////////////////////////////////////////////////////////
/**
 * Isotropic point source
 * As the name suggests this point source emits photons in all direction
 * with equal probability.
 */
template <typename P>
class PointSource {
 public:
  using Photon = P;
  using Vector = typename Photon::Vector;
  using Float = typename Vector::Float;
  using C = Constants<Float>;

  /**
   * Create point source at given position
   * @param position The position at which to place the source. */
  PointSource(Vector position) : position_(position) {}

  template<typename Generator, typename Atmosphere>
      Photon sample_photon(Generator &generator,
                           Atmosphere &atmosphere) {
    auto theta = generator.sample_uniform_angle();
    auto phi = generator.sample_zenith_angle();
    Vector v = Vector{cos(phi) * cos(theta),
                      cos(phi) * sin(theta),
                      sin(phi)};
    return Photon(atmosphere.place_on_grid(position_, v));
  }

 private:
 Vector position_;
};

////////////////////////////////////////////////////////////////////////////////
// Beam source
////////////////////////////////////////////////////////////////////////////////
/**
 * A perfectly collimated source.
 *
 * Beam source that emits a perfectly collimated beam of photons.
 */
template <typename P>
class BeamSource {
 public:
  using Photon = P;
  using Vector = typename Photon::Vector;
  using Float = typename Vector::Float;
  using C = Constants<Float>;

  /**
   * Create beam source.
   * @param position The position at which to place the source.
   * @param direction The direction in which to emits photons.
   */
BeamSource(Vector position, Vector direction) : position_(position), direction_(direction) {}

  /**
   * Emit photon in given direction and place on atmosphere grid.
   * @param generator Not used, process is deterministic
   * @param atmosphere The atmosphere into which to place the photon.
   */
  template<typename Generator, typename Atmosphere>
      Photon sample_photon(Generator /*&generator*/,
                           Atmosphere &atmosphere) {
      return Photon(atmosphere.place_on_grid(position_, direction_));
  }

private:
  Vector position_, direction_;
};

////////////////////////////////////////////////////////////////////////////////
// Random offset
////////////////////////////////////////////////////////////////////////////////
/** RandomOffset
 *
 * A 1D random offset that is applied to the photon position after creation.
 */
template <typename Source>
class RandomOffset : public Source {
public:

  using Float   = typename Source::Float;
  using Vector = typename Source::Vector;
  using Photon = typename Source::Photon;

  /**
   * Create source with given random offset. The random is created
   * by sampling a scaling factor uniformly from the interval [limit_low, limit_high] and
   * adding the scaled direction vector to the photon.
   *
   * @param direction Vector defining direction of offset.
   * @param limit_low Lower bound on
   * @param limit_high
   * @param ts Remaining parameters are passed on to the base class
   * constructor.
   */
  template <typename ... Ts>
  RandomOffset(Vector direction,
               Float limit_low,
               Float limit_high,
               Ts ... ts)
      : Source(std::forward<Ts>(ts) ...),
      direction_(direction),
      limit_low_(limit_low),
      limit_high_(limit_high)
  {}

  /**
   * Create photon from base source and add offset.
   *
   * Note that the photon is currently placed on grid twice,
   * so this is probably not optimal. So this is mainly for testing purposes.
   *
   * @param generator Random number generator to use to generate offset.
   * @param atmosphere The atmosphere into which to place the photon.
   */
  template<typename Generator, typename Atmosphere>
      Photon sample_photon(Generator &generator,
                           Atmosphere &atmosphere) {
      auto p = Source::sample_photon(generator, atmosphere);
      auto s = generator.sample_uniform(limit_low_, limit_high_);
      p.position = p.position + s * direction_;
      p = Photon(atmosphere.place_on_grid(p.position, p.direction));
      return p;
  }

 private:
  Vector direction_;
  Float limit_low_ = 0.0;
  Float limit_high_ = 0.0;
};

/** RandomDirection
 *
 * Randomize direction from source.
 */
template <typename Source>
class RandomDirection : public Source {
public:

    using Float   = typename Source::Float;
    using Vector = typename Source::Vector;
    using Photon = typename Source::Photon;

    /**
     * Create random direction modifier.
     *
     * @param ts Remaining parameters are passed on to the base class
     */
    template <typename ... Ts>
        RandomDirection(Float theta_max,
                        Ts ... ts)
        : Source(ts ...), theta_max_(theta_max)
    {}

    /**
     * Create photon from base source apply random rotation offset.
     *
     * Note that the photon is currently placed on grid twice,
     * so this is probably not optimal. So this is mainly for testing purposes.
     *
     * @param generator Random number generator to use to generate offset.
     * @param atmosphere The atmosphere into which to place the photon.
     */
    template<typename Generator, typename Atmosphere>
        Photon sample_photon(Generator &generator,
                             Atmosphere &atmosphere) {
        auto p = Source::sample_photon(generator, atmosphere);
        auto d = p.direction;
        auto n = geometry::RandomPlane::get_normal(generator, d);
        auto theta = generator.sample_uniform(0, theta_max_);
        auto dn = rotations::rotate(d, n, theta);
        p.set_direction(dn);

        return p;
    }

private:
    Float theta_max_ = 0.0;
};

////////////////////////////////////////////////////////////////////////////////
// Planar source
////////////////////////////////////////////////////////////////////////////////
/** Planar Source
 *
 * This is a planar source, which implements the backwards interface, i.e.
 * it provides a get_itensity method, which can be used to determin source
 * irradiance given an outgoing backwards photon.
 */
template <typename P>
class PlanarSource {
 public:
  using Photon = P;
  using Vector = typename Photon::Vector;
  using Float = typename Vector::Float;
  using C = Constants<Float>;

  /**
   * Create a collimated planar source with given intensity normal vector
   * and beam width.
   * @param intensity Spectral intensity of the source.
   * @param Normal of the plane
   * @param Beam width of source.
   */
  PlanarSource(Float intensity, Vector normal,
               Float dtheta = 2.0 * Constants<Float>::pi / 180)
      : domega_(2.0 * dtheta / C::pi),
        dtheta_(dtheta),
        intensity_(intensity),
        normal_(normal * -1.0) {}

  /**
   * Get intensity for outgoing photon.
   * @param p The outgoing photon.
   */
  Float get_intensity(const Photon &p) {
    if (angle(p.get_direction(), normal_) < dtheta_) {
      return 1.0 / domega_;
    } else {
      return 0.0;
    }
  }

 private:
  Float domega_, dtheta_, intensity_;
  Vector normal_;
};

}  // namespace llrte
#endif
