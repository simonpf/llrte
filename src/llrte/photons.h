#ifndef _LLRTE_PHOTONS_
#define _LLRTE_PHOTONS_

#include "llrte/tracers.h"
#include "math.h"

namespace llrte {

/** Photon
 * This class implements a basic photon template class.
 * @tparam V The vector class to use to represent 3D vectors.
 * @tparam The grid on which the photon will live.
 */
template <
    typename V,
    template <typename> typename GP
    >
class Photon : public GP<V> {
public:

  using Vector = V;
  using GridPosition = GP<V>;
  using Float = typename V::Float;

__DEV__ Photon(GridPosition grid_position)
    : GridPosition(grid_position),
      n_scattered_(0),
      energy_(1.0) {
    // Nothing to do here.
  }

  __DEV__ Photon & operator=(GridPosition gp) {
      GridPosition::operator=(gp);
      return *this;
  }

  /**
   * Scatter photon.
   * @param generate The random number generator to use.
   * @param p Phase function object representing the scattering
   * properties.
   */
  template <typename Generator, typename PhaseFunction>
  __DEV__ void scatter(Generator& generator, const PhaseFunction& p) {
    p.scatter(generator, *this);
    n_scattered_++;
  }

  /** Return number of scattering events the photon has undergone. */
  size_t get_scattering_events() const { return n_scattered_; }

  /** Return current photon energy.
   * @return The energy of the photon.
   **/
  __DEV__ Float get_energy() const { return energy_; }
  /** Set photon energy.
   * @param e Energy to set photon energy to
   */
  void set_energy(Float e) {energy_ = e; }
  /** Scale the photon energy
   * @param f The scaling factor.
   **/
  __DEV__ void scale_energy(Float f) { energy_ *= f;}

 private:
  size_t n_scattered_ = 0;
  Float energy_ = 1.0;
};

template <typename Vector, template <typename> typename GridPosition>
std::ostream& operator<<(std::ostream& os,
                         const Photon<Vector, GridPosition>& p) {
    os << "Photon(e = " << p.get_energy() << ") @ " << p.position << " -> " << p.direction << " (" << p.direction.length() << ")" << std::endl;
    os << (GridPosition<Vector>) p << std::endl;
    return os;
}

}  // namespace llrte
#endif
