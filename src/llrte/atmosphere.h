#ifndef _LLRTE_ATMOSPHERE_H_
#define _LLRTE_ATMOSPHERE_H_

#include "llrte/common.h"
#include "llrte/utils/tuple.h"

namespace llrte {

namespace detail {
/** Boundary acting on photon.
 *
 * Dummy struct that performs the action of a boundary on
 * a photon:
 * - Checks whether photon has crossed boundary
 * - If so, applies action of boundary on photon
 */
struct BoundaryApplier {
  template <typename Boundary, typename Generator, typename Photon>
  void apply(Boundary &b, Generator &g, Photon &p) {
    if (b.has_crossed(p)) {
      b.apply(g, p);
      hit = true;
    }
  }
  bool hit = false;
};
}  // namespace detail

/** Atmosphere
 *
 * This class represents the atmosphere in which an RT calculation is performed.
 * Its main role is to bundle all objects that make up the atmosphere, that
 * means determine its optical properties, and provides an interface for the
 * photons and beams that propagate through it.
 *
 */
template <typename Grid, typename AbsorptionModel, typename ScatteringModel,
          typename Boundaries = std::tuple<> >
class Atmosphere {
 public:
  /** The floating point type used to represent scalars. */
 using Float = typename std::remove_reference<typename std::remove_reference<Grid>::type::Float>::type;
  /** The class representing the phase function. */
  using PhaseFunction = typename ScatteringModel::PhaseFunction;

  Atmosphere(Grid grid,
             AbsorptionModel absorption_model,
             ScatteringModel scattering_model,
             Boundaries boundaries = Boundaries{})
 : grid_(grid),
        absorption_model_(absorption_model),
        scattering_model_(scattering_model),
        boundaries_(boundaries) {
    // Nothing to do here.
  }

  /** Compute absorption coefficient for given position.
   *
   * Calls absorption model to obtain absorption coefficient
   * for the current position.
   * @tparam GridPosition the griposition type representing position
   * in the grid.
   * @param gp The position in the atmospheric grid
   * @return the absorption coefficient at gp
   */
  template <typename GridPosition>
  __DEV__ Float get_absorption_coefficient(GridPosition gp) {
    return absorption_model_.get_absorption_coefficient(grid_, gp);
  }

  /** Compute scattering coefficient for given position.
   *
   * Calls scattering model to obtain scattering coefficient
   * for the current position.
   *
   * @tparam GridPosition the griposition type representing position
   * in the grid.
   * @param gp The position in the atmospheric grid
   * @return the scatttering coefficient at gp
   */
  template <typename GridPosition>
  __DEV__ Float get_scattering_coefficient(GridPosition gp) {
    return scattering_model_.get_scattering_coefficient(grid_, gp);
  }

  /** Return phase function a current position.
   *
   * Return a phase function object representing the scattering phase
   * function a the current position.
   *
   * @tparam GridPosition the griposition type representing position
   * in the grid.
   * @param gp The position in the atmospheric grid
   * @return the phase function at gp
   */
  template <typename GridPosition>
  __DEV__ PhaseFunction get_phase_function(GridPosition gp) {
    return scattering_model_.get_phase_function(grid_, gp);
  }

  /** Locate position in grid
   *
   * Compute the grid position of a given Vector in the atmosphere
   * grid.
   *
   * @tparam Vector The vector type representing a position in R^n
   * @param position The position in R^n
   * @return A gridposition object representing the position in the
   * grid.
   */
  template <typename Vector>
  __DEV__ auto get_grid_position(Vector &position) {
    return grid_.get_grid_position(position);
  }

  /** Determine whether position is in grid.
   *
   * @tparam GridPosition The type used to represent positions in the
   * atmospheric grid.
   * @param position The position
   * @return Boolean value indicating whether the current position is
   * in the grid.
   */
  template <typename GridPosition>
  __DEV__ auto is_inside(GridPosition position) {
    return grid_.is_inside(position);
  }

    template <typename Photon>
    __DEV__ bool is_leaving(Photon &photon) {
        return grid_.is_leaving(photon);
   }
  /** Compute next intersection with grid.
   *
   * Computes the position where the photon will hit the next grid
   * boundary or, if this distance is longer than the given step
   * length, the position that is at a distance of step length
   * from the current position in the given direction.
   *
   * @tparam Vector The type used to represent directions.
   * @tparam GridPosition The type used to represent positions in the
   * @tparam Float The type used to represent floating point numbers.
   * @param gp The current position.
   * @param direction The direction in which to move.
   * @param step_length The maximum step length.
   * @return Pair containing the length of the performed step
   * and the new position.
   */
  template <typename Photon>
  __DEV__ Float step(Photon &photon, Float step_length) {
      return grid_.step(photon, step_length);
  }

  template <typename GridPosition>
  __DEV__ size_t get_boundary_index(GridPosition gp) const {
    return grid_.get_boundary_index(gp);
  }

  template <typename Generator, typename Photon>
  __DEV__ bool apply_boundaries(Generator &g, Photon &p) {
    detail::BoundaryApplier ba{};
    tuple::map(ba, boundaries_, g, p);
    return ba.hit;
  }

  template <size_t i>
  __DEV__ auto get_boundary() -> typename std::tuple_element<i, Boundaries>::type & {
    return std::get<i>(boundaries_);
  }

 template <typename Vector>
 __DEV__ auto place_on_grid(Vector position,
                            Vector direction) {
     return grid_.place_on_grid(position, direction);
 }

 #ifdef CUDA
 void device() {
     grid_.device();
 }
 void host() {
     grid_.host();
 }
 #endif

 private:
  Grid grid_;
  AbsorptionModel absorption_model_;
  ScatteringModel scattering_model_;
  Boundaries boundaries_;
};
}  // namespace llrte

#endif
