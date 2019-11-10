#ifndef _LLRTE_ATMOSPHERE_H_
#define _LLRTE_ATMOSPHERE_H_

#include <llrte/utils/tuple.h>

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
      template <typename Boundary, typename Photon>
      void apply(Boundary b, Photon p) {
        if (b.has_crossed(p)) {
          b.apply(p);
        }
      }
    };
  }

/** Atmosphere
 *
 * This class represents the atmosphere in which an RT calculation is performed.
 * Its main role is to bundle all objects that make up the atmosphere, that means
 * determine its optical properties, and provides an interface for the photons
 * and beams that propagate through it.
 *
 */
template <
    typename Grid,
    typename AbsorptionModel,
    typename ScatteringModel,
    typename ... Boundaries
>
class Atmosphere {
 public:

  /** The floating point type used to represent scalars. */
  using Float = typename Grid::Float;
  /** The class representing the phase function. */
  using PhaseFunction = typename ScatteringModel::PhaseFunction;

  Atmosphere(Grid grid,
             AbsorptionModel absorption_model,
             ScatteringModel scattering_model,
             std::tuple<Boundaries ...> boundaries = std::tuple<Boundaries ...>())
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
  Float get_absorption_coefficient(GridPosition gp) {
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
  Float get_scattering_coefficient(GridPosition gp) {
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
  PhaseFunction get_phase_function(GridPosition gp) {
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
  auto get_grid_position(Vector position) {
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
  auto is_inside(GridPosition position) {
    return grid_.is_inside(position);
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
  template <typename Vector, typename GridPosition, typename Float>
  std::pair<Float, GridPosition> get_intersection(GridPosition gp,
                                                  Vector direction,
                                                  Float step_length) {
    return grid_.get_intersection(gp, direction, step_length);
  }


  template <typename GridPosition>
  size_t get_boundary_index(GridPosition gp) {
    return grid_.get_boundary_index(gp);
  }

  template <typename Photon>
  void apply_boundaries(Photon &p) {
      detail::BoundaryApplier ba{};
      tuple::map(ba, boundaries_);
  }

  template <size_t i>
  auto get_boundary() {
      return std::get<i>(boundaries_);
  }

 private:
  Grid grid_;
  AbsorptionModel absorption_model_;
  ScatteringModel scattering_model_;
  std::tuple<Boundaries ...> boundaries_;

};
}  // namespace llrte

#endif
