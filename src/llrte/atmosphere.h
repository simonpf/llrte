#ifndef _LLRTE_ATMOSPHERE_H_
#define _LLRTE_ATMOSPHERE_H_

namespace llrte {

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

  template <typename GridPosition>
  Float get_absorption_coefficient(GridPosition gp) {
    return absorption_model_.get_absorption_coefficient(grid_, gp);
  }

  template <typename GridPosition>
  Float get_scattering_coefficient(GridPosition gp) {
    return scattering_model_.get_scattering_coefficient(grid_, gp);
  }

  template <typename GridPosition>
  PhaseFunction get_phase_function(GridPosition gp) {
    return scattering_model_.get_phase_function(grid_, gp);
  }

  template <typename Vector>
  auto get_grid_position(Vector position) {
    return grid_.get_grid_position(position);
  }

  template <typename GridPosition>
  auto is_inside(GridPosition position) {
    return grid_.is_inside(position);
  }

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

 private:
  Grid grid_;
  AbsorptionModel absorption_model_;
  ScatteringModel scattering_model_;
  std::tuple<Boundaries ...> boundaries_;

};
}  // namespace llrte

#endif
