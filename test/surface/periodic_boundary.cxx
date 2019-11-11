#include <llrte/surfaces.h>
#include <llrte/types/vector.h>

#include "utils.h"

int main(int args, const char **argv) {
  using V3 = llrte::Vector<3, float>;

  // Setup black boundary.
  auto base_b = V3{};
  base_b[0] = 1.0;
  base_b[1] = 0.0;
  base_b[2] = 0.0;

  auto normal_b = V3{};
  normal_b[0] = -1.0;
  normal_b[1] = 0.0;
  normal_b[2] = 0.0;

  // Setup periodic boundary.
  auto base_1_p = V3{};
  base_1_p[0] = 0.0;
  base_1_p[1] = -0.5;
  base_1_p[2] = 0.0;

  auto base_2_p = V3{};
  base_2_p[0] = 0.0;
  base_2_p[1] = 0.5;
  base_2_p[2] = 0.0;

  auto normal_p = V3{};
  normal_p[0] = 0.0;
  normal_p[1] = -1.0;
  normal_p[2] = 0.0;

  auto surfaces = std::make_tuple(
      llrte::surfaces::BlackPlane<V3>(base_b, normal_b),
      llrte::surfaces::PeriodicBoundary<V3> (base_1_p, base_2_p, normal_p)
      );

  auto test_setup = make_test_atmosphere(surfaces);
  auto atmosphere = std::get<0>(test_setup);

  auto source_direction = V3{};
  source_direction[0] = 1.0;
  source_direction[1] = -10.0;
  source_direction[2] = 0.0;
  auto source = std::get<1>(test_setup);
  source.set_direction(source_direction);

  using Solver =
      llrte::MonteCarloSolver<decltype(atmosphere) &, decltype(source) &>;
  Solver solver(atmosphere, source);
  solver.sample_photon();
  auto &boundary = atmosphere.get_boundary<0>();

  if (boundary.get_absorbed_energy() == 1.0) {
    return 0;
  } else {
    return 1;
  }
}
