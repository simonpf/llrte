#include <llrte/surfaces.h>
#include <llrte/types/vector.h>

#include "utils.h"

int main(int args, const char **argv) {
  using V3 = llrte::Vector<3, float>;
  auto surface_base = V3{};
  surface_base[0] = 1.0;
  surface_base[1] = 0.0;
  surface_base[2] = 0.0;

  auto surface_normal = V3{};
  surface_normal[0] = -1.0;
  surface_normal[1] = 0.0;
  surface_normal[2] = 0.0;

  using Surface = llrte::surfaces::BlackPlane<V3>;
  auto surface = std::make_tuple(Surface(surface_base, surface_normal));
  auto test_setup = make_test_atmosphere(surface);
  auto atmosphere = std::get<0>(test_setup);
  auto source = std::get<1>(test_setup);

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
