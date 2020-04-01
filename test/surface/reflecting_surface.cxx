#include <llrte/surfaces.h>
#include <llrte/types/vector.h>

#include "utils.h"

int main(int /*args*/, const char **/*argv*/) {
    using Float = float;
  using V3 = llrte::Vector3<float>;

  size_t n = 1;
  auto surface_base = V3{1.0, 0.0, 0.0};
  auto surface_normal = V3{-1.0, -1.0, 0.0};

  using Tracer = llrte::tracers::PhotonTracer<Float>;
  using Surface = llrte::surfaces::ReflectingPlane<V3, llrte::surfaces::Lambertian<>>;
  auto surface = std::make_tuple(Surface(surface_base, surface_normal, 0.5));

  Tracer tracer{n};
  auto test_setup = make_test_atmosphere(surface, tracer);
  auto solver = std::get<0>(test_setup);
  auto source = std::get<1>(test_setup);

  for (size_t i = 0; i < n; ++i) {
      solver.forward(source);
  }
  solver.tracer().save("surface_lambertian.nc");
}
