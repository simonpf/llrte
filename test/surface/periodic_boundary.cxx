#include <llrte/surfaces.h>
#include <llrte/types/vector.h>

#include "utils.h"

#include <llrte/surfaces.h>
#include <llrte/types/vector.h>

#include "utils.h"

void periodic_boundary() {

    using Float = float;
    using V3 = llrte::Vector3<float>;

    auto p1_base1 = V3{0.0, -0.5, 0.0};
    auto p1_base2 = V3{0.0, 0.5, 0.0};
    auto p1_normal = V3{0.0, 1.0, 0.0};
    auto p1_normal2 = V3{0.0, -1.0, 0.0};

    auto p2_base1 = V3{0.0, 0.0, -0.5};
    auto p2_base2 = V3{0.0, 0.0, 0.5};
    auto p2_normal = V3{0.0, 0.0, 1.0};
    auto p2_normal2 = V3{0.0, 0.0, -1.0};

    long n = 1000000;

    using Tracer = llrte::tracers::PhotonTracer<Float>;
    using Surface = llrte::surfaces::PeriodicBoundary<V3>;
    auto surface = std::make_tuple(Surface(p1_base1, p1_base2, p1_normal),
                                   Surface(p2_base1, p2_base2, p2_normal));

    Tracer tracer{n};
    auto test_setup = make_test_atmosphere_random_source(surface, tracer);
    auto solver = std::get<0>(test_setup);
    auto source = std::get<1>(test_setup);

    for (long i = 0; i < n; ++i) {
        solver.forward(source);
    }
    solver.tracer().save("surface_periodic.nc");
}

int main(int /*args*/, const char **/*argv*/) {
    periodic_boundary();
}
