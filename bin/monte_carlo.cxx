//#include <llrte/compiler/backends.h>
//#include <llrte/compiler.h>
#include <llrte/solvers/monte_carlo.h>
#include <llrte/grids.h>
#include <memory>

struct Deleter {
    template <typename T>
    void operator()(T *) {
    }
};

template<typename F>
std::shared_ptr<F[]> make_linear_vector(F start,
                                        F stop,
                                        size_t steps) {


    std::shared_ptr<F[]> v{new F[steps], Deleter()};

    F d = (stop - start) / steps;
    F x = start;
    for (size_t i = 0; i < steps;  ++i) {
        v[i] = x;
        x = x + d;
    }
    return v;
}

int main(int /*argc*/, const char **/***argv*/) {


    using V3 = llrte::Vector<3, float>;
    using Float = float;
    using Grid = llrte::RegularGrid<Float>;
    using AbsorptionModel = llrte::ConstantAbsorption<Float>;
    using Atmosphere = llrte::Atmosphere<Grid, AbsorptionModel>;
    using Source = llrte::PointSource<V3>;
    using Results = llrte::Histogram<Grid>;
    using Solver = llrte::MonteCarloSolver<Atmosphere, Source, Results>;

    auto source_position = V3{};
    source_position[0] = 0.0;
    source_position[0] = 0.0;
    source_position[0] = 0.0;

    auto source = llrte::PointSource<V3>(source_position);


    float start = -10.0e3;
    float stop = 10.0e3;
    auto x = make_linear_vector(start, stop, 201);
    auto y = make_linear_vector(start, stop, 201);
    auto z = make_linear_vector(start, stop, 201);
    size_t shape[3] = {201, 201, 201};

    auto grid = Grid{shape, x, y, z};
    auto absorption_model = llrte::ConstantAbsorption<float>(1e-9);
    auto atmosphere = Atmosphere{grid, absorption_model};
    auto results = Results{grid};

    auto solver = Solver(atmosphere, source, results);

    for (size_t i = 0; i < 1000000; i++) {
        solver.sample_photon();
    }
    results.dump("results.bin");
}
