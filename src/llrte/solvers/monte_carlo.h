#ifndef _LLRTE_SOLVERS_MONTE_CARLO_H_
#define _LLRTE_SOLVERS_MONTE_CARLO_H_

#include <math.h>
#include <random>
#include <tuple>
#include <llrte/constants.h>
#include <iostream>
#include <fstream>
#include <memory>

namespace llrte {

template <size_t N, typename F>
class Vector {
public:

    using Float = F;

    Vector () {
        for (size_t i = 0; i < N; ++i) {
            elements_[i] = 0.0;
        }
    }

    template <typename ... Ts>
    Vector(Ts ... ts) {
        auto t = std::make_tuple(ts ...);
        for (size_t i = 0; i < N; ++i) {
            //elements_[i] = std::get<i>(t);
        }
    }

    Float operator[](size_t i) const {
        return elements_[i];
    }

    Float& operator[](size_t i) {
        return elements_[i];
    }

    Vector operator+(const Vector &v) {
        Vector w;
        for (size_t i = 0; i < N; ++i) {
            w[i] = v[i] + elements_[i];
        }
    }

    Vector operator-(const Vector &v) {
        Vector w;
        for (size_t i = 0; i < N; ++i) {
            w[i] = v[i] + elements_[i];
        }
    }

    Float length() {
        Float s = 0.0;
        for (size_t  i = 0; i < N; ++i) {
            s += elements_[i] * elements_[i];
        }
        return sqrt(s);
    }

public:
    Float elements_[N];
};

template <size_t N, typename Real>
    std::ostream& operator<<(std::ostream& os, const Vector<N, Real>& v)
{
    os << "[";
    for (size_t i = 0; i < N - 1; ++i) {
        os << v[i] << ",";
    }
    os << v[N - 1] << "]" << std::endl;
    return os;
}


template <typename Vector>
class Photon {
public:

    Photon(Vector position,
           Vector direction)
        : position_(position), direction_(direction) {
        // Nothing to do here.
    }

    Vector get_position() {
        return position_;
    }

    Vector get_direction() {
        return direction_;
    }


private:
    Vector position_;
    Vector direction_;
};

template <typename V>
class PointSource {

public:

    using Vector = V;
    using Float = typename V::Float;
    using C = Constants<Float>;

    std::default_random_engine generator;
    std::uniform_real_distribution<Float> theta_d{-C::pi, C::pi};
    std::uniform_real_distribution<Float> phi_d{-C::pi / 2.0, C::pi};

    PointSource(Vector position) :
    position_(position) {


    }

    Photon<Vector> sample_photon() {
        auto theta = theta_d(generator);
        auto phi = phi_d(generator);
        Vector v = Vector{};
        v[0] = cos(phi) * cos(theta);
        v[1] = cos(phi) * sin(theta);
        v[2] = sin(phi);
        return Photon<Vector>(position_, v);
    }

public:
    Vector position_;
};

template <typename Float>
class ConstantAbsorption {
public:
    ConstantAbsorption(Float absorption)
        : absorption_(absorption){}

    template <typename ... Ts>
    Float get_absorption(Ts ...) {
        return absorption_;
    }

private:
    Float absorption_;
};

template<typename Grid,
         typename AbsorptionModel>
class Atmosphere {

    public:

    using Float = typename Grid::Float;

    Atmosphere(Grid grid,
               AbsorptionModel absorption_model)
        : grid_(grid), absorption_model_(absorption_model) {
        // Nothing to do here.
    }

    template<typename GridPosition>
    Float get_absorption(GridPosition gp) {
        return absorption_model_.get_absorption(grid_, gp);
    }

    template <typename Vector>
        auto get_grid_position(Vector position)
    {
        return grid_.get_grid_position(position);
    }

    template <typename GridPosition>
    auto is_inside(GridPosition position)
        {
            return grid_.is_inside(position);
        }

    template <typename Vector, typename GridPosition>
    std::pair<Float, GridPosition> get_intersection(GridPosition gp,
                                                    Vector direction)
    {
        return grid_.get_intersection(gp, direction);
    }

    template<typename T>
    void trace(T) {
        
    }


public:
    Grid grid_;
    AbsorptionModel absorption_model_;
};

template <typename Grid>
class Histogram {
public:

    using Index = typename Grid::Index;

    Histogram(const Grid &grid) {
        std::tie(shape_[0], shape_[1], shape_[2]) = grid.get_extent();
        size_t n = shape_[0] * shape_[1] * shape_[2];
        data_ = std::shared_ptr<Index[]>(new Index[n]);
        for (size_t i = 0; i < n; ++i) {
            data_[i] = 0;
        }
    }

    template<typename GridPos>
    void trace(GridPos gp) {
        size_t index = gp.i * gp.j * gp.k;
        data_[index] += 1;
    }

    void dump(std::string filename) {
        std::ofstream file;
        file.open (filename, std::ios::out | std::ios::binary); 
        size_t n = shape_[0] * shape_[1] * shape_[2];
        for (size_t i = 0; i < n; ++i) {
            file << data_[i];
        }
    }

private:

    Index shape_[3] = {0, 0, 0};
    std::shared_ptr<Index[]> data_ = nullptr;
};

template<typename Atmosphere,
         typename Source,
         typename Results>
class MonteCarloSolver {
public:

    using Float = typename Atmosphere::Float;

    MonteCarloSolver(Atmosphere atmosphere,
                     Source source,
                     Results results)
        : atmosphere_(atmosphere), source_(source), results_(results)
    {
        // Nothing to do here.
    }

    void sample_photon() {

        auto photon = source_.sample_photon();

        auto photon_position = photon.get_position();
        auto photon_direction = photon.get_direction();

        auto position = atmosphere_.get_grid_position(photon_position);

        while (true) {
            auto absorption = atmosphere_.get_absorption(position);
            auto intersection = atmosphere_.get_intersection(position,
                                                             photon_direction);

            if (!atmosphere_.is_inside(std::get<1>(intersection))) {
                break;
            }

            auto d = std::get<0>(intersection);
            if (d > sample_path_length(1.0 / absorption)) {
                break;
            }

            results_.trace(position);
            position = std::get<1>(intersection);

        }
        //bool graveyard = false;
    }

    Float sample_path_length(Float m) {
        auto y = distribution(generator);
        return -m * log(y);
    }

private:
    size_t n_photons = 0;
    std::default_random_engine generator;
    std::uniform_real_distribution<Float> distribution{0.0,1.0};

    Atmosphere atmosphere_;
    Source source_;
    Results results_;

};

}

#endif
