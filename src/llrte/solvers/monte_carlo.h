#ifndef _LLRTE_SOLVERS_MONTE_CARLO_H_
#define _LLRTE_SOLVERS_MONTE_CARLO_H_

#include <math.h>
#include <random>
#include <tuple>
#include <llrte/constants.h>
#include <iostream>

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

template<typename Atmosphere,
         typename Source>
class MonteCarloSolver {
public:

    using Float = typename Atmosphere::Float;

    MonteCarloSolver(Atmosphere atmosphere,
                     Source source)
        : atmosphere_(atmosphere), source_(source)
    {
        // Nothing to do here.
    }

    void sample_photon() {

        auto photon = source_.sample_photon();

        auto photon_position = photon.get_position();
        auto photon_direction = photon.get_direction();

        auto position = atmosphere_.get_grid_position(photon_position);
        std::cout << position << std::endl;

        while (true) {
            auto absorption = atmosphere_.get_absorption(position);
            auto intersection = atmosphere_.get_intersection(position,
                                                             photon_direction);

            std::cout << std::get<1>(intersection) << std::endl;
            std::cout << photon_direction << std::endl;

            if (!std::get<1>(intersection).is_inside()) {
                break;
            }

            auto d = std::get<0>(intersection);
            if (d > sample_path_length(absorption)) {
                break;
            }

            atmosphere_.trace(position);
            position = std::get<1>(intersection);

        }
        //bool graveyard = false;
    }

    Float sample_path_length(Float m) {
        auto y = distribution(generator);
        return -m * log(1 - y);
    }

private:
    size_t n_photons = 0;
    std::default_random_engine generator;
    std::uniform_real_distribution<Float> distribution{0.0,1.0};

    Atmosphere atmosphere_;
    Source source_;

};

}

#endif
