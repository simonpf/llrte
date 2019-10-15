#ifndef _LLRTE_GRID_REGULAR_H_
#define _LLRTE_GRID_REGULAR_H_

#include <utility>
#include <limits>
#include <iostream>
#include "llrte/definitions.h"

namespace llrte {

template<typename Float,
    typename Index = size_t>
struct GridPosition {
    size_t i, j, k;
    Float x, y, z;
    Float w_i, w_j, w_k;

    bool is_inside() {
        return (w_i > 0.0) && (w_j > 0.0) && (w_k > 0.0);
    }

};

template <typename Real>
    std::ostream& operator<<(std::ostream& os, const GridPosition<Real>& gp)
{
    os << "[" << gp.x << ", " << gp.y << ", " << gp.z << "]" << std::endl;
    os << "[" << gp.i << ", " << gp.j << ", " << gp.k << "]" << std::endl;
    os << "[" << gp.w_i << ", " << gp.w_j << ", " << gp.w_j << "]" << std::endl;
    return os ;
}

template <typename F>
class RegularGrid {

public:

    using Float = F;


    RegularGrid(size_t shape[3],
                std::shared_ptr<Float[]> x,
                std::shared_ptr<Float[]> y,
                std::shared_ptr<Float[]> z)
        : x_(x), y_(y), z_(z)
    {
        for (size_t i = 0; i < 3; ++i){
            shape_[i] = shape[i];
        }
    }

    template <typename Backend>
    auto abstract() -> typename Backend::ClassType {
        auto cl = typename Backend::ClassType("RegularGrid");

        //auto fl = Float<Backend>{4};
        //auto array_of_float = Array<Float<Backend>, Backend>{fl};

        //cl.add_member(ClassAttribute("x", array_of_float));
        //cl.add_member(ClassAttribute("y", array_of_float));
        //cl.add_member(ClassAttribute("z", array_of_float));

        return cl;
    }

    template<typename Vector>
    GridPosition<Float> get_grid_position(Vector pos) {
        size_t i = 0;
        size_t j = 0;
        size_t k = 0;
        while ((x_[i] <= pos[0]) && (i < shape_[0])) {
            i++;
        }
        Float w_i = -1.0;
        if ((i > 0) && (i < shape_[0])) {
            w_i = (pos[0] - x_[i - 1]) / (x_[i] - x_[i - 1]);}

        while ((y_[j] <= pos[1]) && (j < shape_[1])) {
            j++;
        }
        Float w_j = -1.0;
        if ((j > 0) && (k < shape_[1])) {
            w_j = (pos[1] - y_[j - 1]) / (y_[j] - x_[j - 1]);
        }

        while ((z_[k] <= pos[2]) && (k < shape_[2])) {
            k++;
        }
        Float w_k = -1.0;
        if ((k > 0) && (k < shape_[2])) {
            w_k = (pos[2] - z_[k - 1]) / (z_[k] - z_[k - 1]);
        }
        return GridPosition<Float>{i, j, k, pos[0], pos[1], pos[2], w_i, w_j, w_k};
    }

    //template<typename Vector>
    //GridPosition<FloatType> get_intersection(VectorType pos,
    //                                         VectorType dir) {

    //}

    template <size_t axis>
    size_t first() const {
        if (axis == 0) {
            return x_[0];
        } else if (axis == 1) {
            return y_[1];
        } else if (axis == 2) {
            return z_[0];
        }
    }

    template <size_t axis>
    size_t last() const {
        if (axis == 0) {
            return x_[shape_[0] - 1];
        } else if (axis == 1) {
            return y_[shape_[1] - 1];
        } else if (axis == 2) {
            return z_[shape_[2] - 1];
        }
    }

    template <size_t axis>
    size_t first_index() const {
        return 0;
    }

    template <size_t axis>
    size_t last_index() const {
        if (axis == 0) {
            return shape_[0] - 1;
        } else if (axis == 1) {
            return shape_[1] - 1;
        } else if (axis == 2) {
            return shape_[2] - 1;
        }
    }

    template<size_t axis>
    std::tuple<Float, Float, Float, size_t, size_t> get_pos(GridPosition<Float> gp) {
        if (axis == 0) {
            if (gp.w_i > 0.0) {
                return std::make_tuple(gp.x, x_[gp.i - 1], x_[gp.i], gp.i - 1, gp.i);
            } else {
                return std::make_tuple(gp.x, first<0>(),
                                       last<0>(),
                                       first_index<0>(),
                                       last_index<0>());
            }
        } else if (axis == 1) {
            if (gp.w_j > 0.0) {
                return std::make_tuple(gp.y, y_[gp.j - 1], y_[gp.j], gp.j - 1, gp.j);
            } else {
                return std::make_tuple(gp.y, first<1>(), last<1>(), first_index<1>(), last_index<1>());
            }
        } else if (axis == 2) {
            if (gp.w_k > 0.0) {
                return std::make_tuple(gp.z, z_[gp.k - 1], z_[gp.k], gp.k - 1, gp.k);
            } else {
                return std::make_tuple(gp.z, first<2>(), last<2>(), first_index<2>(), last_index<2>());
            }
        }
    }

    template<typename Vector>
    std::pair<Float, GridPosition<Float>> get_intersection(GridPosition<Float> gp,
                                                           Vector dir) {
        Float d_min = std::numeric_limits<Float>::max();
        size_t axis = 4;
        size_t i_new = 0;

        Float x, x1, x2, d;
        size_t i1, i2;

        d = std::numeric_limits<Float>::max();
        std::tie(x, x1, x2, i1, i2) = get_pos<0>(gp);
        if (dir[0] > 0.0) {
            d = (x2 - x) / dir[0];
        } else if (dir[0] < 0.0) {
            d = (x1 - x) / dir[0];
        }
        if (d < d_min) {
            d_min = d;
            axis = 0;
        }

        d = std::numeric_limits<Float>::max();
        std::tie(x, x1, x2, i1, i2) = get_pos<1>(gp);
        if (dir[1] > 0.0) {
            d = (x2 - x) / dir[1];
        } else if (dir[1] < 0.0) {
            d = (x1 - x) / dir[1];
        }
        if (d < d_min) {
            d_min = d;
            axis = 0;
        }

        d = std::numeric_limits<Float>::max();
        std::tie(x, x1, x2, i1, i2) = get_pos<2>(gp);
        if (dir[2] > 0.0) {
            d = (x2 - x) / dir[2];
        } else if (dir[2] < 0.0) {
            d = (x1 - x) / dir[2];
        }
        if (d < d_min) {
            d_min = d;
            axis = 0;
        }


        GridPosition<Float> gp_new(gp);

        if (axis == 0) {
            gp_new.i = i_new;
            gp_new.w_i = 0.0;
        }

        if (axis == 1) {
            gp_new.j = i_new;
            gp_new.w_j = 0.0;
        }

        if (axis == 2) {
            gp_new.k = i_new;
            gp_new.w_k = 0.0;
        }
        return std::make_pair(d_min * dir.length(), gp_new);
    }

    void set_x(Float *x, size_t n) {
        shape_[0] = n;
        x_ = std::make_shared<Float[]>(n);
        std::copy(x, x + n, x_.get());
    }

    std::pair<Float, size_t> get_x() const {
        return std::make_pair(x_.get(), shape_[0]);
    }

    void set_y(Float *y, size_t n) {
        shape_[1] = n;
        y_ = std::make_shared<Float[]>(n);
        std::copy(y, y + n, y_.get());
    }

    std::pair<Float, size_t> get_y() const {
        return std::make_pair(y_.get(), shape_[1]);
    }

    void set_z(Float *z, size_t n) {
        shape_[2] = n;
        z_ = std::make_shared<Float[]>(n);
        std::copy(z, z + n, z_.get());
    }

    std::pair<Float, size_t> get_z() const {
        return std::make_pair(z_.get(), shape_[2]);
    }

private:

    size_t shape_[3];
    std::shared_ptr<Float[]> x_;
    std::shared_ptr<Float[]> y_;
    std::shared_ptr<Float[]> z_;

};

}
#endif
