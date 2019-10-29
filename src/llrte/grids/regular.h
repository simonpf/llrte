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
};

template <typename Real>
    std::ostream& operator<<(std::ostream& os, const GridPosition<Real>& gp)
{
    os << "[" << gp.x << ", " << gp.y << ", " << gp.z << "]" << std::endl;
    os << "[" << gp.i << ", " << gp.j << ", " << gp.k << "]" << std::endl;
    return os ;
}

template <typename F>
class RegularGrid {

public:

    using Index = size_t;

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

    std::tuple<Index, Index, Index> get_extent() const {
        return std::make_tuple(shape_[0], shape_[1], shape_[2]);
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
        while ((y_[j] <= pos[1]) && (j < shape_[1])) {
            j++;
        }
        while ((z_[k] <= pos[2]) && (k < shape_[2])) {
            k++;
        }
        return GridPosition<Float>{i, j, k, pos[0], pos[1], pos[2]};
    }

    template <typename GridPosition>
    bool is_inside(GridPosition gp) {
        bool inside = ((gp.x > x_[0]) && (gp.x < x_[shape_[0] - 1]));
        inside &= ((gp.y > y_[0]) && (gp.y < y_[shape_[1] - 1]));
        inside &= ((gp.z > z_[0]) && (gp.z < z_[shape_[2] - 1]));
        return inside;
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
    size_t get_lower(GridPosition<Float> gp) {
        if (axis == 0) {
            return std::max<size_t>(gp.i - 1, 0);
        }
        if (axis == 1) {
            return std::max<size_t>(gp.j - 1, 0);
        }
        if (axis == 2) {
            return std::max<size_t>(gp.k - 1, 0);
        }
    }

    template<size_t axis>
    size_t get_higher(GridPosition<Float> gp) {
        if (axis == 0) {
            return std::min(gp.i, shape_[0] - 1);
        }
        if (axis == 1) {
            return std::min(gp.j, shape_[1] - 1);
        }
        if (axis == 2) {
            return std::min(gp.k, shape_[2] - 1);
        }
    }

    template<typename Vector>
    std::pair<Float, GridPosition<Float>> get_intersection(GridPosition<Float> gp,
                                                           Vector dir) {

        float d = std::numeric_limits<Float>::max();
        size_t direction = 0;

        float dx = -1.0;
        size_t i = 0;

        if (dir[0] < 0.0) {
            if (gp.i > 0) {
                i = get_lower<0>(gp);
                dx = (x_[i] - gp.x) / dir[0];
            }
        } else {
            if (gp.i < shape_[0]) {
                i = get_higher<0>(gp);
                dx = (x_[i] - gp.x) / dir[0];
                i++;
            }
        }

        if ((dx >= 0.0) && (dx < d)) {
            direction = 0;
            d = dx;
        }

        float dy = -1.0;
        size_t j = 0;

        if (dir[1] < 0.0) {
            if (gp.j > 0) {
                j = get_lower<1>(gp);
                dy = (y_[j] - gp.y) / dir[1];
            }
        } else {
            if (gp.j < shape_[1]) {
                j = get_higher<1>(gp);
                dy = (y_[j] - gp.y) / dir[1];
                j++;
            }
        }

        if ((dy >= 0.0) && (dy < d)) {
            direction = 1;
            d = dy;
        }

        float dz = -1.0;
        size_t k = 0;

        if (dir[2] < 0.0) {
            if (gp.k > 0) {
                k = get_lower<2>(gp);
                dz = (z_[k] - gp.z) / dir[2];
            }
        } else {
            if (gp.k < shape_[2]) {
                k = get_higher<2>(gp);
                dz = (z_[k] - gp.z) / dir[2];
                k++;
            }
        }

        if ((dz >= 0.0) && (dz < d)) {
            direction = 2;
            d = dz;
        }

        GridPosition<Float> gp_new(gp);

        gp_new.x = gp.x + d * dir[0];
        gp_new.y = gp.y + d * dir[1];
        gp_new.z = gp.z + d * dir[2];

        if (direction == 0) {
            gp_new.i = i;
        }
        if (direction == 1) {
            gp_new.j = j;
        }
        if (direction == 2) {
            gp_new.k = k;
        }
        return std::make_pair(d * dir.length(), gp_new);
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
