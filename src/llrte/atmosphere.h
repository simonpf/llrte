#ifndef _LLRTE_ATMOSPHERE_H_
#define _LLRTE_ATMOSPHERE_H_

namespace llrte {

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


public:
    Grid grid_;
    AbsorptionModel absorption_model_;
};
}

#endif
