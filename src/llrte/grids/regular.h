#ifndef _LLRTE_GRID_REGULAR_H_
#define _LLRTE_GRID_REGULAR_H_

#include <iostream>
#include <limits>
#include <memory>
#include <utility>

#include "llrte/common.h"
#include "llrte/data.h"
#include "llrte/maths.h"

namespace llrte {

//////////////////////////////////////////////////////////////////////////////////////////
// Grid Position
//////////////////////////////////////////////////////////////////////////////////////////

/** GridPosition
 *
 * The grid position struct represents a position in the grid.
 * The i, j, k indices point to the next grid boundary that the photon will
 * traverse given its current moving direction. The x, y, z values hold the
 * current position holds the indices of the grid cell that the photon is in
 * using 1-based indexing.
 *
 * @tparam Vector The vector type to use to represent positions and directions.
 * @tparam Index The integer type to use to represent indices.
 */
template <
    typename Vector,
    typename Index = short unsigned int
    >
struct GridPosition {
  using Float = typename Vector::Float;

  /*! 3D Vector representing the current position. */
  Vector position;
  /*! 3D Vector representing the current direction. */
  Vector direction;
  /*! 1-based index for next grid boundary in x direction. */
  Index i;
  /*! 1-based index for next grid boundary in y direction. */
  Index j;
  /*! 1-based index for next grid boundary in z direction. */
  Index k;

  void change_direction(Vector d) {
      if ((d.x * direction.x) < 0) {
          if (d.x < 0.0) {
              --i;
          } else if (d.x > 0.0) {
              ++i;
          }
      }
      if ((d.y * direction.y) < 0) {
          if (d.y < 0.0) {
              --j;
          } else if (d.y > 0.0) {
              ++j;
          }
      }
      if ((d.z * direction.z) < 0) {
          if (d.z < 0.0) {
              --k;
          } else if (d.z > 0.0) {
              ++k;
          }
      }
      direction = d;
  }

  /*! x coordinate of current position*/
  Float &x() { return position.x; }
  /*! y coordinate of current position*/
  Float &y() { return position.y; }
  /*! z coordinate of current position*/
  Float &z() { return position.z; }
  const Float &x() const { return position.x; }
  const Float &y() const { return position.y; }
  const Float &z() const { return position.z; }
};

template <typename Float, typename Index>
std::ostream &operator<<(std::ostream &os,
                         const GridPosition<Float, Index> &gp) {
  os << "[" << gp.x() << ", " << gp.y() << ", " << gp.z() << "] :: ";
  os << "[" << gp.i << ", " << gp.j << ", " << gp.k << "]" << std::endl;
  return os;
}

//////////////////////////////////////////////////////////////////////////////////////////
// Regular Grid
//////////////////////////////////////////////////////////////////////////////////////////

/** RegularGrid
 *
 * The regular grid class implements a regular, Cartesian grid. It is defined by three
 * arrays, holding the grid boundaries in each dimension. It provides functionality
 * to place objects on the grid and let them move across the grid.
 *
 * @tparam Float The floating point type to use for real numbers.
 * @tparam Index The integer type to use for indices.
 */
template <
    typename F,
    typename I = short unsigned int>
class RegularGrid {
 public:
 using Index = I;
 using Float = F;
  /**
   * Create a regular grid given arrays containing the edges in x-, y-, and z-
   * direction.
   *
   * @param x Array holding cell boundaries along the x-dimension
   * @param y Array holding cell boundaries along the y-dimension
   * @param z Array holding cell boundaries along the z-dimension
   *
   */
  RegularGrid(Array<Float> &&x,
              Array<Float> &&y,
              Array<Float> &&z)
 : x_(x), y_(y), z_(z) {}

 std::array<Index, 3> get_extent() const {
     return std::array{static_cast<Index>(x_.size()),
             static_cast<Index>(y_.size()),
             static_cast<Index>(z_.size())};
  }

  /**
   * Finds index of the grid boundary that a particle placed
   * at a given position would intersect next.
   *
   * If d is zero, movement into left (negative) direction is
   * assumed.
   *
   * @param grid Array holding boundaries
   * @param p The current position on the grid
   * @param d The moving direction
   * @return The 1-based index of the intersected boundary.
   */
  __DEV__ Index find_next_index(const Array<Float> &grid,
                        Float p,
                        Float d) {
    Index i = 0;
    // Find index, which is left of p
    while ((grid[i] < p) && (i < grid.size())) {
      ++i;
    }
    // If moving to right, increase index.
    if (d > 0.0) {
      ++i;
    }
    return i;
  }

  /**
   * Place an object on the grid by determining the next intersection
   * points of its given position and direction.
   *
   * @tparam Vector The vector type used to encode positions and moving
   * and moving directions.
   * @param position The position of the object to place on the grid.
   * @param direction The direction of the object to place on the grid.
   * @return GridPosition object representing the object on the grid.
   */
  template <typename Vector>
  __DEV__ GridPosition<Vector, Index> place_on_grid(const Vector &position,
                                            const Vector &direction) {
    Index i = find_next_index(x_, position.x, direction.x);
    Index j = find_next_index(y_, position.y, direction.y);
    Index k = find_next_index(z_, position.z, direction.z);

    return GridPosition<Vector, Index>{position, direction, i, j, k};
  }

  /**
   * Determine whether given grid position is inside the grid.
   *
   * @tparam Vector The vector type to represent positions and directions.
   * @param gp The grid position representing the object on the grid.
   * @return Boolean indicating whether the object is on the grid.
   */
  template <typename Vector>
  bool is_inside(GridPosition<Vector, Index> gp) {
    bool inside = ((gp.x > x_[0]) && (gp.x < x_[x_.size() - 1]));
    inside &= ((gp.y > y_[0]) && (gp.y < y_[y_.size() - 1]));
    inside &= ((gp.z > z_[0]) && (gp.z < z_[z_.size() - 1]));
    return inside;
  }

  /**
   * Get relative distance to next intersecting boundary in a given dimension.
   *
   * @param grid Array holding cell boundaries for given dimension
   * @param index The index of the next boundary for given dimension
   * @param direction Component of the direction vector in this dimension
   * @return Distance to next intersection given as multiple of direction. -1
   * if particle is leaving the grid.
   */
  __DEV__ Float next_plane(const Array<Float> &grid,
                           Index index,
                           Float position,
                           Float direction) {
    Float d = std::numeric_limits<Float>::max();
    if (!maths::small(direction)) {
      // We're on the left side of domain.
      if (index == 0) {
        d = -1.0;
        // We're on the right side of domain.
      } else if (index == grid.size() + 1) {
        d = -1.0;
        // We're in the middle
      } else {
        d = (grid[index - 1] - position) / direction;
      }
    }
    return d;
  }

 /**
  * Performs a step of a given maximum step length on the grid.
  * The step will move to the next grid boundary unless this
  * would exceed the provided step length. In that case the
  * the particle will move until the maximum length is reached.
  *
  * @param gp On entry: The position of the object to be moved.
  * On exit: The position after a step on the grid.
  * @param step_length: The maximum length of the step.
  * @return The actual step length or -1 if the particle is
  * leaving the grid.
  */
  template <typename Vector>
  __DEV__ Float step(GridPosition<Vector> &gp,
                     Float step_length) {

    Vector &direction = gp.direction;
    Vector &position = gp.position;
    Float d = std::numeric_limits<Float>::max();
    ushort di = 0;

    Float dx = next_plane(x_, gp.i, position.x, direction.x);
    Float dy = next_plane(y_, gp.j, position.y, direction.y);
    Float dz = next_plane(z_, gp.k, position.z, direction.z);

    if (dx < 0.0 || dy < 0.0 || dz < 0.0) {
      return -1.0;
    }

    if (dx < d) {
      di = 0;
      d = dx;
    }
    if (dy < d) {
      di = 1;
      d = dy;
    }
    if (dz < d) {
      di = 2;
      d = dz;
    }

    Float dl = direction.length();
    Float l = d * dl;

    // Compute new position.
    if (l > step_length) {
      l = step_length;
      d = step_length / dl;
      position.x += d * direction.x;
      position.y += d * direction.y;
      position.z += d * direction.z;
    } else {
      position.x += d * direction.x;
      position.y += d * direction.y;
      position.z += d * direction.z;
      if (di == 0) {
        if (direction.x < 0.0)
          --gp.i;
        else
          ++gp.i;
      } else if (di == 1) {
        if (direction.y < 0.0)
          --gp.j;
        else
          ++gp.j;
      } else if (di == 2) {
        if (direction.z < 0.0)
          --gp.k;
        else
          ++gp.k;
      }
    }

    if (direction.x < 0.0) {
      if (maths::small(position.x - x_[gp.i - 1]) && (gp.i > 0)) --gp.i;
    } else if (direction.x > 0.0) {
      if (maths::small(position.x - x_[gp.i - 1]) && (gp.i <= x_.size()))
        ++gp.i;
    }

    if (direction.y < 0.0) {
      if (maths::small(position.y - y_[gp.j - 1]) && (gp.j > 0)) --gp.j;
    } else if (direction.y > 0.0) {
      if (maths::small(position.y - y_[gp.j - 1]) && (gp.j <= y_.size()))
        ++gp.j;
    }

    if (direction.z < 0.0) {
      if (maths::small(position.z - z_[gp.k - 1]) && (gp.k > 0)) --gp.k;
    } else if (direction.z > 0.0) {
      if (maths::small(position.z - z_[gp.k - 1]) && (gp.k <= z_.size()))
        ++gp.k;
    }
    return l;
  }

  template <typename Vector>
 /**
  * @param gp The grid position of the particle
  * @return Integer representation of the boundary through
  * which the particle left the grid.
  */
  __DEV__ size_t get_boundary_index(GridPosition<Vector, Index> gp) const {
    if (gp.x() <= x_[0]) {
      return 0;
    }
    if (gp.x() >= x_[x_.size()]) {
      return 1;
    }
    if (gp.y() <= y_[0]) {
      return 2;
    }
    if (gp.y() >= y_[y_.size()]) {
      return 3;
    }
    if (gp.z() <= z_[0]) {
      return 4;
    }
    if (gp.z() >= z_[z_.size()]) {
      return 5;
    }
    return 999;
  }

#ifdef CUDA
void device() {
    x_.device();
    y_.device();
    z_.device();
}
void host() {
    x_.host();
    y_.host();
    z_.host();
}

#endif

 private:
  Array<Float> x_;
  Array<Float> y_;
  Array<Float> z_;
};

}  // namespace llrte
#endif
