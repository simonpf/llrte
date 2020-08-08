#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "llrte/grids/regular.h"
#include "llrte/data.h"
#include "llrte/maths.h"
#include "llrte/types/vector.h"

TEMPLATE_TEST_CASE("Single cell, x-direction",
                   "[RegularGrid][place_on_grid]", float, double) {
  using llrte::maths::small;
  using Vector = llrte::Vector3<TestType>;
  auto x = llrte::Array<TestType>::range(2);
  auto y = llrte::Array<TestType>::range(2);
  auto z = llrte::Array<TestType>::range(2);
  auto rg = llrte::RegularGrid<TestType>{x, y, z};

  //
  // Moving on lower edge.
  //

  // Left of cell, moving in.
  auto gp = rg.place_on_grid(Vector{-1.0, 0.0, 0.0}, Vector{1.0, 0.0, 0.0});
  REQUIRE(((gp.i == 1) && (gp.j == 0) && (gp.k == 0)));
  // Left boundary.
  auto l = rg.step(gp, 2.0);
  REQUIRE(((gp.i == 2) && (gp.j == 0) && (gp.k == 0)));
  REQUIRE(small(l - 1.0));
  // Right boundary.
  l = rg.step(gp, 2.0);
  REQUIRE(((gp.i == 3) && (gp.j == 0) && (gp.k == 0)));
  REQUIRE(small(l - 1.0));
  // Leaving.
  l = rg.step(gp, 2.0);
  REQUIRE(l < 0.0);
  REQUIRE(((gp.i == 3) && (gp.j == 0) && (gp.k == 0)));
  REQUIRE(small(gp.position.x - 1.0));

  // Right of cell, moving in
  gp = rg.place_on_grid(Vector{2.0, 0.0, 0.0}, Vector{-1.0, 0.0, 0.0});
  REQUIRE(((gp.i == 2) && (gp.j == 0) && (gp.k == 0)));
  // Right boundary.
  l = rg.step(gp, 2.0);
  REQUIRE(((gp.i == 1) && (gp.j == 0) && (gp.k == 0)));
  REQUIRE(small(l - 1.0));
  // Left boundary
  l = rg.step(gp, 2.0);
  REQUIRE(((gp.i == 0) && (gp.j == 0) && (gp.k == 0)));
  REQUIRE(small(l - 1.0));
  // Leaving
  l = rg.step(gp, 2.0);
  REQUIRE(((gp.i == 0) && (gp.j == 0) && (gp.k == 0)));
  REQUIRE(l < 0.0);
  REQUIRE(small(gp.position.x));

  //
  // Moving on upper edge.
  //

  // Left of cell, moving in.
  gp = rg.place_on_grid(Vector{-1.0, 1.0, 1.0}, Vector{1.0, 0.0, 0.0});
  REQUIRE(((gp.i == 1) && (gp.j == 1) && (gp.k == 1)));
  // Left boundary.
  l = rg.step(gp, 2.0);
  REQUIRE(small(l - 1.0));
  REQUIRE(((gp.i == 2) && (gp.j == 1) && (gp.k == 1)));
  // Right boundary.
  l = rg.step(gp, 2.0);
  REQUIRE(small(l - 1.0));
  REQUIRE(((gp.i == 3) && (gp.j == 1) && (gp.k == 1)));
  // Leaving.
  l = rg.step(gp, 2.0);
  REQUIRE(((gp.i == 3) && (gp.j == 1) && (gp.k == 1)));
  REQUIRE(l < 0.0);
  REQUIRE(small(gp.position.x - 1.0));

  // Right of cell, moving in
  gp = rg.place_on_grid(Vector{2.0, 1.0, 1.0}, Vector{-1.0, 0.0, 0.0});
  REQUIRE(((gp.i == 2) && (gp.j == 1) && (gp.k == 1)));
  // Right boundary.
  l = rg.step(gp, 2.0);
  REQUIRE(((gp.i == 1) && (gp.j == 1) && (gp.k == 1)));
  REQUIRE(small(l - 1.0));
  // Left boundary.
  l = rg.step(gp, 2.0);
  REQUIRE(small(l - 1.0));
  REQUIRE(((gp.i == 0) && (gp.j == 1) && (gp.k == 1)));
  // Leaving.
  l = rg.step(gp, 2.0);
  REQUIRE(((gp.i == 0) && (gp.j == 1) && (gp.k == 1)));
  REQUIRE(small(gp.position.x));
  REQUIRE(l < 0.0);
}

TEMPLATE_TEST_CASE("Single cell, y-direction",
                   "[RegularGrid][place_on_grid]", float, double) {
  using llrte::maths::small;
  using Vector = llrte::Vector3<TestType>;
  auto x = llrte::Array<TestType>::range(2);
  auto y = llrte::Array<TestType>::range(2);
  auto z = llrte::Array<TestType>::range(2);
  auto rg = llrte::RegularGrid<TestType>{x, y, z};

  //
  // Moving on lower edge.
  //

  // Left of cell, moving in.
  auto gp = rg.place_on_grid(Vector{0.0, -1.0, 0.0}, Vector{0.0, 1.0, 0.0});
  REQUIRE(((gp.i == 0) && (gp.j == 1) && (gp.k == 0)));
  // Left boundary.
  auto l = rg.step(gp, 2.0);
  REQUIRE(small(l - 1.0));
  REQUIRE(((gp.i == 0) && (gp.j == 2) && (gp.k == 0)));
  // Right boundary.
  l = rg.step(gp, 2.0);
  REQUIRE(small(l - 1.0));
  REQUIRE(((gp.i == 0) && (gp.j == 3) && (gp.k == 0)));
  // Leaving.
  l = rg.step(gp, 2.0);
  REQUIRE(l < 0.0);
  REQUIRE(((gp.i == 0) && (gp.j == 3) && (gp.k == 0)));
  REQUIRE(small(gp.position.y - 1.0));

  // Right of cell, moving in
  gp = rg.place_on_grid(Vector{0.0, 2.0, 0.0}, Vector{0.0, -1.0, 0.0});
  REQUIRE(((gp.i == 0) && (gp.j == 2) && (gp.k == 0)));
  // Right boundary.
  l = rg.step(gp, 2.0);
  REQUIRE(small(l - 1.0));
  REQUIRE(((gp.i == 0) && (gp.j == 1) && (gp.k == 0)));
  // Left boundary.
  l = rg.step(gp, 2.0);
  REQUIRE(small(l - 1.0));
  REQUIRE(((gp.i == 0) && (gp.j == 0) && (gp.k == 0)));
  // Leaving.
  l = rg.step(gp, 2.0);
  REQUIRE(l < 0.0);
  REQUIRE(((gp.i == 0) && (gp.j == 0) && (gp.k == 0)));
  REQUIRE(small(gp.position.y));

  //
  // Moving on upper edge.
  //

  // Left of cell, moving in.
  gp = rg.place_on_grid(Vector{1.0, -1.0, 1.0}, Vector{0.0, 1.0, 0.0});
  REQUIRE(((gp.i == 1) && (gp.j == 1) && (gp.k == 1)));
  // Left boundary.
  l = rg.step(gp, 2.0);
  REQUIRE(small(l - 1.0));
  REQUIRE(((gp.i == 1) && (gp.j == 2) && (gp.k == 1)));
  // Right boundary.
  l = rg.step(gp, 2.0);
  REQUIRE(small(l - 1.0));
  REQUIRE(((gp.i == 1) && (gp.j == 3) && (gp.k == 1)));
  // Leaving.
  l = rg.step(gp, 2.0);
  REQUIRE(l < 0.0);
  REQUIRE(((gp.i == 1) && (gp.j == 3) && (gp.k == 1)));
  REQUIRE(small(gp.position.y - 1.0));

  // Right of cell, moving in
  gp = rg.place_on_grid(Vector{1.0, 2.0, 1.0}, Vector{0.0, -1.0, 0.0});
  REQUIRE(((gp.i == 1) && (gp.j == 2) && (gp.k == 1)));
  // Right boundary.
  l = rg.step(gp, 2.0);
  REQUIRE(small(l - 1.0));
  REQUIRE(((gp.i == 1) && (gp.j == 1) && (gp.k == 1)));
  // Left boundary.
  l = rg.step(gp, 2.0);
  REQUIRE(small(l - 1.0));
  REQUIRE(((gp.i == 1) && (gp.j == 0) && (gp.k == 1)));
  // Leaving.
  l = rg.step(gp, 2.0);
  REQUIRE(l < 0.0);
  REQUIRE(((gp.i == 1) && (gp.j == 0) && (gp.k == 1)));
  REQUIRE(small(gp.position.y));
}

TEMPLATE_TEST_CASE("Single cell, z-direction",
                   "[RegularGrid][place_on_grid]", float, double) {
  using llrte::maths::small;
  using Vector = llrte::Vector3<TestType>;
  auto x = llrte::Array<TestType>::range(2);
  auto y = llrte::Array<TestType>::range(2);
  auto z = llrte::Array<TestType>::range(2);
  auto rg = llrte::RegularGrid<TestType>{x, y, z};

  //
  // Moving on lower edge.
  //

  // Left of cell, moving in.
  auto gp = rg.place_on_grid(Vector{0.0, 0.0, -1.0}, Vector{0.0, 0.0, 1.0});
  REQUIRE(((gp.i == 0) && (gp.j == 0) && (gp.k == 1)));
  // Left boundary.
  auto l = rg.step(gp, 2.0);
  REQUIRE(small(l - 1.0));
  REQUIRE(((gp.i == 0) && (gp.j == 0) && (gp.k == 2)));
  // Right boundary.
  l = rg.step(gp, 2.0);
  REQUIRE(small(l - 1.0));
  REQUIRE(((gp.i == 0) && (gp.j == 0) && (gp.k == 3)));
  // Leaving.
  l = rg.step(gp, 2.0);
  REQUIRE(l < 0.0);
  REQUIRE(((gp.i == 0) && (gp.j == 0) && (gp.k == 3)));
  REQUIRE(small(gp.position.z - 1.0));

  // Right of cell, moving in
  gp = rg.place_on_grid(Vector{0.0, 0.0, 2.0}, Vector{0.0, 0.0, -1.0});
  REQUIRE(((gp.i == 0) && (gp.j == 0) && (gp.k == 2)));
  // Right boundary.
  l = rg.step(gp, 2.0);
  REQUIRE(small(l - 1.0));
  REQUIRE(((gp.i == 0) && (gp.j == 0) && (gp.k == 1)));
  // Left boundary.
  l = rg.step(gp, 2.0);
  REQUIRE(small(l - 1.0));
  REQUIRE(((gp.i == 0) && (gp.j == 0) && (gp.k == 0)));
  // Leaving.
  l = rg.step(gp, 2.0);
  REQUIRE(l < 0.0);
  REQUIRE(((gp.i == 0) && (gp.j == 0) && (gp.k == 0)));
  REQUIRE(small(gp.position.z));

  //
  // Moving on upper edge.
  //

  // Left of cell, moving in.
  gp = rg.place_on_grid(Vector{1.0, 1.0, -1.0}, Vector{0.0, 0.0, 1.0});
  REQUIRE(((gp.i == 1) && (gp.j == 1) && (gp.k == 1)));
  // Left boundary.
  l = rg.step(gp, 2.0);
  REQUIRE(small(l - 1.0));
  REQUIRE(((gp.i == 1) && (gp.j == 1) && (gp.k == 2)));
  // Right boundary.
  l = rg.step(gp, 2.0);
  REQUIRE(small(l - 1.0));
  REQUIRE(((gp.i == 1) && (gp.j == 1) && (gp.k == 3)));
  // Leaving.
  l = rg.step(gp, 2.0);
  REQUIRE(l < 0.0);
  REQUIRE(((gp.i == 1) && (gp.j == 1) && (gp.k == 3)));
  REQUIRE(small(gp.position.z - 1.0));

  // Right of cell, moving in
  gp = rg.place_on_grid(Vector{1.0, 1.0, 2.0}, Vector{0.0, 0.0, -1.0});
  REQUIRE(((gp.i == 1) && (gp.j == 1) && (gp.k == 2)));
  // Right boundary.
  l = rg.step(gp, 2.0);
  REQUIRE(small(l - 1.0));
  REQUIRE(((gp.i == 1) && (gp.j == 1) && (gp.k == 1)));
  // Left boundary.
  l = rg.step(gp, 2.0);
  REQUIRE(small(l - 1.0));
  REQUIRE(((gp.i == 1) && (gp.j == 1) && (gp.k == 0)));
  // Leaving.
  l = rg.step(gp, 2.0);
  REQUIRE(l < 0.0);
  REQUIRE(((gp.i == 1) && (gp.j == 1) && (gp.k == 0)));
  REQUIRE(small(gp.position.z));
}

TEMPLATE_TEST_CASE("Diagonal direction", "[RegularGrid][place_on_grid]",
                   float, double) {
  using llrte::maths::small;
  using Vector = llrte::Vector3<TestType>;
  auto x = llrte::Array<TestType>::range(3);
  auto y = llrte::Array<TestType>::range(3);
  auto z = llrte::Array<TestType>::range(3);
  auto rg = llrte::RegularGrid<TestType>{x, y, z};

  auto gp =
      rg.place_on_grid(Vector{-1.0, -1.0, -1.0}, Vector{1.0, 1.0, 1.0});
  REQUIRE(((gp.i == 1) && (gp.j == 1) && (gp.k == 1)));

  auto l = rg.step(gp, 2.0);
  REQUIRE(((gp.i == 2) && (gp.j == 2) && (gp.k == 2)));

  l = rg.step(gp, 2.0);
  REQUIRE(((gp.i == 3) && (gp.j == 3) && (gp.k == 3)));

  l = rg.step(gp, 2.0);
  l = rg.step(gp, 2.0);
  REQUIRE(l < -0.0);

  // Moving in opposite direction.
  gp = rg.place_on_grid(Vector{3.0, 3.0, 3.0}, Vector{-1.0, -1.0, -1.0});
  REQUIRE(((gp.i == 3) && (gp.j == 3) && (gp.k == 3)));

  l = rg.step(gp, 2.0);
  REQUIRE(((gp.i == 2) && (gp.j == 2) && (gp.k == 2)));

  l = rg.step(gp, 2.0);
  REQUIRE(((gp.i == 1) && (gp.j == 1) && (gp.k == 1)));

  l = rg.step(gp, 2.0);
  REQUIRE(((gp.i == 0) && (gp.j == 0) && (gp.k == 0)));

  l = rg.step(gp, 2.0);
  l = rg.step(gp, 2.0);
  REQUIRE(l < 0.0);
}
