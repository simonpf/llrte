#include <math.h>
#include <iostream>

#include <llrte/constants.h>
#include <llrte/maths/general.h>
#include <llrte/maths/geometry.h>
#include <llrte/random.h>
#include <llrte/rotations.h>
#include <llrte/surfaces.h>
#include <llrte/types/matrix.h>
#include <llrte/types/vector.h>

template <typename Generator>
bool test_geometry(Generator &generator) {
  using Float = typename Generator::Float;
  using Vector = llrte::Vector<3, Float>;
  using Matrix = llrte::Matrix<3, 3, Float>;

  using llrte::maths::geometry::is_close;
  using llrte::maths::geometry::orthonormal_basis;
  using llrte::maths::geometry::unit_vector;
  using llrte::random_direction;
  using llrte::column;

  auto r = random_direction<Vector>(generator);
  auto b = orthonormal_basis<Matrix>(r);
  auto b3 = llrte::cross(column<Vector>(b, 0), column<Vector>(b, 1));

  return is_close(b3 * (1.0 / b3.length()),
                  llrte::column<Vector>(b, 2));
}

template <typename Generator>
bool test_perpendicular(Generator &generator) {
    using Float = typename Generator::Float;
    using Vector = llrte::Vector<3, Float>;

    using llrte::maths::geometry::is_close;
    using llrte::maths::geometry::perpendicular;
    using llrte::maths::geometry::unit_vector;
    using llrte::random_direction;
    using llrte::column;

    auto r = random_direction<Vector>(generator);
    auto pr = perpendicular(r);
    return llrte::maths::small(dot(r, pr));
}

int main(int /*args*/, const char **/*argv*/) {
  auto generator = llrte::Generator<float>();
  for (size_t i = 0; i < 100; ++i) {
    if (!test_geometry(generator)) {
      return 1;
    }
    if (!test_perpendicular(generator)) {
        return 1;
    }
  }
  return 0;
}
