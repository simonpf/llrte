#include <llrte/surfaces.h>
#include <llrte/types/vector.h>
#include <llrte/rotations.h>
#include <llrte/maths/geometry.h>
#include <llrte/types/matrix.h>
#include <llrte/constants.h>

template <typename F>
bool test_geometry() {
  using Float = F;
  using Vector = llrte::Vector<3, Float>;
  using Matrix = llrte::Matrix<3, 3, Float>;

  using llrte::maths::geometry::unit_vector;
  using llrte::maths::geometry::orthonormal_basis;

  Vector x = unit_vector<Vector, 0>();
  Vector y = unit_vector<Vector, 1>();
  Vector z = unit_vector<Vector, 2>();

  std::cout << orthonormal_basis<Matrix>(x) << std::endl;
  std::cout << orthonormal_basis<Matrix>(y) << std::endl;
  std::cout << orthonormal_basis<Matrix>(z) << std::endl;

}

int main(int args, const char **argv) {
    test_geometry<float>();
}
