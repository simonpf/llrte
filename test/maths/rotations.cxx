#include <llrte/surfaces.h>
#include <llrte/types/vector.h>
#include <llrte/rotations.h>
#include <llrte/random.h>
#include <llrte/maths/geometry.h>
#include <llrte/types/matrix.h>
#include <llrte/constants.h>

//bool test_random_rotation

template <typename F>
bool test_rotations() {
  using Float = F;
  using Vector = llrte::Vector<3, Float>;

  using llrte::rotations::rotate_y;
  using llrte::rotations::rotate_z;
  using llrte::rotations::rotate;
  using llrte::maths::geometry::unit_vector;
  using llrte::maths::geometry::is_close;

  Vector x = unit_vector<Vector, 0>();
  Vector y = unit_vector<Vector, 1>();
  Vector z = unit_vector<Vector, 2>();

  // Rotate x-vector
  Float phi = llrte::Constants<Float>::pi / 2.0;
  Float theta = 0;

  if (!is_close(y, rotate(x, phi, theta))) {
    return false;
  }

  phi = llrte::Constants<Float>::pi / 2.0;
  theta = -llrte::Constants<Float>::pi / 2.0;
  if (!is_close(z, rotate(x, phi, theta))) {
    return false;
  }

  // Rotate y-vector
  phi = 2.0 * llrte::Constants<Float>::pi;
  theta = -llrte::Constants<Float>::pi / 2.0;
  if (!is_close(z, rotate(y, phi, theta))) {
    return false;
  }

  phi = -llrte::Constants<Float>::pi / 2.0;
  theta = llrte::Constants<Float>::pi / 2.0;
  if (!is_close(x, rotate(y, phi, theta))) {
    return false;
  }

  // Rotate z-vector
  phi = llrte::Constants<Float>::pi / 2.0;
  theta = llrte::Constants<Float>::pi / 2.0;
  if (!is_close(y, rotate(y, phi, theta))) {
    return false;
  }

  phi = 0.0;
  theta = llrte::Constants<Float>::pi / 2.0;
  if (!is_close(x, rotate(y, phi, theta))) {
    return false;
  }
  return true;
}

template <typename F, typename Generator>
bool test_rotation_matrices(Generator &generator) {
  using Float = F;
  using Vector = llrte::Vector<3, Float>;
  using Matrix = llrte::Matrix<3, 3, Float>;

  using llrte::rotations::make_rotation_x;
  using llrte::rotations::make_rotation_y;
  using llrte::rotations::make_rotation_z;
  using llrte::rotations::make_rotation;
  using llrte::maths::geometry::unit_vector;
  using llrte::maths::geometry::is_close;

  Vector x = unit_vector<Vector, 0>();
  Vector y = unit_vector<Vector, 1>();
  Vector z = unit_vector<Vector, 2>();

  auto theta = generator.sample_angle_uniform();

  //auto rx = make_rotation_x<llrte::Matrix>(theta);
  auto rr = make_rotation<Matrix>(x, theta);

  theta = generator.sample_angle_uniform();

  auto ry = make_rotation_y<llrte::Matrix>(theta);
  rr = make_rotation<Matrix>(y, theta);

  std::cout << ry << std::endl;
  std::cout << rr << std::endl;

  theta = generator.sample_angle_uniform();
  //auto rz = make_rotation_z<llrte::Matrix>(theta);
  rr = make_rotation<Matrix>(z, theta);

  std::cout << theta << std::endl;



  return true;
}

int main(int /*args*/, const char **/*argv*/) {
  std::cout << "rotations: " << test_rotations<float>() << std::endl;

  llrte::Generator<float> g{};
  std::cout << "rotations: " << test_rotation_matrices<float>(g) << std::endl;
}
