#ifndef _LLRTE_SENSOR_H_
#define _LLRTE_SENSOR_H_

#include <string>
#include <algorithm>
#include "llrte/data.h"
#include "llrte/io/netcdf.h"
#include "llrte/rotations.h"

namespace llrte {

template <typename Vector>
class SensorArray {
public:
  using Float = typename Vector::Float;
  SensorArray(Vector position,
              Vector x,
              Vector y,
              Array<Float> dx,
              Array<Float> dy,
              Array<Float> zenith_angles,
              Array<Float> azimuth_angles
      )
      : position_(position),
        x_(x),
        y_(y),
        n_(cross(x, y)),
        dx_(dx),
        dy_(dy),
        zenith_angles_(zenith_angles),
        azimuth_angles_(azimuth_angles),
        data_{{dx.size(),
               dy.size(),
              std::max<size_t>(zenith_angles.size() - 1, 1),
              std::max<size_t>(azimuth_angles.size() - 1, 1)}} {
      data_.fill(0.0);
  }

  template <typename Generator, typename Solver>
  void sample(Generator &generator,  Solver &solver, size_t n) {
    for (size_t i_x = 0; i_x < dx_.size(); ++i_x) {
        std::cout << i_x << " // " << dx_.size() << std::endl;
      for (size_t i_y = 0; i_y < dy_.size(); ++i_y) {
          size_t n_a = std::max<size_t>(azimuth_angles_.size() - 1, 1);
          for (size_t i_a = 0; i_a < n_a; ++i_a) {
              size_t n_z = std::max<size_t>(zenith_angles_.size() - 1, 1);
              for (size_t i_z = 0; i_z < n_z; ++i_z) {
                  for (size_t i_s = 0; i_s < n; ++ i_s) {
                      Float phi, theta;

                      if (azimuth_angles_.size() > 1) {
                          phi = generator.sample_uniform(azimuth_angles_[i_a],
                                                         azimuth_angles_[i_a + 1]);
                          if (zenith_angles_.size() > 1) {
                              theta = generator.sample_zenith_angle(zenith_angles_[i_z],
                                                                    zenith_angles_[i_z + 1]);
                          } else {
                              theta = zenith_angles_[0];
                          }
                      } else {
                          phi = azimuth_angles_[0];
                          if (zenith_angles_.size() > 1) {
                              theta = generator.sample_uniform(zenith_angles_[i_z],
                                                               zenith_angles_[i_z + 1]);
                          } else {
                              theta = zenith_angles_[0];
                          }
                      }


                      auto y = rotations::rotate(y_, n_, phi);
                      auto d = rotations::rotate(n_, y, theta);
                      auto p = position_ + dx_[i_x] * x_ + dy_[i_y] * y_;

                      Photon<Vector> photon{p, d};
                      photon = solver.propagate_photon(generator, photon);
                      data_(i_x, i_y, i_z, i_a) += photon.get_energy();
                  }
                  data_(i_x, i_y, i_z, i_a) /= static_cast<Float>(n);
              }
          }
        }
      }
    }

  void dump(std::string filename) {
      io::NetCDFFile netcdf_file(filename, true);
      netcdf_file.add_dimension("x", dx_.size());
      netcdf_file.add_dimension("y", dy_.size());
      netcdf_file.add_dimension("zenith_angle", std::max<size_t>(zenith_angles_.size() - 1, 1));
      netcdf_file.add_dimension("azimuth_angle", std::max<size_t>(azimuth_angles_.size() - 1, 1));
      netcdf_file.store_variable(data_,
                                 "data",
                                 {"x", "y", "zenith_angle", "azimuth_angle"});
  }

 private:

  Vector position_, x_, y_, n_;
  Array<Float> dx_, dy_, zenith_angles_, azimuth_angles_;
  Tensor<Float, 4> data_;
};

}  // namespace llrte
#endif
