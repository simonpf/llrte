#ifndef _LLRTE_SENSOR_H_
#define _LLRTE_SENSOR_H_

#include <algorithm>
#include <string>

#include "llrte/eigen.h"
#include "llrte/io/netcdf.h"
#include "llrte/rotations.h"

namespace llrte {

using eigen::Vector;
using eigen::Tensor;

template <typename Photon>
class SensorArray {
 public:
  using V3 = typename Photon::Vector;
  using Float = typename V3::Float;
  SensorArray(V3 position, V3 x, V3 y, Vector<Float> dx, Vector<Float> dy,
              Vector<Float> zenith_angles, Vector<Float> azimuth_angles)
      : position_(position),
        x_(x),
        y_(y),
        n_(cross(x, y)),
        dx_(dx),
        dy_(dy),
        zenith_angles_(zenith_angles),
        azimuth_angles_(azimuth_angles),
        data_{{dx.size(), dy.size(),
               std::max<size_t>(zenith_angles.size() - 1, 1),
               std::max<size_t>(azimuth_angles.size() - 1, 1)}} {
    data_.fill(0.0);
  }

  template <typename Solver, typename Source>
  void sample(Solver &solver, Source &source, size_t n) {
    auto generator = solver.generator();
    for (size_t i_x = 0; i_x < dx_.size(); ++i_x) {
      std::cout << i_x << " // " << dx_.size() << std::endl;
      for (size_t i_y = 0; i_y < dy_.size(); ++i_y) {
        size_t n_a = std::max<size_t>(azimuth_angles_.size() - 1, 1);
        for (size_t i_a = 0; i_a < n_a; ++i_a) {
          size_t n_z = std::max<size_t>(zenith_angles_.size() - 1, 1);
          for (size_t i_z = 0; i_z < n_z; ++i_z) {
            for (size_t i_s = 0; i_s < n; ++i_s) {
              Float phi, theta;

              if (azimuth_angles_.size() > 1) {
                phi = generator.sample_uniform(azimuth_angles_[i_a],
                                               azimuth_angles_[i_a + 1]);
                if (zenith_angles_.size() > 1) {
                  theta = generator.sample_zenith_angle(
                      zenith_angles_[i_z], zenith_angles_[i_z + 1]);
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

              Photon photon = solver.template backward<Photon>(p, d, source);
              data_(i_x, i_y, i_z, i_a) += photon.get_energy();
            }
          }
        }
      }
    }
  }

  void dump(std::string filename) {
    io::NetCDFFile netcdf_file(filename, true);
    netcdf_file.add_dimension("x", dx_.size());
    netcdf_file.add_dimension("y", dy_.size());
    netcdf_file.add_dimension("zenith_angle",
                              std::max<size_t>(zenith_angles_.size() - 1, 1));
    netcdf_file.add_dimension("azimuth_angle",
                              std::max<size_t>(azimuth_angles_.size() - 1, 1));
    netcdf_file.store_variable(data_, "data",
                               {"x", "y", "zenith_angle", "azimuth_angle"});
  }

 private:
  V3 position_, x_, y_, n_;
  Vector<Float> dx_, dy_, zenith_angles_, azimuth_angles_;
  Tensor<Float, 4> data_;
};

template <typename Photon>
class SphericalSensor {
 public:
  using V3 = typename Photon::Vector;
  using Float = typename V3::Float;
  SphericalSensor(V3 position, V3 direction, V3 axis_1, V3 axis_2,
                  Vector<Float> d1, Vector<Float> d2)
      : position_(position),
        direction_(direction),
        axis_1_(axis_1),
        axis_2_(axis_2),
        d1_(d1),
        d2_(d2),
        data_{Tensor<Float, 2>{d1.size(), d2.size()}.setZero()} {}

  template <typename Solver, typename Source>
  void sample(Solver &solver, Source &source, size_t n) {
    auto generator = solver.generator();
    for (size_t i_1 = 0; i_1 < d1_.size(); ++i_1) {
      for (size_t i_2 = 0; i_2 < d2_.size(); ++i_2) {
          std::cout << i_1 << " / " << i_2 << std::endl;
        auto d = rotations::rotate(direction_, axis_1_, d1_[i_1]);
        d = rotations::rotate(d, axis_2_, d2_[i_2]);
        for (size_t i_s = 0; i_s < n; ++i_s) {
          Photon photon =
              solver.template backward<Photon>(position_, d, source);
              auto e = photon.get_energy();
              data_.coeffRef(i_1, i_2) += e;
        }
      }
    }
  }

  void reset() {
      data_ = 0.0 * data_;
  }

  void dump(std::string filename) {
    io::NetCDFFile netcdf_file(filename, true);
    netcdf_file.add_dimension("d1", d1_.size());
    netcdf_file.add_dimension("d2", d2_.size());
    netcdf_file.store_variable(data_, "data", {"d1", "d2"});
  }

  const Tensor<Float, 2> & get_data() const {return data_;}
  Tensor<Float, 2> & get_data() {return data_;}

 private:
  V3 position_, direction_, axis_1_, axis_2_;
  Vector<Float> d1_, d2_;
  Tensor<Float, 2> data_;
};

}  // namespace llrte
#endif
