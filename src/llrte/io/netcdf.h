#ifndef _LLRTE_IO_NETCDF_H_
#define _LLRTE_IO_NETCDF_H_

extern "C" {
#include <netcdf.h>
}

#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <iostream>
#include <sys/stat.h>

namespace llrte {
template <typename T, size_t rank>
class Tensor;
}

#define __NC_ERROR__(a, b)                                    \
  ({                                                          \
    int retval = b;                                           \
    if (retval) {                                             \
      std::string msg = std::string(a) + nc_strerror(retval); \
      throw std::runtime_error(msg);                          \
    }                                                         \
  })

namespace llrte::io {

/* Alias for unlimited dimension in NetCDFFile. */
size_t unlimited = NC_UNLIMITED;
/* Map used to store dimensions, their ids and sizes. */
using DimensionMap = std::map<std::string, std::pair<int, size_t>>;

//******************************************************************************
// NetCDF traits
//******************************************************************************

template <typename T>
struct NetCDFType;

template <>
struct NetCDFType<double> {
  constexpr static int id = NC_DOUBLE;
  static int put(int ncid, int varid, const double *ptr) {
    return nc_put_var_double(ncid, varid, ptr);
  }
};

template <>
struct NetCDFType<float> {
  constexpr static int id = NC_FLOAT;
  static int put(int ncid, int varid, const float *ptr) {
    return nc_put_var_float(ncid, varid, ptr);
  }
};

enum class Mode {read = NC_SHARE, write = NC_WRITE};

//*****************************************************************************
// NetCDFFile
//*****************************************************************************
/**
 * \brief NetCDF 4 file IO
 *
 * This class implements an interface for the writing and reading of NetCDF 4
 * files.
 *
 */
class NetCDFFile {

 public:

  /**
   * \brief Open a NetCDF 4 file.
   *
   * Opens and existing file (if clobber is set to false) or creates a
   * new file with the given name.
   *
   * @param filename The name of the file to open.
   * @param mode The mode (read or write) in which to open the file.
   * @param clobber If true an existing file will be overwritten. Otherwise
   *    the file is opened in read mode.
   */
  NetCDFFile(const std::string &filename,
             Mode mode = Moder::read,
             bool clobber = false) {
    // Check if file exists.
    struct stat buffer;
    bool exists = stat(name.c_str(), &buffer) == 0;

    if (exists && !clobber) {
        // Open existing file.
      __NC_ERROR__("Error opening NetCDF File: ",
                   nc_open(filename.c_str(), mode, &ncid_));

    } else {
        // Create new file.
      __NC_ERROR__("Error creating NetCDF File: ",
                   nc_create(filename.c_str(), clobber, &ncid_));
    }
  }

  void add_dimension(std::string name,
                     size_t size = unlimited) {
    auto it = dimensions_.find(name);
    if (it != dimensions_.end()) {
      int id2;
      size_t size2;
      std::tie(id2, size2) = std::get<1>((*it));
      if (size != size2) {
        std::stringstream ss{};
        ss << "Dimension " << name << "added to file with different size "
           << size2 << " (instead of " << size << ").";
        throw std::runtime_error(ss.str());
      }
    } else {

      int dim_id;
      int retval = nc_def_dim(ncid_, name.c_str(), size, &dim_id);
      if (retval) {
        std::string msg =
            std::string("Error adding dimension: ") + nc_strerror(retval);
        throw std::runtime_error(msg);
      }

      dimensions_.emplace(name, std::make_pair(dim_id, size));
    }
  }

  template <typename T, size_t rank>
  void store_variable(const Tensor<T, rank> &t, const std::string &name,
                      const std::array<std::string, rank> &dimensions) {
    std::array<int, rank> dim_ids;
    std::array<size_t, rank> sizes;
    for (size_t i = 0; i < rank; ++i) {
      const std::string &dim_name = dimensions[i];
      auto it = dimensions_.find(dim_name);
      if (it == dimensions_.end()) {
        std::stringstream ss{};
        ss << "Dimension " << dim_name << " not found.";
        throw std::runtime_error(ss.str());
      } else {
        std::tie(dim_ids[i], sizes[i]) = std::get<1>((*it));
      }
    }

    if (sizes != t.shape()) {
      std::stringstream ss{};
      ss << "Dimensions " << t.shape()
         << "of tensor do not match "
            "diemsions "
         << sizes << "in NetCDF file.";
      throw std::runtime_error(ss.str());
    }

    // Enter define mode.
    int retval = nc_redef(ncid_);
    if ((retval) && (retval != NC_EINDEFINE)) {
      std::string msg = std::string("Error adding variable to NetCDF File: ") +
                        nc_strerror(retval);
      throw std::runtime_error(msg);
    }

    // Add the variable.
    int varid;
    retval = nc_def_var(ncid_, name.c_str(), NetCDFType<T>::id, rank,
                        dim_ids.data(), &varid);
    if (retval) {
      std::string msg = std::string("Error adding variable to NetCDF File: ") +
                        nc_strerror(retval);
      throw std::runtime_error(msg);
    }
    retval = nc_enddef(ncid_);
    if (retval) {
      std::string msg = std::string("Error adding variable to NetCDF File: ") +
                        nc_strerror(retval);
      throw std::runtime_error(msg);
    }

    retval = NetCDFType<T>::put(ncid_, varid, t.get_data_pointer());
    if (retval) {
      std::string msg = std::string("Error writing variable to NetCDF File: ") +
                        nc_strerror(retval);
      throw std::runtime_error(msg);
    }
  }

  ~NetCDFFile() {
    try {
      nc_close(ncid_);
    } catch (...) {
    }
  }

 private:
  int ncid_;
  DimensionMap dimensions_ = {};
};

}  // namespace llrte::io

#endif
