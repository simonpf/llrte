#ifndef _LLRTE_IO_NETCDF_H_
#define _LLRTE_IO_NETCDF_H_

extern "C" {
#include <netcdf.h>
}


#include <llrte/utils/array.h>
#include <llrte/data.h>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <sys/stat.h>

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
             bool clobber = false,
             Mode mode = Mode::read) {
    // Check if file exists.
    struct stat buffer;
    bool exists = stat(filename.c_str(), &buffer) == 0;

    if (exists && !clobber) {
        // Open existing file.
      __NC_ERROR__("Error opening NetCDF File: ",
                   nc_open(filename.c_str(),
                           static_cast<int>(mode),
                           &ncid_));

    } else {
        // Create new file.
      __NC_ERROR__("Error creating NetCDF File: ",
                   nc_create(filename.c_str(), clobber, &ncid_));
    }
  }

  NetCDFFile(const NetCDFFile &other) {
      try {
          nc_close(ncid_);
      } catch (...) {}
      ncid_ = other.ncid_;
      dimensions_ = other.dimensions_;
  }

  NetCDFFile(NetCDFFile &&other) {
      try {
          nc_close(ncid_);
      } catch (...) {}
      ncid_ = other.ncid_;
      other.ncid_=-1;
      dimensions_ = other.dimensions_;
  }

  NetCDFFile& operator=(const NetCDFFile &other) {
      try {
          nc_close(ncid_);
      } catch (...) {
      }
      ncid_ = other.ncid_;
      dimensions_ = other.dimensions_;
      return *this;
  }

  NetCDFFile& operator=(NetCDFFile &&other) {
      try {
          nc_close(ncid_);
      } catch (...) {
      }
      ncid_ = other.ncid_;
      other.ncid_=-1;
      dimensions_ = other.dimensions_;
      return *this;
  }

  /**
   * \brief Close the file.
   */
  ~NetCDFFile() {
      try {
          nc_close(ncid_);
      } catch (...) {
      }
  }
  /**
   * \brief Add a dimension to a NetCDF file.
   *
   * @param name The name of the dimension to add.
   * @param size The size of the dimension.
   */
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

  /**
   * \brief Store a variable in a NetCDF file.
   *
   * @param name The tensor variable to store.
   * @param size The dimensions of the file to associate
   *    with the axes of the tensor.
   */
  template <typename T, size_t rank>
  void store_variable(const Tensor<T, rank> &t,
                      const std::string &name,
                      const std::array<std::string, rank> &dimensions) {

    // Resolve dimensions.
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
    __NC_ERROR__("Error adding variable to NetCDF File: ",
                 nc_def_var(ncid_, name.c_str(), NetCDFType<T>::id, rank,
                            dim_ids.data(), &varid));

    __NC_ERROR__("Error adding variable to NetCDF File:",
                 nc_enddef(ncid_));

    // Write variable to file.
    std::array<size_t, rank> size = t.shape();
    std::array<size_t, rank> start{0};
    __NC_ERROR__("Error writing variable to NetCDF File:",
                 nc_put_vara(ncid_, varid, start.data(), size.data(), t.get_data_pointer()));
  }

  template <typename T, size_t rank>
  Tensor<T, rank> load_variable(const std::string &name) {
      int varid;
      __NC_ERROR__("Error finding variable in file:",
                   nc_inq_varid(ncid_, name.c_str(), &varid));

      int ndims;
      __NC_ERROR__("Error getting number of dimensions:",
                   nc_inq_varndims(ncid_, varid, &ndims));

      if (ndims != rank) {
          throw std::runtime_error("Number of dimensions does not match rank "
                                   " of tensor to return.");
      }

      std::array<int, rank> dimids;
      __NC_ERROR__("Error retrieving dimension ids of variable:",
                   nc_inq_vardimid(ncid_, varid, dimids.data()));

      std::array<size_t, rank> shape;
      for (size_t i = 0; i < static_cast<size_t>(ndims); ++i) {
          __NC_ERROR__("Error retrieving dimensions of variable.",
                       nc_inq_dim(ncid_, dimids[i], nullptr, shape.data() + i)) ;
      }

      Tensor<T, rank> data_out{shape};

      /* Read the data. */
      std::array<size_t, rank> start;
      for (size_t i = 0; i < rank; ++i) start[i] = 0.0;

      __NC_ERROR__("Error retrieving variable data.",
                   retval = nc_get_vara(ncid_,
                                        varid,
                                        start.data(),
                                        shape.data(),
                                        data_out.get_data_pointer()));
      return data_out;
  }

  void close() {
      try {
          nc_close(ncid_);
      } catch (...) {}
  }


 private:
  int ncid_;
  DimensionMap dimensions_ = {};
};

}  // namespace llrte::io

#endif
