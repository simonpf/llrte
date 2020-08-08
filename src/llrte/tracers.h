#ifndef _LLRTE_TRACERS_
#define _LLRTE_TRACERS_

#include "llrte/common.h"
#include "llrte/definitions.h"
#include "llrte/io/netcdf.h"

namespace llrte::tracers {

using llrte::eigen::Index;
using llrte::eigen::Tensor;
template<typename Scalar>
using Array = Tensor<Scalar, 1>;

////////////////////////////////////////////////////////////////////////////////
// No Trace
////////////////////////////////////////////////////////////////////////////////

/**
 * NoTrace
 *
 * Don't trace any results in Monte Carlo simulation. The NoTrace class
 * implements all function of the tracer interface as NOOPS. It's useful for
 * backward simulations as well as base class for custom tracers.
 */
struct NoTrace {
  /**
   * Photon creation. This function is called when a photon is created.
   * @param Reference to the photon
   */
  template <typename Photon>
  __DEV__ void created(const Photon & /*photon*/) {}
  /**
   * Trace absorption. This method is called for every step of the MC algorithm
   * and can be used to trace photon position as well as absorbed energy.
   * @param Reference to the photon
   */
  template <typename Photon, typename Float>
  __DEV__ void absorption(const Photon & /*photon*/, Float /*energy*/) {}

  /**
   * Trace out of energy. This method is called when a photon's energy falls
   * below it's minimum energy threshold.
   * @param Reference to the photon
   */
  template <typename Photon>
  __DEV__ void out_of_energy(const Photon & /*photon*/) {}

  /**
   * Trace scattering. This method is called when a scattering event occurs.
   * @param Reference to the photon
   */
  template <typename Photon>
  __DEV__ void scattering(const Photon & /*photon*/) {}

  /**
   * Left atmosphere. This method is called when a photon leaves the atmosphere.
   * @param Reference to the photon
   */
  template <typename Photon, typename Atmosphere>
  __DEV__ void left_atmosphere(const Photon & /*photon*/,
                               const Atmosphere & /*atmosphere*/) {}
};

////////////////////////////////////////////////////////////////////////////////
// Histogram
////////////////////////////////////////////////////////////////////////////////
/** Histogram
 * Histogram tracer that only counts of photons in different grid cells.
 * @tparam Grid The 3D grid type which described the atmosphere.
 */
template <typename Grid>
class Histogram : public NoTrace {
 public:
  using Index = typename Grid::Index;
  using Float = typename Grid::Float;

  /**
   * Create histogram to trace counts in given grid.
   * @param grid The atmosphere grid
   */
  Histogram(const Grid &grid) : i_(0), j_(0), k_(0) {
    std::array<size_t, 3> extent;
    extent[0] = grid.get_extent()[0] + 1;
    extent[1] = grid.get_extent()[1] + 1;
    extent[2] = grid.get_extent()[2] + 1;
    irradiance_ = std::move(extent);
  }

  /**
   * Photon creation is traced in order to store cell indices.
   */
  template <typename Photon>
  __DEV__ void created(const Photon &photon) {
    i_ = photon.i;
    j_ = photon.j;
    k_ = photon.k;
  }
  /**
   * The histogram only traces which cells a photons enters and
   * therefore only implements the step method.
   */
  template <typename Photon>
  __DEV__ void absorption(const Photon &photon, Float /*absorbed_energy*/) {
    atomic_add(irradiance_, photon.get_energy(), i_, j_, k_);
    i_ = photon.i;
    j_ = photon.j;
    k_ = photon.k;
  }

  /**
   * Save histogram results to file.
   * @param filename Name of the file to which to store results.
   */
  void save(std::string filename) {
    llrte::io::NetCDFFile file(filename, true);
    file.add_dimension("x", irradiance_.dimension(0));
    file.add_dimension("y", irradiance_.dimension(1));
    file.add_dimension("z", irradiance_.dimension(2));
    file.store_variable(irradiance_, "irradiance", {"x", "y", "z"});
  }

  #ifdef CUDA
  void device() {irradiance_.device();}
  void host() {irradiance_.host();}
  #endif

 public:
  Tensor<Float, 3> irradiance_{{1, 1, 1}};
  size_t i_, j_, k_;
};

////////////////////////////////////////////////////////////////////////////////
// Absorption Tracer
////////////////////////////////////////////////////////////////////////////////
/** AbsorptionTracer
 *
 * The absorption tracer class traces absorption across the atmospheric grid
 * as well as the different atmospheric boundaries.
 * @tparam Grid The 3D grid type which describes the atmosphere.
 */
template <typename Grid>
class AbsorptionTracer : public NoTrace {
 public:
  using Index = typename Grid::Index;
  using Float = typename Grid::Float;

  AbsorptionTracer(const Grid &grid) {
    std::array<size_t, 3> extent;
    extent[0] += grid.get_extent()[0] + 1;
    extent[1] += grid.get_extent()[1] + 1;
    extent[2] += grid.get_extent()[2] + 1;

    absorbed_energy_ = std::move(Tensor<Float, 3>{extent});
    scattered_energy_ = std::move(Tensor<Float, 3>{extent});

    leaving_photons_ = std::move(Tensor<Float, 1>(6));
    scattering_frequencies_ = std::move(Array<int>(11));
  }

  /**
   * Photon creation is traced in order to store cell indices.
   */
  template <typename Photon>
  __DEV__ void created(const Photon &photon) {
    i_ = photon.i;
    j_ = photon.j;
    k_ = photon.k;
  }

  /**
   * Records absorbed energy.
   */
  template <typename Photon>
  __DEV__ void absorption(const Photon &photon, Float absorbed_energy) {
    absorbed_energy_(i_, j_, k_) += absorbed_energy;
    i_ = photon.i;
    j_ = photon.j;
    k_ = photon.k;
  }

  template <typename Photon>
  void scattering(const Photon &photon) {
    scattered_energy_(i_, j_, k_) += photon.get_energy();
    i_ = photon.i;
    j_ = photon.j;
    k_ = photon.k;
  }

  /**
   * Trace out of energy. This method is called when a photon's energy falls
   * below it's minimum energy threshold.
   * @param Reference to the photon
   */
  template <typename Photon>
  void out_of_energy(const Photon &photon) {
    auto i = photon.get_scattering_events();
    if (i < 10) {
      ++scattering_frequencies_[i];
    } else {
      ++scattering_frequencies_[10];
    }
  }

  /**
   * Left atmosphere. This method is called when a photon leaves the atmosphere.
   * @param Reference to the photon
   */
  template <typename Photon, typename Atmosphere>
  __DEV__ void left_atmosphere(const Photon &photon, const Atmosphere &atmosphere) {
    auto i = atmosphere.get_boundary_index(photon);
    leaving_photons_[i] += photon.get_energy();
  }

  void save(std::string filename) {
    llrte::io::NetCDFFile file(filename, true);
    file.add_dimension("x", absorbed_energy_.dimensions(0));
    file.add_dimension("y", absorbed_energy_.dimensions(1));
    file.add_dimension("z", absorbed_energy_.dimensions(2));
    file.add_dimension("boundaries", 6);
    file.add_dimension("events", 11);
    file.store_variable(absorbed_energy_, "absorbed_energy", {"x", "y", "z"});
    file.store_variable(scattered_energy_, "scattered_energy", {"x", "y", "z"});
    file.store_variable(leaving_photons_, "leaving_photons_", {"boundaries"});
    file.store_variable(scattering_frequencies_, "scattering_frequencies",
                        {"events"});
  }

  Float get_total_absorption_counts() { return absorbed_energy_.sum(); }

  Float get_total_leaving_counts(size_t i) { return leaving_photons_[i]; }

#ifdef CUDA
  void device() {
    absorbed_energy_.device();
    scattered_energy_.device();
    scattering_frequencies_.device();
    leaving_photons_.device();
  }
  void host() {
    absorbed_energy_.host();
    scattered_energy_.host();
    scattering_frequencies_.host();
    leaving_photons_.host();
  }
#endif

   private:
    size_t i_, j_, k_;
    Tensor<Float, 3> absorbed_energy_{{1, 1, 1}};
    Tensor<Float, 3> scattered_energy_{{1, 1, 1}};
    Array<int> scattering_frequencies_{1};
    Array<Float> leaving_photons_{1};
  };

  ////////////////////////////////////////////////////////////////////////////////
  // Photon Tracer
  ////////////////////////////////////////////////////////////////////////////////

  /** Photon Tracer
   * Stores positions and direction of all photons that enter or leave the
   * domain.
   */
  template<typename Float>
  class PhotonTracer : public NoTrace {

   public:

  /** Create PhotonTracer for given number of photons
   *@param n: The number of photons
   */
    PhotonTracer(Index n)
        : n_(n),
          i_(0),
          incoming_positions_{{n, 3}},
          incoming_directions_{{n, 3}},
          outgoing_positions_{{n, 3}},
              outgoing_directions_{{n, 3}} {}

    template <typename Photon>
    void created(const Photon &photon) {
      incoming_positions_.coeffRef(i_, 0) = photon.position.x;
      incoming_positions_.coeffRef(i_, 1) = photon.position.y;
      incoming_positions_.coeffRef(i_, 2) = photon.position.z;
      incoming_directions_.coeffRef(i_, 0) = photon.direction.x;
      incoming_directions_.coeffRef(i_, 1) = photon.direction.y;
      incoming_directions_.coeffRef(i_, 2) = photon.direction.z;
    }

    template <typename Photon>
    void out_of_energy(const Photon &photon) {
      outgoing_positions_.coeffRef(i_, 0) = photon.position.x;
      outgoing_positions_.coeffRef(i_, 1) = photon.position.y;
      outgoing_positions_.coeffRef(i_, 2) = photon.position.z;
      outgoing_directions_.coeffRef(i_, 0) = photon.direction.x;
      outgoing_directions_.coeffRef(i_, 1) = photon.direction.y;
      outgoing_directions_.coeffRef(i_, 2) = photon.direction.z;
      i_ = (i_ + 1) % n_;
    }

    template <typename Photon, typename Atmosphere>
    void left_atmosphere(const Photon &photon,
                         const Atmosphere &) {
      outgoing_positions_.coeffRef(i_, 0) = photon.position.x;
      outgoing_positions_.coeffRef(i_, 1) = photon.position.y;
      outgoing_positions_.coeffRef(i_, 2) = photon.position.z;
      outgoing_directions_.coeffRef(i_, 0) = photon.direction.x;
      outgoing_directions_.coeffRef(i_, 1) = photon.direction.y;
      outgoing_directions_.coeffRef(i_, 2) = photon.direction.z;
      i_ = (i_ + 1) % n_;
    }

    void save(std::string filename) {
      llrte::io::NetCDFFile file(filename, true);
      file.add_dimension("photons", n_);
      file.add_dimension("coordinates", 3);
      file.store_variable(incoming_positions_, "incoming_positions",
                          {"photons", "coordinates"});
      file.store_variable(incoming_directions_, "incoming_directions",
                          {"photons", "coordinates"});
      file.store_variable(outgoing_positions_, "outgoing_positions",
                          {"photons", "coordinates"});
      file.store_variable(outgoing_directions_, "outgoing_directions",
                          {"photons", "coordinates"});
    }

   private:
    size_t n_, i_;
    Tensor<Float, 2> incoming_positions_;
    Tensor<Float, 2> incoming_directions_;
    Tensor<Float, 2> outgoing_positions_;
    Tensor<Float, 2> outgoing_directions_;
  };

  }  // namespace llrte::tracers
#endif
