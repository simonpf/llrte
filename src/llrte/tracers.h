#ifndef _LLRTE_TRACERS_
#define _LLRTE_TRACERS_

#include "llrte/common.h"
#include "llrte/definitions.h"
#include "llrte/io/netcdf.h"

namespace llrte::tracers {

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
  template <typename Photon, typename Float>
  __DEV__ void created(const Photon & /*photon*/, Float /*energy*/) {}
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
    irradiance_ = std::move(Tensor<Float, 3>(extent));
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
    irradiance_.atomic_add(photon.get_energy(), i_, j_, k_);
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
    file.add_dimension("x", irradiance_.shape()[0]);
    file.add_dimension("y", irradiance_.shape()[1]);
    file.add_dimension("z", irradiance_.shape()[2]);
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

    leaving_photons_ = std::move(Array<int>(6));
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
    ++leaving_photons_[i];
  }

  void save(std::string filename) {
    llrte::io::NetCDFFile file(filename, true);
    file.add_dimension("x", absorbed_energy_.shape()[0]);
    file.add_dimension("y", absorbed_energy_.shape()[1]);
    file.add_dimension("z", absorbed_energy_.shape()[2]);
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
    Array<int> leaving_photons_{1};
  };

  ////////////////////////////////////////////////////////////////////////////////
  // Photon Tracer
  ////////////////////////////////////////////////////////////////////////////////

  /** Photon Tracer
   * Stores all leaving photons in a vector.
   */
  template <typename Photon>
  class PhotonTracer : public NoTrace {
   public:
    void out_of_energy(const Photon &photon) { photons_.push_back(photon); }

    template <typename Atmosphere>
    __DEV__ void left_atmosphere(const Photon &photon,
                                 const Atmosphere /*&atmosphere*/) {
      photons_.push_back(photon);
    }

    static std::vector<Photon> &get_photons() { return photons_; }

   private:
    static std::vector<Photon> photons_;
  };

}  // namespace llrte::tracers
#endif
