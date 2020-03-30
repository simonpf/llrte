#ifndef _LLRTE_TRACERS_
#define _LLRTE_TRACERS_

#include <fstream>

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
  void created(const Photon &photon, Float energy) {}
  /**
   * Trace absorption. This method is called for every step of the MC algorithm
   * and can be used to trace photon position as well as absorbed energy.
   * @param Reference to the photon
   */
  template <typename Photon, typename Float>
  void absorption(const Photon &photon, Float energy) {}

  /**
   * Trace out of energy. This method is called when a photon's energy falls
   * below it's minimum energy threshold.
   * @param Reference to the photon
   */
  template <typename Photon>
  void out_of_energy(const Photon &photon) {}

  /**
   * Trace scattering. This method is called when a scattering event occurs.
   * @param Reference to the photon
   */
  template <typename Photon>
  void scattering(const Photon &photon) {}

  /**
   * Left atmosphere. This method is called when a photon leaves the atmosphere.
   * @param Reference to the photon
   */
  template <typename Photon>
  void left_atmosphere(const Photon &phton) {}
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
      extent[0] += grid.get_extent()[0] + 1;
      extent[1] += grid.get_extent()[1] + 1;
      extent[2] += grid.get_extent()[2] + 1;
    counts_ = Tensor<int, 3>(extent);
  }

  /**
   * Photon creation is traced in order to store cell indices.
   */
  template <typename Photon>
  void created(const Photon &photon) {
    i_ = photon.i;
    j_ = photon.j;
    k_ = photon.k;
  }
  /**
   * The histogram only traces which cells a photons enters and
   * therefore only implements the step method.
   */
  template <typename Photon>
  void absorption(const Photon &photon, Float /*absorbed_energy*/) {
    counts_(i_, j_, k_) += 1;
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
    file.add_dimension("x", counts_.shape()[0]);
    file.add_dimension("y", counts_.shape()[1]);
    file.add_dimension("z", counts_.shape()[2]);
    file.store_variable(counts_, "counts", {"x", "y", "z"});
  }

 private:
  Tensor<int, 3> counts_{{0, 0, 0}};
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
    auto extent = grid.get_extent();
    extent[0] += 1;
    extent[1] += 1;
    extent[2] += 1;

    absorbed_energy_ = Tensor<Float, 3>{extent};
    scattered_energy_ = Tensor<Float, 3>{extent};

    leaving_photons_ = Array<int>(6);
    scattering_frequencies_ = Array<int>(11);
  }

  /**
   * Photon creation is traced in order to store cell indices.
   */
  template <typename Photon>
  void created(const Photon &photon) {
    i_ = photon.i;
    j_ = photon.j;
    k_ = photon.k;
  }

  /**
   * Records absorbed energy.
   */
  template <typename Photon>
  void absorption(const Photon &photon, Float absorbed_energy) {
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
  void left_atmosphere(const Photon &photon, const Atmosphere &atmosphere) {
    auto i = atmosphere.get_boundary_index(atmosphere);
    ++leaving_photons_[i];
  }

  void save(std::string filename) {
    llrte::io::NetCDFFile file(filename, true);
    file.add_dimension("x", absorbed_energy_.shape()[0]);
    file.add_dimension("y", absorbed_energy_.shape()[1]);
    file.add_dimension("z", absorbed_energy_.shape()[2]);
    file.store_variable(absorbed_energy_, "absorbed_energy", {"x", "y", "z"});
    file.store_variable(scattered_energy_, "scattered_energy", {"x", "y", "z"});
    file.add_dimension("boundaries", 6);
    file.store_variable(leaving_photons_, "leaving_photons_", {"boundaries"});
    file.add_dimension("events", 11);
    file.store_variable(scattering_frequencies_, "scattering_frequencies", {"events"});
  }

  Float get_total_absorption_counts() { return absorbed_energy_.sum(); }

  Float get_total_leaving_counts(size_t i) {
    return leaving_photons_[i];
  }

 private:
  size_t i_, j_, k_;
  Tensor<Float, 3> absorbed_energy_{{0, 0, 0}};
  Tensor<Float, 3> scattered_energy_{{0, 0, 0}};
  Array<int> scattering_frequencies_{0};
  Array<int> leaving_photons_{0};
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
  void left_atmosphere(const Photon &photon, const Atmosphere /*&atmosphere*/) {
      photons_.push_back(photon);
  }

  static std::vector<Photon> &get_photons() { return photons_; }

 private:
  static std::vector<Photon> photons_;
};

}  // namespace llrte::tracers
#endif
