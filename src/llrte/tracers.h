#ifndef _LLRTE_TRACERS_
#define _LLRTE_TRACERS_

#include <fstream>

#include "llrte/definitions.h"

namespace llrte {

////////////////////////////////////////////////////////////////////////////////
// No Trace
////////////////////////////////////////////////////////////////////////////////

struct NoTrace {
  template <typename... T>
  static void trace(T...) {
    // Nothing.
  }
};

////////////////////////////////////////////////////////////////////////////////
// Histogram
////////////////////////////////////////////////////////////////////////////////

template <typename Grid>
class Histogram {
 public:
  using Index = typename Grid::Index;
  using Float = typename Grid::Float;

  static void initialize(const Grid &grid) {
    std::tie(shape_[0], shape_[1], shape_[2]) = grid.get_extent();
    size_t n = (shape_[0] - 1) * (shape_[1] - 1) * (shape_[2] - 1);
    data_ = std::make_unique<Float[]>(n);
    for (size_t i = 0; i < n; ++i) {
      data_[i] = 0.0;
    }
  }

  template <typename Photon, typename GridPos, typename... Ts>
  static void trace(const Photon & /*p*/, GridPos gp, Event e, Ts...) {
    if (e == Event::step) {
      if ((gp.i < shape_[0]) && (gp.j < shape_[1]) && (gp.k < shape_[2])) {
        size_t index = gp.k - 1;
        index *= (shape_[2] - 1);
        index += gp.j - 1;
        index *= (shape_[1] - 1);
        index += gp.i - 1;
        data_[index] += 1.0;
      }
    }
  }

  static void dump(std::string filename) {
    std::ofstream file;
    file.open(filename, std::ios::out | std::ios::binary);
    size_t n = (shape_[0] - 1) * (shape_[1] - 1) * (shape_[2] - 1);
    file.write((char *)data_.get(), n * sizeof(Float));
    file.close();
  }

 private:
  static Float shape_[3];
  static std::unique_ptr<Float[]> data_;
};

template <typename Grid>
typename Histogram<Grid>::Float Histogram<Grid>::shape_[3] = {0, 0, 0};
template <typename Grid>
std::unique_ptr<typename Histogram<Grid>::Float[]> Histogram<Grid>::data_ =
    nullptr;

////////////////////////////////////////////////////////////////////////////////
// Absorption Tracer
////////////////////////////////////////////////////////////////////////////////

template <typename Grid>
class AbsorptionTracer {
 public:
  using Index = typename Grid::Index;
  using Float = typename Grid::Float;

  static void initialize(const Grid &grid) {
    grid_ = &grid;
    std::tie(shape_[0], shape_[1], shape_[2]) = grid.get_extent();
    size_t n = (shape_[0] - 1) * (shape_[1] - 1) * (shape_[2] - 1);

    absorption_counts_ = std::make_unique<Float[]>(n);
    scattering_counts_ = std::make_unique<Float[]>(n);

    for (size_t i = 0; i < n; ++i) {
      absorption_counts_[i] = 0.0;
      scattering_counts_[i] = 0.0;
    }

    n = 3 * 2;
    leaving_photons_ = std::make_unique<Float[]>(n);
    for (size_t i = 0; i < n; ++i) {
      leaving_photons_[i] = 0;
    }

    n = 11;
    scattering_frequencies_ = std::make_unique<Float[]>(n);
    for (size_t i = 0; i < n; ++i) {
      scattering_frequencies_[i] = 0;
    }
  }

  template <typename GridPos>
  static size_t get_index(const GridPos &gp) {
      size_t index = (gp.k - 1) % (shape_[2] - 1);
      index *= (shape_[1] - 1);

      index += (gp.j - 1) % (shape_[1] - 1);
      index *= (shape_[0] - 1);

      index += (gp.i - 1) % (shape_[0] - 1);
      return index;
  }

  template <typename Photon, typename GridPos>
  static void trace(const Photon &p, GridPos gp, Event e) {
    if (e == Event::absorption) {
      size_t index = get_index(gp);
      absorption_counts_[index] += 1.0;
    } else if (e == Event::scattering) {
      size_t index = get_index(gp);
      scattering_counts_[index] += 1.0;
    } else if (e == Event::left_domain) {
      auto i = grid_->get_boundary_index(gp);
      leaving_photons_[i] += 1.0;
      i = p.get_scattering_events();
      if (i < 10) {
        scattering_frequencies_[i] += 1.0;
      } else {
        scattering_frequencies_[10] += 1.0;
      }
    } else if (e == Event::out_of_energy) {
      auto i = p.get_scattering_events();
      if (i < 10) {
        scattering_frequencies_[i] += 1.0;
      } else {
        scattering_frequencies_[10] += 1.0;
      }
    }
  }

  template <typename Photon, typename GridPos, typename... Ts>
  static void trace(const Photon &p, GridPos gp, Event e,
                    typename Photon::Float value, Ts...) {
    if (e == Event::absorption) {
      size_t index = get_index(gp);
      absorption_counts_[index] += value;
    } else if (e == Event::scattering) {
      size_t index = get_index(gp);
      scattering_counts_[index] += value;
    } else if (e == Event::left_domain) {
      auto i = grid_->get_boundary_index(gp);
      leaving_photons_[i] += value;
      i = p.get_scattering_events();
      if (i < 10) {
        scattering_frequencies_[i] += 1.0;
      } else {
        scattering_frequencies_[10] += 1.0;
      }
    } else if (e == Event::out_of_energy) {
      auto i = p.get_scattering_events();
      if (i < 10) {
        scattering_frequencies_[i] += 1.0;
      } else {
        scattering_frequencies_[10] += 1.0;
      }
    }
  }

  static void dump(std::string filename) {
    std::ofstream file;
    file.open(filename, std::ios::out | std::ios::binary);
    size_t n = (shape_[0] - 1) * (shape_[1] - 1) * (shape_[2] - 1);
    file.write((char *)absorption_counts_.get(), n * sizeof(Float));
    file.write((char *)scattering_counts_.get(), n * sizeof(Float));
    file.write((char *)leaving_photons_.get(), 6 * sizeof(Float));
    file.write((char *)scattering_frequencies_.get(), 11 * sizeof(Float));
    file.close();
  }

  static Float get_total_absorption_counts() {
    Float sum = 0.0;
    size_t n = (shape_[0] - 1) * (shape_[1] - 1) * (shape_[2] - 1);
    for (size_t i = 0; i < n; ++i) {
      sum += absorption_counts_[i];
    }
    return sum;
  }

  static Float get_total_leaving_counts(size_t i) {
    return leaving_photons_[i];
  }

 private:
  static size_t shape_[3];
  static std::unique_ptr<Float[]> absorption_counts_;
  static std::unique_ptr<Float[]> scattering_counts_;
  static std::unique_ptr<Float[]> scattering_frequencies_;
  static std::unique_ptr<Float[]> leaving_photons_;
  static const Grid *grid_;
};

template <typename Grid>
size_t AbsorptionTracer<Grid>::shape_[3] = {0, 0, 0};
template <typename Grid>
std::unique_ptr<typename Grid::Float[]>
    AbsorptionTracer<Grid>::absorption_counts_{nullptr};
template <typename Grid>
std::unique_ptr<typename Grid::Float[]>
    AbsorptionTracer<Grid>::scattering_counts_{nullptr};
template <typename Grid>
std::unique_ptr<typename Grid::Float[]>
    AbsorptionTracer<Grid>::scattering_frequencies_{nullptr};
template <typename Grid>
std::unique_ptr<typename Grid::Float[]>
    AbsorptionTracer<Grid>::leaving_photons_{nullptr};
template <typename Grid>
const Grid *AbsorptionTracer<Grid>::grid_{nullptr};

////////////////////////////////////////////////////////////////////////////////
// Photon Tracer
////////////////////////////////////////////////////////////////////////////////

template <typename Photon>
class PhotonTracer {
 public:
  using Float = typename Photon::Float;

  template <typename GridPos, typename... Ts>
  static void trace(const Photon &/*p*/, GridPos /*gp*/, Event /*e*/, Ts...) {}

  template <typename GridPos, typename... Ts>
  static void trace(const Photon &p, GridPos /*gp*/, Event e,
                    typename Photon::Float /*value*/, Ts...) {
    if (e == Event::out_of_energy) {
      photons_.push_back(p);
    } else if (e == Event::left_domain) {
      photons_.push_back(p);
    }
  }

  static std::vector<Photon> & get_photons() { return photons_; }

 private:
  static std::vector<Photon> photons_;
};

template <typename Photon>
std::vector<Photon> PhotonTracer<Photon>::photons_{};

}  // namespace llrte
#endif
