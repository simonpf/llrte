#ifndef _LLRTE_SOLVERS_MONTE_CARLO_H_
#define _LLRTE_SOLVERS_MONTE_CARLO_H_

#include <llrte/constants.h>
#include <llrte/random.h>
#include <math.h>

#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <tuple>

namespace llrte {

template <typename Atmosphere, typename Source>
class MonteCarloSolver {
 public:
  using Vector = typename std::remove_reference_t<Source>::Vector;
  using Float = typename std::remove_reference_t<Source>::Float;

  MonteCarloSolver(Atmosphere atmosphere, Source source)
      : atmosphere_(atmosphere), source_(source), generator_() {
    // Nothing to do here.
  }

  void sample_photon() {
    auto photon = source_.sample_photon();
    photon.propagate(atmosphere_, generator_);
  }

 private:
  size_t n_photons = 0;

  Atmosphere atmosphere_;
  Source source_;
  Generator<Vector> generator_;
};

}  // namespace llrte

#endif
