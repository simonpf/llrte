#ifndef _LLRTE_CONFIGURATIONS_H_
#define _LLRTE_CONFIGURATIONS_H_

namespace llrte {

  template<typename F>
  struct Limits;

  template<>
  struct Limits<float> {
    static constexpr float eps = 1e-2;
  };

  template<>
  struct Limits<double> {
    static constexpr double eps = 1e-6;
  };
} // namspace llrte
#endif

