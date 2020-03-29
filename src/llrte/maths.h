#ifndef _LLRTE_MATHS_H_
#define _LLRTE_MATHS_H_

#include "llrte/configurations.h"
#include "cmath"

namespace llrte::maths {

  template<typename Float>
  bool small(Float f) {
    auto eps = Limits<Float>::eps;
    return fabs(f) < eps;
  }

}  // namespace llrte::maths

#endif
