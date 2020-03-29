#ifndef _LLRTE_MATHS_H_
#define _LLRTE_MATHS_H_

#include "llrte/configurations.h"

namespace llrte::maths {

  template<typename Float>
  bool small(Float f) {
    auto eps = Limits<Float>::eps;
    return abs(f) < eps;
  }

}  // namespace llrte::maths

#endif
