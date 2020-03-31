#ifndef _LLRTE_MATHS_H_
#define _LLRTE_MATHS_H_

#include "cmath"

#include "llrte/common.h"
#include "llrte/configurations.h"

namespace llrte::maths {

  template<typename Float>
  __DEV__ bool small(Float f) {
    auto eps = Limits<Float>::eps;
    return fabs(f) < eps;
  }

}  // namespace llrte::maths

#endif
