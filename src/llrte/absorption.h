#ifndef _LLRTE_ABSORPTION_H_
#define _LLRTE_ABSORPTION_H_

namespace llrte {

template <typename Float>
class NoAbsorption {
 public:
  NoAbsorption() {}

  template <typename... Ts>
  constexpr Float get_absorption_coefficient(Ts...) {
    return 0.0;
  }
};

template <typename Float>
class ConstantAbsorption {
 public:
  ConstantAbsorption(Float absorption) : absorption_coefficient_(absorption) {}

  template <typename... Ts>
  Float get_absorption_coefficient(Ts...) {
    return absorption_coefficient_;
  }

 private:
  Float absorption_coefficient_;
};

}  // namespace llrte
#endif
