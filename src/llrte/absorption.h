#ifndef _LLRTE_ABSORPTION_H_

namespace llrte {
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
