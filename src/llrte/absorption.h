#ifndef _LLRTE_ABSORPTION_H_

namespace llrte {
template <typename Float>
class ConstantAbsorption {
public:
ConstantAbsorption(Float absorption)
    : absorption_(absorption){}

    template <typename ... Ts>
        Float get_absorption(Ts ...) {
        return absorption_;
    }

private:
    Float absorption_;
};

}
#endif
