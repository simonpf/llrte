#include "llrte/data.h"

int main(int argc, char **argv) {
    auto t = llrte::Tensor<double, 2>({3, 3});
    t.fill(1.0);
    std::cout << t << std::endl;
}
