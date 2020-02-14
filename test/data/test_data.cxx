#include "llrte/data.h"

int main(int argc, char **argv) {
    auto t = llrte::Tensor<double, 3>({3, 3, 3});
    std::cout << t << std::endl;
}
