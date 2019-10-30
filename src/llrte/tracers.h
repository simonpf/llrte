#ifndef _LLRTE_TRACERS_
#define _LLRTE_TRACERS_

#include <fstream>

namespace llrte {

struct NoTrace {
    template <typename ... T>
    static void trace(T ... ) {
        // Nothing.
    }
};

template <typename Grid>
class Histogram {
public:

    using Index = typename Grid::Index;
    using Float = typename Grid::Float;

    static void initialize(const Grid &grid) {
        std::tie(shape_[0], shape_[1], shape_[2]) = grid.get_extent();
        size_t n = (shape_[0] - 1) * (shape_[1] - 1) * (shape_[2] - 1);
        data_ = std::shared_ptr<Float[]>(new Float[n]);
        for (size_t i = 0; i < n; ++i) {
            data_[i] = 0.0;
        }
    }


    template<typename GridPos, typename ... Ts>
    static void trace(GridPos gp, Event /*e*/, Ts ...) {
        size_t index = gp.k - 1;
        index *= (shape_[2] - 1);
        index += gp.j - 1;
        index *= (shape_[1] - 1);
        index += gp.i - 1;
        data_[index] += 1.0;
    }


    static void dump(std::string filename) {
        std::ofstream file;
        file.open (filename, std::ios::out | std::ios::binary); 
        size_t n = (shape_[0] - 1) * (shape_[1] - 1) * (shape_[2] - 1);
        file.write((char*) data_.get(), n * sizeof(Float));
        file.close();
    }

private:

    static Float shape_[3];
    static std::shared_ptr<Float[]> data_;
};

template<typename Grid>
typename Histogram<Grid>::Float Histogram<Grid>::shape_[3] = {0, 0, 0};
template<typename Grid>
std::shared_ptr<typename Histogram<Grid>::Float[]> Histogram<Grid>::data_ = nullptr;

}
#endif

