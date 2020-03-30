#ifndef _LLRTE_COMMON_H_
#define _LLRTE_COMMON_H_

#ifdef __CUDA__

#include <cuda_runtime_api.h>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <iostream>

//
// Handle CUDA error.
//

#define CUDA_CALL(x) print_cuda_error(x, __FILE__, __LINE__, false);
inline void print_cuda_error(cudaError_t code,
                             const char* file,
                             int line,
                             bool abort = false) {
    if (code != cudaSuccess) {
        std::cerr << "CUDA Error in " << file << ", l.: " << line << ":" << std::endl;
        std::cerr << cudaGetErrorString(code) << std::endl;
    }
}

#define CURAND_CALL(x) print_curand_error(x, __FILE__, __LINE__, false);
inline void print_curand_error(int code,
                               const char *file,
                               int line,
                               bool abort = false) {
    if (code != CURAND_STATUS_SUCCESS) {
        std::cerr << "CURAND Error in " << file << ", l.: " << line << std::endl;
    }
}

#define __DEV__ __device__ __host__
#define __HOST__ __host__

#else

#define __DEV__
#define __HOST__

#endif
#endif
