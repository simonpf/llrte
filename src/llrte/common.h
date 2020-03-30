#ifndef _LLRTE_COMMON_H_
#define _LLRTE_COMMON_H_

#ifdef __CUDA__

#include <cuda_runtime_api.h>
#include <cuda.h>
#define CUDAERROR(X) if (X != cudaSuccess) { std::cout << "CUDA ERROR: " << cudaGetErrorString(X);}

#define __DEV__ __device__ __host__

#else

#define __DEV__

#endif
#endif
