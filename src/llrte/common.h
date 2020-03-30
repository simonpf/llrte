#ifndef _LLRTE_COMMON_H_
#define _LLRTE_COMMON_H_

#ifdef __CUDA__

#include <cuda_runtime_api.h>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include <iostream>

#define CUDA_CALL(x) if((x) != cudaSuccess) {             \
            printf("Error at %s:%d\n",__FILE__,__LINE__); \
            }
#define CURAND_CALL(x) if((x)!=CURAND_STATUS_SUCCESS) {   \
            printf("Error at %s:%d\n",__FILE__,__LINE__); \
            }

#define __DEV__ __device__ __host__
#define __HOST__ __host__

#else

#define __DEV__
#define __HOST__

#endif
#endif
