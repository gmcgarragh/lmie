/******************************************************************************%
**
**    Copyright (C) 2008-2020 Greg McGarragh <gregm@atmos.colostate.edu>
**
**    This source code is licensed under the GNU General Public License (GPL),
**    Version 2.  See the file COPYING for more details.
**
*******************************************************************************/

#include <stdio.h>

#include <cuda.h>
#include <cuda_runtime_api.h>

#include "cudacomplex.h"


#define MAX_DERIVS 7

#ifdef UNIX
#define PI	(M_PI)
#else
#define PI	3.14159265358979323846
#endif
/*
#define USE_DOUBLE_PRECISION
*/
#ifndef USE_DOUBLE_PRECISION

#define FLOAT		float
#define atomicAdd_x	atomicAdd_f
#define COMPLEX		singlecomplex
#define MAKE_COMPLEX	make_singlecomplex

#else

#define FLOAT		double
#define atomicAdd_x	atomicAdd_d
#define COMPLEX		doublecomplex
#define MAKE_COMPLEX	make_doublecomplex

#endif

/*******************************************************************************
 *
 ******************************************************************************/
#define CUDA_CALL_CHECK(ERROR) do {			\
     if (cuda_call_check(ERROR, __FILE__, __LINE__))	\
          return -1;					\
} while (0)


int cuda_call_check(cudaError_t error, char *file, int line) {

     if (error != cudaSuccess) {
          fprintf(stderr, "ERROR, file %s, line %d: %s\n", file, line, cudaGetErrorString(error));
          return -1;
     }

     return 0;
}



/*******************************************************************************
 *
 ******************************************************************************/
__device__ __host__ float realsc(singlecomplex a) {

      return a.value.x;
}

__device__  __host__ double realdc(doublecomplex a) {

      return a.value.x;
}



__device__ __host__ float imagsc(singlecomplex a) {

      return a.value.y;
}

__device__ __host__ double imagdc(doublecomplex a) {

      return a.value.y;
}



__device__ __host__ float abssc(singlecomplex a) {

      return sqrt(a.value.x*a.value.x + a.value.y*a.value.y);
}

__device__ __host__ double absdc(doublecomplex a) {

      return sqrt(a.value.x*a.value.x + a.value.y*a.value.y);
}


#ifndef USE_DOUBLE_PRECISION

#define realxc	realsc
#define imagxc	imagsc
#define absxc	abssc

#else

#define realxc	realdc
#define imagxc	imagdc
#define absxc	absdc

#endif


/*******************************************************************************
 *
 ******************************************************************************/
#if __CUDA_ARCH__ < 200

static __inline__ __device__ float atomicAdd_f(float *address, float val)
{
    float old= * address, assumed;

    do {
        assumed = old;
        old = int_as_float(atomicCAS((int *) address, float_as_int(assumed), float_as_int(val+assumed)));
    } while(float_as_int(assumed) != float_as_int(old));

    return old;
}

#endif // __CUDA_ARCH__ < 200

#if __CUDA_ARCH__ >= 130

static __inline__ __device__ double atomicAdd_d(double *address, double val)
{
    double old = *address, assumed;

    do {
        assumed = old;
        old = __longlong_as_double(atomicCAS((unsigned long long int *) address, __double_as_longlong(assumed), __double_as_longlong(val+assumed)));
    } while(__double_as_longlong(assumed) != __double_as_longlong(old));

    return old;
}

#endif // __CUDA_ARCH__ >= 130


/*******************************************************************************
 *
 ******************************************************************************/
#define MIN(a,b) (a < b ? a : b)
#define MAX(a,b) (a < b ? b : a)

#define A2(i,j) ((i) * n_derivs + (j))

#define USE_ATOMIC

/*
#define THREADSPERBLOCK 128
#include "lmie_cuda1.h"

#define THREADSPERBLOCK 1
#include "lmie_cuda2.h"
*/

#define THREADSPERBLOCK 64 /* hertz2, r2 = 32. */
#ifdef CRAP
#define THREADSPERBLOCK 320 /* amazon, r2 = 100. */
#endif
#include "lmie_cuda3.h"
