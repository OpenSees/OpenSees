/*! @file
 * <pre>
 * -- Distributed SuperLU routine (version 4.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * October 1, 2014
 * </pre>
 */

#ifndef CUBLAS_UTILS_H
#define CUBLAS_UTILS_H

#include <cublas_v2.h>
#include "cuda.h"
#include "cuda_runtime_api.h"
#include "cuda_runtime.h"

extern void DisplayHeader();
extern const char* cublasGetErrorString(cublasStatus_t status);
extern cudaError_t checkCuda(cudaError_t);
extern cublasStatus_t checkCublas(cublasStatus_t);
extern cublasHandle_t create_handle ();
extern void destroy_handle (cublasHandle_t handle);

#endif 
