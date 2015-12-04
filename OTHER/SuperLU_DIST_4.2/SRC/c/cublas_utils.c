#include <stdio.h>
#include "cublas_utils.h"

 void DisplayHeader()
{
    const int kb = 1024;
    const int mb = kb * kb;
    // cout << "NBody.GPU" << endl << "=========" << endl << endl;

    printf("CUDA version:   v %d\n",CUDART_VERSION);
    //cout << "Thrust version: v" << THRUST_MAJOR_VERSION << "." << THRUST_MINOR_VERSION << endl << endl; 

    int devCount;
    cudaGetDeviceCount(&devCount);
    printf( "CUDA Devices: \n \n"); 

    for(int i = 0; i < devCount; ++i)
    {
        struct cudaDeviceProp props;       
        cudaGetDeviceProperties(&props, i);
        printf("%d : %s %d %d\n",i, props.name,props.major,props.minor );
        // cout << i << ": " << props.name << ": " << props.major << "." << props.minor << endl;
        printf("  Global memory:   %ld mb \n", props.totalGlobalMem / mb);
        // cout << "  Global memory:   " << props.totalGlobalMem / mb << "mb" << endl;
        printf("  Shared memory:   %ld kb \n", props.sharedMemPerBlock / kb ); //<<  << "kb" << endl;
        printf("  Constant memory: %ld kb \n", props.totalConstMem / kb );
        printf("  Block registers: %d \n\n", props.regsPerBlock );

        // to do these later
        // printf("  Warp size:         %d" << props.warpSize << endl;
        // printf("  Threads per block: %d" << props.maxThreadsPerBlock << endl;
        // printf("  Max block dimensions: [ %d" << props.maxThreadsDim[0] << ", " << props.maxThreadsDim[1]  << ", " << props.maxThreadsDim[2] << " ]" << endl;
        // printf("  Max grid dimensions:  [ %d" << props.maxGridSize[0] << ", " << props.maxGridSize[1]  << ", " << props.maxGridSize[2] << " ]" << endl;

        // cout << "  Shared memory:   " << props.sharedMemPerBlock / kb << "kb" << endl;
        // cout << "  Constant memory: " << props.totalConstMem / kb << "kb" << endl;
        // cout << "  Block registers: " << props.regsPerBlock << endl << endl;

        // cout << "  Warp size:         " << props.warpSize << endl;
        // cout << "  Threads per block: " << props.maxThreadsPerBlock << endl;
        // cout << "  Max block dimensions: [ " << props.maxThreadsDim[0] << ", " << props.maxThreadsDim[1]  << ", " << props.maxThreadsDim[2] << " ]" << endl;
        // cout << "  Max grid dimensions:  [ " << props.maxGridSize[0] << ", " << props.maxGridSize[1]  << ", " << props.maxGridSize[2] << " ]" << endl;
        // cout << endl;
    }
}


const char* cublasGetErrorString(cublasStatus_t status)
{
    switch(status)
    {
        case CUBLAS_STATUS_SUCCESS: return "CUBLAS_STATUS_SUCCESS";
        case CUBLAS_STATUS_NOT_INITIALIZED: return "CUBLAS_STATUS_NOT_INITIALIZED";
        case CUBLAS_STATUS_ALLOC_FAILED: return "CUBLAS_STATUS_ALLOC_FAILED";
        case CUBLAS_STATUS_INVALID_VALUE: return "CUBLAS_STATUS_INVALID_VALUE"; 
        case CUBLAS_STATUS_ARCH_MISMATCH: return "CUBLAS_STATUS_ARCH_MISMATCH"; 
        case CUBLAS_STATUS_MAPPING_ERROR: return "CUBLAS_STATUS_MAPPING_ERROR";
        case CUBLAS_STATUS_EXECUTION_FAILED: return "CUBLAS_STATUS_EXECUTION_FAILED"; 
        case CUBLAS_STATUS_INTERNAL_ERROR: return "CUBLAS_STATUS_INTERNAL_ERROR"; 
    }
    return "unknown error";
}

inline
cudaError_t checkCuda(cudaError_t result)
{
#if defined(DEBUG) || defined(_DEBUG)
    if (result != cudaSuccess) {
        fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
        assert(result == cudaSuccess);
    }
#endif
    return result;
}

cublasStatus_t checkCublas(cublasStatus_t result)
{
#if defined(DEBUG) || defined(_DEBUG)
  if (result != CUBLAS_STATUS_SUCCESS) {
    fprintf(stderr, "CUDA Blas Runtime Error: %s\n", cublasGetErrorString(result));
    assert(result == CUBLAS_STATUS_SUCCESS);
  }
#endif
  return result;
}


cublasHandle_t create_handle ()
{
       cublasHandle_t handle;
       checkCublas(cublasCreate(&handle));
       return handle;
 }

 void destroy_handle (cublasHandle_t handle)
 {
      checkCublas(cublasDestroy(handle));
 }

