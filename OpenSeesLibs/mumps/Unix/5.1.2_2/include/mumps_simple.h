#include "mpi.h"
#include "smumps_c.h"  // Real, single precision
#include "dmumps_c.h"  // Real, double precision
#include "cmumps_c.h"  // Complex, single precision
#include "zmumps_c.h"  // Complex, double precision
#include "complex.h"

void* mumps_initialize_float(int, int*, float*);
void  mumps_associate_matrix_float(void*, int, int, float*, int*, int*);
void  mumps_factorize_float(void*);
void  mumps_associate_rhs_float(void*, int, float*);
void  mumps_solve_float(void*, int*);
int   mumps_get_nrhs_float(void*);
void  mumps_get_solution_float(void*, float*);
void  mumps_get_info_float(void*, int*, float*);
void  mumps_finalize_float(void*);
void  mumps_alloc_float(SMUMPS_STRUC_C**);
void  mumps_free_float(SMUMPS_STRUC_C**);

void* mumps_initialize_double(int, int*, double*);
void  mumps_associate_matrix_double(void*, int, int, double*, int*, int*);
void  mumps_factorize_double(void*);
void  mumps_associate_rhs_double(void*, int, double*);
void  mumps_solve_double(void*, int*);
int   mumps_get_nrhs_double(void*);
void  mumps_get_solution_double(void*, double*);
void  mumps_get_info_double(void*, int*, double*);
void  mumps_finalize_double(void*);
void  mumps_alloc_double(DMUMPS_STRUC_C**);
void  mumps_free_double(DMUMPS_STRUC_C**);

void* mumps_initialize_float_complex(int, int*, float*);
void  mumps_associate_matrix_float_complex(void*, int, int, float complex*, int*, int*);
void  mumps_factorize_float_complex(void*);
void  mumps_associate_rhs_float_complex(void*, int, float complex*);
void  mumps_solve_float_complex(void*, int*);
int   mumps_get_nrhs_float_complex(void*);
void  mumps_get_solution_float_complex(void*, float complex*);
void  mumps_get_info_float_complex(void*, int*, float*);
void  mumps_finalize_float_complex(void*);
void  mumps_alloc_float_complex(CMUMPS_STRUC_C**);
void  mumps_free_float_complex(CMUMPS_STRUC_C**);

void* mumps_initialize_double_complex(int, int*, double*);
void  mumps_associate_matrix_double_complex(void*, int, int, double complex*, int*, int*);
void  mumps_factorize_double_complex(void*);
void  mumps_associate_rhs_double_complex(void*, int, double complex*);
void  mumps_solve_double_complex(void*, int*);
int   mumps_get_nrhs_double_complex(void*);
void  mumps_get_solution_double_complex(void*, double complex*);
void  mumps_get_info_double_complex(void*, int*, double*);
void  mumps_finalize_double_complex(void*);
void  mumps_alloc_double_complex(ZMUMPS_STRUC_C**);
void  mumps_free_double_complex(ZMUMPS_STRUC_C**);
