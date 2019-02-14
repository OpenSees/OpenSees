/*
* Copyright (c) 1997-1999, 2003 Massachusetts Institute of Technology
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
*
* 1. Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
*
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
* OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
* GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
* WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef FFTW_MPI_H
#define FFTW_MPI_H

#include <mpi.h> /* need access to the MPI type definitions */
#include "fftw.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef MKL_INT
#define MKL_INT int
#endif

/***********************************************************************/

typedef fftw_real TRANSPOSE_EL_TYPE;

typedef struct {
     int block_num, dest_pe, send_size, recv_size;
} transpose_mpi_exchange;

typedef struct {
     MPI_Comm comm;
     int n_pes, my_pe;

     int nx,ny,local_nx,local_ny;

     transpose_mpi_exchange *exchange;
     int num_steps, send_block_size, recv_block_size;

     MPI_Datatype el_type;

     MPI_Request request[2];

     int *perm_block_dest;
     int num_perm_blocks, perm_block_size;

     int all_blocks_equal;
     int *send_block_sizes, *send_block_offsets;
     int *recv_block_sizes, *recv_block_offsets;

     char *move;
     int move_size;
} transpose_mpi_plan_struct;

typedef transpose_mpi_plan_struct *transpose_mpi_plan;

extern void transpose_mpi_get_local_size(int n, int my_pe, int n_pes,
					 int *local_n, int *local_start);
extern int transpose_mpi_get_local_storage_size(int nx, int ny,
						int my_pe, int n_pes);

extern transpose_mpi_plan transpose_mpi_create_plan(int nx, int ny,
						    MPI_Comm comm);
extern void transpose_mpi_destroy_plan(transpose_mpi_plan p);

extern void transpose_mpi(transpose_mpi_plan p, int el_size,
			  TRANSPOSE_EL_TYPE *local_data,
			  TRANSPOSE_EL_TYPE *work);

typedef enum { BEFORE_TRANSPOSE, AFTER_TRANSPOSE } transpose_in_place_which;

typedef enum { TRANSPOSE_SYNC, TRANSPOSE_ASYNC } transpose_sync_type;

extern void transpose_in_place_local(transpose_mpi_plan p,
                              int el_size, TRANSPOSE_EL_TYPE *local_data,
                              transpose_in_place_which which);

extern TRANSPOSE_EL_TYPE *transpose_allocate_send_buf(transpose_mpi_plan p,
						      int el_size);
extern void transpose_get_send_block(transpose_mpi_plan p, int step,
				     int *block_y_start, int *block_ny);
extern void transpose_start_exchange_step(transpose_mpi_plan p,
					  int el_size,
					  TRANSPOSE_EL_TYPE *local_data,
					  TRANSPOSE_EL_TYPE *send_buf,
					  int step,
					  transpose_sync_type sync_type);
extern void transpose_finish_exchange_step(transpose_mpi_plan p, int step);

/***********************************************************************/

typedef struct {
     fftw_plan p_fft_x;  /* plan for first dimension */
     fftwnd_plan p_fft;  /* plan for subsequent dimensions */
     transpose_mpi_plan p_transpose, p_transpose_inv;
     fftw_complex *work; /* extra workspace, if needed */
} fftwnd_mpi_plan_data;

typedef fftwnd_mpi_plan_data *fftwnd_mpi_plan;

typedef enum {
    FFTW_NORMAL_ORDER,
    FFTW_TRANSPOSED_ORDER
} fftwnd_mpi_output_order;

extern fftwnd_mpi_plan fftwnd_mpi_create_plan(MPI_Comm comm,
					      int rank, const MKL_INT *n,
					      fftw_direction dir,
					      int flags);
extern fftwnd_mpi_plan fftw2d_mpi_create_plan(MPI_Comm comm,
					      MKL_INT nx, MKL_INT ny,
					  fftw_direction dir, int flags);
extern fftwnd_mpi_plan fftw3d_mpi_create_plan(MPI_Comm comm,
					      MKL_INT nx, MKL_INT ny, MKL_INT nz,
					  fftw_direction dir, int flags);

extern void fftwnd_mpi_destroy_plan(fftwnd_mpi_plan p);

extern void fftwnd_mpi_local_sizes(fftwnd_mpi_plan p,
				   MKL_INT *local_nx,
				   MKL_INT *local_x_start,
				   MKL_INT *local_ny_after_transpose,
				   MKL_INT *local_y_start_after_transpose,
				   MKL_INT *total_local_size);

extern void fftwnd_mpi(fftwnd_mpi_plan p,
		       int n_fields,
		       fftw_complex *local_data, fftw_complex *work,
		       fftwnd_mpi_output_order output_order);

extern void fftw_mpi_die(const char *error_string);

/***********************************************************************/

typedef struct fftw_mpi_twiddle_struct {
     int rows, rowstart, cols, n;
     fftw_complex *W;
     int refcount;
     struct fftw_mpi_twiddle_struct *next;
} fftw_mpi_twiddle;

typedef struct fftw_mpi_plan_struct {
     int n, m, r, local_m, local_m_start, local_r, local_r_start;
     fftw_complex *fft_work;
     fftw_mpi_twiddle *tw;
     transpose_mpi_plan p_transpose, p_transpose_inv;
     fftw_plan pm, pr;
     int flags;
     fftw_direction dir;
} *fftw_mpi_plan;

/* new flags for the MPI planner: */
#define FFTW_SCRAMBLED_INPUT (8192)
#define FFTW_SCRAMBLED_OUTPUT (16384)

extern void fftw_mpi_local_sizes(fftw_mpi_plan p,
				 MKL_INT *local_n,
				 MKL_INT *local_start,
				 MKL_INT *local_n_after_transform,
				 MKL_INT *local_start_after_transform,
				 MKL_INT *total_local_size);

extern fftw_mpi_plan fftw_mpi_create_plan(MPI_Comm comm,
					  MKL_INT n,
					  fftw_direction dir, int flags);

extern void fftw_mpi_destroy_plan(fftw_mpi_plan p);

extern void fftw_mpi(fftw_mpi_plan p, int n_fields,
		     fftw_complex *local_data, fftw_complex *work);

extern void fftw_mpi_print_plan(fftw_mpi_plan p);

/***********************************************************************/

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* FFTW_MPI_H */
