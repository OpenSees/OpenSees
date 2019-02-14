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

#ifndef RFFTW_MPI_H
#define RFFTW_MPI_H

#include "fftw_mpi.h"
#include "rfftw.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef MKL_INT
#define MKL_INT int
#endif

/***********************************************************************/

typedef struct {
     fftw_plan p_fft_x;  /* plan for first dimension */
     rfftwnd_plan p_fft;  /* plan for subsequent dimensions */
     transpose_mpi_plan p_transpose, p_transpose_inv;
     fftw_complex *work; /* extra workspace, if needed */
} rfftwnd_mpi_plan_data;

typedef rfftwnd_mpi_plan_data *rfftwnd_mpi_plan;

extern rfftwnd_mpi_plan rfftwnd_mpi_create_plan(MPI_Comm comm,
					      int rank, const MKL_INT *n,
					      fftw_direction dir,
					      int flags);
extern rfftwnd_mpi_plan rfftw2d_mpi_create_plan(MPI_Comm comm,
					      MKL_INT nx, MKL_INT ny,
					  fftw_direction dir, int flags);
extern rfftwnd_mpi_plan rfftw3d_mpi_create_plan(MPI_Comm comm,
					      MKL_INT nx, MKL_INT ny, MKL_INT nz,
					  fftw_direction dir, int flags);

extern void rfftwnd_mpi_destroy_plan(rfftwnd_mpi_plan p);

extern void rfftwnd_mpi_local_sizes(rfftwnd_mpi_plan p,
				   MKL_INT *local_nx,
				   MKL_INT *local_x_start,
				   MKL_INT *local_ny_after_transpose,
				   MKL_INT *local_y_start_after_transpose,
				   MKL_INT *total_local_size);

extern void rfftwnd_mpi(rfftwnd_mpi_plan p,
		       int n_fields,
		       fftw_real *local_data, fftw_real *work,
		       fftwnd_mpi_output_order output_order);

/***********************************************************************/

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* RFFTW_MPI_H */
