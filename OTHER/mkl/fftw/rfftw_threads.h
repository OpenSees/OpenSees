/*******************************************************************************
!   Copyright (c) 2003 Matteo Frigo
!   Copyright (c) 2003 Massachusetts Institute of Technology
!
!   This program is distributed with permission
!
!*******************************************************************************/

#ifndef RFFTW_THREADS_H
#define RFFTW_THREADS_H

#include "rfftw.h"
#include "fftw_threads.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/******************** User Interface *********************/

extern void rfftw_threads(int nthreads,
                   fftw_plan plan, int howmany, fftw_real *in, int istride,
                   int idist, fftw_real *out, int ostride, int odist);
extern void rfftw_threads_one(int nthread, fftw_plan plan,
			      fftw_real *in, fftw_real *out);

extern void rfftwnd_threads_real_to_complex(int nthreads, fftwnd_plan p,
					    int howmany,
					    fftw_real *in,
					    int istride, int idist,
					    fftw_complex *out,
					    int ostride, int odist);
extern void rfftwnd_threads_complex_to_real(int nthreads, fftwnd_plan p,
					    int howmany,
					    fftw_complex *in,
					    int istride, int idist,
					    fftw_real *out,
					    int ostride, int odist);
extern void rfftwnd_threads_one_real_to_complex(int nthreads, fftwnd_plan p,
						fftw_real *in,
						fftw_complex *out);
extern void rfftwnd_threads_one_complex_to_real(int nthreads, fftwnd_plan p,
						fftw_complex *in,
						fftw_real *out);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* RFFTW_THREADS_H */
