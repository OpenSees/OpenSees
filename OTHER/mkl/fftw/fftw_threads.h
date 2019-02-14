/*******************************************************************************
!   Copyright (c) 2003 Matteo Frigo
!   Copyright (c) 2003 Massachusetts Institute of Technology
!
!   This program is distributed with permission
!
!*******************************************************************************/

#ifndef FFTW_THREADS_H
#define FFTW_THREADS_H

#include "fftw.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/******************** User Interface *********************/

extern void fftw_threads(int nthreads,
		  fftw_plan plan, int howmany, fftw_complex *in, int istride,
		  int idist, fftw_complex *out, int ostride, int odist);
extern void fftwnd_threads(int nthreads,
			   fftwnd_plan plan, int howmany,
			   fftw_complex *in, int istride, int idist,
			   fftw_complex *out, int ostride, int odist);

extern void fftw_threads_one(int nthreads,
			     fftw_plan plan,
			     fftw_complex *in, fftw_complex *out);
extern void fftwnd_threads_one(int nthreads,
			       fftwnd_plan plan,
			       fftw_complex *in, fftw_complex *out);

extern int fftw_threads_init(void);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* FFTW_THREADS_H */
