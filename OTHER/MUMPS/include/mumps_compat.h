/*
 *
 *  This file is part of MUMPS 5.5.1, released
 *  on Tue Jul 12 13:17:24 UTC 2022
 *
 *
 *  Copyright 1991-2022 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
 *  Mumps Technologies, University of Bordeaux.
 *
 *  This version of MUMPS is provided to you free of charge. It is
 *  released under the CeCILL-C license 
 *  (see doc/CeCILL-C_V1-en.txt, doc/CeCILL-C_V1-fr.txt, and
 *  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html)
 *
 */

/* Compatibility issues between various Windows versions */
#ifndef MUMPS_COMPAT_H
#define MUMPS_COMPAT_H


#if defined(_WIN32) && ! defined(__MINGW32__)
# define MUMPS_WIN32 1
#endif

#ifndef MUMPS_CALL
# ifdef MUMPS_WIN32
/* Choose between next lines or modify according
 * to your Windows calling conventions:       */
/*   #define MUMPS_CALL                       */
/*   #define MUMPS_CALL __stdcall             */
/*   #define MUMPS_CALL __declspec(dllexport) */
#  define MUMPS_CALL
# else
#  define MUMPS_CALL
# endif
#endif

#if (__STDC_VERSION__ >= 199901L)
# define MUMPS_INLINE static inline
#else
# define MUMPS_INLINE
#endif


#endif /* MUMPS_COMPAT_H */
