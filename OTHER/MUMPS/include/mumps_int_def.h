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
#if ! defined(MUMPS_INT_H)
#   define MUMPS_INT_H
/*  MUMPS has been configured without -DINTSIZE64:
 *  both 32-bit and 64-bit integers are used, depending on usage
 *  (e.g., the order of a matrix, N, is a 32-bit integer, and the
 *  number of nonzeros, NNZ, is a 64-bit integers) */
#   define MUMPS_INTSIZE32
#endif
