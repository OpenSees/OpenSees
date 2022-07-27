/*
 *
 *  This file is part of MUMPS 5.4.1, released
 *  on Tue Aug  3 09:49:43 UTC 2021
 *
 *
 *  Copyright 1991-2021 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
 *  Mumps Technologies, University of Bordeaux.
 *
 *  This version of MUMPS is provided to you free of charge. It is
 *  released under the CeCILL-C license 
 *  (see doc/CeCILL-C_V1-en.txt, doc/CeCILL-C_V1-fr.txt, and
 *  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html)
 *
 */
#ifndef MUMPS_SAVE_RESTORE_C_H
#define MUMPS_SAVE_RESTORE_C_H
#include "mumps_common.h"
#define MUMPS_GET_SAVE_DIR_C \
    F_SYMBOL(get_save_dir_c,GET_SAVE_DIR_C)
void MUMPS_CALL
MUMPS_GET_SAVE_DIR_C(MUMPS_INT *len_save_dir, char* save_dir, mumps_ftnlen l1);
#define MUMPS_GET_SAVE_PREFIX_C \
    F_SYMBOL(get_save_prefix_c,GET_SAVE_PREFIX_C)
void MUMPS_CALL
MUMPS_GET_SAVE_PREFIX_C(MUMPS_INT *len_save_prefix, char* save_prefix, mumps_ftnlen l1);
#define MUMPS_SAVE_RESTORE_RETURN_C \
    F_SYMBOL(save_restore_return_c,SAVE_RESTORE_RETURN_C)
void MUMPS_CALL
MUMPS_SAVE_RESTORE_RETURN_C();
#endif /* MUMPS_SAVE_RESTORE_C_H */
