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
#include "mumps_common.h"
/* Special case of mapping and pivnul_list -- allocated from MUMPS */
static MUMPS_INT * MUMPS_MAPPING;
static MUMPS_INT * MUMPS_PIVNUL_LIST;
/* as uns_perm and sym_perm */
static MUMPS_INT * MUMPS_SYM_PERM;
static MUMPS_INT * MUMPS_UNS_PERM;
MUMPS_INT*
mumps_get_mapping()
{
    return MUMPS_MAPPING;
}
void MUMPS_CALL
MUMPS_ASSIGN_MAPPING(MUMPS_INT * f77mapping)
{
    MUMPS_MAPPING = f77mapping;
}
void MUMPS_CALL
MUMPS_NULLIFY_C_MAPPING()
{
    MUMPS_MAPPING = 0;
}
MUMPS_INT*
mumps_get_pivnul_list()
{
    return MUMPS_PIVNUL_LIST;
}
void MUMPS_CALL
MUMPS_ASSIGN_PIVNUL_LIST(MUMPS_INT * f77pivnul_list)
{
    MUMPS_PIVNUL_LIST = f77pivnul_list;
}
void MUMPS_CALL
MUMPS_NULLIFY_C_PIVNUL_LIST()
{
    MUMPS_PIVNUL_LIST = 0;
}
MUMPS_INT*
mumps_get_sym_perm()
{
    return MUMPS_SYM_PERM;
}
void MUMPS_CALL
MUMPS_ASSIGN_SYM_PERM(MUMPS_INT * f77sym_perm)
{
    MUMPS_SYM_PERM = f77sym_perm;
}
void MUMPS_CALL
MUMPS_NULLIFY_C_SYM_PERM()
{
    MUMPS_SYM_PERM = 0;
}
MUMPS_INT*
mumps_get_uns_perm()
{
    return MUMPS_UNS_PERM;
}
void MUMPS_CALL
MUMPS_ASSIGN_UNS_PERM(MUMPS_INT * f77uns_perm)
{
    MUMPS_UNS_PERM = f77uns_perm;
}
void MUMPS_CALL
MUMPS_NULLIFY_C_UNS_PERM()
{
    MUMPS_UNS_PERM = 0;
}
void MUMPS_CALL
MUMPS_ICOPY_32TO64_64C_IP_C(MUMPS_INT *inouttab, MUMPS_INT8 *sizetab)
/* Copies in-place *sizetab int values starting at address inouttab
   into *sizetab int64_t values starting at the same address.
*/
{
   MUMPS_INT8 i8; /* signed integer needed for reversed loop below */
   for (i8=*sizetab-1; i8 >=0; i8--)
     {
       /* outtab8[i8]=(MUMPS_INT8)intab4[i8]; */
       ((MUMPS_INT8 *)inouttab)[i8]=(MUMPS_INT8)inouttab[i8];
     }
}
void MUMPS_CALL
MUMPS_ICOPY_64TO32_64C_IP_C(MUMPS_INT8 *inouttab, MUMPS_INT8 *sizetab)
/* Copies in-place *sizetab int64_t values starting at address inouttab
   into *sizetab int values starting at the same address */
{
   MUMPS_INT8 i8;
   for (i8=0; i8 < *sizetab; i8++)
     {
       /*       outtab4[i8]=(MUMPS_INT)intab8[i8]; */
       ((MUMPS_INT *)inouttab)[i8]=(MUMPS_INT)inouttab[i8];
     }
}
