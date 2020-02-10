#if !defined(_LABELIMPL_H)
#define _LABELIMPL_H

#include <petscdmlabel.h>
#include <petscbt.h>
#include <petscistypes.h>
#include <petsc/private/hashseti.h>

/* This is an integer map, in addition it is also a container class
   Design points:
     - Low storage is the most important design point
     - We want flexible insertion and deletion
     - We can live with O(log) query, but we need O(1) iteration over strata
*/
struct _n_DMLabel {
  PetscInt         refct;
  PetscObjectState state;
  char       *name;           /* Label name */
  PetscInt    numStrata;      /* Number of integer values */
  PetscInt    defaultValue;   /* Background value when no value explicitly given */
  PetscInt   *stratumValues;  /* Value of each stratum */
  /* Basic IS storage */
  PetscBool  *validIS;        /* The IS is valid (no additions need to be merged in) */
  PetscInt   *stratumSizes;   /* Size of each stratum */
  IS         *points;         /* Points for each stratum, always sorted */
  /* Hashtable for fast insertion */
  PetscHSetI *ht;             /* Hash table for fast insertion */
  /* Index for fast search */
  PetscInt    pStart, pEnd;   /* Bounds for index lookup */
  PetscBT     bt;             /* A bit-wise index */
};

PETSC_INTERN PetscErrorCode PetscSectionSymCreate_Label(PetscSectionSym);
#endif
