#if !defined(_PETSCISTYPES_H)
#define _PETSCISTYPES_H

/*S
     IS - Abstract PETSc object that allows indexing.

   Level: beginner

  Concepts: indexing, stride

.seealso:  ISCreateGeneral(), ISCreateBlock(), ISCreateStride(), ISGetIndices(), ISDestroy()
S*/
typedef struct _p_IS* IS;

/*S
   ISLocalToGlobalMapping - mappings from an arbitrary
      local ordering from 0 to n-1 to a global PETSc ordering
      used by a vector or matrix.

   Level: intermediate

   Note: mapping from Local to Global is scalable; but Global
  to Local may not be if the range of global values represented locally
  is very large.

   Note: the ISLocalToGlobalMapping is actually a private object; it is included
  here for the inline function ISLocalToGlobalMappingApply() to allow it to be inlined since
  it is used so often.

.seealso:  ISLocalToGlobalMappingCreate()
S*/
typedef struct _p_ISLocalToGlobalMapping* ISLocalToGlobalMapping;

/*S
     ISColoring - sets of IS's that define a coloring
              of the underlying indices

   Level: intermediate

    Notes:
        One should not access the *is records below directly because they may not yet
    have been created. One should use ISColoringGetIS() to make sure they are
    created when needed.

    Developer Note: this is not a PetscObject

.seealso:  ISColoringCreate(), ISColoringGetIS(), ISColoringView(), ISColoringGetIS()
S*/
typedef struct _n_ISColoring* ISColoring;

/*S
     PetscLayout - defines layout of vectors and matrices across processes (which rows are owned by which processes)

   Level: developer


.seealso:  PetscLayoutCreate(), PetscLayoutDestroy()
S*/
typedef struct _n_PetscLayout* PetscLayout;

/*S
  PetscSection - Mapping from integers in a designated range to contiguous sets of integers.

  In contrast to IS, which maps from integers to single integers, the range of a PetscSection is in the space of
  contiguous sets of integers. These ranges are frequently interpreted as domains of other array-like objects,
  especially other PetscSections, Vecs, and ISs. The domain is set with PetscSectionSetChart() and does not need to
  start at 0. For each point in the domain of a PetscSection, the output set is represented through an offset and a
  count, which are set using PetscSectionSetOffset() and PetscSectionSetDof() respectively. Lookup is typically using
  accessors or routines like VecGetValuesSection().

  Level: developer

.seealso:  PetscSectionCreate(), PetscSectionDestroy()
S*/
typedef struct _p_PetscSection *PetscSection;

/*S
  PetscSectionSym - Symmetries of the data referenced by a PetscSection.

  Often the order of data index by a PetscSection is meaningful, and describes additional structure, such as points on a
  line, grid, or lattice.  If the data is accessed from a different "orientation", then the image of the data under
  access then undergoes a symmetry transformation.  A PetscSectionSym specifies these symmetries.  The types of
  symmetries that can be specified are of the form R * P, where R is a diagonal matrix of scalars, and P is a permutation.

  Level: developer

.seealso: PetscSectionSymCreate(), PetscSectionSymDestroy(), PetscSectionSetSym(), PetscSectionGetSym(), PetscSectionSetFieldSym(), PetscSectionGetFieldSym(), PetscSectionGetSymPoints()
S*/
typedef struct _p_PetscSectionSym *PetscSectionSym;

#endif
