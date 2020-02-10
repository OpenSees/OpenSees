#if !defined(_DMMBIMPL_H)
#define _DMMBIMPL_H

#include <petscdmmoab.h>    /*I      "petscdmmoab.h"    I*/
#include "petsc/private/dmimpl.h"

/* This is an integer map, in addition it is also a container class
   Design points:
     - Low storage is the most important design point
     - We want flexible insertion and deletion
     - We can live with O(log) query, but we need O(1) iteration over strata
*/
typedef struct {
  moab::Interface    *mbiface;
#ifdef MOAB_HAVE_MPI
  moab::ParallelComm *pcomm;
#endif
  moab::Range        *tag_range; /* entities to which this tag applies */
  moab::Tag           tag;
  PetscInt            tag_size;
  PetscBool           new_tag;
  PetscBool           is_global_vec;
  PetscBool           is_native_vec;
  Vec                 local;
} Vec_MOAB;

namespace moab {
class NestedRefine;
}

typedef struct {
  /* common data */
  PetscInt                dim;                            /* Current topological dimension handled by DMMoab */
  PetscInt                bs;                             /* Block size that controls the strided vs interlaced configuration in discrete systems -
                                                             This affects the layout and hence the degree-of-freedom of the different fields (components) */

  PetscInt               *dfill, *ofill;                  /* The diagonal and off-diagonal block-fill to indicate coupling between components */
  PetscInt               *materials;                      /* The array that caches the material data for each element */

  PetscInt                numFields;
  const char              **fieldNames;

  /* level specific data */
  PetscInt                n, nloc, nghost;                /* Number of global, local only and shared vertices for current partition */
  PetscInt                nele, neleloc, neleghost;       /* Number of global, local only and shared elements for current partition */
  PetscInt                *gsindices;                     /* Global ID for all local+ghosted vertices */
  PetscInt                *gidmap, *lidmap;               /* Global ID indices, Local ID indices, field-based local map, field-based global map */
  PetscInt                seqstart, seqend;               /* Local start and end entity IDs for vertices */
  PetscInt                vstart, vend;                   /* Global start and end index for distributed Vec */
  PetscInt                nghostrings;                    /* Number of ghost ring layers */
  PetscInt                gminmax[2], lminmax[2];         /* Local and global min/max in the ID sequence */
  PetscInt                refct;

  /* store the mapping information */
  ISLocalToGlobalMapping  ltog_map;
  VecScatter              ltog_sendrecv;

  /* MOAB objects cached internally in DMMoab */

  /* common data */
  moab::Interface         *mbiface;                       /* MOAB Interface/Core reference */
#ifdef MOAB_HAVE_MPI
  moab::ParallelComm      *pcomm;                         /* MOAB ParallelComm reference */
#endif
  moab::Tag               ltog_tag;                       /* MOAB supports "global id" tags */
  moab::Tag               material_tag;                   /* MOAB supports "material_set" tags */
  PetscBool               icreatedinstance;               /* true if DM created moab instance internally, will destroy instance in DMDestroy */

  /* store options to customize DMMoab I/O */
  PetscInt                rw_dbglevel;
  PetscBool               partition_by_rank;
  char                    extra_read_options[PETSC_MAX_PATH_LEN];
  char                    extra_write_options[PETSC_MAX_PATH_LEN];
  MoabReadMode            read_mode;
  MoabWriteMode           write_mode;

  /* level specific data */
  moab::Range             *vowned, *vghost, *vlocal;      /* Vertex entities: strictly owned, strictly ghosted, owned+ghosted */
  moab::Range             *elocal, *eghost;               /* Topological dimensional entities: strictly owned, strictly ghosted */
  moab::Range             *bndyvtx, *bndyfaces, *bndyelems; /* Boundary entities: skin vertices, skin faces and elements on the outer skin */
  moab::EntityHandle      fileset;                        /* The Global set to which all local entities belong */

  /* level hierarchy in MOAB */
  moab::NestedRefine     *hierarchy;
  PetscInt                nhlevels, hlevel;
  moab::EntityHandle     *hsets;

  /* Sub-mesh level data-strucuture */
  DM                     *parent;

} DM_Moab;

typedef struct {
  DM_Moab             *pdmmoab;
  moab::NestedRefine  *hierarchy;
  PetscInt             nhlevels, hlevel;
  moab::EntityHandle  *hsets;
} SubDM_MOAB;

#endif /* _DMMBIMPL_H */

