#if !defined(__PETSCMOAB_H)
#define __PETSCMOAB_H

#include <petscvec.h>    /*I      "petscvec.h"    I*/
#include <petscmat.h>    /*I      "petscmat.h"    I*/
#include <petscdm.h>     /*I      "petscdm.h"    I*/
#include <petscdt.h>            /*I "petscdt.h" I*/

#include <string>
#include <moab/Core.hpp> /*I      "moab/Core.hpp"    I*/
#ifdef MOAB_HAVE_MPI
#include <moab/ParallelComm.hpp> /*I      "moab/ParallelComm.hpp"    I*/
#endif

/* The MBERR macro is used to save typing. It checks a MOAB error code
 * (rval) and calls SETERRQ if not MB_SUCCESS. A message (msg) can
 * also be passed in. */
#define MBERR(msg,rval) do{if(rval != moab::MB_SUCCESS) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_LIB,"MOAB ERROR (%i): %s",(PetscErrorCode)rval,msg);} while(0)
#define MBERRNM(rval) do{if(rval != moab::MB_SUCCESS) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_LIB,"MOAB ERROR (%i)",rval);} while(0)
#define MBERRV(mbif,rval) do{if(rval != moab::MB_SUCCESS) { std::string emsg; mbif->get_last_error(emsg); SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_LIB,"MOAB ERROR (%i): %s",(PetscErrorCode)rval,emsg.c_str());} } while(0)
#define MBERRVM(mbif,msg,rval) do{if(rval != moab::MB_SUCCESS) { std::string emsg; mbif->get_last_error(emsg); SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_LIB,"MOAB ERROR (%i): %s :: %s",(PetscErrorCode)rval,msg,emsg.c_str());} } while(0)


/* define enums for options to read and write MOAB files in parallel */
typedef enum {READ_PART, READ_DELETE, BCAST_DELETE} MoabReadMode;
static const char *const MoabReadModes[] = {"READ_PART", "READ_DELETE", "BCAST_DELETE", "MoabReadMode", "", 0};
typedef enum {WRITE_PART, FORMAT} MoabWriteMode;
static const char *const MoabWriteModes[] = {"WRITE_PART", "FORMAT", "MoabWriteMode", "", 0};

PETSC_EXTERN PetscErrorCode DMMoabCreate(MPI_Comm comm, DM *moab);
PETSC_EXTERN PetscErrorCode DMMoabCreateMoab(MPI_Comm comm, moab::Interface *mbiface, moab::Tag *ltog_tag, moab::Range *range, DM *moab);
PETSC_EXTERN PetscErrorCode DMMoabOutput(DM, const char*, const char*);

PETSC_EXTERN PetscErrorCode DMMoabSetInterface(DM, moab::Interface *);
PETSC_EXTERN PetscErrorCode DMMoabGetInterface(DM, moab::Interface **);
#ifdef MOAB_HAVE_MPI
PETSC_EXTERN PetscErrorCode DMMoabGetParallelComm(DM, moab::ParallelComm**);
#endif

PETSC_EXTERN PetscErrorCode DMMoabSetLocalVertices(DM, moab::Range *);
PETSC_EXTERN PetscErrorCode DMMoabGetAllVertices(DM, moab::Range *local);
PETSC_EXTERN PetscErrorCode DMMoabGetLocalVertices(DM, const moab::Range**, const moab::Range**);
PETSC_EXTERN PetscErrorCode DMMoabSetLocalElements(DM, moab::Range*);
PETSC_EXTERN PetscErrorCode DMMoabGetLocalElements(DM, const moab::Range**);
PETSC_EXTERN PetscErrorCode DMMoabSetLocalToGlobalTag(DM, moab::Tag);
PETSC_EXTERN PetscErrorCode DMMoabGetLocalToGlobalTag(DM, moab::Tag*);
PETSC_EXTERN PetscErrorCode DMMoabSetBlockSize(DM, PetscInt bs);
PETSC_EXTERN PetscErrorCode DMMoabGetBlockSize(DM, PetscInt *bs);
PETSC_EXTERN PetscErrorCode DMMoabSetBlockFills(DM, const PetscInt*, const PetscInt*);
PETSC_EXTERN PetscErrorCode DMMoabGetHierarchyLevel(DM, PetscInt *);

PETSC_EXTERN PetscErrorCode DMMoabGetDimension(DM dm, PetscInt *dim);
PETSC_EXTERN PetscErrorCode DMMoabGetBoundaryEntities(DM dm, moab::Range*, moab::Range*, moab::Range*);
PETSC_EXTERN PetscErrorCode DMMoabGetMaterialBlock(DM dm, const moab::EntityHandle, PetscInt*);

PETSC_EXTERN PetscErrorCode DMMoabGetSize(DM dm, PetscInt*, PetscInt*);
PETSC_EXTERN PetscErrorCode DMMoabGetLocalSize(DM dm, PetscInt*, PetscInt*, PetscInt*, PetscInt*);
PETSC_EXTERN PetscErrorCode DMMoabGetOffset(DM dm, PetscInt*);

PETSC_EXTERN PetscErrorCode DMMoabVecGetArrayRead(DM, Vec, void*);
PETSC_EXTERN PetscErrorCode DMMoabVecRestoreArrayRead(DM, Vec, void*);
PETSC_EXTERN PetscErrorCode DMMoabVecGetArray(DM, Vec, void*);
PETSC_EXTERN PetscErrorCode DMMoabVecRestoreArray(DM, Vec, void*);

PETSC_EXTERN PetscErrorCode DMMoabCreateVector(DM dm, moab::Tag tag, const moab::Range *range, PetscBool serial, PetscBool destroy_tag, Vec *X);
PETSC_EXTERN PetscErrorCode DMMoabGetVecTag(Vec vec, moab::Tag *tag);
PETSC_EXTERN PetscErrorCode DMMoabGetVecRange(Vec vec, moab::Range *range);

PETSC_EXTERN PetscErrorCode DMMoabSetFieldVector(DM, PetscInt, Vec);
PETSC_EXTERN PetscErrorCode DMMoabSetGlobalFieldVector(DM, Vec);

PETSC_EXTERN PetscErrorCode DMMoabCreateVertices(DM, const PetscReal* coords, PetscInt nverts, moab::Range*);
PETSC_EXTERN PetscErrorCode DMMoabCreateElement(DM, const moab::EntityType type, const moab::EntityHandle* conn, PetscInt nverts, moab::EntityHandle* elem);
PETSC_EXTERN PetscErrorCode DMMoabCreateSubmesh(DM dm, DM *newdm);
PETSC_EXTERN PetscErrorCode DMMoabRenumberMeshEntities(DM dm);

PETSC_EXTERN PetscErrorCode DMMoabGetFieldName(DM dm, PetscInt field, const char **fieldName);
PETSC_EXTERN PetscErrorCode DMMoabSetFieldName(DM dm, PetscInt field, const char *fieldName);
PETSC_EXTERN PetscErrorCode DMMoabSetFieldNames(DM dm, PetscInt nfields, const char* fields[]);
PETSC_EXTERN PetscErrorCode DMMoabGetFieldDof(DM dm, moab::EntityHandle point, PetscInt field, PetscInt* dof);
PETSC_EXTERN PetscErrorCode DMMoabGetFieldDofs(DM dm, PetscInt npoints, const moab::EntityHandle* points, PetscInt field, PetscInt* dof);
PETSC_EXTERN PetscErrorCode DMMoabGetFieldDofsLocal(DM dm, PetscInt npoints, const moab::EntityHandle* points, PetscInt field, PetscInt* dof);
PETSC_EXTERN PetscErrorCode DMMoabGetDofs(DM dm, PetscInt npoints, const moab::EntityHandle* points, PetscInt* dof);
PETSC_EXTERN PetscErrorCode DMMoabGetDofsLocal(DM dm, PetscInt npoints, const moab::EntityHandle* points, PetscInt* dof);
PETSC_EXTERN PetscErrorCode DMMoabGetDofsBlocked(DM dm, PetscInt npoints, const moab::EntityHandle* points, PetscInt* dof);
PETSC_EXTERN PetscErrorCode DMMoabGetDofsBlockedLocal(DM dm, PetscInt npoints, const moab::EntityHandle* points, PetscInt* dof);

PETSC_EXTERN PetscErrorCode DMMoabGetVertexDofsBlocked(DM dm, PetscInt** dof);
PETSC_EXTERN PetscErrorCode DMMoabGetVertexDofsBlockedLocal(DM dm, PetscInt** dof);

/* discretization and assembly specific DMMoab interface functions */
PETSC_EXTERN PetscErrorCode DMMoabGetElementConnectivity(DM dm, moab::EntityHandle ehandle, PetscInt* nconn, const moab::EntityHandle **conn);
PETSC_EXTERN PetscErrorCode DMMoabGetVertexConnectivity(DM dm, moab::EntityHandle ehandle, PetscInt* nconn, moab::EntityHandle **conn);
PETSC_EXTERN PetscErrorCode DMMoabRestoreVertexConnectivity(DM dm, moab::EntityHandle ehandle, PetscInt* nconn, moab::EntityHandle **conn);
PETSC_EXTERN PetscErrorCode DMMoabGetVertexCoordinates(DM dm, PetscInt nconn, const moab::EntityHandle *conn, PetscReal *vpos);
PETSC_EXTERN PetscErrorCode DMMoabIsEntityOnBoundary(DM dm, const moab::EntityHandle ent, PetscBool* ent_on_boundary);
PETSC_EXTERN PetscErrorCode DMMoabCheckBoundaryVertices(DM dm, PetscInt nconn, const moab::EntityHandle *cnt, PetscBool* isbdvtx);
PETSC_EXTERN PetscErrorCode DMMoabGetBoundaryMarkers(DM dm, const moab::Range **bdvtx, const moab::Range** bdelems, const moab::Range** bdfaces);

/* TODO: Replace nverts/coords with just moab::EntityHandle -- can also eliminate dim */
/* TODO: Replace quad/npts with PetscDT */
PETSC_EXTERN PetscErrorCode DMMoabFEMCreateQuadratureDefault ( const PetscInt dim, const PetscInt nverts, PetscQuadrature *quadrature );
PETSC_EXTERN PetscErrorCode DMMoabFEMComputeBasis ( const PetscInt dim, const PetscInt nverts, const PetscReal *coordinates, const PetscQuadrature quadrature, PetscReal *phypts, PetscReal *jxw, PetscReal *phi, PetscReal **dphi);
PETSC_EXTERN PetscErrorCode DMMoabPToRMapping ( const PetscInt dim, const PetscInt nverts, const PetscReal *coordinates, const PetscReal* xphy, PetscReal* natparam, PetscReal* phi);

/* DM utility creation interface */
PETSC_EXTERN PetscErrorCode DMMoabCreateBoxMesh(MPI_Comm, PetscInt, PetscBool, const PetscReal*, PetscInt, PetscInt, DM*);
PETSC_EXTERN PetscErrorCode DMMoabLoadFromFile(MPI_Comm, PetscInt, PetscInt, const char*, const char*, DM*);

/* Uniform refinement hierarchy interface */
PETSC_EXTERN PetscErrorCode DMMoabGenerateHierarchy(DM dm, PetscInt nlevels, PetscInt *ldegrees);

#endif
