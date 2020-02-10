
#if !defined(__PETSCVIEWERHDF5_H)
#define __PETSCVIEWERHDF5_H

#include <petscviewer.h>

#if defined(PETSC_HAVE_HDF5)

#include <hdf5.h>
PETSC_EXTERN PetscErrorCode PetscViewerHDF5GetFileId(PetscViewer,hid_t*);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5OpenGroup(PetscViewer, hid_t *, hid_t *);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5ReadSizes(PetscViewer, const char[], PetscInt *, PetscInt *);

/* On 32 bit systems HDF5 is limited by size of integer, because hsize_t is defined as size_t */
#define PETSC_HDF5_INT_MAX  2147483647
#define PETSC_HDF5_INT_MIN -2147483647

PETSC_STATIC_INLINE PetscErrorCode PetscHDF5IntCast(PetscInt a,hsize_t *b)
{
  PetscFunctionBegin;
#if defined(PETSC_USE_64BIT_INDICES) && (PETSC_SIZEOF_SIZE_T == 4)
  if ((a) > PETSC_HDF5_INT_MAX) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Array too long for HDF5");
#endif
  *b =  (hsize_t)(a);
  PetscFunctionReturn(0);
}
PETSC_EXTERN PetscErrorCode PetscDataTypeToHDF5DataType(PetscDataType,hid_t*);
PETSC_EXTERN PetscErrorCode PetscHDF5DataTypeToPetscDataType(hid_t,PetscDataType*);

#define PetscStackCallHDF5(func,args) do {                        \
    herr_t _status;                                               \
    PetscStackPush(#func);_status = func args;PetscStackPop; if (_status) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_LIB,"Error in HDF5 call %s() Status %d",#func,(int)_status); \
  } while (0)

#define PetscStackCallHDF5Return(ret,func,args) do {              \
    PetscStackPush(#func);ret = func args;PetscStackPop; if (ret < 0) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_LIB,"Error in HDF5 call %s() Status %d",#func,(int)ret); \
  } while (0)

#endif  /* defined(PETSC_HAVE_HDF5) */

PETSC_EXTERN PetscErrorCode PetscViewerHDF5ReadAttribute(PetscViewer,const char[],const char[],PetscDataType,void*);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5WriteAttribute(PetscViewer,const char[],const char[],PetscDataType,const void*);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5HasAttribute(PetscViewer, const char[], const char[], PetscBool *);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5WriteSDS(PetscViewer,float *,int,int *,int);

PETSC_EXTERN PetscErrorCode PetscViewerHDF5Open(MPI_Comm,const char[],PetscFileMode,PetscViewer*);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5PushGroup(PetscViewer,const char *);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5PopGroup(PetscViewer);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5GetGroup(PetscViewer, const char **);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5IncrementTimestep(PetscViewer);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5SetTimestep(PetscViewer,PetscInt);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5GetTimestep(PetscViewer,PetscInt*);

PETSC_EXTERN PetscErrorCode PetscViewerHDF5SetBaseDimension2(PetscViewer,PetscBool);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5GetBaseDimension2(PetscViewer,PetscBool*);

PETSC_EXTERN PetscErrorCode PetscViewerHDF5SetSPOutput(PetscViewer,PetscBool);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5GetSPOutput(PetscViewer,PetscBool*);

#endif
