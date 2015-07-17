#include <elementAPI.h>
#include <OPS_Stream.h>

//#include <Domain.h>
//#include <TclModelBuilder.h>

#include <windows.h>

#include <SimulationInformation.h>
//SimulationInformation simulationInfo;

#define DllExport _declspec(dllexport)

BOOL APIENTRY DllMain(HANDLE hModule, 
                      DWORD  ul_reason_for_call, 
                      LPVOID lpReserved)
{
    return TRUE;
}

OPS_Stream *opserrPtr = 0;
SimulationInformation *theSimulationInfo = 0;
//Domain *ops_TheActiveDomain = 0;

//double ops_Dt = 0;

typedef int (*OPS_ErrorPtrType)(char *, int);
typedef int (*OPS_GetNumRemainingInputArgsType)();
typedef int (*OPS_GetIntInputPtrType)(int *, int *);
typedef int (*OPS_GetDoubleInputPtrType)(int *, double *);
typedef int (*OPS_GetStringType)(char *, int);
typedef int (*OPS_GetStringCopyType)(char **);
typedef int (*OPS_AllocateElementPtrType)(eleObj *, int *matTags, int *maType);
typedef int (*OPS_AllocateMaterialPtrType)(matObj *);
typedef UniaxialMaterial *(*OPS_GetUniaxialMaterialPtrType)(int matTag);
typedef NDMaterial *(*OPS_GetNDMaterialPtrType)(int matTag);
typedef SectionForceDeformation *(*OPS_GetSectionForceDeformationPtrType)(int matTag);
typedef CrdTransf *(*OPS_GetCrdTransfPtrType)(int tag);
typedef int (*OPS_GetNodeInfoPtrType)(int *, int *, double *);
typedef int (*OPS_InvokeMaterialDirectlyPtrType)(matObject **, modelState *, double *, double *, double *, int *);
typedef int (*OPS_GetIntPtrType)();
typedef FE_Datastore *(*OPS_GetFEDatastorePtrType)();
typedef const char *(_cdecl *OPS_GetInterpPWD_PtrType)();


//int    OPS_InvokeMaterial(struct eleObj *, int *,modelState *, double *, double *, double *, int *);

OPS_ErrorPtrType OPS_ErrorPtr = 0;
OPS_GetIntInputPtrType OPS_GetIntInputPtr = 0;
OPS_GetDoubleInputPtrType OPS_GetDoubleInputPtr = 0;
OPS_AllocateElementPtrType OPS_AllocateElementPtr = 0;
OPS_AllocateMaterialPtrType OPS_AllocateMaterialPtr = 0;
OPS_GetUniaxialMaterialPtrType OPS_GetUniaxialMaterialPtr = 0;
OPS_GetNDMaterialPtrType OPS_GetNDMaterialPtr = 0;
OPS_GetSectionForceDeformationPtrType OPS_GetSectionForceDeformationPtr = 0;
OPS_GetCrdTransfPtrType OPS_GetCrdTransfPtrFunc = 0;
OPS_GetNodeInfoPtrType OPS_GetNodeCrdPtr = 0;
OPS_GetNodeInfoPtrType OPS_GetNodeDispPtr = 0;
OPS_GetNodeInfoPtrType OPS_GetNodeVelPtr = 0;
OPS_GetNodeInfoPtrType OPS_GetNodeAccelPtr = 0;
OPS_GetNodeInfoPtrType OPS_GetNodeIncrDispPtr = 0;
OPS_GetNodeInfoPtrType OPS_GetNodeIncrDeltaDispPtr = 0;
OPS_InvokeMaterialDirectlyPtrType OPS_InvokeMaterialDirectlyPtr = 0;
OPS_GetNumRemainingInputArgsType OPS_GetNumRemainingInputArgsPtr = 0;
OPS_GetStringType OPS_GetStringPtr = 0;
OPS_GetStringCopyType OPS_GetStringCopyPtr = 0;
OPS_GetIntPtrType OPS_GetNDM_Ptr = 0;
OPS_GetIntPtrType OPS_GetNDF_Ptr = 0;
OPS_GetFEDatastorePtrType OPS_GetFEDatastorePtr = 0;
OPS_GetInterpPWD_PtrType OPS_GetInterpPWD_Ptr = 0;


extern "C" DllExport
void setGlobalPointers(OPS_Stream *theErrorStreamPtr,
		       Domain *theDomain,
                       SimulationInformation *theSimulationInfoPtr,
                       OPS_ErrorPtrType errorFunct,
                       OPS_GetIntInputPtrType getIntInputFunct,
                       OPS_GetDoubleInputPtrType getDoubleInputFunct,
                       OPS_AllocateElementPtrType allocateElementFunct,
                       OPS_AllocateMaterialPtrType allocateMaterialFunct,
                       OPS_GetUniaxialMaterialPtrType OPS_GetUniaxialMaterialFunct,
                       OPS_GetNDMaterialPtrType OPS_GetNDMaterialFunct,
		       OPS_GetSectionForceDeformationPtrType OPS_GetSectionForceDeformationFunct,
                       OPS_InvokeMaterialDirectlyPtrType OPS_InvokeMaterialDirectlyFunct,
                       OPS_GetNodeInfoPtrType OPS_GetNodeCrdFunct,
                       OPS_GetNodeInfoPtrType OPS_GetNodeDispFunct,
                       OPS_GetNodeInfoPtrType OPS_GetNodeVelFunct,
                       OPS_GetNodeInfoPtrType OPS_GetNodeAccelFunct,
                       OPS_GetNodeInfoPtrType OPS_GetNodeIncrDispFunct,
                       OPS_GetNodeInfoPtrType OPS_GetNodeIncrDeltaDispFunct,
                       OPS_GetNumRemainingInputArgsType OPS_GetNumRemainingArgsFunct,
                       OPS_GetStringType OPS_GetStringFunct,
                       OPS_GetStringCopyType OPS_GetStringCopyFunct,
                       OPS_GetCrdTransfPtrType OPS_GetCrdTransfFunct,
                       OPS_GetIntPtrType OPS_GetNDM_Funct,
                       OPS_GetIntPtrType OPS_GetNDF_Funct,
                       OPS_GetFEDatastorePtrType OPS_GetFEDatastoreFunct,
                       OPS_GetInterpPWD_PtrType OPS_GetInterpPWD_Funct)
{
    opserrPtr = theErrorStreamPtr;
    ops_TheActiveDomain = theDomain;
    theSimulationInfo = theSimulationInfoPtr;
    OPS_ErrorPtr = errorFunct;
    OPS_GetIntInputPtr = getIntInputFunct;
    OPS_GetDoubleInputPtr = getDoubleInputFunct;
    OPS_AllocateElementPtr = allocateElementFunct;
    OPS_AllocateMaterialPtr = allocateMaterialFunct;
    OPS_GetUniaxialMaterialPtr = OPS_GetUniaxialMaterialFunct;
    OPS_GetNDMaterialPtr = OPS_GetNDMaterialFunct;
	OPS_GetSectionForceDeformationPtr = OPS_GetSectionForceDeformationFunct;
    OPS_GetNodeCrdPtr = OPS_GetNodeCrdFunct;
    OPS_GetNodeDispPtr = OPS_GetNodeDispFunct;
    OPS_GetNodeVelPtr = OPS_GetNodeVelFunct;
    OPS_GetNodeAccelPtr = OPS_GetNodeAccelFunct;
    OPS_GetNodeIncrDispPtr = OPS_GetNodeIncrDispFunct;
    OPS_GetNodeIncrDeltaDispPtr = OPS_GetNodeIncrDeltaDispFunct;
    OPS_InvokeMaterialDirectlyPtr = OPS_InvokeMaterialDirectlyFunct;
    OPS_GetNumRemainingInputArgsPtr = OPS_GetNumRemainingArgsFunct;
    OPS_GetStringPtr = OPS_GetStringFunct;
    OPS_GetStringCopyPtr = OPS_GetStringCopyFunct;
    OPS_GetCrdTransfPtrFunc = OPS_GetCrdTransfFunct;
    OPS_GetNDM_Ptr = OPS_GetNDM_Funct;
    OPS_GetNDF_Ptr = OPS_GetNDF_Funct;
    OPS_GetFEDatastorePtr = OPS_GetFEDatastoreFunct;
    OPS_GetInterpPWD_Ptr = OPS_GetInterpPWD_Funct;
}


UniaxialMaterial *
OPS_GetUniaxialMaterial(int matTag)
{
    return (*OPS_GetUniaxialMaterialPtr)(matTag);
}

NDMaterial *
OPS_GetNDMaterial(int matTag)
{
    return (*OPS_GetNDMaterialPtr)(matTag);
}

SectionForceDeformation *
OPS_GetSectionForceDeformation(int matTag)
{
return (*OPS_GetSectionForceDeformationPtr)(matTag);
}

CrdTransf *
OPS_GetCrdTransfPtr(int tag)
{
    return (*OPS_GetCrdTransfPtrFunc)(tag);
}

int OPS_Error(char *data, int length)
{
    return (*OPS_ErrorPtr)(data, length);
}

extern "C" int OPS_GetIntInput(int *numData, int*data)
{
    return (*OPS_GetIntInputPtr)(numData, data);
}

extern "C" int OPS_GetDoubleInput(int *numData, double *data)
{
    return (*OPS_GetDoubleInputPtr)(numData, data);
}

extern "C" int OPS_AllocateMaterial(matObj *mat)
{
    return (*OPS_AllocateMaterialPtr)(mat);
}

extern "C" int OPS_AllocateElement(eleObj *ele, int *matTags, int *matType)
{
    return (*OPS_AllocateElementPtr)(ele, matTags, matType);
}

extern "C" int OPS_GetNodeCrd(int *nodeTag, int *sizeData, double *data)
{
    return (*OPS_GetNodeCrdPtr)(nodeTag, sizeData, data);
}

extern "C" int OPS_GetNodeDisp(int *nodeTag, int *sizeData, double *data)
{
    return (*OPS_GetNodeDispPtr)(nodeTag, sizeData, data);
}

extern "C" int OPS_GetNodeVel(int *nodeTag, int *sizeData, double *data)
{
    return (*OPS_GetNodeVelPtr)(nodeTag, sizeData, data);
}

extern "C" int OPS_GetNodeAccel(int *nodeTag, int *sizeData, double *data)
{
    return (*OPS_GetNodeAccelPtr)(nodeTag, sizeData, data);
}

extern "C" int OPS_GetNodeIncrDisp(int *nodeTag, int *sizeData, double *data)
{
    return (*OPS_GetNodeIncrDispPtr)(nodeTag, sizeData, data);
}

extern "C" int OPS_GetNodeIncrDeltaDisp(int *nodeTag, int *sizeData, double *data)
{
    return (*OPS_GetNodeIncrDeltaDispPtr)(nodeTag, sizeData, data);
}

extern "C" int OPS_InvokeMaterialDirectly(matObject **theMat, modelState *model,
                                          double *strain, double *stress, double *tang, int *isw)
{
    return (*OPS_InvokeMaterialDirectlyPtr)(theMat, model, strain, stress, tang, isw);
}

extern "C" int OPS_GetString(char *cArray, int sizeArray)
{
    return (*OPS_GetStringPtr)(cArray, sizeArray);  
}

extern "C" int OPS_GetStringCopy(char **cArray)
{
    return (*OPS_GetStringCopyPtr)(cArray);  
}

extern "C" int OPS_GetNumRemainingInputArgs()
{
    return (*OPS_GetNumRemainingInputArgsPtr)();  
}

extern "C" int OPS_GetNDM()
{
    return (*OPS_GetNDM_Ptr)();
}

extern "C" int OPS_GetNDF()
{
    return (*OPS_GetNDF_Ptr)();
}

FE_Datastore *
OPS_GetFEDatastore()
{
    return (*OPS_GetFEDatastorePtr)();
}

extern "C" const char *OPS_GetInterpPWD() 
{
    return (*OPS_GetInterpPWD_Ptr)();
}
