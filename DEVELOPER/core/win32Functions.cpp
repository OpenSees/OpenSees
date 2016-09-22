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
typedef int (*OPS_ResetCurrentInputArgType)(int);
typedef int (*OPS_GetIntInputPtrType)(int *, int *);
typedef int (*OPS_GetDoubleInputPtrType)(int *, double *);
typedef const char *(*OPS_GetStringType)();
typedef int (*OPS_GetStringCopyType)(char **);
typedef int (*OPS_AllocateElementPtrType)(eleObj *, int *matTags, int *maType);
typedef int (*OPS_AllocateMaterialPtrType)(matObj *);
typedef UniaxialMaterial *(*OPS_GetUniaxialMaterialPtrType)(int matTag);
typedef NDMaterial *(*OPS_GetNDMaterialPtrType)(int matTag);
typedef SectionForceDeformation *(*OPS_GetSectionForceDeformationPtrType)(int secTag);
typedef CrdTransf *(*OPS_GetCrdTransfPtrType)(int crdTag);
typedef FrictionModel *(*OPS_GetFrictionModelPtrType)(int frnTag);
typedef int (*OPS_GetNodeInfoPtrType)(int *, int *, double *);
typedef int (*OPS_InvokeMaterialDirectlyPtrType)(matObject **, modelState *, double *, double *, double *, int *);
typedef int (*OPS_GetIntPtrType)();
typedef FE_Datastore *(*OPS_GetFEDatastorePtrType)();
typedef const char *(_cdecl *OPS_GetInterpPWD_PtrType)();
typedef Domain *(*OPS_GetDomainPointer)(void);
typedef AnalysisModel **(*OPS_GetAnalysisModelPtrType)(void);
typedef EquiSolnAlgo **(*OPS_GetAlgorithmPtrType)(void);
typedef ConstraintHandler **(*OPS_GetHandlerPtrType)(void);
typedef DOF_Numberer **(*OPS_GetNumbererPtrType)(void);
typedef LinearSOE **(*OPS_GetSOEPtrType)(void);
typedef EigenSOE **(*OPS_GetEigenSOEPtrType)(void);
typedef StaticAnalysis **(*OPS_GetStaticAnalysisPtrType)(void);
typedef DirectIntegrationAnalysis **(*OPS_GetTransientAnalysisPtrType)(void);
typedef VariableTimeStepDirectIntegrationAnalysis **(*OPS_GetVariableTimeStepTransientAnalysisPtrType)(void);
typedef int *(*OPS_GetNumEigenPtrType)(void);
typedef StaticIntegrator **(*OPS_GetStaticIntegratorPtrType)(void);
typedef TransientIntegrator **(*OPS_GetTransientIntegratorPtrType)(void);
typedef ConvergenceTest **(*OPS_GetTestPtrType)(void);
typedef bool *(*OPS_builtModelPtrType)(void);

//int    OPS_InvokeMaterial(struct eleObj *, int *,modelState *, double *, double *, double *, int *);

OPS_ErrorPtrType OPS_ErrorPtr = 0;
OPS_GetIntInputPtrType OPS_GetIntInputPtr = 0;
OPS_GetDoubleInputPtrType OPS_GetDoubleInputPtr = 0;
OPS_AllocateElementPtrType OPS_AllocateElementPtr = 0;
OPS_AllocateMaterialPtrType OPS_AllocateMaterialPtr = 0;
OPS_GetUniaxialMaterialPtrType OPS_GetUniaxialMaterialPtr = 0;
OPS_GetNDMaterialPtrType OPS_GetNDMaterialPtr = 0;
OPS_GetSectionForceDeformationPtrType OPS_GetSectionForceDeformationPtr = 0;
OPS_GetCrdTransfPtrType OPS_GetCrdTransfPtr = 0;
OPS_GetFrictionModelPtrType OPS_GetFrictionModelPtr = 0;
OPS_GetNodeInfoPtrType OPS_GetNodeCrdPtr = 0;
OPS_GetNodeInfoPtrType OPS_GetNodeDispPtr = 0;
OPS_GetNodeInfoPtrType OPS_GetNodeVelPtr = 0;
OPS_GetNodeInfoPtrType OPS_GetNodeAccelPtr = 0;
OPS_GetNodeInfoPtrType OPS_GetNodeIncrDispPtr = 0;
OPS_GetNodeInfoPtrType OPS_GetNodeIncrDeltaDispPtr = 0;
OPS_InvokeMaterialDirectlyPtrType OPS_InvokeMaterialDirectlyPtr = 0;
OPS_GetNumRemainingInputArgsType OPS_GetNumRemainingInputArgsPtr = 0;
OPS_ResetCurrentInputArgType OPS_ResetCurrentInputArgPtr = 0;
OPS_GetStringType OPS_GetStringPtr = 0;
OPS_GetStringCopyType OPS_GetStringCopyPtr = 0;
OPS_GetIntPtrType OPS_GetNDM_Ptr = 0;
OPS_GetIntPtrType OPS_GetNDF_Ptr = 0;
OPS_GetFEDatastorePtrType OPS_GetFEDatastorePtr = 0;
OPS_GetInterpPWD_PtrType OPS_GetInterpPWD_Ptr = 0;
OPS_GetDomainPointer OPS_GetDomainPtr  = 0;
OPS_GetAnalysisModelPtrType OPS_GetAnalysisModelPtr  = 0;
OPS_GetAlgorithmPtrType OPS_GetAlgorithmPtr  = 0;
OPS_GetHandlerPtrType OPS_GetHandlerPtr  = 0;
OPS_GetNumbererPtrType OPS_GetNumbererPtr  = 0;
OPS_GetSOEPtrType OPS_GetSOEPtr  = 0;
OPS_GetEigenSOEPtrType OPS_GetEigenSOEPtr  = 0;
OPS_GetStaticAnalysisPtrType OPS_GetStaticAnalysisPtr  = 0;
OPS_GetTransientAnalysisPtrType OPS_GetTransientAnalysisPtr  = 0;
OPS_GetVariableTimeStepTransientAnalysisPtrType OPS_GetVariableTimeStepTransientAnalysisPtr  = 0;
OPS_GetNumEigenPtrType OPS_GetNumEigenPtr  = 0;
OPS_GetStaticIntegratorPtrType OPS_GetStaticIntegratorPtr  = 0;
OPS_GetTransientIntegratorPtrType OPS_GetTransientIntegratorPtr  = 0;
OPS_GetTestPtrType OPS_GetTestPtr  = 0;
OPS_builtModelPtrType OPS_builtModelPtr = 0;


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
                       OPS_GetCrdTransfPtrType OPS_GetCrdTransfFunct,
                       OPS_GetFrictionModelPtrType OPS_GetFrictionModelFunct,
                       OPS_InvokeMaterialDirectlyPtrType OPS_InvokeMaterialDirectlyFunct,
                       OPS_GetNodeInfoPtrType OPS_GetNodeCrdFunct,
                       OPS_GetNodeInfoPtrType OPS_GetNodeDispFunct,
                       OPS_GetNodeInfoPtrType OPS_GetNodeVelFunct,
                       OPS_GetNodeInfoPtrType OPS_GetNodeAccelFunct,
                       OPS_GetNodeInfoPtrType OPS_GetNodeIncrDispFunct,
                       OPS_GetNodeInfoPtrType OPS_GetNodeIncrDeltaDispFunct,
                       OPS_GetNumRemainingInputArgsType OPS_GetNumRemainingArgsFunct,
                       OPS_ResetCurrentInputArgType OPS_ResetCurrentInputArgFunct,
                       OPS_GetStringType OPS_GetStringFunct,
                       OPS_GetStringCopyType OPS_GetStringCopyFunct,
                       OPS_GetIntPtrType OPS_GetNDM_Funct,
                       OPS_GetIntPtrType OPS_GetNDF_Funct,
                       OPS_GetFEDatastorePtrType OPS_GetFEDatastoreFunct,
                       OPS_GetInterpPWD_PtrType OPS_GetInterpPWD_Funct,
		               OPS_GetAnalysisModelPtrType OPS_GetAnalysisModelFunct,
		               OPS_GetAlgorithmPtrType OPS_GetAlgorithmFunct,
		               OPS_GetHandlerPtrType OPS_GetHandlerFunct,
		               OPS_GetNumbererPtrType OPS_GetNumbererFunct,
		               OPS_GetSOEPtrType OPS_GetSOEFunct,
		               OPS_GetEigenSOEPtrType OPS_GetEigenSOEFunct,
		               OPS_GetStaticAnalysisPtrType OPS_GetStaticAnalysisFunct,
		               OPS_GetTransientAnalysisPtrType OPS_GetTransientAnalysisFunct,
		               OPS_GetVariableTimeStepTransientAnalysisPtrType OPS_GetVariableTimeStepTransientAnalysisFunct,
		               OPS_GetNumEigenPtrType OPS_GetNumEigenFunct,
		               OPS_GetStaticIntegratorPtrType OPS_GetStaticIntegratorFunct,
		               OPS_GetTransientIntegratorPtrType OPS_GetTransientIntegratorFunct,
		               OPS_GetTestPtrType OPS_GetTestFunct,
		               OPS_builtModelPtrType OPS_builtModelFunct,
		               OPS_GetDomainPointer OPS_getDomainPtr)
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
    OPS_GetCrdTransfPtr = OPS_GetCrdTransfFunct;
    OPS_GetFrictionModelPtr = OPS_GetFrictionModelFunct;
    OPS_GetNodeCrdPtr = OPS_GetNodeCrdFunct;
    OPS_GetNodeDispPtr = OPS_GetNodeDispFunct;
    OPS_GetNodeVelPtr = OPS_GetNodeVelFunct;
    OPS_GetNodeAccelPtr = OPS_GetNodeAccelFunct;
    OPS_GetNodeIncrDispPtr = OPS_GetNodeIncrDispFunct;
    OPS_GetNodeIncrDeltaDispPtr = OPS_GetNodeIncrDeltaDispFunct;
    OPS_InvokeMaterialDirectlyPtr = OPS_InvokeMaterialDirectlyFunct;
    OPS_GetNumRemainingInputArgsPtr = OPS_GetNumRemainingArgsFunct;
    OPS_ResetCurrentInputArgPtr = OPS_ResetCurrentInputArgFunct;
    OPS_GetStringPtr = OPS_GetStringFunct;
    OPS_GetStringCopyPtr = OPS_GetStringCopyFunct;
    OPS_GetNDM_Ptr = OPS_GetNDM_Funct;
    OPS_GetNDF_Ptr = OPS_GetNDF_Funct;
    OPS_GetFEDatastorePtr = OPS_GetFEDatastoreFunct;
    OPS_GetInterpPWD_Ptr = OPS_GetInterpPWD_Funct;
    OPS_GetAnalysisModelPtr = OPS_GetAnalysisModelFunct;
    OPS_GetAlgorithmPtr = OPS_GetAlgorithmFunct;
    OPS_GetHandlerPtr = OPS_GetHandlerFunct;
    OPS_GetNumbererPtr = OPS_GetNumbererFunct;
    OPS_GetSOEPtr = OPS_GetSOEFunct;
    OPS_GetEigenSOEPtr = OPS_GetEigenSOEFunct;
    OPS_GetStaticAnalysisPtr = OPS_GetStaticAnalysisFunct;
    OPS_GetTransientAnalysisPtr = OPS_GetTransientAnalysisFunct;
    OPS_GetVariableTimeStepTransientAnalysisPtr = OPS_GetVariableTimeStepTransientAnalysisFunct;
    OPS_GetNumEigenPtr = OPS_GetNumEigenFunct;
    OPS_GetStaticIntegratorPtr = OPS_GetStaticIntegratorFunct;
    OPS_GetTransientIntegratorPtr = OPS_GetTransientIntegratorFunct;
    OPS_GetTestPtr = OPS_GetTestFunct;
    OPS_builtModelPtr = OPS_builtModelFunct;
    OPS_GetDomainPtr = OPS_getDomainPtr;
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
OPS_GetSectionForceDeformation(int secTag)
{
    return (*OPS_GetSectionForceDeformationPtr)(secTag);
}

CrdTransf *
OPS_GetCrdTransf(int crdTag)
{
    return (*OPS_GetCrdTransfPtr)(crdTag);
}

FrictionModel *
OPS_GetFrictionModel(int frnTag)
{
    return (*OPS_GetFrictionModelPtr)(frnTag);
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

extern "C" const char *OPS_GetString()
{
    return (*OPS_GetStringPtr)();
}

extern "C" int OPS_GetStringCopy(char **cArray)
{
    return (*OPS_GetStringCopyPtr)(cArray);  
}

extern "C" int OPS_GetNumRemainingInputArgs()
{
    return (*OPS_GetNumRemainingInputArgsPtr)();  
}

extern "C" int OPS_ResetCurrentInputArg(int cArg)
{
    return (*OPS_ResetCurrentInputArgPtr)(cArg);  
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

extern "C" AnalysisModel **OPS_GetAnalysisModel(void)
{
	return (*OPS_GetAnalysisModelPtr)();
}

extern "C" EquiSolnAlgo **OPS_GetAlgorithm(void)
{
	return (*OPS_GetAlgorithmPtr)();
}

extern "C" ConstraintHandler **OPS_GetHandler(void)
{
	return (*OPS_GetHandlerPtr)();
}

extern "C" DOF_Numberer **OPS_GetNumberer(void)
{
	return (*OPS_GetNumbererPtr)();
}

extern "C" LinearSOE **OPS_GetSOE(void)
{
	return (*OPS_GetSOEPtr)();
}

extern "C" EigenSOE **OPS_GetEigenSOE(void)
{
	return (*OPS_GetEigenSOEPtr)();
}

extern "C" StaticAnalysis **OPS_GetStaticAnalysis(void)
{
	return (*OPS_GetStaticAnalysisPtr)();
}

extern "C" DirectIntegrationAnalysis **OPS_GetTransientAnalysis(void)
{
	return (*OPS_GetTransientAnalysisPtr)();
}

extern "C" VariableTimeStepDirectIntegrationAnalysis **OPS_GetVariableTimeStepTransientAnalysis(void)
{
	return (*OPS_GetVariableTimeStepTransientAnalysisPtr)();
}

extern "C" int *OPS_GetNumEigen(void)
{
	return (*OPS_GetNumEigenPtr)();
}

extern "C" StaticIntegrator **OPS_GetStaticIntegrator(void)
{
	return (*OPS_GetStaticIntegratorPtr)();
}

extern "C" TransientIntegrator **OPS_GetTransientIntegrator(void)
{
	return (*OPS_GetTransientIntegratorPtr)();
}

extern "C" ConvergenceTest **OPS_GetTest(void)
{
	return (*OPS_GetTestPtr)();
}

extern "C" bool *OPS_builtModel(void)
{
	return (*OPS_builtModelPtr)();
}

Domain *OPS_GetDomain(void)
{
    return (*OPS_GetDomainPtr)();
}
