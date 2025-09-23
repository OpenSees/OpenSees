//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Written: fmk
//
#include <stdlib.h>
#include <string.h>
#include <OPS_Globals.h>
#include <tcl.h>
#include <sys/stat.h>

#ifdef _WIN32
#  define byte win_byte_override
#  include <windows.h>
#  include <elementAPI.h>
#  include <OpenSeesFFI.h>
   class VariableTimeStepDirectIntegrationAnalysis;
#else
#  include <dlfcn.h>
#endif

int
getLibraryFunction(const char *libName, const char *funcName, void **libHandle,
                   void **funcHandle)
{

  int result = 0;

  *libHandle = nullptr;
  *funcHandle = nullptr;

#ifdef _WIN32

    //
    // first try and open dll
    //

    int libNameLength = (int)strlen(libName);
    char* localLibName = new char[libNameLength + 5];
    strcpy(localLibName, libName);
    strcpy(&localLibName[libNameLength], ".dll");

    HINSTANCE hLib = LoadLibrary(localLibName);

    delete[] localLibName;

    if (hLib != nullptr) {

        char mod[124];
        GetModuleFileName((HMODULE)hLib, (LPTSTR)mod, 124);

        //
        // Now look for function with funcName
        //

        (*funcHandle) = (void*)GetProcAddress((HMODULE)hLib, funcName);

        if (*funcHandle == nullptr) {
            char* underscoreFunctionName = new char[strlen(funcName) + 2];
            strcpy(underscoreFunctionName, funcName);
            strcpy(&underscoreFunctionName[strlen(funcName)], "_");
            (*funcHandle) = (void*)GetProcAddress((HMODULE)hLib, underscoreFunctionName);
            delete[] underscoreFunctionName;
        }


        if (*funcHandle == nullptr) {
            FreeLibrary((HMODULE)hLib);
            return -2;
        }

        //
        // we need to set the OpenSees pointer global variables if function there
        //

        typedef int(_cdecl* LocalInitPtrType)();
        typedef int(_cdecl* OPS_ErrorPtrType)(char*, int);
        typedef int(_cdecl* OPS_GetNumRemainingInputArgsType)();
        typedef int(_cdecl* OPS_ResetCurrentInputArgType)(int);
        //typedef int(_cdecl* OPS_ResetInputType)(ClientData, Tcl_Interp*, int, int, TCL_Char**, Domain*, TclModelBuilder*);
        typedef int(_cdecl* OPS_ResetInputNoBuilderType)(ClientData, Tcl_Interp*, int, int, TCL_Char**, Domain*);
        typedef int(_cdecl* OPS_GetIntInputPtrType)(int*, int*);
        typedef int(_cdecl* OPS_GetDoubleInputPtrType)(int*, double*);
        typedef const char* (_cdecl* OPS_GetStringType)();
        typedef int(_cdecl* OPS_GetStringCopyType)(char**);
        typedef int(_cdecl* OPS_AllocateElementPtrType)(eleObj*, int*, int*);
        typedef int(_cdecl* OPS_AllocateMaterialPtrType)(matObj*);
        typedef UniaxialMaterial* (*OPS_GetUniaxialMaterialPtrType)(int);
        typedef NDMaterial* (*OPS_GetNDMaterialPtrType)(int);
        typedef SectionForceDeformation* (*OPS_GetSectionForceDeformationPtrType)(int);
        typedef CrdTransf* (*OPS_GetCrdTransfPtrType)(int);
        typedef FrictionModel* (*OPS_GetFrictionModelPtrType)(int);
        typedef int(_cdecl* OPS_GetNodeInfoPtrType)(int*, int*, double*);
        typedef int(_cdecl* OPS_InvokeMaterialDirectlyPtrType)(matObject**, modelState*, double*, double*, double*, int*);
        typedef int(_cdecl* OPS_GetIntPtrType)();

        typedef FE_Datastore* (*OPS_GetFEDatastorePtrType)();
        typedef const char* (_cdecl* OPS_GetInterpPWD_PtrType)();

        typedef Domain* (*OPS_GetDomainPointerType)();
        typedef AnalysisModel** (*OPS_GetAnalysisModelPtrType)();
        typedef EquiSolnAlgo** (*OPS_GetAlgorithmPtrType)();
        typedef ConstraintHandler** (*OPS_GetHandlerPtrType)();
        typedef DOF_Numberer** (*OPS_GetNumbererPtrType)();
        typedef LinearSOE** (*OPS_GetSOEPtrType)();
        typedef EigenSOE** (*OPS_GetEigenSOEPtrType)();
        typedef StaticAnalysis** (*OPS_GetStaticAnalysisPtrType)();
        typedef DirectIntegrationAnalysis** (*OPS_GetTransientAnalysisPtrType)();
        typedef VariableTimeStepDirectIntegrationAnalysis** (*OPS_GetVariableTimeStepTransientAnalysisPtrType)();
        typedef int* (*OPS_GetNumEigenPtrType)();
        typedef StaticIntegrator** (*OPS_GetStaticIntegratorPtrType)();
        typedef TransientIntegrator** (*OPS_GetTransientIntegratorPtrType)();
        typedef ConvergenceTest** (*OPS_GetTestPtrType)();
        typedef bool* (*OPS_builtModelPtrType)();
#if 0
        typedef void(_cdecl* setGlobalPointersFunction)(
            OPS_Stream*,
            Domain*,
            SimulationInformation*,
            OPS_ErrorPtrType,
            OPS_GetIntInputPtrType,
            OPS_GetDoubleInputPtrType,
            OPS_AllocateElementPtrType,
            OPS_AllocateMaterialPtrType,
            OPS_GetUniaxialMaterialPtrType,
            OPS_GetNDMaterialPtrType,
            OPS_GetSectionForceDeformationPtrType,
            OPS_GetCrdTransfPtrType,
            OPS_GetFrictionModelPtrType,
            OPS_InvokeMaterialDirectlyPtrType,
            OPS_GetNodeInfoPtrType,
            OPS_GetNodeInfoPtrType,
            OPS_GetNodeInfoPtrType,
            OPS_GetNodeInfoPtrType,
            OPS_GetNodeInfoPtrType,
            OPS_GetNodeInfoPtrType,
            OPS_GetNumRemainingInputArgsType,
            OPS_ResetCurrentInputArgType,
            //OPS_ResetInputType,
            OPS_ResetInputNoBuilderType,
            OPS_GetStringType,
            OPS_GetStringCopyType,
            OPS_GetIntPtrType,
            OPS_GetIntPtrType,
            OPS_GetFEDatastorePtrType,
            OPS_GetInterpPWD_PtrType,
            OPS_GetAnalysisModelPtrType,
            OPS_GetAlgorithmPtrType,
            OPS_GetHandlerPtrType,
            OPS_GetNumbererPtrType,
            OPS_GetSOEPtrType,
            OPS_GetEigenSOEPtrType,
            OPS_GetStaticAnalysisPtrType,
            OPS_GetTransientAnalysisPtrType,
            OPS_GetVariableTimeStepTransientAnalysisPtrType,
            OPS_GetNumEigenPtrType,
            OPS_GetStaticIntegratorPtrType,
            OPS_GetTransientIntegratorPtrType,
            OPS_GetTestPtrType,
            OPS_builtModelPtrType,
            OPS_GetDomainPointerType);

        setGlobalPointersFunction funcPtr;
        // look for pointer function
        funcPtr = (setGlobalPointersFunction)GetProcAddress((HMODULE)hLib, "setGlobalPointers");
        if (funcPtr == 0) {
            FreeLibrary((HMODULE)hLib);
            return -2;
        }
#endif
#if 0
        // invoke pointer function
        (funcPtr)(opserrPtr,
            ops_TheActiveDomain,
            nullptr, // theSimulationInfoPtr,
            OPS_Error,
            OPS_GetIntInput,
            OPS_GetDoubleInput,
            OPS_AllocateElement,
            OPS_AllocateMaterial,
            OPS_GetUniaxialMaterial,
            OPS_GetNDMaterial,
            OPS_GetSectionForceDeformation,
            OPS_GetCrdTransf,
            OPS_GetFrictionModel,
            OPS_InvokeMaterialDirectly,
            OPS_GetNodeCrd,
            OPS_GetNodeDisp,
            OPS_GetNodeVel,
            OPS_GetNodeAccel,
            OPS_GetNodeIncrDisp,
            OPS_GetNodeIncrDeltaDisp,
            OPS_GetNumRemainingInputArgs,
            OPS_ResetCurrentInputArg,
            //OPS_ResetInput,
            OPS_ResetInputNoBuilder,
            OPS_GetString,
            OPS_GetStringCopy,
            OPS_GetNDM,
            OPS_GetNDF,
            OPS_GetFEDatastore,
            OPS_GetInterpPWD,
            OPS_GetAnalysisModel,
            OPS_GetAlgorithm,
            OPS_GetHandler,
            OPS_GetNumberer,
            OPS_GetSOE,
            OPS_GetEigenSOE,
            OPS_GetStaticAnalysis,
            OPS_GetTransientAnalysis,
            OPS_GetVariableTimeStepTransientAnalysis,
            OPS_GetNumEigen,
            OPS_GetStaticIntegrator,
            OPS_GetTransientIntegrator,
            OPS_GetTest,
            OPS_builtModel,
            OPS_GetDomain);
#  endif
        LocalInitPtrType initPtr;
        initPtr = (LocalInitPtrType)GetProcAddress((HMODULE)hLib, "localInit");
        if (initPtr != 0) {
            initPtr();
        }
        else {
            initPtr = (LocalInitPtrType)GetProcAddress((HMODULE)hLib, "localinit_");
            if (initPtr != 0) {
                initPtr();
            }
        }

    }
    else // no lib exists
        return -1;

    libHandle = (void**)&hLib;

#else

  int libNameLength = strlen(libName);
  char *localLibName = new char[libNameLength + 10];
  strcpy(localLibName, libName);

#ifdef _MACOSX
  strcpy(&localLibName[libNameLength], ".dylib");
#else
  strcpy(&localLibName[libNameLength], ".so");
#endif

  // Attempt to get the file attributes
  // intintStat = stat(localLibName, &stFileInfo);
  /* get library
  if(intStat != 0) {
    opserr << "packages.cpp - NO FILE EXISTS: - trying OpenSees" << localLibName
  << endln; int res = httpGET_File("opensees.berkeley.edu", localLibName, 80,
  localLibName); if (res != 0) { opserr << "packages.cpp - NO FILE EXISTS: " <<
  localLibName << endln; return -1;
    }
  }
  */

  *libHandle = dlopen(localLibName, RTLD_NOW);

  if (*libHandle == nullptr) {
    delete[] localLibName;
    return -1; // no lib exists
  }

  void *funcPtr = dlsym(*libHandle, funcName);

  //
  // look for fortran procedure, trailing underscore
  //

  if (funcPtr == nullptr) {
    int funcNameLength = strlen(funcName);
    char *underscoreFunctionName = new char[funcNameLength + 2];
    strcpy(underscoreFunctionName, funcName);
    strcpy(&underscoreFunctionName[funcNameLength], "_");
    strcpy(&underscoreFunctionName[funcNameLength + 1], "");
    funcPtr = dlsym(*libHandle, underscoreFunctionName);
    delete[] underscoreFunctionName;
  }

  if (funcPtr == nullptr) {
    dlclose(*libHandle);
    delete[] localLibName;
    return -1;
  }

  *funcHandle = funcPtr;

  typedef int (*localInitPtrType)();
  localInitPtrType initFunct;
  funcPtr = dlsym(*libHandle, "localInit");

  if (funcPtr != nullptr) {
    initFunct = (localInitPtrType)funcPtr;
    initFunct();
  } else {
    funcPtr = dlsym(*libHandle, "localinit_");
    if (funcPtr != nullptr) {
      initFunct = (localInitPtrType)funcPtr;
      initFunct();
    }
  }

  delete[] localLibName;

#endif

  return result;
}
