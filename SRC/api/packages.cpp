/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */

/*                                                                        
** $Revision: 1.13 $
** $Date: 2009-10-02 22:20:35 $
** $Source: /usr/local/cvs/OpenSees/SRC/api/packages.cpp,v $
                                                                        
** Written: fmk 
*/                                                                        

#include <stdlib.h>
#include <string.h>
#include <OPS_Globals.h>
#include <sys/stat.h>
#include <SimulationInformation.h>



extern
#ifdef _WIN32
int __cdecl
#else
int
#endif
httpGET_File(char const *URL, char const *page, unsigned int port, const char *filename);



#ifdef _WIN32

#include <windows.h>
#include <elementAPI.h>
extern SimulationInformation *theSimulationInfoPtr;

#else
#include <dlfcn.h>
#endif

int 
getLibraryFunction(const char *libName, const char *funcName, void **libHandle, void **funcHandle) {

  int result = 0;
  
  *libHandle = NULL;
  *funcHandle = NULL;

  struct stat stFileInfo;
  bool blnReturn;
  int intStat;

  
#ifdef _WIN32
  
  //
  // first try and open dll
  //
  
  int libNameLength = strlen(libName);
  char *localLibName = new char[libNameLength+5];
  strcpy(localLibName, libName);
  strcpy(&localLibName[libNameLength], ".dll");
  
  HINSTANCE hLib = LoadLibrary(localLibName);
  
  delete [] localLibName;

  if (hLib != NULL) {

    char mod[124];
    GetModuleFileName((HMODULE)hLib, (LPTSTR)mod, 124);
    
    //
    // Now look for function with funcName
    //
    
    (*funcHandle) = (void *)GetProcAddress((HMODULE)hLib, funcName);    
    
    if (*funcHandle == NULL) {
      char *underscoreFunctionName = new char[strlen(funcName)+2];
      strcpy(underscoreFunctionName, funcName);
      strcpy(&underscoreFunctionName[strlen(funcName)], "_");   
      (*funcHandle) = (void *)GetProcAddress((HMODULE)hLib, underscoreFunctionName);
      delete [] underscoreFunctionName;
    }

    
    if (*funcHandle == NULL) {
      FreeLibrary((HMODULE)hLib);
      return -2;
    } 

    //
    // we need to set the OpenSees pointer global variables if function there
    //

    typedef int (_cdecl *LocalInitPtrType)();
    typedef int (_cdecl *OPS_ErrorPtrType)(char *, int);
    typedef int (_cdecl *OPS_GetNumRemainingInputArgsType)();
    typedef int (_cdecl *OPS_GetIntInputPtrType)(int *, int *);
    typedef int (_cdecl *OPS_GetDoubleInputPtrType)(int *, double *);
    typedef int (_cdecl *OPS_GetStringType)(char *, int);
    typedef int (_cdecl *OPS_GetStringCopyType)(char **);
    typedef int (_cdecl *OPS_AllocateElementPtrType)(eleObj *, int *matTags, int *maType);
    typedef int (_cdecl *OPS_AllocateMaterialPtrType)(matObj *);
    typedef UniaxialMaterial *(*OPS_GetUniaxialMaterialPtrType)(int matTag);
    typedef NDMaterial *(*OPS_GetNDMaterialPtrType)(int matTag);
    typedef CrdTransf *(*OPS_GetCrdTransfPtrType)(int matTag);
    typedef int (_cdecl *OPS_GetNodeInfoPtrType)(int *, int *, double *);
    typedef int (_cdecl *OPS_InvokeMaterialDirectlyPtrType)(matObject **, modelState *, double *, double *, double *, int *);
    typedef int (_cdecl *OPS_GetIntPtrType)();

    typedef FE_Datastore *(*OPS_GetFEDatastorePtrType)();
    typedef const char * (_cdecl *OPS_GetInterpPWD_PtrType)();
    
    typedef void (_cdecl *setGlobalPointersFunction)(OPS_Stream *,
						     SimulationInformation *,
						     OPS_ErrorPtrType,
						     OPS_GetIntInputPtrType,
						     OPS_GetDoubleInputPtrType,
						     OPS_AllocateElementPtrType,
						     OPS_AllocateMaterialPtrType,
						     OPS_GetUniaxialMaterialPtrType,
						     OPS_GetNDMaterialPtrType,
						     OPS_InvokeMaterialDirectlyPtrType,
						     OPS_GetNodeInfoPtrType,
						     OPS_GetNodeInfoPtrType,
						     OPS_GetNodeInfoPtrType,
						     OPS_GetNodeInfoPtrType,
						     OPS_GetNodeInfoPtrType,
						     OPS_GetNodeInfoPtrType,
						     OPS_GetNumRemainingInputArgsType,
						     OPS_GetStringType,
						     OPS_GetStringCopyType,
						     OPS_GetCrdTransfPtrType,
						     OPS_GetIntPtrType,
						     OPS_GetIntPtrType,
						     OPS_GetFEDatastorePtrType,
						     OPS_GetInterpPWD_PtrType);

    setGlobalPointersFunction funcPtr;
    
    // look for pointer function
    funcPtr = (setGlobalPointersFunction)GetProcAddress((HMODULE)hLib,"setGlobalPointers");
    if (funcPtr == 0) {
      FreeLibrary((HMODULE)hLib);
      return -2;
    }
  
    // invoke pointer function
    (funcPtr)(opserrPtr, theSimulationInfoPtr, OPS_Error, OPS_GetIntInput, OPS_GetDoubleInput,
	      OPS_AllocateElement, OPS_AllocateMaterial, OPS_GetUniaxialMaterial, 
	      OPS_GetNDMaterial, OPS_InvokeMaterialDirectly, OPS_GetNodeCrd, 
	      OPS_GetNodeDisp, OPS_GetNodeVel, OPS_GetNodeAcc, 
	      OPS_GetNodeIncrDisp, OPS_GetNodeIncrDeltaDisp,
	      OPS_GetNumRemainingInputArgs, OPS_GetString, OPS_GetStringCopy, OPS_GetCrdTransfPtr, OPS_GetNDM, OPS_GetNDF,
	      OPS_GetFEDatastore, OPS_GetInterpPWD);

   LocalInitPtrType initPtr;
   initPtr = (LocalInitPtrType)GetProcAddress((HMODULE)hLib,"localInit");
   if (initPtr !=0) {
     initPtr();
   } else {
	  initPtr = (LocalInitPtrType)GetProcAddress((HMODULE)hLib,"localinit_");
	  if (initPtr !=0) {
	    initPtr();
	  }
   }
    
  } else // no lib exists
    return -1;
  
  libHandle =  (void **)&hLib;

#else

  int libNameLength = strlen(libName);
  char *localLibName = new char[libNameLength+10];
  strcpy(localLibName, libName);

#ifdef _MACOSX
  strcpy(&localLibName[libNameLength], ".dylib");
#else
  strcpy(&localLibName[libNameLength], ".so");
#endif

  // Attempt to get the file attributes
  intStat = stat(localLibName, &stFileInfo);
  /* get library
  if(intStat != 0) {
    opserr << "packages.cpp - NO FILE EXISTS: - trying OpenSees" << localLibName << endln;
    int res = httpGET_File("opensees.berkeley.edu", localLibName, 80, localLibName);
    if (res != 0) {
      opserr << "packages.cpp - NO FILE EXISTS: " << localLibName << endln;
      return -1;
    }
  } 
  */
  char *error;

  *libHandle = dlopen (localLibName, RTLD_NOW);
  
  if (*libHandle == NULL)  {
    delete [] localLibName;
    return -1; // no lib exists
  }

  void *funcPtr = dlsym(*libHandle, funcName);
  
  error = dlerror();
  
  //
  // look for fortran procedure, trailing underscore
  //
  
  if (funcPtr == NULL ) {
    int funcNameLength  =strlen(funcName);
    char *underscoreFunctionName = new char[funcNameLength+2];
    strcpy(underscoreFunctionName, funcName);
    strcpy(&underscoreFunctionName[funcNameLength], "_");   
    strcpy(&underscoreFunctionName[funcNameLength+1], "");    
    funcPtr = dlsym(*libHandle, underscoreFunctionName);
    delete [] underscoreFunctionName;
  } 
  
  if (funcPtr == NULL)  {
    dlclose(*libHandle);
    delete [] localLibName;
    return -1;
  }
  
  
  *funcHandle = funcPtr;
  
  typedef int (*localInitPtrType)();
  localInitPtrType initFunct;
  funcPtr = dlsym(*libHandle, "localInit");
  
  if (funcPtr != NULL ) {
    initFunct = (localInitPtrType)funcPtr;
    initFunct();
  } else {
    funcPtr = dlsym(*libHandle, "localinit_");
    if (funcPtr != NULL ) {
      initFunct = (localInitPtrType)funcPtr;
      initFunct();
    }
  }
  
  delete [] localLibName;

#endif

  return result;
}
