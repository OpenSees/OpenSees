// NewCommand.cpp : Defines the entry point for the DLL application and
//    a function that can be called to set the global pointer variables in 
//    the dll to be the same as those in the existing process address space.


#include <stdio.h>
#include <stdlib.h>

/*
#include <OPS_ProceduralAPI.h>
typedef int(*OPS_ErrorPtrType)(char *data, int length);
typedef int(*OPS_GetIntInputPtrType)(int *numData, int*data);
typedef int(*OPS_GetDoubleInputPtrType)(int *numData, double *data);
OPS_ErrorPtrType OPS_ErrorPtr = NULL;
OPS_GetIntInputPtrType OPS_GetIntInputPtr = NULL;
OPS_GetDoubleInputPtrType OPS_GetDoubleInputPtr = NULL;
*/


#ifdef _Win32
#include <OPS_Stream.h>
#include <Domain.h>
#include <TclModelBuilder.h>

#include <windows.h>
#endif

//#include <SimulationInformation.h>
//SimulationInformation simulationInfo;

#ifdef _USRDLL
#define DllExport _declspec(dllexport)
#elif _MACOSX
#define DllExport __attribute__((visibility("default")))
#else
#define DllExport
#endif

#ifdef _USERDLL

BOOL APIENTRY DllMain( HANDLE hModule, 
                       DWORD  ul_reason_for_call, 
                       LPVOID lpReserved
					 )
{
    return TRUE;
}

int __cdecl OpenSeesExit(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{

  Tcl_Exit(0);
  return 0;
}

#elif _MACOSX

// Initializer.
__attribute__((constructor))
static void initializer(void) {                             // 2
  //  printf("[%s] initializer()\n", __FILE__);
}

__attribute__((destructor))
static void finalizer(void) {                               // 3
  //  printf("[%s] finalizer()\n", __FILE__);
}

#endif

/*
extern "C" DllExport
void setGlobalPointers(OPS_Stream *theErrorStreamPtr,
		       OPS_ErrorPtrType          errorFunct,
		       OPS_GetIntInputPtrType    getIntInputFunct,
		       OPS_GetDoubleInputPtrType getDoubleInputFunct)
{
  opserrPtr = theErrorStreamPtr;
  OPS_ErrorPtr = errorFunct;
  OPS_GetIntInputPtr = getIntInputFunct;
  OPS_GetDoubleInputPtr =getDoubleInputFunct;;
}


int        OPS_Error(char *data, int length)
{
  return (*OPS_ErrorPtr)(data, length);
}


int        OPS_GetIntInput(int *numData, int*data)
{
  opserr << "int OPS_GetIntInput(int *numData, int*data)\n";
  return -1;
  //  return (*OPS_GetIntInputPtr)(numData, data);
}


int        OPS_GetDoubleInput(int *numData, double *data)
{
  opserr << "int        OPS_GetDoubleInput(int *numData, int*data)\n";
  return (*OPS_GetDoubleInputPtr)(numData, data);
}
*/
