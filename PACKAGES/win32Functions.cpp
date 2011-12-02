// NewCommand.cpp : Defines the entry point for the DLL application and
//    a function that can be called to set the global pointer variables in 
//    the dll to be the same as those in the existing process address space.

#include <OPS_Stream.h>
#include <Domain.h>
#include <TclModelBuilder.h>

#include <windows.h>

#define DllExport _declspec(dllexport)

BOOL APIENTRY DllMain( HANDLE hModule, 
                       DWORD  ul_reason_for_call, 
                       LPVOID lpReserved
					 )
{
    return TRUE;
}


extern "C" DllExport
void setGlobalPointers(OPS_Stream *theErrorStreamPtr)
{
	opserrPtr = theErrorStreamPtr;
}