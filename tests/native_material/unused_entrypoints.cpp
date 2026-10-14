// Standalone material tests do not invoke interpreter parsing or dense solves.
// Fail loudly if either becomes necessary; all constitutive functions are real.
#include <cstdlib>
#include <map>
#include <vector>
#include <elementAPI.h>
#include <DummyStream.h>
namespace { DummyStream sink; }
OPS_Stream *opserrPtr = &sink;
extern "C" int OPS_GetNumRemainingInputArgs() {std::abort();}
extern "C" int OPS_GetIntInput(int*,int*) {std::abort();}
extern "C" int OPS_GetDoubleInput(int*,double*) {std::abort();}
extern "C" int OPS_Error(char*,int) {std::abort();}
extern "C" int OPS_ResetCurrentInputArg(int) {std::abort();}
extern "C" const char* OPS_GetString() {std::abort();}
extern "C" int OPS_SetDoubleDictListOutput(std::map<const char*,std::vector<double>>&) {std::abort();}
#ifdef _WIN32
extern "C" int DGESV(int*,int*,double*,int*,int*,double*,int*,int*) {std::abort();}
extern "C" int DGETRF(int*,int*,double*,int*,int*,int*) {std::abort();}
extern "C" int DGETRI(int*,double*,int*,int*,double*,int*,int*) {std::abort();}
#else
extern "C" int dgesv_(int*,int*,double*,int*,int*,double*,int*,int*) {std::abort();}
extern "C" int dgetrf_(int*,int*,double*,int*,int*,int*) {std::abort();}
extern "C" int dgetri_(int*,double*,int*,int*,double*,int*,int*) {std::abort();}
#endif
