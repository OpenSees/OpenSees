

#include "OpenSeesEvaluator.h"

#include <string.h>
#include <iostream>
using std::cerr;

extern "C" {
#include <tcl.h>
#include <tk.h>
EXTERN int		TclFormatInt _ANSI_ARGS_((char *buffer, long n));
EXTERN int		TclObjCommandComplete _ANSI_ARGS_((Tcl_Obj *cmdPtr));
}


#include <XmlFileStream.h>
#include <SimulationInformation.h>
SimulationInformation simulationInfo;
char *simulationInfoOutputFilename =0;

#include "commands.h"

/*
 * The following code ensures that tclLink.c is linked whenever
 * Tcl is linked.  Without this code there's no reference to the
 * code in that file from anywhere in Tcl, so it may not be
 * linked into the application.
 */
#ifdef _WIN32

#else
#ifdef _MAC
EXTERN int Tcl_LinkVar(Tcl_Interp*, char*, char*, int);
#else
//EXTERN int Tcl_LinkVar();

#ifdef _TCL84
int (*tclDummyLinkVarPtr)(Tcl_Interp *interp, const char *a,
			  char *b, int c) = Tcl_LinkVar;
#else
int (*tclDummyLinkVarPtr)(Tcl_Interp *interp, const char *a,
			  char *b, int c) = Tcl_LinkVar;
#endif
//int (*tclDummyLinkVarPtr)() = Tcl_LinkVar;
#endif
#endif

# undef TCL_STORAGE_CLASS
# define TCL_STORAGE_CLASS DLLEXPORT

/*
 * Declarations for various library procedures and variables (don't want
 * to include tclPort.h here, because people might copy this file out of
 * the Tcl source directory to make their own modified versions).
 * Note:  "exit" should really be declared here, but there's no way to
 * declare it without causing conflicts with other definitions elsewher
 * on some systems, so it's better just to leave it out.
 */

extern "C" int	isatty _ANSI_ARGS_((int fd));
extern "C" char * strcpy _ANSI_ARGS_((char *dst, CONST char *src)) throw();
static char *tclStartupScriptFileName = NULL;

static Tcl_Channel inChannel;
static Tcl_Channel outChannel;
static Tcl_Channel errChannel;
static Tcl_Interp *interp =0;

// the following is a little kludgy but it works!
#ifdef _USING_STL_STREAMS

#include <iomanip>
using std::ios;
#include <iostream>
using std::ofstream;

#else

/*
#include <FileStream.h>
FileStream sserr;
OPS_Stream *opserrPtr = &sserr;
*/

#endif


/*
 *----------------------------------------------------------------------
 *
 * TclSetStartupScriptFileName --
 *
 *	Primes the startup script file name, used to override the
 *      command line processing.
 *
 * Results:
 *	None. 
 *
 * Side effects:
 *	This procedure initializes the file name of the Tcl script to
 *      run at startup.
 *
 *----------------------------------------------------------------------
 */
void TclSetStartupScriptFileName(char *fileName)
{
    tclStartupScriptFileName = fileName;
}


/*
 *----------------------------------------------------------------------
 *
 * TclGetStartupScriptFileName --
 *
 *	Gets the startup script file name, used to override the
 *      command line processing.
 *
 * Results:
 *	The startup script file name, NULL if none has been set.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */
char *TclGetStartupScriptFileName()
{
    return tclStartupScriptFileName;
}

void g3TclMain(int argc, char **argv, Tcl_AppInitProc *appInitProc);
int		Tcl_AppInit _ANSI_ARGS_((Tcl_Interp *interp));



JNIEXPORT jint JNICALL Java_OpenSeesEvaluator_openSeesInit(JNIEnv *env, jobject obj)
{

#ifndef TCL_LOCAL_APPINIT
#define TCL_LOCAL_APPINIT Tcl_AppInit    
#endif

  Tcl_Obj *commandPtr = NULL;
  int gotPartial, tty;
  
  // fmk    Tcl_FindExecutable(argv[0]);
  interp = Tcl_CreateInterp();

  /*  sserr.setFile("opensees.err", OVERWRITE); */
  
  tty = isatty(0);
  char one[2] = "1";
  char zero[2] = "0";
  
  Tcl_SetVar(interp, "tcl_interactive",
	     ((tclStartupScriptFileName == NULL) && tty) ? one : zero,
	     TCL_GLOBAL_ONLY);
  
  /*
   * Invoke application-specific initialization.
   */
  
  if ((*TCL_LOCAL_APPINIT)(interp) != TCL_OK) {
    errChannel = Tcl_GetStdChannel(TCL_STDERR);
    if (errChannel) {
      Tcl_WriteChars(errChannel,
		     "application-specific initialization failed: ", -1);
      Tcl_WriteObj(errChannel, Tcl_GetObjResult(interp));
      Tcl_WriteChars(errChannel, "\n", 1);
    }
  }
  
  /*
   * Process commands from stdin until there's an end-of-file.  Note
   * that we need to fetch the standard channels again after every
   * eval, since they may have been changed.
   */
  
  commandPtr = Tcl_NewObj();
  Tcl_IncrRefCount(commandPtr);
  
  inChannel = Tcl_GetStdChannel(TCL_STDIN);
  outChannel = Tcl_GetStdChannel(TCL_STDOUT);
  gotPartial = 0;
  
  return 0;
}
  
JNIEXPORT jstring JNICALL Java_OpenSeesEvaluator_openSeesEval(JNIEnv *env, 
							      jobject obj, 
							      jstring expression, 
							      jint error)
{
  const char *exp = env->GetStringUTFChars(expression, 0);
  char *copyExp = new char[strlen(exp) +1];
  strcpy(copyExp, exp);
  error = Tcl_Eval(interp, copyExp);
  /* cerr << error << "\n"; */
  const char *res = Tcl_GetStringResult(interp);

  env->ReleaseStringUTFChars(expression, exp);
  // result = NewString(env, res, strlen(res));
  return  env->NewStringUTF(res);
}

JNIEXPORT jint JNICALL Java_OpenSeesEvaluator_openSeesQuit(JNIEnv *env, jobject obj)
{
  Tcl_Eval(interp, "quit");

  /*  sserr.close(); */

  return 0;
}


int Tcl_AppInit(Tcl_Interp *interp)
{
    if (Tcl_Init(interp) == TCL_ERROR) {
	return TCL_ERROR;
    }

    /*
     * Call the init procedures for included packages.  Each call should
     * look like this:
     *
     * if (Mod_Init(interp) == TCL_ERROR) {
     *     return TCL_ERROR;
     * }
     *
     * where "Mod" is the name of the module.
     */

    /*
     * Call Tcl_CreateCommand for application-specific commands, if
     * they weren't already created by the init procedures called above.
     */

    if (g3AppInit(interp) < 0)
	  return TCL_ERROR; 

    /*
     * Specify a user-specific startup file to invoke if the application
     * is run interactively.  Typically the startup file is "~/.apprc"
     * where "app" is the name of the application.  If this line is deleted
     * then no user-specific startup file will be run under any conditions.
     */

    Tcl_SetVar(interp, "tcl_rcFileName", "~/.tclshrc", TCL_GLOBAL_ONLY);
    return TCL_OK;
}

int OpenSeesExit(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	
  if (simulationInfoOutputFilename != 0) {
    simulationInfo.end();
    XmlFileStream simulationInfoOutputFile;
    simulationInfoOutputFile.setFile(simulationInfoOutputFilename);
    simulationInfoOutputFile.open();
    simulationInfoOutputFile << simulationInfo;
    simulationInfoOutputFile.close();
  }

  Tcl_Exit(0);
  return 0;
}
