//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
#include <string.h>
#include <stdlib.h>
#include <string>

#ifndef _WIN32
#  include <unistd.h>
#endif

#ifdef USE_TCL_STUBS
#undef USE_TCL_STUBS
#endif

extern "C" {
#include <tcl.h>
#include <tclDecls.h>
 #ifdef _TCL85
 #define TclFormatInt(buf, n)   sprintf((buf),"%ld", (long)(n))
 #else
 EXTERN int  TclFormatInt _ANSI_ARGS_((char *buffer, long n));
 #endif
 EXTERN int  TclObjCommandComplete _ANSI_ARGS_((Tcl_Obj *cmdPtr));
}

extern "C" int 
#ifdef _WIN32
__declspec(dllexport)
#endif
Openseesrt_Init(Tcl_Interp *interp);

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

#ifdef _WIN32
extern "C" int          isatty _ANSI_ARGS_((int fd));
//extern "C" char * strcpy _ANSI_ARGS_((char *dst, CONST char *src)) throw();
#endif

static char *tclStartupScriptFileName = NULL;

void 
TclSetStartupScriptFileName(char *fileName)
{
  tclStartupScriptFileName = fileName;
}

char *
TclGetStartupScriptFileName()
{
  return tclStartupScriptFileName;
}

/*
 *----------------------------------------------------------------------
 *
 * Main program for the xara interactive shell.
 * 
 *----------------------------------------------------------------------
 */


int
main(int argc, char **argv)
{
  Tcl_Obj *resultPtr;
  Tcl_Obj *commandPtr = NULL;
  char buffer[1000], *args;
  int code, gotPartial, tty, length;
  int exitCode = 0;
  Tcl_Channel inChannel, outChannel, errChannel;
  Tcl_DString argString;
  static bool OPS_showHeader = false;

  Tcl_Interp *interp = Tcl_CreateInterp();
  if (Tcl_InitStubs(interp, "8.5-10", 0) == NULL) {
    fprintf(stderr, "Tcl_InitStubs failed: %s\n", Tcl_GetStringResult(interp));
    exit(1);
  }


  if (OPS_showHeader) {
    fprintf(stderr,"\n\n");
    fprintf(stderr,"         OpenSees -- Open System For Earthquake Engineering Simulation\n");
    fprintf(stderr,"                 Pacific Earthquake Engineering Research Center\n");
    
    fprintf(stderr,"      (c) Copyright 1999-2025 The Regents of the University of California\n");
    fprintf(stderr,"                              All Rights Reserved\n");
    fprintf(stderr,"  (Copyright and Disclaimer @ http://www.berkeley.edu/OpenSees/copyright.html)\n\n\n");
  }

  Tcl_Eval(interp, "rename load import;");
  Tcl_Eval(interp, "interp alias {} load {} import;");


#ifdef TCL_MEM_DEBUG
  Tcl_InitMemory(interp);
#endif

  /*
    * Make command-line arguments available in the Tcl variables "argc"
    * and "argv".  If the first argument doesn't start with a "-" then
    * strip it off and use it as the name of a script file to process.
    */
  tclStartupScriptFileName = argv[1];
  if (tclStartupScriptFileName == NULL) {
    if ((argc > 1) && (argv[1][0] != '-')) {
      tclStartupScriptFileName = argv[1];
      argc--;
      argv++;
    }
  }

  args = Tcl_Merge(argc-1, argv+1);
  Tcl_ExternalToUtfDString(NULL, args, -1, &argString);
  Tcl_SetVar(interp, "argv", Tcl_DStringValue(&argString), TCL_GLOBAL_ONLY);
  Tcl_DStringFree(&argString);
  ckfree(args);


  if (tclStartupScriptFileName == NULL) {
    Tcl_ExternalToUtfDString(NULL, argv[0], -1, &argString);
  } else {
    tclStartupScriptFileName = Tcl_ExternalToUtfDString(NULL,
              tclStartupScriptFileName, -1, &argString);
  }

  /*
    * Set the "tcl_interactive" variable.
    */

  tty = isatty(0);
  char one[2] = "1";
  char zero[2] = "0";

  Tcl_SetVar(interp, "tcl_interactive",
            ((tclStartupScriptFileName == NULL) && tty) ? one : zero,
            TCL_GLOBAL_ONLY);
  
  //
  // Load the OpenSeesRT library
  //
  if (Openseesrt_Init(interp) != TCL_OK) {
      fprintf(stderr, "Error loading OpenSeesRT library: %s\n",
              Tcl_GetStringResult(interp));
  }

  /*
    * If a script file was specified then just source that file
    * and quit.
    */

  if (tclStartupScriptFileName != NULL) {

    code = Tcl_EvalFile(interp, tclStartupScriptFileName);
    
    if (code != TCL_OK) {
      errChannel = Tcl_GetStdChannel(TCL_STDERR);
      if (errChannel) {
        //
        // The following statement guarantees that the errorInfo
        // variable is set properly.
        //
        Tcl_AddErrorInfo(interp, "");
        Tcl_WriteObj(errChannel, Tcl_GetVar2Ex(interp, "errorInfo",
                                                          NULL, TCL_GLOBAL_ONLY));
        Tcl_WriteChars(errChannel, "\n", 1);
      }
      exitCode = 1;
    }
    goto done;
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

  while (true) {
    if (tty) {
      Tcl_Obj *promptCmdPtr;
      
      char one[12] = "tcl_prompt1";
      char two[12] = "tcl_prompt2";
      promptCmdPtr = Tcl_GetVar2Ex(interp,
                                          (gotPartial ? one : two),
                                          NULL, TCL_GLOBAL_ONLY);
      if (promptCmdPtr == NULL) {
      defaultPrompt:
        if (!gotPartial && outChannel) {
          Tcl_WriteChars(outChannel, "OpenSees > ", 11);
        }
      } else {
        
        code = Tcl_EvalObjEx(interp, promptCmdPtr, 0);
        
        inChannel = Tcl_GetStdChannel(TCL_STDIN);
        outChannel = Tcl_GetStdChannel(TCL_STDOUT);
        errChannel = Tcl_GetStdChannel(TCL_STDERR);
        if (code != TCL_OK) {
          if (errChannel) {
              Tcl_WriteObj(errChannel, Tcl_GetObjResult(interp));
              Tcl_WriteChars(errChannel, "\n", 1);
          }
          Tcl_AddErrorInfo(interp,
                                "\n    (script that generates prompt)");
          goto defaultPrompt;
        }
      }
      if (outChannel) {
        Tcl_Flush(outChannel);
        }
    }
    if (!inChannel) {
      goto done;
    }
    length = Tcl_GetsObj(inChannel, commandPtr);
    if (length < 0) {
      goto done;
    }
    if ((length == 0) && Tcl_Eof(inChannel) && (!gotPartial)) {
      goto done;
    }
          
    /*
      * Add the newline removed by Tcl_GetsObj back to the string.
      */
    
    Tcl_AppendToObj(commandPtr, "\n", 1);

    gotPartial = 0;
    code = Tcl_RecordAndEvalObj(interp, commandPtr, 0);
    inChannel = Tcl_GetStdChannel(TCL_STDIN);
    outChannel = Tcl_GetStdChannel(TCL_STDOUT);
    errChannel = Tcl_GetStdChannel(TCL_STDERR);
    Tcl_DecrRefCount(commandPtr);
    commandPtr = Tcl_NewObj();
    Tcl_IncrRefCount(commandPtr);
    if (code != TCL_OK) {
      if (errChannel) {
        Tcl_WriteObj(errChannel, Tcl_GetObjResult(interp));
        Tcl_WriteChars(errChannel, "\n", 1);
      }
    } else if (tty) {
      resultPtr = Tcl_GetObjResult(interp);
      Tcl_GetStringFromObj(resultPtr, &length);
      if ((length > 0) && outChannel) {
        Tcl_WriteObj(outChannel, resultPtr);
        Tcl_WriteChars(outChannel, "\n", 1);
      }
    }
  }

done:

  if (commandPtr != NULL)
    Tcl_DecrRefCount(commandPtr);

#if defined(_PARALLEL_PROCESSING) || defined( _PARALLEL_INTERPRETERS)
  return;
#endif

  /*
    * Rather than calling exit, invoke the "exit" command so that
    * users can replace "exit" with some other command to do additional
    * cleanup on exit.  The Tcl_Eval call should never return.
    */
  Tcl_Eval(interp, buffer);

  Tcl_Eval(interp, "quit"); 

  return exitCode;
}
