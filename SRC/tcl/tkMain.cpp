/* 
 * 
 * TkMain.c --
 *
 *	This file contains a generic main program for Tk-based applications.
 *	It can be used as-is for many applications, just by supplying a
 *	different appInitProc procedure for each specific application.
 *	Or, it can be used as a template for creating new main programs
 *	for Tk applications.
 *
 * Copyright (c) 1990-1994 The Regents of the University of California.
 * Copyright (c) 1994-1997 Sun Microsystems, Inc.
 *
 * See Tcl/Tk License Terms for information on usage and redistribution
 * of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 *
 * RCS: @(#) $Id: tkMain.cpp,v 1.23 2010-04-23 23:01:15 fmk Exp $

Tcl/Tk License Terms:
This software is copyrighted by the Regents of the University of
California, Sun Microsystems, Inc., Scriptics Corporation, ActiveState
Corporation and other parties.  The following terms apply to all files
associated with the software unless explicitly disclaimed in
individual files.

The authors hereby grant permission to use, copy, modify, distribute,
and license this software and its documentation for any purpose, provided
that existing copyright notices are retained in all copies and that this
notice is included verbatim in any distributions. No written agreement,
license, or royalty fee is required for any of the authorized uses.
Modifications to this software may be copyrighted by their authors
and need not follow the licensing terms described here, provided that
the new terms are clearly indicated on the first page of each file where
they apply.

IN NO EVENT SHALL THE AUTHORS OR DISTRIBUTORS BE LIABLE TO ANY PARTY
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
ARISING OUT OF THE USE OF THIS SOFTWARE, ITS DOCUMENTATION, OR ANY
DERIVATIVES THEREOF, EVEN IF THE AUTHORS HAVE BEEN ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

THE AUTHORS AND DISTRIBUTORS SPECIFICALLY DISCLAIM ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT.  THIS SOFTWARE
IS PROVIDED ON AN "AS IS" BASIS, AND THE AUTHORS AND DISTRIBUTORS HAVE
NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
MODIFICATIONS.

GOVERNMENT USE: If you are acquiring this software on behalf of the
U.S. government, the Government shall have only "Restricted Rights"
in the software and related documentation as defined in the Federal 
Acquisition Regulations (FARs) in Clause 52.227.19 (c) (2).  If you
are acquiring the software on behalf of the Department of Defense, the
software shall be classified as "Commercial Computer Software" and the
Government shall have only "Restricted Rights" as defined in Clause
252.227-7013 (c) (1) of DFARs.  Notwithstanding the foregoing, the
authors grant the U.S. Government and others acting in its behalf
permission to use and distribute the software in accordance with the
terms specified in this license. 

 */

/*                       MODIFIED   FOR                              */

/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */


extern "C" {
#define WIN32_LEAN_AND_MEAN
#ifdef _WIN32
#include <windows.h>
#endif
#undef WIN32_LEAN_AND_MEAN
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <tcl.h>
#include <tk.h>
#ifdef NO_STDLIB_H
#   include "../compat/stdlib.h"
#else
#   include <stdlib.h>
#endif
}
#ifdef __WIN32__
#   include <tkWin.h>
#endif

#include <OPS_Globals.h>
extern "C" int OpenSeesParseArgv(int argc, char **argv);

typedef struct ThreadSpecificData {
    Tcl_Interp *interp;         /* Interpreter for this thread. */
    Tcl_DString command;        /* Used to assemble lines of terminal input
				 * into Tcl commands. */
    Tcl_DString line;           /* Used to read the next line from the
				 * terminal input. */
    int tty;                    /* Non-zero means standard input is a 
				 * terminal-like device.  Zero means it's
				 * a file. */
} ThreadSpecificData;
Tcl_ThreadDataKey dataKey;

/*
 * Declarations for various library procedures and variables (don't want
 * to include tkInt.h or tkPort.h here, because people might copy this
 * file out of the Tk source directory to make their own modified versions).
 * Note: don't declare "exit" here even though a declaration is really
 * needed, because it will conflict with a declaration elsewhere on
 * some systems.
 */

#if !defined(__WIN32__) && !defined(_WIN32) && !defined(_KAI)
extern "C" int		isatty _ANSI_ARGS_((int fd));
//extern "C" char *	strrchr _ANSI_ARGS_((CONST char *string, int c)) throw();
#endif

#if !defined(__WIN32__) && !defined(_WIN32) && defined(_KAI)
extern "C" int		isatty _ANSI_ARGS_((int fd));
//extern "C" char *	strrchr _ANSI_ARGS_((CONST char *string, int c));
#endif

#ifdef _TCL84
extern "C" void  TkpDisplayWarning _ANSI_ARGS_((const char *msg, char *title));
#elif _TCL85
extern "C" void  TkpDisplayWarning _ANSI_ARGS_((const char *msg, char *title));
#else
extern "C" void  TkpDisplayWarning _ANSI_ARGS_((char *msg, char *title));
#endif

/*
 * Forward declarations for procedures defined later in this file.
 */

static void		Prompt _ANSI_ARGS_((Tcl_Interp *interp, int partial));
static void		StdinProc _ANSI_ARGS_((ClientData clientData,
			    int mask));


static char *tclStartupScriptFileName = NULL;

void TclSetStartupScriptFileName(char *fileName)
{
    tclStartupScriptFileName = fileName;
}


char *TclGetStartupScriptFileName()
{
    return tclStartupScriptFileName;
}




/*
 *----------------------------------------------------------------------
 *
 * Tk_MainOpenSees --
 *
 *	Main program for Wish and most other Tk-based applications.
 *
 * Results:
 *	None. This procedure never returns (it exits the process when
 *	it's done.
 *
 * Side effects:
 *	This procedure initializes the Tk world and then starts
 *	interpreting commands;  almost anything could happen, depending
 *	on the script being interpreted.
 *
 *----------------------------------------------------------------------
 */
void
Tk_MainOpenSees(int argc, char **argv, Tcl_AppInitProc *appInitProc, Tcl_Interp *interp)
{
    char *args, *fileName;
    char buf[TCL_INTEGER_SPACE];
    int code;
    size_t length;
    Tcl_Channel inChannel, outChannel;
    Tcl_DString argString;
    ThreadSpecificData *tsdPtr;

#ifdef __WIN32__
    HANDLE handle;
#endif

    /* fmk - beginning of modifications for OpenSees */
    fprintf(stderr,"\n\n");
	fprintf(stderr,"         OpenSees -- Open System For Earthquake Engineering Simulation\n");
	fprintf(stderr,"                 Pacific Earthquake Engineering Research Center\n");
	fprintf(stderr,"                        Version %s %s\n\n", OPS_VERSION, WIN_ARCH);
	
	fprintf(stderr,"      (c) Copyright 1999-2016 The Regents of the University of California\n");
	fprintf(stderr,"                              All Rights Reserved\n");
	fprintf(stderr,"  (Copyright and Disclaimer @ http://www.berkeley.edu/OpenSees/copyright.html)\n\n\n");
    /* fmk - end of modifications for OpenSees */

    /*
     * Ensure that we are getting the matching version of Tcl.  This is
     * really only an issue when Tk is loaded dynamically.
     */

    if (Tcl_InitStubs(interp, TCL_VERSION, 0) == NULL) {
	abort();
    }

    tsdPtr = (ThreadSpecificData *) 
	Tcl_GetThreadData(&dataKey, sizeof(ThreadSpecificData));
    
    Tcl_FindExecutable(argv[0]);
    tsdPtr->interp = interp;

#if (defined(__WIN32__) || defined(MAC_TCL))
    Tk_InitConsoleChannels(interp);
#endif
    
#ifdef TCL_MEM_DEBUG
    Tcl_InitMemory(interp);
#endif

    /*
     * Parse command-line arguments.  A leading "-file" argument is
     * ignored (a historical relic from the distant past).  If the
     * next argument doesn't start with a "-" then strip it off and
     * use it as the name of a script file to process.
     */

    fileName = TclGetStartupScriptFileName();

    if (argc > 1) {
	length = strlen(argv[1]);
	if ((length >= 2) && (strncmp(argv[1], "-file", length) == 0)) {
	    argc--;
	    argv++;
	}
    }
    if (fileName == NULL) {
	if ((argc > 1) && (argv[1][0] != '-')) {
	    fileName = argv[1];
	    argc--;
	    argv++;
	}
    }
    
	OpenSeesParseArgv(argc, argv);

    /*
     * Make command-line arguments available in the Tcl variables "argc"
     * and "argv".
     */

    args = Tcl_Merge(argc-1, argv+1);
    Tcl_ExternalToUtfDString(NULL, args, -1, &argString);
    Tcl_SetVar(interp, "argv", Tcl_DStringValue(&argString), TCL_GLOBAL_ONLY);
    Tcl_DStringFree(&argString);
    ckfree(args);
    sprintf(buf, "%d", argc-1);

    if (fileName == NULL) {
	Tcl_ExternalToUtfDString(NULL, argv[0], -1, &argString);
    } else {
	fileName = Tcl_ExternalToUtfDString(NULL, fileName, -1, &argString);
    }
    Tcl_SetVar(interp, "argc", buf, TCL_GLOBAL_ONLY);
    Tcl_SetVar(interp, "argv0", Tcl_DStringValue(&argString), TCL_GLOBAL_ONLY);

    /*
     * Set the "tcl_interactive" variable.
     */

    /*
     * For now, under Windows, we assume we are not running as a console mode
     * app, so we need to use the GUI console.  In order to enable this, we
     * always claim to be running on a tty.  This probably isn't the right
     * way to do it.
     */

#ifdef __WIN32__
    handle = GetStdHandle(STD_INPUT_HANDLE);

    if ((handle == INVALID_HANDLE_VALUE) || (handle == 0) 
	     || (GetFileType(handle) == FILE_TYPE_UNKNOWN)) {
	/*
	 * If it's a bad or closed handle, then it's been connected
	 * to a wish console window.
	 */

	tsdPtr->tty = 1;
    } else if (GetFileType(handle) == FILE_TYPE_CHAR) {
	/*
	 * A character file handle is a tty by definition.
	 */

	tsdPtr->tty = 1;
    } else {
	tsdPtr->tty = 0;
    }

#else
    tsdPtr->tty = isatty(0);
#endif
    char one[2] = "1";
    char zero[2] = "0";
    Tcl_SetVar(interp, "tcl_interactive",
	    ((fileName == NULL) && tsdPtr->tty) ? one : zero, TCL_GLOBAL_ONLY);

    /*
     * Invoke application-specific initialization.
     */
	if ((*appInitProc)(interp) != TCL_OK) {
      TkpDisplayWarning(Tcl_GetStringResult(interp), "Application Initialization Failed");
	}
    /*
     * Invoke the script specified on the command line, if any.
     */

    if (fileName != NULL) {
	Tcl_ResetResult(interp);
	code = Tcl_EvalFile(interp, fileName);
	if (code != TCL_OK) {
	    /*
	     * The following statement guarantees that the errorInfo
	     * variable is set properly.
	     */

	    Tcl_AddErrorInfo(interp, "");
	    TkpDisplayWarning(Tcl_GetVar(interp, "error Info",
					 TCL_GLOBAL_ONLY), "Error in startup script");
	    Tcl_DeleteInterp(interp);
	    Tcl_Exit(1);
	}
	tsdPtr->tty = 0;
    } else {

	/*
	 * Evaluate the .rc file, if one has been specified.
	 */

	Tcl_SourceRCFile(interp);

	/*
	 * Establish a channel handler for stdin.
	 */

	inChannel = Tcl_GetStdChannel(TCL_STDIN);
	if (inChannel) {
	    Tcl_CreateChannelHandler(inChannel, TCL_READABLE, StdinProc,
		    (ClientData) inChannel);
	}
	if (tsdPtr->tty) {
	    Prompt(interp, 0);
	}
    }
    Tcl_DStringFree(&argString);

    outChannel = Tcl_GetStdChannel(TCL_STDOUT);
    if (outChannel) {
	Tcl_Flush(outChannel);
    }
    Tcl_DStringInit(&tsdPtr->command);
    Tcl_DStringInit(&tsdPtr->line);
    Tcl_ResetResult(interp);

    /*
     * Loop infinitely, waiting for commands to execute.  When there
     * are no windows left, Tk_MainLoop returns and we exit.
     */

    Tk_MainLoop();
    Tcl_DeleteInterp(interp);
    Tcl_Exit(0);
}

/*
 *----------------------------------------------------------------------
 *
 * StdinProc --
 *
 *	This procedure is invoked by the event dispatcher whenever
 *	standard input becomes readable.  It grabs the next line of
 *	input characters, adds them to a command being assembled, and
 *	executes the command if it's complete.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	Could be almost arbitrary, depending on the command that's
 *	typed.
 *
 *----------------------------------------------------------------------
 */

    /* ARGSUSED */
static void
StdinProc(ClientData clientData, int mask)
{
    static int gotPartial = 0;
    char *cmd;
    int code, count;
    Tcl_Channel chan = (Tcl_Channel) clientData;
    ThreadSpecificData *tsdPtr = (ThreadSpecificData *) 
            Tcl_GetThreadData(&dataKey, sizeof(ThreadSpecificData));
    Tcl_Interp *interp = tsdPtr->interp;

    count = Tcl_Gets(chan, &tsdPtr->line);

    if (count < 0) {
	if (!gotPartial) {
	    if (tsdPtr->tty) {
		Tcl_Exit(0);
	    } else {
		Tcl_DeleteChannelHandler(chan, StdinProc, (ClientData) chan);
	    }
	    return;
	} 
    }

    (void) Tcl_DStringAppend(&tsdPtr->command, Tcl_DStringValue(
            &tsdPtr->line), -1);
    cmd = Tcl_DStringAppend(&tsdPtr->command, "\n", -1);
    Tcl_DStringFree(&tsdPtr->line);
    if (!Tcl_CommandComplete(cmd)) {
        gotPartial = 1;
        goto prompt;
    }
    gotPartial = 0;

    /*
     * Disable the stdin channel handler while evaluating the command;
     * otherwise if the command re-enters the event loop we might
     * process commands from stdin before the current command is
     * finished.  Among other things, this will trash the text of the
     * command being evaluated.
     */

    Tcl_CreateChannelHandler(chan, 0, StdinProc, (ClientData) chan);
    code = Tcl_RecordAndEval(interp, cmd, TCL_EVAL_GLOBAL);
    
    chan = Tcl_GetStdChannel(TCL_STDIN);
    if (chan) {
	Tcl_CreateChannelHandler(chan, TCL_READABLE, StdinProc,
		(ClientData) chan);
    }
    Tcl_DStringFree(&tsdPtr->command);
    if (Tcl_GetStringResult(interp)[0] != '\0') {
	if ((code != TCL_OK) || (tsdPtr->tty)) {
	    chan = Tcl_GetStdChannel(TCL_STDOUT);
	    if (chan) {
		Tcl_WriteObj(chan, Tcl_GetObjResult(interp));
		Tcl_WriteChars(chan, "\n", 1);
	    }
	}
    }

    /*
     * Output a prompt.
     */

    prompt:
    if (tsdPtr->tty) {
	Prompt(interp, gotPartial);
    }
    Tcl_ResetResult(interp);
}

/*
 *----------------------------------------------------------------------
 *
 * Prompt --
 *
 *	Issue a prompt on standard output, or invoke a script
 *	to issue the prompt.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	A prompt gets output, and a Tcl script may be evaluated
 *	in interp.
 *
 *----------------------------------------------------------------------
 */

static void
Prompt(Tcl_Interp *interp, int partial)
{

#ifdef _TCL84
  const char *promptCmd;
  const char one[12] = "tcl_prompt1";
  const char two[12] = "tcl_prompt2";
#elif _TCL85
  const char *promptCmd;
  const char one[12] = "tcl_prompt1";
  const char two[12] = "tcl_prompt2";
#else
  char *promptCmd;
  char one[12] = "tcl_prompt1";
  char two[12] = "tcl_prompt2";
#endif

  int code;
  Tcl_Channel outChannel, errChannel;

  promptCmd = Tcl_GetVar(interp, partial ? two : one, TCL_GLOBAL_ONLY);
			   
  if (promptCmd == NULL) {
  defaultPrompt:
	if (!partial) {
	  
	  /*
	   * We must check that outChannel is a real channel - it
	   * is possible that someone has transferred stdout out of
	   * this interpreter with "interp transfer".
	   */
	  
	  outChannel = Tcl_GetChannel(interp, "stdout", NULL);
	  if (outChannel != (Tcl_Channel) NULL) {
	    Tcl_WriteChars(outChannel, "OpenSees > ", 11);
	  }
	}
    } else {
      code = Tcl_Eval(interp, promptCmd);
      if (code != TCL_OK) {
	Tcl_AddErrorInfo(interp,
			 "\n    (script that generates prompt)");
	/*
	 * We must check that errChannel is a real channel - it
	 * is possible that someone has transferred stderr out of
	 * this interpreter with "interp transfer".
	 */
	
	errChannel = Tcl_GetChannel(interp, "stderr", NULL);
	if (errChannel != (Tcl_Channel) NULL) {
	  Tcl_WriteObj(errChannel, Tcl_GetObjResult(interp));
	  Tcl_WriteChars(errChannel, "\n", 1);
	}
	goto defaultPrompt;
      }
    }
  outChannel = Tcl_GetChannel(interp, "stdout", NULL);
  if (outChannel != (Tcl_Channel) NULL) {
    Tcl_Flush(outChannel);
  }
}



