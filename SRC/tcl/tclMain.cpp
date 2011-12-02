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
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.2 $
// $Date: 2000-12-12 06:00:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/tcl/tclMain.cpp,v $
                                                                        
                                                                        
/* 
 * tclMain.c --
 *
 *	Main program for Tcl shells and other Tcl-based applications.
 *
 * Copyright (c) 1988-1994 The Regents of the University of California.
 * Copyright (c) 1994-1996 Sun Microsystems, Inc.
 *
 * See the file "license.terms" for information on usage and redistribution
 * of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 *
 * SCCS: @(#) tclMain.c 1.54 97/08/07 19:04:43
 *
 * revision: modified somewhat for g3 and C++ by fmk.
 */

#include <iostream.h>

extern "C" {
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <tcl.h>
}


# undef TCL_STORAGE_CLASS
# define TCL_STORAGE_CLASS DLLEXPORT

/*
 * The following code ensures that tclLink.c is linked whenever
 * Tcl is linked.  Without this code there's no reference to the
 * code in that file from anywhere in Tcl, so it may not be
 * linked into the application.
 */

// EXTERN int Tcl_LinkVar();
int (*tclDummyLinkVarPtr)(Tcl_Interp *interp, char *a,
			  char *b, int c) = Tcl_LinkVar;

/*
 * Declarations for various library procedures and variables (don't want
 * to include tclPort.h here, because people might copy this file out of
 * the Tcl source directory to make their own modified versions).
 * Note:  "exit" should really be declared here, but there's no way to
 * declare it without causing conflicts with other definitions elsewher
 * on some systems, so it's better just to leave it out.
 */

extern "C" int	isatty _ANSI_ARGS_((int fd));
extern "C" char * strcpy _ANSI_ARGS_((char *dst, CONST char *src));

extern "C" int TclFormatInt _ANSI_ARGS_((char *buffer, long n));
extern "C" int TclObjCommandComplete _ANSI_ARGS_((Tcl_Obj *cmdPtr));

static Tcl_Interp *interp;	/* Interpreter for application. */

#ifdef TCL_MEM_DEBUG
static char dumpFile[100];	/* Records where to dump memory allocation
				 * information. */
static int quitFlag = 0;	/* 1 means "checkmem" command was called,
				 * so the application should quit and dump
				 * memory allocation information. */
#endif

/*
 * Forward references for procedures defined later in this file:
 */

#ifdef TCL_MEM_DEBUG
static "C" int		CheckmemCmd _ANSI_ARGS_((ClientData clientData,
			    Tcl_Interp *interp, int argc, char *argv[]));
#endif

/*
 *----------------------------------------------------------------------
 *
 * g3TclMain --
 *
 *	Main program for tclsh and most other Tcl-based applications.
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
g3TclMain(int argc, char **argv, Tcl_AppInitProc *appInitProc)
{
    Tcl_Obj *prompt1NamePtr = NULL;
    Tcl_Obj *prompt2NamePtr = NULL;
    Tcl_Obj *resultPtr;
    Tcl_Obj *commandPtr = NULL;
    char buffer[1000], *args, *fileName, *bytes;
    int code, gotPartial, tty, length;
    int exitCode = 0;
    Tcl_Channel inChannel, outChannel, errChannel;

    fprintf(stderr,"\n\n\t    OpenSees -- Open System For Earthquake Engineering Simulation");
    fprintf(stderr,"\n\t\t    Pacific Earthquake Engineering Research Center\n\n");
    
    fprintf(stderr,"\t    (c) Copyright 1999 The Regents of the University of California");
    fprintf(stderr,"\n\t\t\t\t All Rights Reserved \n\n\n");    

    Tcl_FindExecutable(argv[0]);
    interp = Tcl_CreateInterp();
#ifdef TCL_MEM_DEBUG
    Tcl_InitMemory(interp);
    Tcl_CreateCommand(interp, "checkmem", CheckmemCmd, (ClientData) 0,
	    (Tcl_CmdDeleteProc *) NULL);
#endif

    /*
     * Make command-line arguments available in the Tcl variables "argc"
     * and "argv".  If the first argument doesn't start with a "-" then
     * strip it off and use it as the name of a script file to process.
     */

    fileName = NULL;
    if ((argc > 1) && (argv[1][0] != '-')) {
	fileName = argv[1];
	argc--;
	argv++;
    }
    args = Tcl_Merge(argc-1, argv+1);
    Tcl_SetVar(interp, "argv", args, TCL_GLOBAL_ONLY);
    Tcl_Free(args);
    TclFormatInt(buffer, argc-1);
    Tcl_SetVar(interp, "argc", buffer, TCL_GLOBAL_ONLY);
    Tcl_SetVar(interp, "argv0", (fileName != NULL) ? fileName : argv[0],
	    TCL_GLOBAL_ONLY);

    /*
     * Set the "tcl_interactive" variable.
     */

    tty = isatty(0);
    char one[2] = "1";
    char zero[2] = "0";
    Tcl_SetVar(interp, "tcl_interactive",
	    ((fileName == NULL) && tty) ? one : zero, TCL_GLOBAL_ONLY);
    
    /*
     * Invoke application-specific initialization.
     */

    if ((*appInitProc)(interp) != TCL_OK) {
	errChannel = Tcl_GetStdChannel(TCL_STDERR);
	if (errChannel) {
	    Tcl_Write(errChannel,
		    "application-specific initialization failed: ", -1);
	    Tcl_Write(errChannel, interp->result, -1);
	    Tcl_Write(errChannel, "\n", 1);
	}
    }

    /*
     * If a script file was specified then just source that file
     * and quit.
     */

    if (fileName != NULL) {
	code = Tcl_EvalFile(interp, fileName);
	if (code != TCL_OK) {
	    errChannel = Tcl_GetStdChannel(TCL_STDERR);
	    if (errChannel) {
		/*
		 * The following statement guarantees that the errorInfo
		 * variable is set properly.
		 */

		Tcl_AddErrorInfo(interp, "");
		Tcl_Write(errChannel,
			Tcl_GetVar(interp, "errorInfo", TCL_GLOBAL_ONLY), -1);
		Tcl_Write(errChannel, "\n", 1);
	    }
	    exitCode = 1;
	}
	goto done;
    }

    /*
     * We're running interactively.  Source a user-specific startup
     * file if the application specified one and if the file exists.
     */

    Tcl_SourceRCFile(interp);

    /*
     * Process commands from stdin until there's an end-of-file.  Note
     * that we need to fetch the standard channels again after every
     * eval, since they may have been changed.
     */

    commandPtr = Tcl_NewObj();
    Tcl_IncrRefCount(commandPtr);
    prompt1NamePtr = Tcl_NewStringObj("tcl_prompt1", -1);
    Tcl_IncrRefCount(prompt1NamePtr);
    prompt2NamePtr = Tcl_NewStringObj("tcl_prompt2", -1);
    Tcl_IncrRefCount(prompt2NamePtr);
    
    inChannel = Tcl_GetStdChannel(TCL_STDIN);
    outChannel = Tcl_GetStdChannel(TCL_STDOUT);
    gotPartial = 0;
    while (1) {
	if (tty) {
	    Tcl_Obj *promptCmdPtr;

	    promptCmdPtr = Tcl_ObjGetVar2(interp,
		    (gotPartial? prompt2NamePtr : prompt1NamePtr),
		    (Tcl_Obj *) NULL, TCL_GLOBAL_ONLY);
	    if (promptCmdPtr == NULL) {
                defaultPrompt:
		if (!gotPartial && outChannel) {
		    Tcl_Write(outChannel, "OpenSees > ", 11);
		}
	    } else {
		code = Tcl_EvalObj(interp, promptCmdPtr);
		inChannel = Tcl_GetStdChannel(TCL_STDIN);
		outChannel = Tcl_GetStdChannel(TCL_STDOUT);
		errChannel = Tcl_GetStdChannel(TCL_STDERR);
		if (code != TCL_OK) {
		    if (errChannel) {
			resultPtr = Tcl_GetObjResult(interp);
			bytes = Tcl_GetStringFromObj(resultPtr, &length);
			Tcl_Write(errChannel, bytes, length);
			Tcl_Write(errChannel, "\n", 1);
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
	if (!TclObjCommandComplete(commandPtr)) {
	    gotPartial = 1;
	    continue;
	}


	gotPartial = 0;
	code = Tcl_RecordAndEvalObj(interp, commandPtr, 0);
	inChannel = Tcl_GetStdChannel(TCL_STDIN);
	outChannel = Tcl_GetStdChannel(TCL_STDOUT);
	errChannel = Tcl_GetStdChannel(TCL_STDERR);
	Tcl_SetObjLength(commandPtr, 0);
	if (code != TCL_OK) {
	    if (errChannel) {
		resultPtr = Tcl_GetObjResult(interp);
		bytes = Tcl_GetStringFromObj(resultPtr, &length);
		Tcl_Write(errChannel, bytes, length);
		Tcl_Write(errChannel, "\n", 1);
	    }
	} else if (tty) {
	    resultPtr = Tcl_GetObjResult(interp);
	    bytes = Tcl_GetStringFromObj(resultPtr, &length);
	    if ((length > 0) && outChannel) {
		Tcl_Write(outChannel, bytes, length);
		Tcl_Write(outChannel, "\n", 1);
	    }
	}
#ifdef TCL_MEM_DEBUG
	if (quitFlag) {
	    Tcl_DecrRefCount(commandPtr);
	    Tcl_DecrRefCount(prompt1NamePtr);
	    Tcl_DecrRefCount(prompt2NamePtr);
	    Tcl_DeleteInterp(interp);
	    Tcl_Exit(0);
	}
#endif
    }

    /*
     * Rather than calling exit, invoke the "exit" command so that
     * users can replace "exit" with some other command to do additional
     * cleanup on exit.  The Tcl_Eval call should never return.
     */

    done:
    if (commandPtr != NULL) {
	Tcl_DecrRefCount(commandPtr);
    }
    if (prompt1NamePtr != NULL) {
	Tcl_DecrRefCount(prompt1NamePtr);
    }
    if (prompt2NamePtr != NULL) {
	Tcl_DecrRefCount(prompt2NamePtr);
    }
    sprintf(buffer, "exit %d", exitCode);
    Tcl_Eval(interp, buffer);
}

/*
 *----------------------------------------------------------------------
 *
 * CheckmemCmd --
 *
 *	This is the command procedure for the "checkmem" command, which
 *	causes the application to exit after printing information about
 *	memory usage to the file passed to this command as its first
 *	argument.
 *
 * Results:
 *	Returns a standard Tcl completion code.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */
#ifdef TCL_MEM_DEBUG

	/* ARGSUSED */
static int
CheckmemCmd(clientData, interp, argc, argv)
    ClientData clientData;		/* Not used. */
    Tcl_Interp *interp;			/* Interpreter for evaluation. */
    int argc;				/* Number of arguments. */
    char *argv[];			/* String values of arguments. */
{
    extern char *tclMemDumpFileName;
    if (argc != 2) {
	Tcl_AppendResult(interp, "wrong # args: should be \"", argv[0],
		" fileName\"", (char *) NULL);
	return TCL_ERROR;
    }
    strcpy(dumpFile, argv[1]);
    tclMemDumpFileName = dumpFile;
    quitFlag = 1;
    return TCL_OK;
}
#endif
