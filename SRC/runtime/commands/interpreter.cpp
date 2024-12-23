//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains basic commands that enhance the
// experience of the interpreter. This file should not reference
// any analysis or modeling classes.
//
#include <tcl.h>
#include <assert.h>
#include <runtimeAPI.h>
#include <G3_Runtime.h>
#include <OPS_Globals.h>
#include <Timer.h>

static Tcl_ObjCmdProc *Tcl_putsCommand = nullptr;
static Timer *theTimer = nullptr;

Tcl_CmdProc TclCommand_wipeModel;
Tcl_CmdProc TclCommand_clearAnalysis;
Tcl_CmdProc TclCommand_specifyModel;

class ProgressBar;
Tcl_ObjCmdProc TclObjCommand_progress;
extern ProgressBar* progress_bar_ptr;


const char *getInterpPWD(Tcl_Interp *interp);

int TclObjCommand_pragma([[maybe_unused]] ClientData clientData, 
                     Tcl_Interp *interp, int objc, Tcl_Obj *const objv[]);

// formats.cpp
Tcl_CmdProc convertBinaryToText;
Tcl_CmdProc convertTextToBinary;
Tcl_CmdProc stripOpenSeesXML;

//
// Consider reimplmenting to use Tcl built-ins; see
// https://wiki.tcl-lang.org/page/timers
//
static int
startTimer(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  if (theTimer == nullptr)
    theTimer = new Timer();

  theTimer->start();
  return TCL_OK;
}

static int
stopTimer(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  if (theTimer == nullptr)
    return TCL_OK;

  theTimer->pause();
  opserr << *theTimer;
  return TCL_OK;
}

static int
timer(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  if ((argc == 1) || (strcmp(argv[1], "start")==0)) {
    stopTimer(clientData, interp, argc, argv);
    return startTimer(clientData, interp, argc, argv);
  } else if (strcmp(argv[1], "stop")==0) {
    return stopTimer(clientData, interp, argc, argv);
  }
  opserr << "Unknown argument '" << argv[1] << "'\n";
  return TCL_ERROR;
}

//
// revised puts command to send to stderr
//
static int
OpenSees_putsCommand(ClientData dummy, Tcl_Interp *interp, int objc,
                     Tcl_Obj *const objv[])
{
  // Tcl_Channel chan;           // The channel to puts on.
  Tcl_Obj *string;            /* String to write. */
  Tcl_Obj *chanObjPtr = NULL; /* channel object. */
  int newline;                /* Add a newline at end? */

  switch (objc) {
  case 2: /* [puts $x] */
    string = objv[1];
    newline = 1;
    break;

  case 3: /* [puts -nonewline $x] or [puts $chan $x] */
    if (strcmp(Tcl_GetString(objv[1]), "-nonewline") == 0) {
      newline = 0;
    } else {
      newline = 1;
      chanObjPtr = objv[1];
    }
    string = objv[2];
    break;

  case 4: /* [puts -nonewline $chan $x] or [puts $chan $x nonewline] */
    newline = 0;
    if (strcmp(Tcl_GetString(objv[1]), "-nonewline") == 0) {
      chanObjPtr = objv[2];
      string = objv[3];
      break;

    } else if (strcmp(Tcl_GetString(objv[3]), "nonewline") == 0) {
      /*
       * The code below provides backwards compatibility with an old
       * form of the command that is no longer recommended or
       * documented. See also [Bug #3151675]. Will be removed in Tcl 9,
       * maybe even earlier.
       */
      chanObjPtr = objv[1];
      string = objv[2];
      break;
    }
    /* Fall through */
  default:
    /* [puts] or [puts some bad number of arguments...] */
    Tcl_WrongNumArgs(interp, 1, objv, "?-nonewline? ?channelId? string");
    return TCL_ERROR;
  }

  if (chanObjPtr == NULL) {
    G3_Runtime *rt;
    if ((rt = G3_getRuntime(interp))) {
      if (newline == 0)
        fprintf(rt->streams[1], "%s", Tcl_GetString(string));
      else
        fprintf(rt->streams[1], "%s\n", Tcl_GetString(string));

    } else {
      if (newline == 0)
        opserr << Tcl_GetString(string);
      else
        opserr << Tcl_GetString(string) << "\n";
    }
    return TCL_OK;
  } else {
    assert(Tcl_putsCommand != nullptr);
    return Tcl_putsCommand(dummy, interp, objc, objv);

    return TCL_OK;
  }
}


static int
OPS_SetObjCmd(ClientData clientData, Tcl_Interp *interp, int objc,
              Tcl_Obj *const objv[])
{

  Tcl_Obj *varValueObj;

  if (objc == 2) {
    varValueObj = Tcl_ObjGetVar2(interp, objv[1], NULL, TCL_LEAVE_ERR_MSG);
    if (varValueObj == NULL) {
      return TCL_ERROR;
    }
    Tcl_SetObjResult(interp, varValueObj);
    return TCL_OK;
  } else if (objc == 3) {
    varValueObj =
        Tcl_ObjSetVar2(interp, objv[1], NULL, objv[2], TCL_LEAVE_ERR_MSG);
    if (varValueObj == NULL) {
      return TCL_ERROR;
    }
    Tcl_SetObjResult(interp, varValueObj);
    return TCL_OK;
  } else {
    Tcl_WrongNumArgs(interp, 1, objv, "varName ?newValue?");
    return TCL_ERROR;
  }

  return 0;
}


static int
logFile(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{

  if (argc < 2) {
    opserr << "WARNING logFile fileName? - no filename supplied\n";
    return TCL_ERROR;
  }
  bool echo = true;
  openMode mode = openMode::OVERWRITE;

  int cArg = 2;
  while (cArg < argc) {
    if (strcmp(argv[cArg], "-append") == 0)
      mode = openMode::APPEND;

    if (strcmp(argv[cArg], "-noEcho") == 0)
      echo = false;
    cArg++;
  }

  if (opserr.setFile(argv[1], mode, echo) < 0)
    opserr << "WARNING logFile " << argv[1] << " failed to set the file\n";


  return TCL_OK;
}

static int
setPrecision(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char ** const argv)
{

  if (argc < 2) {
    opserr << "WARNING setPrecision precision? - no precision value supplied\n";
    return TCL_ERROR;
  }
  int precision;
  if (Tcl_GetInt(interp, argv[1], &precision) != TCL_OK) {
    opserr << "WARNING setPrecision precision? - error reading precision value "
              "supplied\n";
    return TCL_ERROR;
  }
  opserr.setPrecision(precision);

  return TCL_OK;
}

const char *
getInterpPWD(Tcl_Interp *interp)
{
  static char *pwd = 0;

  if (pwd != 0)
    delete[] pwd;

#ifdef _TCL84
  Tcl_Obj *cwd = Tcl_FSGetCwd(interp);
  if (cwd != NULL) {
    int length;
    const char *objPWD = Tcl_GetStringFromObj(cwd, &length);
    pwd = new char[length + 1];
    strcpy(pwd, objPWD);
    Tcl_DecrRefCount(cwd);
  }
#else
  Tcl_DString buf;
  const char *objPWD = Tcl_GetCwd(interp, &buf);

  pwd = new char[strlen(objPWD) + 1];
  strcpy(pwd, objPWD);

  Tcl_DStringFree(&buf);
#endif
  return pwd;
}

int
OPS_SourceCmd(ClientData dummy,      /* Not used. */
              Tcl_Interp *interp,    /* Current interpreter. */
              int objc,              /* Number of arguments. */
              Tcl_Obj *CONST objv[]) /* Argument objects. */
{
  Tcl_Obj *fileName;

  if (objc != 2 && objc != 4) {
    Tcl_WrongNumArgs(interp, 1, objv, "?-encoding name? fileName");
    return TCL_ERROR;
  }

  fileName = objv[objc - 1];

  if (objc == 4) {
    static CONST char *options[] = {"-encoding", NULL};
    int index;

    if (TCL_ERROR == Tcl_GetIndexFromObj(interp, objv[1], options, "option",
                                         TCL_EXACT, &index)) {
      return TCL_ERROR;
    }
  }

  const char *fileN = Tcl_GetString(fileName);

  return Tcl_EvalFile(interp, fileN);
}

static int
OpenSeesExit(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char ** const argv)
{
  //  theDomain.clearAll();

#ifdef _PARALLEL_PROCESSING
  // mpi clean up
  if (theMachineBroker != 0) {
    theMachineBroker->shutdown();
    fprintf(stderr, "Process Terminating\n");
  }
  MPI_Finalize();
#endif

#ifdef _PARALLEL_INTERPRETERS
  // mpi clean up
  if (theMachineBroker != 0) {
    theMachineBroker->shutdown();
    fprintf(stderr, "Process Terminating\n");
  }
  MPI_Finalize();
#endif

  int returnCode = 0;
  if (argc > 1) {
    if (Tcl_GetInt(interp, argv[1], &returnCode) != TCL_OK)
      opserr << "WARNING: OpenSeesExit - failed to read return code\n";
  }
  Tcl_Exit(returnCode);

  return 0;
}

static int
maxOpenFiles(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char ** const argv)
{
  int maxOpenFiles;

  if (Tcl_GetInt(interp, argv[1], &maxOpenFiles) != TCL_OK) {
    return TCL_ERROR;
  }

#ifdef _WIN32
  int newMax = _setmaxstdio(maxOpenFiles);
  if (maxOpenFiles > 2048) {
    opserr << "setMaxOpenFiles: too many files specified (2048 max)\n";
  } else {
    if (newMax != maxOpenFiles) {
      opserr << "setMaxOpenFiles FAILED: max allowed files: " << newMax;
      return TCL_ERROR;
    }
  }
  return TCL_OK;
#endif

  opserr << "setMaxOpenFiles FAILED: - command not available on this machine\n";
  return TCL_OK;
}

int
OpenSeesAppInit(Tcl_Interp *interp)
{

  // redo puts command so we can capture puts into std:cerr
  Tcl_CmdInfo putsCommandInfo;
  Tcl_GetCommandInfo(interp, "puts", &putsCommandInfo);
  Tcl_putsCommand = putsCommandInfo.objProc;

  // if handle, use our procedure as opposed to theirs
  if (Tcl_putsCommand != nullptr) {
    Tcl_CreateObjCommand(interp, "oldputs", Tcl_putsCommand,      nullptr, nullptr);
    Tcl_CreateObjCommand(interp, "puts",    OpenSees_putsCommand, nullptr, nullptr);
  }


#ifndef _LINUX
  opserr.setFloatField(OPS_Stream::Float::Scientific);
  opserr.setFloatField(OPS_Stream::Float::Fixed);
#endif
  Tcl_Eval(interp, "rename load opensees::import;");
  Tcl_Eval(interp, "interp alias {} import {} opensees::import");


  Tcl_CreateCommand(interp, "logFile",             logFile,      nullptr, nullptr);
  Tcl_CreateCommand(interp, "setPrecision",        setPrecision, nullptr, nullptr);
  Tcl_CreateCommand(interp, "exit",                OpenSeesExit, nullptr, nullptr);
  Tcl_CreateCommand(interp, "quit",                OpenSeesExit, nullptr, nullptr);
  Tcl_CreateCommand(interp, "fault", 
    [](ClientData, Tcl_Interp*, int, G3_Char**)->int{throw 20; return 0;}, nullptr, nullptr);

  // Timer
  Tcl_CreateCommand(interp, "start",               startTimer,   nullptr, nullptr);
  Tcl_CreateCommand(interp, "stop",                stopTimer,    nullptr, nullptr);
  Tcl_CreateCommand(interp, "timer",               timer,        nullptr, nullptr);

  // File utilities
  Tcl_CreateCommand(interp, "stripXML",            stripOpenSeesXML,    nullptr, NULL);
  Tcl_CreateCommand(interp, "convertBinaryToText", convertBinaryToText, nullptr, NULL);
  Tcl_CreateCommand(interp, "convertTextToBinary", convertTextToBinary, nullptr, NULL);
  Tcl_CreateCommand(interp, "setMaxOpenFiles",     maxOpenFiles,        nullptr, nullptr);

  // Some entry points
  Tcl_CreateCommand(interp, "model",               TclCommand_specifyModel,   nullptr, nullptr);
  Tcl_CreateCommand(interp, "opensees::model",     TclCommand_specifyModel,   nullptr, nullptr);
  Tcl_CreateCommand(interp, "wipe",                TclCommand_wipeModel,      nullptr, nullptr);
  Tcl_CreateCommand(interp, "_clearAnalysis",      TclCommand_clearAnalysis,  nullptr, nullptr);

  Tcl_CreateObjCommand(interp, "pset",             OPS_SetObjCmd, nullptr, nullptr);
  Tcl_CreateObjCommand(interp, "source",           OPS_SourceCmd, nullptr, nullptr);
  Tcl_CreateObjCommand(interp, "pragma",           TclObjCommand_pragma, nullptr, nullptr);
  Tcl_CreateObjCommand(interp, "progress",         TclObjCommand_progress, (ClientData)&progress_bar_ptr, nullptr);

  // Tcl_CreateCommand(interp, "searchPeerNGA", &peerNGA, nullptr, nullptr);
  // Tcl_CreateCommand(interp, "defaultUnits",        &defaultUnits, nullptr, NULL);

  return TCL_OK;
}

