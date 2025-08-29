#include <tcl.h>
#include <ItpackLinSOE.h>
#include <ItpackLinSolver.h>

LinearSOE*
TclDispatch_newItpackLinearSOE(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
    // now must determine the type of solver to create 
    // from rest of args
    int method = 1;
    if (argc == 3) {
      if (Tcl_GetInt(interp, argv[2], &method) != TCL_OK)
        return nullptr;
    }
    ItpackLinSolver *theSolver = new ItpackLinSolver(method);
    return  new ItpackLinSOE(*theSolver);
}



