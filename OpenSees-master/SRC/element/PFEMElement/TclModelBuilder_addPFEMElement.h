#ifndef TclModelBuilder_addPFEMElement_h
#define TclModelBuilder_addPFEMElement_h

#include <tcl.h>
#include <TclModelBuilder.h>
#include <Domain.h>


int
TclModelBuilder_addPFEMElement2D(ClientData clientData, Tcl_Interp *interp,
                                 int argc, TCL_Char **argv,
                                 Domain *theDomain,
                                 TclModelBuilder *theBuilder) ;
int
TclModelBuilder_addPFEMElement3D(ClientData clientData, Tcl_Interp *interp,
                                 int argc, TCL_Char **argv,
                                 Domain *theDomain,
                                 TclModelBuilder *theBuilder) ;

#endif 
