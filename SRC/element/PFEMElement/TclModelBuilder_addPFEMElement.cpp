#include "TclModelBuilder_addPFEMElement.h"
#include "PFEMElement2D.h"
#include <cstring>
#include <Pressure_Constraint.h>

int
TclModelBuilder_addPFEMElement2D(ClientData clientData, Tcl_Interp *interp,
                                 int argc, TCL_Char **argv,
                                 Domain *theDomain,
                                 TclModelBuilder *theBuilder) 
{
    // define PFEM elements
    if(argc < 10) {
        opserr << "Invalid #args: want element PFEMElement2D ";
        opserr << "tag nd1 nd2 nd3 rho mu b1 b2 \n";
        return TCL_ERROR;
    }
    int tag;
    if(Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
        opserr<< "WARNING invalid tag "<< argv[2] << ": element PFEMElement2D\n";
        return TCL_ERROR;
    }
    int nd1;
    if(Tcl_GetInt(interp, argv[3], &nd1) != TCL_OK) {
        opserr<< "WARNING invalid nd1 "<< argv[3] << ": element PFEMElement2D\n";
        return TCL_ERROR;
    }
    int nd2;
    if(Tcl_GetInt(interp, argv[4], &nd2) != TCL_OK) {
        opserr<< "WARNING invalid nd2 "<< argv[4] << ": element PFEMElement2D\n";
        return TCL_ERROR;
    }
    int nd3;
    if(Tcl_GetInt(interp, argv[5], &nd3) != TCL_OK) {
        opserr<< "WARNING invalid nd3 "<< argv[5] << ": element PFEMElement2D\n";
        return TCL_ERROR;
    }
    double rho;
    if(Tcl_GetDouble(interp, argv[6], &rho) != TCL_OK) {
        opserr<< "WARNING invalid rho "<< argv[6] << ": element PFEMElement2D\n";
        return TCL_ERROR;
    }
    double mu;
    if(Tcl_GetDouble(interp, argv[7], &mu) != TCL_OK) {
        opserr<< "WARNING invalid mu "<< argv[7] << ": element PFEMElement2D\n";
        return TCL_ERROR;
    }
    double b1;
    if(Tcl_GetDouble(interp, argv[8], &b1) != TCL_OK) {
        opserr<< "WARNING invalid b1 "<< argv[8] << ": element PFEMElement2D\n";
        return TCL_ERROR;
    }
    double b2;
    if(Tcl_GetDouble(interp, argv[9], &b2) != TCL_OK) {
        opserr<< "WARNING invalid b2 "<< argv[9] << ": element PFEMElement2D\n";
        return TCL_ERROR;
    }

    // regular element
    PFEMElement2D* ele = new PFEMElement2D(tag,nd1,nd2,nd3,rho,mu,b1,b2);
    if (ele == 0) {
        opserr << "WARNING ran out of memory creating element\n";
        opserr << "element: " << tag << endln;
        return TCL_ERROR;
    }

    if (theDomain->addElement(ele) == false) {
        opserr << "WARNING failed to add element to the domain\n";
        opserr << "element: " << tag << endln;
        delete ele; // otherwise memory leak
        return TCL_ERROR;
    }


    return TCL_OK;
}
