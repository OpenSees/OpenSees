#include "TclModelBuilder_addPFEMElement.h"
#include "PFEMElement2D.h"
#include "PFEMElement2DCompressible.h"
#include "PFEMElement2DBubble.h"
#include "PFEMElement3D.h"
#include <cstring>

int
TclModelBuilder_addPFEMElement2D(ClientData clientData, Tcl_Interp *interp,
                                 int argc, TCL_Char **argv,
                                 Domain *theDomain,
                                 TclModelBuilder *theBuilder) 
{
    // define PFEM elements
    if(argc < 11) {
        opserr << "Invalid #args: want element PFEMElement2D type";
        opserr << "tag nd1 nd2 nd3 rho mu b1 b2 <thickness kappa>\n";
        return TCL_ERROR;
    }
    int loc = 2;
    int type = -1;
    if(strcmp(argv[loc], "pressureGradient")==0 || strcmp(argv[loc], "pressuregradient")==0
       || strcmp(argv[loc], "PressureGradient")==0) {
        type = 0;
    } else if(strcmp(argv[loc], "quasi-incompressible")==0 
              || strcmp(argv[loc], "Quasi-Incompressible")==0) {
        type = 1;
    } else if(strcmp(argv[loc], "bubble")==0 
              || strcmp(argv[loc], "Bubble")==0) {
        type = 2;
    } else {
        opserr<<"WARNNG: unknown type for PFEMElement2D \n";
        return TCL_ERROR;
    }
    loc++;
    int tag;
    if(Tcl_GetInt(interp, argv[loc], &tag) != TCL_OK) {
        opserr<< "WARNING invalid tag "<< argv[loc] << ": element PFEMElement2D\n";
        return TCL_ERROR;
    }
    loc++;
    int nd1;
    if(Tcl_GetInt(interp, argv[loc], &nd1) != TCL_OK) {
        opserr<< "WARNING invalid nd1 "<< argv[loc] << ": element PFEMElement2D\n";
        return TCL_ERROR;
    }
    loc++;
    int nd2;
    if(Tcl_GetInt(interp, argv[loc], &nd2) != TCL_OK) {
        opserr<< "WARNING invalid nd2 "<< argv[loc] << ": element PFEMElement2D\n";
        return TCL_ERROR;
    }
    loc++;
    int nd3;
    if(Tcl_GetInt(interp, argv[loc], &nd3) != TCL_OK) {
        opserr<< "WARNING invalid nd3 "<< argv[loc] << ": element PFEMElement2D\n";
        return TCL_ERROR;
    }
    loc++;
    double rho;
    if(Tcl_GetDouble(interp, argv[loc], &rho) != TCL_OK) {
        opserr<< "WARNING invalid rho "<< argv[loc] << ": element PFEMElement2D\n";
        return TCL_ERROR;
    }
    loc++;
    double mu;
    if(Tcl_GetDouble(interp, argv[loc], &mu) != TCL_OK) {
        opserr<< "WARNING invalid mu "<< argv[loc] << ": element PFEMElement2D\n";
        return TCL_ERROR;
    }
    loc++;
    double b1;
    if(Tcl_GetDouble(interp, argv[loc], &b1) != TCL_OK) {
        opserr<< "WARNING invalid b1 "<< argv[loc] << ": element PFEMElement2D\n";
        return TCL_ERROR;
    }
    loc++;
    double b2;
    if(Tcl_GetDouble(interp, argv[loc], &b2) != TCL_OK) {
        opserr<< "WARNING invalid b2 "<< argv[loc] << ": element PFEMElement2D\n";
        return TCL_ERROR;
    }
    loc++;
    double thickness = 1.0;
    if(argc > loc) {
        if(Tcl_GetDouble(interp, argv[loc], &thickness) != TCL_OK) {
            opserr<< "WARNING invalid thickness "<< argv[loc] << ": element PFEMElement2D\n";
            return TCL_ERROR;
        }
    }
    loc++;
    double kappa = 2.2e9;
    if(argc > loc) {
        if(Tcl_GetDouble(interp, argv[loc], &kappa) != TCL_OK) {
            opserr<< "WARNING invalid kappa "<< argv[loc] << ": element PFEMElement2D\n";
            return TCL_ERROR;
        }
    }

    // regular element
    Element* ele = 0;
    if(type == 0) {
        ele = new PFEMElement2D(tag,nd1,nd2,nd3,rho,mu,b1,b2,thickness);
    } else if(type == 1) {
        ele = new PFEMElement2DCompressible(tag,nd1,nd2,nd3,rho,mu,b1,b2,thickness,kappa);
    } else if(type == 2) {
        ele = new PFEMElement2DBubble(tag,nd1,nd2,nd3,rho,mu,b1,b2,thickness,kappa);
    }
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

int
TclModelBuilder_addPFEMElement3D(ClientData clientData, Tcl_Interp *interp,
                                 int argc, TCL_Char **argv,
                                 Domain *theDomain,
                                 TclModelBuilder *theBuilder) 
{
    // define PFEM elements
    if(argc < 12) {
        opserr << "Invalid #args: want element PFEMElement3D ";
        opserr << "tag nd1 nd2 nd3 nd4 rho mu b1 b2 b3\n";
        return TCL_ERROR;
    }
    int tag;
    if(Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
        opserr<< "WARNING invalid tag "<< argv[2] << ": element PFEMElement3D\n";
        return TCL_ERROR;
    }
    int nd1;
    if(Tcl_GetInt(interp, argv[3], &nd1) != TCL_OK) {
        opserr<< "WARNING invalid nd1 "<< argv[3] << ": element PFEMElement3D\n";
        return TCL_ERROR;
    }
    int nd2;
    if(Tcl_GetInt(interp, argv[4], &nd2) != TCL_OK) {
        opserr<< "WARNING invalid nd2 "<< argv[4] << ": element PFEMElement3D\n";
        return TCL_ERROR;
    }
    int nd3;
    if(Tcl_GetInt(interp, argv[5], &nd3) != TCL_OK) {
        opserr<< "WARNING invalid nd3 "<< argv[5] << ": element PFEMElement3D\n";
        return TCL_ERROR;
    }
    int nd4;
    if(Tcl_GetInt(interp, argv[6], &nd4) != TCL_OK) {
        opserr<< "WARNING invalid nd4 "<< argv[6] << ": element PFEMElement3D\n";
        return TCL_ERROR;
    }
    double rho;
    if(Tcl_GetDouble(interp, argv[6], &rho) != TCL_OK) {
        opserr<< "WARNING invalid rho "<< argv[7] << ": element PFEMElement3D\n";
        return TCL_ERROR;
    }
    double mu;
    if(Tcl_GetDouble(interp, argv[8], &mu) != TCL_OK) {
        opserr<< "WARNING invalid mu "<< argv[8] << ": element PFEMElement3D\n";
        return TCL_ERROR;
    }
    double b1;
    if(Tcl_GetDouble(interp, argv[9], &b1) != TCL_OK) {
        opserr<< "WARNING invalid b1 "<< argv[9] << ": element PFEMElement3D\n";
        return TCL_ERROR;
    }
    double b2;
    if(Tcl_GetDouble(interp, argv[10], &b2) != TCL_OK) {
        opserr<< "WARNING invalid b2 "<< argv[10] << ": element PFEMElement3D\n";
        return TCL_ERROR;
    }
    double b3;
    if(Tcl_GetDouble(interp, argv[11], &b3) != TCL_OK) {
        opserr<< "WARNING invalid b3 "<< argv[11] << ": element PFEMElement3D\n";
        return TCL_ERROR;
    }

    // regular element
    PFEMElement3D* ele = new PFEMElement3D(tag,nd1,nd2,nd3,nd4,rho,mu,b1,b2,b3);
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
