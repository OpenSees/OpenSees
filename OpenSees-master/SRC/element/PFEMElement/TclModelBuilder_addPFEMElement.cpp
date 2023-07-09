#include "TclModelBuilder_addPFEMElement.h"
#include "PFEMElement2D.h"
#include "PFEMElement2DFIC.h"
#include "PFEMElement2DCompressible.h"
#include "PFEMElement2DBubble.h"
#include "PFEMElement2Dmini.h"
#include "PFEMElement3D.h"
#include <cstring>
#include <string>

int
TclModelBuilder_addPFEMElement2D(ClientData clientData, Tcl_Interp *interp,
                                 int argc, TCL_Char **argv,
                                 Domain *theDomain,
                                 TclModelBuilder *theBuilder)
{
    // define PFEM elements
    if(argc < 11) {
        opserr << "Invalid #args: want element PFEMElement2D tag";
        opserr << "nd1 nd2 nd3 type rho mu b1 b2 <thickness kappa>\n";
        return TCL_ERROR;
    }
    int loc = 2;
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
    int type = -1;
    if(strcmp(argv[loc], "FIC")==0 || strcmp(argv[loc], "fic")==0
       || strcmp(argv[loc], "Fic")==0) {
        type = 0;
    } else if(strcmp(argv[loc], "quasi-incompressible")==0
              || strcmp(argv[loc], "Quasi-Incompressible")==0) {
        type = 1;
    } else if(strcmp(argv[loc], "bubble")==0
              || strcmp(argv[loc], "Bubble")==0) {
        type = 2;
    } else if(strcmp(argv[loc], "mini")==0
              || strcmp(argv[loc], "Mini")==0) {
        type = 3;
    } else {
        opserr<<"WARNNG: unknown type for PFEMElement2D \n";
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
    loc++;
    int lumped = 0;
    if(argc > loc) {
        if(Tcl_GetInt(interp, argv[loc], &lumped) != TCL_OK) {
            opserr<< "WARNING invalid lumped "<< argv[loc] << ": element PFEMElement2D\n";
            return TCL_ERROR;
        }
    }
    loc++;
    int checkJ = 0;
    if(argc > loc) {
        if(Tcl_GetInt(interp, argv[loc], &checkJ) != TCL_OK) {
            opserr<< "WARNING invalid checkJ "<< argv[loc] << ": element PFEMElement2D\n";
            return TCL_ERROR;
        }
    }
    loc++;

    // regular element
    Element* ele = 0;
    if(type == 0) {
        ele = new PFEMElement2DFIC(tag,nd1,nd2,nd3,rho,mu,b1,b2,thickness);
    } else if(type == 1) {
        //ele = new PFEMElement2DCompressible(tag,nd1,nd2,nd3,rho,mu,b1,b2,thickness,kappa);
    } else if(type == 2) {
        //ele = new PFEMElement2DBubble(tag,nd1,nd2,nd3,rho,mu,b1,b2,thickness,kappa);
    } else if(type == 3) {
        //ele = new PFEMElement2Dmini(tag,nd1,nd2,nd3,rho,mu,b1,b2,thickness,kappa,lumped,checkJ);
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
    if(argc < 13) {
        opserr << "Invalid #args: want element PFEMElement3D ";
        opserr << "tag nd1 nd2 nd3 nd4 type rho mu b1 b2 b3 <kappa lumped checkJ>\n";
        return TCL_ERROR;
    }

    int loc = 2;

    int tag;
    if(Tcl_GetInt(interp, argv[loc], &tag) != TCL_OK) {
        opserr<< "WARNING invalid tag "<< argv[loc] << ": element PFEMElement3D\n";
        return TCL_ERROR;
    }
    loc++;
    int nd1;
    if(Tcl_GetInt(interp, argv[loc], &nd1) != TCL_OK) {
        opserr<< "WARNING invalid nd1 "<< argv[loc] << ": element PFEMElement3D\n";
        return TCL_ERROR;
    }
    loc++;
    int nd2;
    if(Tcl_GetInt(interp, argv[loc], &nd2) != TCL_OK) {
        opserr<< "WARNING invalid nd2 "<< argv[loc] << ": element PFEMElement3D\n";
        return TCL_ERROR;
    }
    loc++;
    int nd3;
    if(Tcl_GetInt(interp, argv[loc], &nd3) != TCL_OK) {
        opserr<< "WARNING invalid nd3 "<< argv[loc] << ": element PFEMElement3D\n";
        return TCL_ERROR;
    }
    loc++;
    int nd4;
    if(Tcl_GetInt(interp, argv[loc], &nd4) != TCL_OK) {
        opserr<< "WARNING invalid nd4 "<< argv[loc] << ": element PFEMElement3D\n";
        return TCL_ERROR;
    }
    loc++;
    std::string type = argv[loc];
    loc++;
    double rho;
    if(Tcl_GetDouble(interp, argv[loc], &rho) != TCL_OK) {
        opserr<< "WARNING invalid rho "<< argv[loc] << ": element PFEMElement3D\n";
        return TCL_ERROR;
    }
    loc++;
    double mu;
    if(Tcl_GetDouble(interp, argv[loc], &mu) != TCL_OK) {
        opserr<< "WARNING invalid mu "<< argv[loc] << ": element PFEMElement3D\n";
        return TCL_ERROR;
    }
    loc++;
    double b1;
    if(Tcl_GetDouble(interp, argv[loc], &b1) != TCL_OK) {
        opserr<< "WARNING invalid b1 "<< argv[loc] << ": element PFEMElement3D\n";
        return TCL_ERROR;
    }
    loc++;
    double b2;
    if(Tcl_GetDouble(interp, argv[loc], &b2) != TCL_OK) {
        opserr<< "WARNING invalid b2 "<< argv[loc] << ": element PFEMElement3D\n";
        return TCL_ERROR;
    }
    loc++;
    double b3;
    if(Tcl_GetDouble(interp, argv[loc], &b3) != TCL_OK) {
        opserr<< "WARNING invalid b3 "<< argv[loc] << ": element PFEMElement3D\n";
        return TCL_ERROR;
    }
    loc++;
    double kappa = 2.2e9;
    if(argc > loc) {
	if(Tcl_GetDouble(interp, argv[loc], &kappa) != TCL_OK) {
	    opserr<<"WARNING invalid kappa "<<argv[loc]<<": element PFEMElement3D\n";
	    return TCL_ERROR;
	}
	loc++;
    }
    int lumped = 0;
    if(argc > loc) {
	if(Tcl_GetInt(interp, argv[loc], &lumped) != TCL_OK) {
	    opserr<<"WARNING invalid lumped "<<argv[loc]<<": element PFEMElement3D\n";
	    return TCL_ERROR;
	}
	loc++;
    }
    int check = 0;
    if(argc > loc) {
	if(Tcl_GetInt(interp, argv[loc], &check) != TCL_OK) {
	    opserr<<"WARNING invalid check "<<argv[loc]<<": element PFEMElement3D\n";
	    return TCL_ERROR;
	}
	loc++;
    }

    // regular element
    Element* ele = 0;
    if(type == "pressuregradient" || type == "PressureGradient") {
        ele = new PFEMElement3D(tag,nd1,nd2,nd3,nd4,rho,mu,b1,b2,b3);
    } else {
	opserr<<"element PFEMElement3D type "<<type.c_str()<<" is not known\n";
	return TCL_ERROR;
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
