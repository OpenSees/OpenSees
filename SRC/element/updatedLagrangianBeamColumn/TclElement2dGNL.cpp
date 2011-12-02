#include <stdlib.h>
#include <string.h>
#include <OPS_Stream.h>

#include <Domain.h>
#include <Node.h>
#include <Matrix.h>

#include <Elastic2DGNL.h>

#include <TclModelBuilder.h>

#define  tcl_debug 1

// Elastic2DGNL(int tag, double A, double E, double I, int Nd1, int Nd2,
//             double rho = 0.0, bool islinear = false);

int
TclModelBuilder_addElastic2dGNL (ClientData clientData, Tcl_Interp *interp,
				 int argc, TCL_Char **argv,
				 Domain *theDomain, TclModelBuilder *theBuilder)
{
	//cerr << "Press key to continue...\n";
	//cin.get();

    if(tcl_debug)
        opserr << " TclModelBuilder_addElastic2dGNL \n";

	if (argc < 8)
	{
		opserr << "WARNING insufficient arguments\n";
		opserr << "element element2dGNL int tag, int Nd1, int Nd2, double A, double E, double Iz, <int linear>\n";

		return TCL_ERROR;
	}

	int tag, ndI, ndJ;
	double E, A, I;
	double massDens = 0.0;
	bool   linear = false;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
	{
		opserr << "WARNING invalid Elastic2dGNL tag" << endln;
		return TCL_ERROR;
	}
    if(tcl_debug) opserr << "\tElement tag = " << tag << "\n";

    if (Tcl_GetInt (interp, argv[3], &ndI) != TCL_OK)
	{
		opserr << "WARNING invalid node I\n";
		opserr << "Elastic2dGNL: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[4], &ndJ) != TCL_OK)
	{
		opserr << "WARNING invalid node J\n";
		opserr << "Elastic2dGNL: " << tag << endln;
		return TCL_ERROR;
	}


	if (Tcl_GetDouble(interp, argv[5], &A) != TCL_OK)
	{
		opserr << "WARNING invalid A\n";
		opserr << "Elastic2dGNL: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[6], &E) != TCL_OK)
	{
		opserr << "WARNING invalid E\n";
		opserr << "Elastic2dGNL: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[7], &I) != TCL_OK)
	{
		opserr << "WARNING invalid I\n";
		opserr << "Elastic2dGNL: " << tag << endln;
		return TCL_ERROR;
	}
	
	if(argc == 9)
	{
		int lin = 0;
		if(Tcl_GetInt(interp, argv[8], &lin) != TCL_OK)
		{
			opserr << "WARNING invalid Linear Flag\n";
			opserr << "Elastic2dGNL: " << tag << endln;
			return TCL_ERROR;
		}
		
		if(lin == 1)
			linear = true;
		
		if(tcl_debug)
			opserr << " 9 arguments - " << lin << endln;
	}
	
	// if(tcl_debug) opserr << "\tAdded upto mass - input parameters\n";


	Element *theElement = new Elastic2dGNL(tag, A, E, I, ndI, ndJ, linear);//, false, massDens);

	if(tcl_debug) opserr << "\tElement created\n";

	// Ensure we have created the element, out of memory if got here and no element
	if (theElement == 0)
	{
		opserr << "WARNING ran out of memory creating element\n";
		opserr << "Elastic2dGNL: " << tag << endln;
		opserr << "\a";
		return TCL_ERROR;
	}

	if (theDomain->addElement(theElement) == false)
	{
		opserr << "WARNING TclElmtBuilder - addElastic2dGNL - could not add element to domain ";
		opserr << tag << endln;
		opserr << "\a";
		return TCL_ERROR;
	}

	if(tcl_debug) opserr << "\tElement number " << tag << " added to domain - returning\n";

	return TCL_OK;
}
