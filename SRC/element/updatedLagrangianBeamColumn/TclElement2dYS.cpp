#include <stdlib.h>
#include <string.h>
#include <OPS_Stream.h>

#include <Domain.h>
#include <Node.h>
#include <Matrix.h>

#include <CyclicModel.h>
#include <Inelastic2DYS01.h>
#include <Inelastic2DYS02.h>
#include <Inelastic2DYS03.h>
//#include <Inelastic2DYS04.h>
//#include <Inelastic2DYS05.h>

#include <YieldSurface_BC.h>
#include <TclModelBuilder.h>

#define  tcl_debug 0

// Element2dGNL(int tag, double A, double E, double I, int Nd1, int Nd2,
//             double rho = 0.0, bool islinear = false);

int
TclModelBuilder_addElement2dYS01 (ClientData clientData, Tcl_Interp *interp,
								   int argc, TCL_Char **argv,
								   Domain *theDomain, TclModelBuilder *theBuilder)
{
	//cerr << "Press key to continue...\n";
	//cin.get();

    if(tcl_debug)
        opserr << " TclModelBuilder_addElement2dGNL \n";

	if (argc < 11)
	{
		opserr << "WARNING insufficient arguments\n";
		opserr << "element element2dYS tag? Nd1? Nd2? A? E? Iz? ysID1? ysID2? algo?";

		return TCL_ERROR;
	}

	int tag, ndI, ndJ;
	double E, A, I;
//	double massDens = 0.0;
	int ysID1, ysID2;
	int rf_algo;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
	{
		opserr << "WARNING invalid element2dYS tag" << endln;
		return TCL_ERROR;
	}
    if(tcl_debug) opserr << "\tElement tag = " << tag << "\n";

    if (Tcl_GetInt (interp, argv[3], &ndI) != TCL_OK)
	{
		opserr << "WARNING invalid node I\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[4], &ndJ) != TCL_OK)
	{
		opserr << "WARNING invalid node J\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}


	if (Tcl_GetDouble(interp, argv[5], &A) != TCL_OK)
	{
		opserr << "WARNING invalid A\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[6], &E) != TCL_OK)
	{
		opserr << "WARNING invalid E\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[7], &I) != TCL_OK)
	{
		opserr << "WARNING invalid I\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[8], &ysID1) != TCL_OK)
	{
		opserr << "WARNING invalid ysID1\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[9], &ysID2) != TCL_OK)
	{
		opserr << "WARNING invalid ysID2\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[10], &rf_algo) != TCL_OK)
	{
		opserr << "WARNING invalid ysID1\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}

	YieldSurface_BC *theYS1 = theBuilder->getYieldSurface_BC(ysID1);
	if(theYS1 == 0)
	{
		opserr << "WARNING element2dYS: " << tag << "\n";
		opserr <<  " no yield surface exists with tag: " << ysID1 << endln;
		return TCL_ERROR;
	}

	YieldSurface_BC *theYS2 = theBuilder->getYieldSurface_BC(ysID2);
	if(theYS2 == 0)
	{
		opserr << "WARNING element2dYS: " << tag << "\n";
		opserr <<  " no yield surface exists with tag: " << ysID2 << endln;
		return TCL_ERROR;
	}

// 		Inelastic2DYS(	int tag, double A, double E, double I, int Nd1, int Nd2,
// 						YieldSurface_BC *ysEnd1, YieldSurface_BC *ysEnd2,
// 						int rf_algo = -1, // updated
	Element *theElement = new Inelastic2DYS01(tag, A, E, I, ndI, ndJ, theYS1, theYS2, rf_algo);

	if(tcl_debug) opserr << "\tElement created\n";

	// Ensure we have created the element, out of memory if got here and no element
	if (theElement == 0)
	{
		opserr << "WARNING ran out of memory creating element\n";
		opserr << "element2dYS: " << tag << endln;
		opserr << "\a";
		return TCL_ERROR;
	}

	if (theDomain->addElement(theElement) == false)
	{
		opserr << "WARNING TclElmtBuilder - addelement2dYS - could not add element to domain ";
		opserr << tag << endln;
		opserr << "\a";
		return TCL_ERROR;
	}

	if(tcl_debug) opserr << "\tElement number " << tag << " added to domain - returning\n";

	return TCL_OK;
}


int
TclModelBuilder_addElement2dYS02 (ClientData clientData, Tcl_Interp *interp,
								   int argc, TCL_Char **argv,
								   Domain *theDomain, TclModelBuilder *theBuilder)
{
	//cerr << "Press key to continue...\n";
	//cin.get();

    if(tcl_debug)
        opserr << " TclModelBuilder_addElement2dGNL \n";

	if (argc < 14)
	{
		opserr << "WARNING insufficient arguments\n";
		opserr << "element element2dYS tag? Nd1? Nd2? A? E? Iz? ysID1? ysID2? cycType? wt? power? algo?";

		return TCL_ERROR;
	}

	int tag, ndI, ndJ;
	double E, A, I;
//	double massDens = 0.0;
	int ysID1, ysID2;
	int cyc_type;
//	double wt;

	int rf_algo=-1;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
	{
		opserr << "WARNING invalid element2dYS tag" << endln;
		return TCL_ERROR;
	}
    if(tcl_debug) opserr << "\tElement tag = " << tag << "\n";

    if (Tcl_GetInt (interp, argv[3], &ndI) != TCL_OK)
	{
		opserr << "WARNING invalid node I\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[4], &ndJ) != TCL_OK)
	{
		opserr << "WARNING invalid node J\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}


	if (Tcl_GetDouble(interp, argv[5], &A) != TCL_OK)
	{
		opserr << "WARNING invalid A\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[6], &E) != TCL_OK)
	{
		opserr << "WARNING invalid E\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[7], &I) != TCL_OK)
	{
		opserr << "WARNING invalid I\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[8], &ysID1) != TCL_OK)
	{
		opserr << "WARNING invalid ysID1\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[9], &ysID2) != TCL_OK)
	{
		opserr << "WARNING invalid ysID2\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}

		if (Tcl_GetInt(interp, argv[10], &cyc_type) != TCL_OK)
	{
		opserr << "WARNING invalid cyc_type\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}

	/*if (Tcl_GetDouble (interp, argv[11], &wt) != TCL_OK)
	{
		opserr << "WARNING invalid wt\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}*/

	double delpmax, alfa, beta;
	if (Tcl_GetDouble (interp, argv[11], &delpmax) != TCL_OK)
	{
		opserr << "WARNING invalid power\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}
	
	if (Tcl_GetDouble (interp, argv[12], &alfa) != TCL_OK)
	{
		opserr << "WARNING invalid power\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble (interp, argv[13], &beta) != TCL_OK)
	{
		opserr << "WARNING invalid rfalgo\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}

	YieldSurface_BC *theYS1 = theBuilder->getYieldSurface_BC(ysID1);
	if(theYS1 == 0)
	{
		opserr << "WARNING element2dYS: " << tag << "\n";
		opserr <<  " no yield surface exists with tag: " << ysID1 << endln;
		return TCL_ERROR;
	}

	YieldSurface_BC *theYS2 = theBuilder->getYieldSurface_BC(ysID2);
	if(theYS2 == 0)
	{
		opserr << "WARNING element2dYS: " << tag << "\n";
		opserr <<  " no yield surface exists with tag: " << ysID2 << endln;
		return TCL_ERROR;
	}

//Inelastic2DYS02(int tag, double a, double e, double i, int Nd1, int Nd2,
//				YieldSurface_BC *ysEnd1,  YieldSurface_BC *ysEnd2,
//				int rf_algo, bool islinear, double rho)

	CyclicModel *theModel = theBuilder->getCyclicModel(cyc_type);
//Element *theElement = new Inelastic2DYS02(tag, A, E, I, ndI, ndJ, theYS1, theYS2, cyc_type, wt, delpmax, alfa, beta, rf_algo);
Element *theElement = new Inelastic2DYS02(tag, A, E, I, ndI, ndJ, theYS1, theYS2, theModel, delpmax, alfa, beta, rf_algo);
    opserr << "Inelastic2DYS02 created\n";

	if(tcl_debug) opserr << "\tElement created\n";

	// Ensure we have created the element, out of memory if got here and no element
	if (theElement == 0)
	{
		opserr << "WARNING ran out of memory creating element\n";
		opserr << "element2dYS: " << tag << endln;
		opserr << "\a";
		return TCL_ERROR;
	}

	opserr << "Inelastic2DYS02 adding to domain\n";

	if (theDomain->addElement(theElement) == false)
	{
		opserr << "WARNING TclElmtBuilder - addelement2dYS - could not add element to domain ";
		opserr << tag << endln;
		opserr << "\a";
		return TCL_ERROR;
	}

	//if(tcl_debug)
		opserr << "Inelastic2DYS02 #" << tag << " added to domain - returning\n";

	return TCL_OK;
}


int
TclModelBuilder_addElement2dYS03 (ClientData clientData, Tcl_Interp *interp,
								   int argc, TCL_Char **argv,
								   Domain *theDomain, TclModelBuilder *theBuilder)
{
	//cerr << "Press key to continue...\n";
	//cin.get();

    if(tcl_debug)
        opserr << " TclModelBuilder_addElement2dGNL \n";

	if (argc < 11)
	{
		opserr << "WARNING insufficient arguments\n";
		opserr << "element element2dYS03 tag? Nd1? Nd2? A_ten? A_com? E? IzPos? IzNeg? ysID1? ysID2? algo?";

		return TCL_ERROR;
	}

	int tag, ndI, ndJ;
	double E, aTens, aComp, Ipos, Ineg;
//	double massDens = 0.0;
	int ysID1, ysID2;

	int rf_algo;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
	{
		opserr << "WARNING invalid element2dYS tag" << endln;
		return TCL_ERROR;
	}
    if(tcl_debug) opserr << "\tElement tag = " << tag << "\n";

    if (Tcl_GetInt (interp, argv[3], &ndI) != TCL_OK)
	{
		opserr << "WARNING invalid node I\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[4], &ndJ) != TCL_OK)
	{
		opserr << "WARNING invalid node J\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}


	if (Tcl_GetDouble(interp, argv[5], &aTens) != TCL_OK)
	{
		opserr << "WARNING invalid A\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[6], &aComp) != TCL_OK)
	{
		opserr << "WARNING invalid A\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[7], &E) != TCL_OK)
	{
		opserr << "WARNING invalid E\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[8], &Ipos) != TCL_OK)
	{
		opserr << "WARNING invalid I\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[9], &Ineg) != TCL_OK)
	{
		opserr << "WARNING invalid I\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[10], &ysID1) != TCL_OK)
	{
		opserr << "WARNING invalid ysID1\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[11], &ysID2) != TCL_OK)
	{
		opserr << "WARNING invalid ysID2\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[12], &rf_algo) != TCL_OK)
	{
		opserr << "WARNING invalid ysID1\n";
		opserr << "element2dYS: " << tag << endln;
		return TCL_ERROR;
	}

	YieldSurface_BC *theYS1 = theBuilder->getYieldSurface_BC(ysID1);
	if(theYS1 == 0)
	{
		opserr << "WARNING element2dYS: " << tag << "\n";
		opserr <<  " no yield surface exists with tag: " << ysID1 << endln;
		return TCL_ERROR;
	}

	YieldSurface_BC *theYS2 = theBuilder->getYieldSurface_BC(ysID2);
	if(theYS2 == 0)
	{
		opserr << "WARNING element2dYS: " << tag << "\n";
		opserr <<  " no yield surface exists with tag: " << ysID2 << endln;
		return TCL_ERROR;
	}

//	Inelastic2DYS03(int tag, double a_ten, double a_com, double e,
//	                double iz_pos, double iz_neg, int Nd1, int Nd2,
//                    YieldSurface_BC *ysEnd1,  YieldSurface_BC *ysEnd2,
//                    int rf_algo, bool islinear, double rho);

Element *theElement = new Inelastic2DYS03(tag, aTens, aComp, E,
                                          Ipos, Ineg, ndI, ndJ,
                                          theYS1, theYS2, rf_algo);

    opserr << "Inelastic2DYS03 created\n";

	if(tcl_debug) opserr << "\tElement created\n";

	// Ensure we have created the element, out of memory if got here and no element
	if (theElement == 0)
	{
		opserr << "WARNING ran out of memory creating element\n";
		opserr << "element2dYS: " << tag << endln;
		opserr << "\a";
		return TCL_ERROR;
	}

	opserr << "Inelastic2DYS03 adding to domain\n";

	if (theDomain->addElement(theElement) == false)
	{
		opserr << "WARNING TclElmtBuilder - addelement2dYS - could not add element to domain ";
		opserr << tag << endln;
		opserr << "\a";
		return TCL_ERROR;
	}

	if(tcl_debug)
		opserr << "Inelastic2DYS03 #" << tag << " added to domain - returning\n";

	return TCL_OK;
}


/*
int
TclModelBuilder_addElement2dYS04 (ClientData clientData, Tcl_Interp *interp,
								   int argc, char **argv,
								   Domain *theDomain, TclModelBuilder *theBuilder)
{
	//cerr << "Press key to continue...\n";
	//cin.get();

	if (argc < 11)
	{
		opserr << "WARNING insufficient arguments\n";
		opserr << "element element2dYS04 tag? Nd1? Nd2? A? E? Iz? ysID1? ysID2? algo?";

		return TCL_ERROR;
	}

	int tag, ndI, ndJ;
	double E, A, I;
//	double massDens = 0.0;
	int ysID1, ysID2;
	int rf_algo;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
	{
		opserr << "WARNING invalid element2dYS04 tag" << endln;
		return TCL_ERROR;
	}
    if(tcl_debug) opserr << "\tElement tag = " << tag << "\n";

    if (Tcl_GetInt (interp, argv[3], &ndI) != TCL_OK)
	{
		opserr << "WARNING invalid node I\n";
		opserr << "element2dYS04: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[4], &ndJ) != TCL_OK)
	{
		opserr << "WARNING invalid node J\n";
		opserr << "element2dYS04: " << tag << endln;
		return TCL_ERROR;
	}


	if (Tcl_GetDouble(interp, argv[5], &A) != TCL_OK)
	{
		opserr << "WARNING invalid A\n";
		opserr << "element2dYS04: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[6], &E) != TCL_OK)
	{
		opserr << "WARNING invalid E\n";
		opserr << "element2dYS04: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[7], &I) != TCL_OK)
	{
		opserr << "WARNING invalid I\n";
		opserr << "element2dYS04: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[8], &ysID1) != TCL_OK)
	{
		opserr << "WARNING invalid ysID1\n";
		opserr << "element2dYS04: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[9], &ysID2) != TCL_OK)
	{
		opserr << "WARNING invalid ysID2\n";
		opserr << "element2dYS04: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[10], &rf_algo) != TCL_OK)
	{
		opserr << "WARNING invalid ysID1\n";
		opserr << "element2dYS04: " << tag << endln;
		return TCL_ERROR;
	}

		YieldSurface_BC *theYS1 = theBuilder->getYieldSurface_BC(ysID1);
	if(theYS1 == 0)
	{
		opserr << "WARNING element2dYS: " << tag << "\n";
		opserr <<  " no yield surface exists with tag: " << ysID1 << endln;
		return TCL_ERROR;
	}

	YieldSurface_BC *theYS2 = theBuilder->getYieldSurface_BC(ysID2);
	if(theYS2 == 0)
	{
		opserr << "WARNING element2dYS: " << tag << "\n";
		opserr <<  " no yield surface exists with tag: " << ysID2 << endln;
		return TCL_ERROR;
	}

// 		Inelastic2DYS(	int tag, double A, double E, double I, int Nd1, int Nd2,
// 						YieldSurface_BC *ysEnd1, YieldSurface_BC *ysEnd2,
// 						int rf_algo = -1, // updated

	Element *theElement = new Inelastic2DYS04(tag, A, E, I, ndI, ndJ, theYS1, theYS2, rf_algo);

	if(tcl_debug) opserr << "\tElement created\n";

	// Ensure we have created the element, out of memory if got here and no element
	if (theElement == 0)
	{
		opserr << "WARNING ran out of memory creating element\n";
		opserr << "element2dYS04: " << tag << endln;
		opserr << "\a";
		return TCL_ERROR;
	}

	if (theDomain->addElement(theElement) == false)
	{
		opserr << "WARNING TclElmtBuilder - addelement2dYS - could not add element to domain ";
		opserr << tag << endln;
		opserr << "\a";
		return TCL_ERROR;
	}

	if(tcl_debug) opserr << "\tElement number " << tag << " added to domain - returning\n";

	return TCL_OK;
}
*/

/*
int
TclModelBuilder_addElement2dYS05 (ClientData clientData, Tcl_Interp *interp,
								   int argc, char **argv,
								   Domain *theDomain, TclModelBuilder *theBuilder)
{
	//cerr << "Press key to continue...\n";
	//cin.get();

	if (argc < 11)
	{
		opserr << "WARNING insufficient arguments\n";
		opserr << "element element2dYS04 tag? Nd1? Nd2? A? E? Iz? ysID1? ysID2? algo?";

		return TCL_ERROR;
	}

	int tag, ndI, ndJ;
	double E, A, I;
//	double massDens = 0.0;
	int ysID1, ysID2;
	int rf_algo;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
	{
		opserr << "WARNING invalid element2dYS05 tag" << endln;
		return TCL_ERROR;
	}
    if(tcl_debug) opserr << "\tElement tag = " << tag << "\n";

    if (Tcl_GetInt (interp, argv[3], &ndI) != TCL_OK)
	{
		opserr << "WARNING invalid node I\n";
		opserr << "element2dYS05: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[4], &ndJ) != TCL_OK)
	{
		opserr << "WARNING invalid node J\n";
		opserr << "element2dYS05: " << tag << endln;
		return TCL_ERROR;
	}


	if (Tcl_GetDouble(interp, argv[5], &A) != TCL_OK)
	{
		opserr << "WARNING invalid A\n";
		opserr << "element2dYS05: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[6], &E) != TCL_OK)
	{
		opserr << "WARNING invalid E\n";
		opserr << "element2dYS05: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[7], &I) != TCL_OK)
	{
		opserr << "WARNING invalid I\n";
		opserr << "element2dYS05: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[8], &ysID1) != TCL_OK)
	{
		opserr << "WARNING invalid ysID1\n";
		opserr << "element2dYS05: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[9], &ysID2) != TCL_OK)
	{
		opserr << "WARNING invalid ysID2\n";
		opserr << "element2dYS05: " << tag << endln;
		return TCL_ERROR;
	}

	if (Tcl_GetInt (interp, argv[10], &rf_algo) != TCL_OK)
	{
		opserr << "WARNING invalid ysID1\n";
		opserr << "element2dYS05: " << tag << endln;
		return TCL_ERROR;
	}

	YieldSurface_BC *theYS1 = theBuilder->getYieldSurface_BC(ysID1);
	if(theYS1 == 0)
	{
		opserr << "WARNING element2dYS05: " << tag << "\n";
		opserr <<  " no yield surface exists with tag: " << ysID1 << endln;
		return TCL_ERROR;
	}

	YieldSurface_BC *theYS2 = theBuilder->getYieldSurface_BC(ysID2);
	if(theYS2 == 0)
	{
		opserr << "WARNING element2dYS05: " << tag << "\n";
		opserr <<  " no yield surface exists with tag: " << ysID2 << endln;
		return TCL_ERROR;
	}

// 		Inelastic2DYS(	int tag, double A, double E, double I, int Nd1, int Nd2,
// 						YieldSurface_BC *ysEnd1, YieldSurface_BC *ysEnd2,
// 						int rf_algo = -1, // updated

	Element *theElement = new Inelastic2DYS05(tag, A, E, I, ndI, ndJ, theYS1, theYS2, rf_algo);

	if(tcl_debug) opserr << "\tElement created\n";

	// Ensure we have created the element, out of memory if got here and no element
	if (theElement == 0)
	{
		opserr << "WARNING ran out of memory creating element\n";
		opserr << "element2dYS05: " << tag << endln;
		opserr << "\a";
		return TCL_ERROR;
	}

	if (theDomain->addElement(theElement) == false)
	{
		opserr << "WARNING TclElmtBuilder - addelement2dYS - could not add element to domain ";
		opserr << tag << endln;
		opserr << "\a";
		return TCL_ERROR;
	}

	if(tcl_debug) opserr << "\tElement number " << tag << " added to domain - returning\n";

	return TCL_OK;
}
*/

/*******************************************************************************************/
int
TclModelBuilder_addElement2dYS (ClientData clientData, Tcl_Interp *interp,
								   int argc, TCL_Char **argv,
								   Domain *theTclDomain, TclModelBuilder *theTclBuilder)
{

  if (strcmp(argv[1],"inelastic2dYS01") == 0) {
	  int result = TclModelBuilder_addElement2dYS01(clientData, interp, argc, argv,
						 theTclDomain, theTclBuilder);
    return result;
  }
  else if (strcmp(argv[1],"inelastic2dYS02") == 0) {
	  int result = TclModelBuilder_addElement2dYS02(clientData, interp, argc, argv,
						 theTclDomain, theTclBuilder);
    return result;
  }
  
   else if (strcmp(argv[1],"inelastic2dYS03") == 0) {
	  int result = TclModelBuilder_addElement2dYS03(clientData, interp, argc, argv,
						 theTclDomain, theTclBuilder);
    return result;
  }

/*	else if (strcmp(argv[1],"inelastic2dYS04") == 0) {
	  int result = TclModelBuilder_addElement2dYS04
	  (clientData, interp, argc, argv,
						 theTclDomain, theTclBuilder);
    return result;
  }*/
  /*else if (strcmp(argv[1],"inelastic2dYS05") == 0) {
	  int result = TclModelBuilder_addElement2dYS05
	  (clientData, interp, argc, argv,
						 theTclDomain, theTclBuilder);
    return result;
  }*/
  else

	return TCL_ERROR;

}
