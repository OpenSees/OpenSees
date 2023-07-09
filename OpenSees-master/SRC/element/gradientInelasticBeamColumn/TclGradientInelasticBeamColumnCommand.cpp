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

/* Written by: Mohammad Salehi (mohammad.salehi@tamu.edu)
** Created: 07/19
** Description: The source code for the Tcl commands of gradient inelastic (GI) force-based beam-column element formulation
*/

#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <TclModelBuilder.h>

#include "GradientInelasticBeamColumn2d.h"
#include "GradientInelasticBeamColumn3d.h"

#include <NewtonCotesBeamIntegration.h>
#include <TrapezoidalBeamIntegration.h>
#include <SimpsonBeamIntegration.h>
#include <LobattoBeamIntegration.h>
#include <LegendreBeamIntegration.h>

extern void printCommand(int argc, TCL_Char** argv);

int
TclModelBuilder_addGradientInelasticBeamColumn(ClientData clientData, Tcl_Interp* interp,
	int argc,
	TCL_Char** argv,
	Domain* theTclDomain,
	TclModelBuilder* theTclBuilder)
{
	// ensure the destructor has not been called
	if (theTclBuilder == 0) {
		opserr << "WARNING builder has been destroyed - gradientInelasticBeamColumn\n";
		return TCL_ERROR;
	}

	Element* theElement = 0;
	int ndm = theTclBuilder->getNDM();
	int ndf = theTclBuilder->getNDF();
	int tag;
	int eleArgStart = 1;

	if (ndm == 2) {
		// check plane frame problem has 3 dof per node
		if (ndf != 3) {
			opserr << "WARNING invalid ndf: " << ndf;
			opserr << ", for plane problem need 3 - gradientInelasticBeamColumn\n";
			return TCL_ERROR;
		}

		// check the number of arguments is correct
		if ((argc - eleArgStart) < 10) {
			opserr << "WARNING insufficient arguments\n";
			printCommand(argc, argv);
			opserr << "Want: gradientInelasticBeamColumn eleTag? iNode? jNode? numIntgrPts? endSecTag1? intSecTag? endSecTag2? secLR1? secLR2? lc? transfTag? <-constH> <-integration integrType?> <-iter maxIter? minTol? maxTol?> <-corControl auto/maxEpsInc? maxPhiInc?>\n";
			return TCL_ERROR;
		}

		// get the id and end nodes 
		int iNode, jNode, endSecTag1, intSecTag, endSecTag2, numIntegrPts, transfTag;
		double lc, secLR1, secLR2;
		int maxIter = 50;
		double minTol = 1E-10, maxTol = 1E-8;
		bool correctionControl = false;
		bool constH = false;
		double maxEpsInc = 0.0, maxPhiInc = 0.0;

		if (Tcl_GetInt(interp, argv[1 + eleArgStart], &tag) != TCL_OK) {
			opserr << "WARNING invalid gradientInelasticBeamColumn eleTag\n";
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[2 + eleArgStart], &iNode) != TCL_OK) {
			opserr << "WARNING invalid iNode";
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[3 + eleArgStart], &jNode) != TCL_OK) {
			opserr << "WARNING invalid jNode";
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[4 + eleArgStart], &numIntegrPts) != TCL_OK) {
			opserr << "WARNING invalid numIntegrPts";
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[5 + eleArgStart], &endSecTag1) != TCL_OK) {
			opserr << "WARNING invalid firstSecTag";
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[6 + eleArgStart], &intSecTag) != TCL_OK) {
			opserr << "WARNING invalid intSecTag";
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[7 + eleArgStart], &endSecTag2) != TCL_OK) {
			opserr << "WARNING invalid lastSecTag";
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		SectionForceDeformation* endSection1 = theTclBuilder->getSection(endSecTag1);

		if (!endSection1) {
			opserr << "WARNING end section not found";
			opserr << " - section: " << endSecTag1;
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		SectionForceDeformation* intSection = theTclBuilder->getSection(intSecTag);

		if (!intSection) {
			opserr << "WARNING intermediate section not found";
			opserr << " - section: " << intSecTag;
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		SectionForceDeformation* endSection2 = theTclBuilder->getSection(endSecTag2);

		if (!endSection2) {
			opserr << "WARNING end section not found";
			opserr << " - section: " << endSecTag2;
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[8 + eleArgStart], &secLR1) != TCL_OK) {
			opserr << "WARNING invalid secLR1";
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[9 + eleArgStart], &secLR2) != TCL_OK) {
			opserr << "WARNING invalid secLR2";
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[10 + eleArgStart], &lc) != TCL_OK) {
			opserr << "WARNING invalid lc";
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[11 + eleArgStart], &transfTag) != TCL_OK) {
			opserr << "WARNING invalid transfTag";
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		CrdTransf* theTransf2d = OPS_getCrdTransf(transfTag);	// was: theTclBuilder->getCrdTransf(transfTag);

		if (!theTransf2d) {
			opserr << "WARNING transformation not found";
			opserr << " - transformation: " << transfTag;
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		// check for optional arguments
		BeamIntegration* beamIntegr = 0;

		for (int i = 12 + eleArgStart; i < argc; i++) {
			if (i + 1 < argc && strcmp(argv[i], "-integration") == 0) {
				if (strcmp(argv[i + 1], "NewtonCotes") == 0) {
					if (numIntegrPts > 20) {
						opserr << "WARNING number of integration points must be less than 20 for Newton-Cotes integration method";
						opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
						return TCL_ERROR;
					}

					beamIntegr = new NewtonCotesBeamIntegration();
				}
				else if (strcmp(argv[i + 1], "Simpson") == 0) {
					if (((numIntegrPts - 1) % 2) != 0) {
						opserr << "WARNING number of integration points must be odd for Simpson's integration method";
						opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
						return TCL_ERROR;
					}

					beamIntegr = new SimpsonBeamIntegration();
				}
				else if (strcmp(argv[i + 1], "Trapezoidal") == 0)
					beamIntegr = new TrapezoidalBeamIntegration();
				else if (strcmp(argv[i + 1], "Lobatto") == 0)
					beamIntegr = new LobattoBeamIntegration();
				else if (strcmp(argv[i + 1], "Legendre") == 0)
					beamIntegr = new LegendreBeamIntegration();

				if (beamIntegr == 0) {
					opserr << "WARNING invalid integration type";
					opserr << " - gradientInelasticBeamColumn element: " << tag;
					opserr << " - Simpson's integration method is used\n";
				}
			}
		}

		if (!beamIntegr) {
			if (((numIntegrPts - 1) % 2) != 0) {
				opserr << "WARNING number of integration points must be odd for Simpson's integration method";
				opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
				return TCL_ERROR;
			}

			beamIntegr = new SimpsonBeamIntegration();
		}

		for (int i = 12 + eleArgStart; i < argc; i++) {
			if (i + 3 < argc && strcmp(argv[i], "-iter") == 0) {
				if (Tcl_GetInt(interp, argv[i + 1], &maxIter) != TCL_OK) {
					opserr << "WARNING invalid maxIter";
					opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
					return TCL_ERROR;
				}

				if (Tcl_GetDouble(interp, argv[i + 2], &minTol) != TCL_OK) {
					opserr << "WARNING invalid minTol";
					opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
					return TCL_ERROR;
				}

				if (Tcl_GetDouble(interp, argv[i + 3], &maxTol) != TCL_OK) {
					opserr << "WARNING invalid maxTol";
					opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
					return TCL_ERROR;
				}
			}
		}

		for (int i = 12 + eleArgStart; i < argc; i++) {
			if (i + 1 < argc && strcmp(argv[i], "-corControl") == 0) {
				correctionControl = true;

				if (i + 2 < argc && strcmp(argv[i + 1], "auto") != 0) {

					if (Tcl_GetDouble(interp, argv[i + 1], &maxEpsInc) != TCL_OK) {
						opserr << "WARNING invalid maxEpsInc";
						opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
						return TCL_ERROR;
					}

					if (Tcl_GetDouble(interp, argv[i + 2], &maxPhiInc) != TCL_OK) {
						opserr << "WARNING invalid maxPhiInc";
						opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
						return TCL_ERROR;
					}
				}
			}
		}

		for (int i = 12 + eleArgStart; i < argc; i++)
			if (strcmp(argv[i], "-constH") == 0)
				constH = true;

		// now create the NonlocalBeamColumn
		theElement = new GradientInelasticBeamColumn2d(tag, iNode, jNode, numIntegrPts, &endSection1, &intSection, &endSection2, secLR1, secLR2, *beamIntegr, *theTransf2d, lc, minTol, maxTol, maxIter, constH, correctionControl, maxEpsInc, maxPhiInc);

		if (!theElement) {
			opserr << "WARNING ran out of memory creating element";
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		// then add the NonlocalBeamColumn to the domain
		if (theTclDomain->addElement(theElement) == false) {
			opserr << "WARNING could not add element to the domain";
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			delete theElement;
			return TCL_ERROR;
		}
	}
	else if (ndm == 3) {
		// check plane frame problem has 6 dof per node
		if (ndf != 6) {
			opserr << "WARNING invalid ndf: " << ndf;
			opserr << ", for plane problem need 6 - gradientInelasticBeamColumn\n";
			return TCL_ERROR;
		}

		// check the number of arguments is correct
		if ((argc - eleArgStart) < 10) {
			opserr << "WARNING insufficient arguments\n";
			printCommand(argc, argv);
			opserr << "Want: gradientInelasticBeamColumn eleTag? iNode? jNode? numIntgrPts? endSecTag1? intSecTag? endSecTag2? secLR1? secLR2? lc? transfTag?  <-constH> <-integration integrType?> <-iter maxIter? minTol? maxTol?> <-corControl auto/maxEpsInc? maxPhiInc?>\n";
			return TCL_ERROR;
		}

		// get the id and end nodes 
		int iNode, jNode, endSecTag1, intSecTag, endSecTag2, secOrder, numIntegrPts, transfTag;
		double lc, secLR1, secLR2;
		int maxIter = 50;
		double minTol = 1E-10, maxTol = 1E-8;
		bool correctionControl = false;
		bool constH = false;
		double maxEpsInc = 0.0, maxPhiInc = 0.0;

		if (Tcl_GetInt(interp, argv[1 + eleArgStart], &tag) != TCL_OK) {
			opserr << "WARNING invalid gradientInelasticBeamColumn eleTag\n";
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[2 + eleArgStart], &iNode) != TCL_OK) {
			opserr << "WARNING invalid iNode";
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[3 + eleArgStart], &jNode) != TCL_OK) {
			opserr << "WARNING invalid jNode";
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[4 + eleArgStart], &numIntegrPts) != TCL_OK) {
			opserr << "WARNING invalid numIntegrPts";
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[5 + eleArgStart], &endSecTag1) != TCL_OK) {
			opserr << "WARNING invalid firstSecTag";
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[6 + eleArgStart], &intSecTag) != TCL_OK) {
			opserr << "WARNING invalid intSecTag";
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[7 + eleArgStart], &endSecTag2) != TCL_OK) {
			opserr << "WARNING invalid lastSecTag";
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		SectionForceDeformation* endSection1 = theTclBuilder->getSection(endSecTag1);

		if (!endSection1) {
			opserr << "WARNING end section not found";
			opserr << " - section: " << endSecTag1;
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		SectionForceDeformation* intSection = theTclBuilder->getSection(intSecTag);

		if (!intSection) {
			opserr << "WARNING intermediate section not found";
			opserr << " - section: " << intSecTag;
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		SectionForceDeformation* endSection2 = theTclBuilder->getSection(endSecTag2);

		if (!endSection2) {
			opserr << "WARNING end section not found";
			opserr << " - section: " << endSecTag2;
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[8 + eleArgStart], &secLR1) != TCL_OK) {
			opserr << "WARNING invalid secLR1";
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[9 + eleArgStart], &secLR2) != TCL_OK) {
			opserr << "WARNING invalid secLR2";
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[10 + eleArgStart], &lc) != TCL_OK) {
			opserr << "WARNING invalid lc";
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[11 + eleArgStart], &transfTag) != TCL_OK) {
			opserr << "WARNING invalid transfTag";
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		CrdTransf* theTransf3d = OPS_getCrdTransf(transfTag);	// was: theTclBuilder->getCrdTransf3d(transfTag);

		if (!theTransf3d) {
			opserr << "WARNING transformation not found";
			opserr << " - transformation: " << transfTag;
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		// check for optional arguments
		BeamIntegration* beamIntegr = 0;

		for (int i = 12 + eleArgStart; i < argc; i++) {
			if (i + 1 < argc && strcmp(argv[i], "-integration") == 0) {
				if (strcmp(argv[i + 1], "NewtonCotes") == 0) {
					if (numIntegrPts > 20) {
						opserr << "WARNING number of integration points must be less than 20 for Newton-Cotes integration method";
						opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
						return TCL_ERROR;
					}

					beamIntegr = new NewtonCotesBeamIntegration();
				}
				else if (strcmp(argv[i + 1], "Simpson") == 0) {
					if (((numIntegrPts - 1) % 2) != 0) {
						opserr << "WARNING number of integration points must be odd for Simpson's integration method";
						opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
						return TCL_ERROR;
					}

					beamIntegr = new SimpsonBeamIntegration();
				}
				else if (strcmp(argv[i + 1], "Trapezoidal") == 0)
					beamIntegr = new TrapezoidalBeamIntegration();
				else if (strcmp(argv[i + 1], "Lobatto") == 0)
					beamIntegr = new LobattoBeamIntegration();
				else if (strcmp(argv[i + 1], "Legendre") == 0)
					beamIntegr = new LegendreBeamIntegration();

				if (beamIntegr == 0) {
					opserr << "WARNING invalid integration type";
					opserr << " - gradientInelasticBeamColumn element: " << tag;
					opserr << " - Simpson's integration method is used\n";
				}
			}
		}

		if (!beamIntegr) {
			if (((numIntegrPts - 1) % 2) != 0) {
				opserr << "WARNING number of integration points must be odd for Simpson's integration method";
				opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
				return TCL_ERROR;
			}

			beamIntegr = new SimpsonBeamIntegration();
		}

		for (int i = 12 + eleArgStart; i < argc; i++) {
			if (i + 3 < argc && strcmp(argv[i], "-iter") == 0) {
				if (Tcl_GetInt(interp, argv[i + 1], &maxIter) != TCL_OK) {
					opserr << "WARNING invalid maxIter";
					opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
					return TCL_ERROR;
				}

				if (Tcl_GetDouble(interp, argv[i + 2], &minTol) != TCL_OK) {
					opserr << "WARNING invalid minTol";
					opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
					return TCL_ERROR;
				}

				if (Tcl_GetDouble(interp, argv[i + 3], &maxTol) != TCL_OK) {
					opserr << "WARNING invalid maxTol";
					opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
					return TCL_ERROR;
				}
			}
		}

		for (int i = 12 + eleArgStart; i < argc; i++) {
			if (i + 1 < argc && strcmp(argv[i], "-corControl") == 0) {
				correctionControl = true;

				if (i + 2 < argc && strcmp(argv[i + 1], "auto") != 0) {

					if (Tcl_GetDouble(interp, argv[i + 1], &maxEpsInc) != TCL_OK) {
						opserr << "WARNING invalid maxEpsInc";
						opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
						return TCL_ERROR;
					}

					if (Tcl_GetDouble(interp, argv[i + 2], &maxPhiInc) != TCL_OK) {
						opserr << "WARNING invalid maxPhiInc";
						opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
						return TCL_ERROR;
					}
				}
			}
		}

		for (int i = 12 + eleArgStart; i < argc; i++)
			if (strcmp(argv[i], "-constH") == 0)
				constH = true;

		// now create the NonlocalBeamColumn
		theElement = new GradientInelasticBeamColumn3d(tag, iNode, jNode, numIntegrPts, &endSection1, &intSection, &endSection2, secLR1, secLR2, *beamIntegr, *theTransf3d, lc, minTol, maxTol, maxIter, constH, correctionControl, maxEpsInc, maxPhiInc);

		if (!theElement) {
			opserr << "WARNING ran out of memory creating element";
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			return TCL_ERROR;
		}

		// then add the NonlocalBeamColumn to the domain
		if (theTclDomain->addElement(theElement) == false) {
			opserr << "WARNING could not add element to the domain";
			opserr << " - gradientInelasticBeamColumn element: " << tag << endln;
			delete theElement;
			return TCL_ERROR;
		}
	}
	else {
		opserr << "WARNING gradientInelasticBeamColumn command only works when ndm is 2 or 3, ndm: ";
		opserr << ndm << endln;
		return TCL_ERROR;
	}

	// if get here we have successfully created the NonlocalBeamColumn and added it to the domain
	return TCL_OK;
}

