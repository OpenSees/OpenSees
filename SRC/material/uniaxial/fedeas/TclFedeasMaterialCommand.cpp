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

// $Revision: 1.6 $
// $Date: 2004-07-15 21:36:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/fedeas/TclFedeasMaterialCommand.cpp,v $

// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the implementation of the
// TclModelBuilder_addFedeasMaterial() function. 

#include <FedeasHardeningMaterial.h>
#include <FedeasBond1Material.h>
#include <FedeasBond2Material.h>
#include <FedeasConcr1Material.h>
#include <FedeasConcr2Material.h>
#include <FedeasConcr3Material.h>
#include <FedeasHyster1Material.h>
#include <FedeasHyster2Material.h>
#include <FedeasSteel2Material.h>
#include <FedeasSteel1Material.h>

#include <TclModelBuilder.h>
#include <Vector.h>
#include <string.h>

static void printCommand(int argc, TCL_Char **argv)
{
    opserr << "Input command: ";
    for (int i=0; i<argc; i++)
	opserr << argv[i] << " ";
    opserr << endln;
} 

UniaxialMaterial *
TclModelBuilder_addFedeasMaterial(ClientData clientData, Tcl_Interp *interp, int argc, 
				  TCL_Char **argv, TclModelBuilder *theTclBuilder)
{
	if (argc < 3) {
		opserr << "WARNING insufficient number of arguments\n";
		printCommand(argc, argv);
		return 0;
	}

	int tag;
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid uniaxialMaterial tag\n";
		printCommand(argc, argv);
	    return 0;
	}

	UniaxialMaterial *theMaterial = 0;

	if (strcmp(argv[1],"Hardening1") == 0 || strcmp(argv[1],"Hardening01") == 0) {
		if (argc < 7) {
			opserr << "WARNING invalid number of arguments\n";
			printCommand(argc,argv);
			opserr << "Want: uniaxialMaterial Hardening01 tag? E? sigY? Hiso? Hkin?" << endln;
			return 0;
		}
		
		double E, sigY, Hiso, Hkin;

		if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
			opserr << "WARNING invalid E\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[4], &sigY) != TCL_OK) {
			opserr << "WARNING invalid sigY\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[5], &Hiso) != TCL_OK) {
			opserr << "WARNING invalid Hiso\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[6], &Hkin) != TCL_OK) {
			opserr << "WARNING invalid Hkin\n";
			printCommand(argc, argv);
			return 0;	
		}

		theMaterial = new FedeasHardeningMaterial(tag, E, sigY, Hiso, Hkin);
	}

	else if (strcmp(argv[1],"Bond1") == 0 || strcmp(argv[1],"Bond01") == 0) {
		if (argc < 15) {
			opserr << "WARNING invalid number of arguments\n";
			printCommand(argc,argv);
			opserr << "Want: uniaxialMaterial Bond01 tag? u1p? q1p? u2p? u3p? q3p? u1n? q1n? u2n? u3n? q3n? s0? bb?" << endln;
			return 0;
		}

		double u1p, q1p, u2p, u3p, q3p;
		double u1n, q1n, u2n, u3n, q3n;
		double s0, bb;

		if (Tcl_GetDouble(interp, argv[3], &u1p) != TCL_OK) {
			opserr << "WARNING invalid u1p\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[4], &q1p) != TCL_OK) {
			opserr << "WARNING invalid q1p\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[5], &u2p) != TCL_OK) {
			opserr << "WARNING invalid u2p\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[6], &u3p) != TCL_OK) {
			opserr << "WARNING invalid u3p\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[7], &q3p) != TCL_OK) {
			opserr << "WARNING invalid q3p\n";
			printCommand(argc, argv);
			return 0;	
		}

		if (Tcl_GetDouble(interp, argv[8], &u1n) != TCL_OK) {
			opserr << "WARNING invalid u1n\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[9], &q1n) != TCL_OK) {
			opserr << "WARNING invalid q1n\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[10], &u2n) != TCL_OK) {
			opserr << "WARNING invalid u2n\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[11], &u3n) != TCL_OK) {
			opserr << "WARNING invalid u3n\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[12], &q3n) != TCL_OK) {
			opserr << "WARNING invalid q3n\n";
			printCommand(argc, argv);
			return 0;	
		}

		if (Tcl_GetDouble(interp, argv[13], &s0) != TCL_OK) {
			opserr << "WARNING invalid s0\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[14], &bb) != TCL_OK) {
			opserr << "WARNING invalid bb\n";
			printCommand(argc, argv);
			return 0;	
		}

		theMaterial = new FedeasBond1Material(tag, u1p, q1p, u2p, u3p, q3p,	u1n, q1n, u2n, u3n, q3n, s0, bb);
	}

	else if (strcmp(argv[1],"Bond2") == 0 || strcmp(argv[1],"Bond02") == 0) {
		if (argc < 17) {
			opserr << "WARNING invalid number of arguments\n";
			printCommand(argc,argv);
			opserr << "Want: uniaxialMaterial Bond02 tag? u1p? q1p? u2p? u3p? q3p? u1n? q1n? u2n? u3n? q3n? s0? bb? alp? aln?" << endln;
			return 0;
		}

		double u1p, q1p, u2p, u3p, q3p;
		double u1n, q1n, u2n, u3n, q3n;
		double s0, bb, alp, aln;

		if (Tcl_GetDouble(interp, argv[3], &u1p) != TCL_OK) {
			opserr << "WARNING invalid u1p\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[4], &q1p) != TCL_OK) {
			opserr << "WARNING invalid q1p\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[5], &u2p) != TCL_OK) {
			opserr << "WARNING invalid u2p\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[6], &u3p) != TCL_OK) {
			opserr << "WARNING invalid u3p\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[7], &q3p) != TCL_OK) {
			opserr << "WARNING invalid q3p\n";
			printCommand(argc, argv);
			return 0;	
		}

		if (Tcl_GetDouble(interp, argv[8], &u1n) != TCL_OK) {
			opserr << "WARNING invalid u1n\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[9], &q1n) != TCL_OK) {
			opserr << "WARNING invalid q1n\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[10], &u2n) != TCL_OK) {
			opserr << "WARNING invalid u2n\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[11], &u3n) != TCL_OK) {
			opserr << "WARNING invalid u3n\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[12], &q3n) != TCL_OK) {
			opserr << "WARNING invalid q3n\n";
			printCommand(argc, argv);
			return 0;	
		}

		if (Tcl_GetDouble(interp, argv[13], &s0) != TCL_OK) {
			opserr << "WARNING invalid s0\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[14], &bb) != TCL_OK) {
			opserr << "WARNING invalid bb\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[15], &alp) != TCL_OK) {
			opserr << "WARNING invalid alp\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[16], &aln) != TCL_OK) {
			opserr << "WARNING invalid aln\n";
			printCommand(argc, argv);
			return 0;	
		}

		theMaterial = new FedeasBond2Material(tag, u1p, q1p, u2p, u3p, q3p,	u1n, q1n, u2n, u3n, q3n, s0, bb, alp, aln);
	}

	else if (strcmp(argv[1],"Concrete1") == 0 || strcmp(argv[1],"concrete01") == 0) {
		if (argc < 7) {
			opserr << "WARNING invalid number of arguments\n";
			printCommand(argc,argv);
			opserr << "Want: uniaxialMaterial Concrete01 tag? fpc? epsc0? fpcu? epscu?" << endln;
			return 0;
		}    

		double fpc, epsc0, fpcu, epscu;

		if (Tcl_GetDouble(interp, argv[3], &fpc) != TCL_OK) {
			opserr << "WARNING invalid fpc\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[4], &epsc0) != TCL_OK) {
			opserr << "WARNING invalid epsc0\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[5], &fpcu) != TCL_OK) {
			opserr << "WARNING invalid fpcu\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[6], &epscu) != TCL_OK) {
			opserr << "WARNING invalid epscu\n";
			printCommand(argc, argv);
			return 0;	
		}

		theMaterial = new FedeasConcr1Material(tag, fpc, epsc0, fpcu, epscu);
	}

	else if (strcmp(argv[1],"Concrete2") == 0 || strcmp(argv[1],"Concrete02") == 0) {
		if (argc < 10) {
			opserr << "WARNING invalid number of arguments\n";
			printCommand(argc,argv);
			opserr << "Want: uniaxialMaterial Concrete02 tag? fpc? epsc0? fpcu? epscu? rat? ft? Ets?" << endln;
			return 0;
		}    

		double fpc, epsc0, fpcu, epscu;
		double rat, ft, Ets;

		if (Tcl_GetDouble(interp, argv[3], &fpc) != TCL_OK) {
			opserr << "WARNING invalid fpc\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[4], &epsc0) != TCL_OK) {
			opserr << "WARNING invalid epsc0\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[5], &fpcu) != TCL_OK) {
			opserr << "WARNING invalid fpcu\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[6], &epscu) != TCL_OK) {
			opserr << "WARNING invalid epscu\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[7], &rat) != TCL_OK) {
			opserr << "WARNING invalid rat\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[8], &ft) != TCL_OK) {
			opserr << "WARNING invalid ft\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[9], &Ets) != TCL_OK) {
			opserr << "WARNING invalid Ets\n";
			printCommand(argc, argv);
			return 0;	
		}

		theMaterial = new FedeasConcr2Material(tag, fpc, epsc0, fpcu, epscu, rat, ft, Ets);
	}

	else if (strcmp(argv[1],"Concrete3") == 0 || strcmp(argv[1],"Concrete03") == 0) {
		if (argc < 13) {
			opserr << "WARNING invalid number of arguments\n";
			printCommand(argc,argv);
			opserr << "Want: uniaxialMaterial Concrete03 tag? fpc? epsc0? fpcu? epscu? rat? ft? epst0? ft0? beta? epstu?" << endln;
			return 0;
		}    

		double fpc, epsc0, fpcu, epscu;
		double rat, ft, epst0, ft0, beta, epstu;

		if (Tcl_GetDouble(interp, argv[3], &fpc) != TCL_OK) {
			opserr << "WARNING invalid fpc\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[4], &epsc0) != TCL_OK) {
			opserr << "WARNING invalid epsc0\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[5], &fpcu) != TCL_OK) {
			opserr << "WARNING invalid fpcu\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[6], &epscu) != TCL_OK) {
			opserr << "WARNING invalid epscu\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[7], &rat) != TCL_OK) {
			opserr << "WARNING invalid rat\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[8], &ft) != TCL_OK) {
			opserr << "WARNING invalid ft\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[9], &epst0) != TCL_OK) {
			opserr << "WARNING invalid epst0\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[10], &ft0) != TCL_OK) {
			opserr << "WARNING invalid ft0\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[11], &beta) != TCL_OK) {
			opserr << "WARNING invalid beta\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[12], &epstu) != TCL_OK) {
			opserr << "WARNING invalid epstu\n";
			printCommand(argc, argv);
			return 0;	
		}

		theMaterial = new FedeasConcr3Material(tag, fpc, epsc0, fpcu, epscu, rat, ft, epst0, ft0, beta, epstu);
	}

	else if (strcmp(argv[1],"Hysteretic1") == 0 || strcmp(argv[1],"Hysteretic01") == 0) {
		if (argc < 15) {
			opserr << "WARNING invalid number of arguments\n";
			printCommand(argc,argv);
			opserr << "Want: uniaxialMaterial Hysteretic01 tag? s1p? e1p? s2p? e2p? s1n? e1n? s2n? e1n? px? py? d1? d2?" << endln;
			return 0;
		}    

		double s1p, e1p, s2p, e2p;
		double s1n, e1n, s2n, e2n;
		double px, py, d1, d2;

		if (Tcl_GetDouble(interp, argv[3], &s1p) != TCL_OK) {
			opserr << "WARNING invalid s1p\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[4], &e1p) != TCL_OK) {
			opserr << "WARNING invalid e1p\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[5], &s2p) != TCL_OK) {
			opserr << "WARNING invalid s2p\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[6], &e2p) != TCL_OK) {
			opserr << "WARNING invalid e2p\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[7], &s1n) != TCL_OK) {
			opserr << "WARNING invalid s1n\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[8], &e1n) != TCL_OK) {
			opserr << "WARNING invalid e1n\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[9], &s2n) != TCL_OK) {
			opserr << "WARNING invalid s2n\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[10], &e2n) != TCL_OK) {
			opserr << "WARNING invalid e2n\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[11], &px) != TCL_OK) {
			opserr << "WARNING invalid px\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[12], &py) != TCL_OK) {
			opserr << "WARNING invalid py\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[13], &d1) != TCL_OK) {
			opserr << "WARNING invalid d1\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[14], &d2) != TCL_OK) {
			opserr << "WARNING invalid d2\n";
			printCommand(argc, argv);
			return 0;	
		}

		theMaterial = new FedeasHyster1Material(tag, s1p, e1p, s2p, e2p, s1n, e1n, s2n, e2n, px, py, d1, d2);
	}

	else if (strcmp(argv[1],"Hysteretic2") == 0 || strcmp(argv[1],"Hysteretic02") == 0) {
		if (argc < 19) {
			opserr << "WARNING invalid number of arguments\n";
			printCommand(argc,argv);
			opserr << "Want: uniaxialMaterial Hysteretic02 tag? s1p? e1p? s2p? e2p? s3p? e3p? s1n? e1n? s2n? e1n? s3n? e3n? px? py? d1? d2?" << endln;
			return 0;
		}    

		double s1p, e1p, s2p, e2p, s3p, e3p;
		double s1n, e1n, s2n, e2n, s3n, e3n;
		double px, py, d1, d2;

		if (Tcl_GetDouble(interp, argv[3], &s1p) != TCL_OK) {
			opserr << "WARNING invalid s1p\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[4], &e1p) != TCL_OK) {
			opserr << "WARNING invalid e1p\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[5], &s2p) != TCL_OK) {
			opserr << "WARNING invalid s2p\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[6], &e2p) != TCL_OK) {
			opserr << "WARNING invalid e2p\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[7], &s3p) != TCL_OK) {
			opserr << "WARNING invalid s2p\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[8], &e3p) != TCL_OK) {
			opserr << "WARNING invalid e2p\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[9], &s1n) != TCL_OK) {
			opserr << "WARNING invalid s1n\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[10], &e1n) != TCL_OK) {
			opserr << "WARNING invalid e1n\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[11], &s2n) != TCL_OK) {
			opserr << "WARNING invalid s2n\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[12], &e2n) != TCL_OK) {
			opserr << "WARNING invalid e2n\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[13], &s3n) != TCL_OK) {
			opserr << "WARNING invalid s2n\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[14], &e3n) != TCL_OK) {
			opserr << "WARNING invalid e2n\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[15], &px) != TCL_OK) {
			opserr << "WARNING invalid px\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[16], &py) != TCL_OK) {
			opserr << "WARNING invalid py\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[17], &d1) != TCL_OK) {
			opserr << "WARNING invalid d1\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[18], &d2) != TCL_OK) {
			opserr << "WARNING invalid d2\n";
			printCommand(argc, argv);
			return 0;	
		}

		theMaterial = new FedeasHyster2Material(tag, s1p, e1p, s2p, e2p, s3p, e3p, s1n, e1n, s2n, e2n, s3n, e3n, px, py, d1, d2);
	}

	else if (strcmp(argv[1],"Steel1") == 0 || strcmp(argv[1],"Steel01") == 0) {
		if (argc < 6) {
			opserr << "WARNING invalid number of arguments\n";
			printCommand(argc,argv);
			opserr << "Want: uniaxialMaterial Steel01 tag? fy? E? b? <a1? a2? a3? a4?>" << endln;
			return 0;
		}

		double fy, E, b;
		double a1, a2, a3, a4;

		if (Tcl_GetDouble(interp, argv[3], &fy) != TCL_OK) {
			opserr << "WARNING invalid fy\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[4], &E) != TCL_OK) {
			opserr << "WARNING invalid E\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[5], &b) != TCL_OK) {
			opserr << "WARNING invalid b\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (argc > 9) {
			if (Tcl_GetDouble(interp, argv[6], &a1) != TCL_OK) {
				opserr << "WARNING invalid a1\n";
				printCommand(argc, argv);
				return 0;	
			}
			if (Tcl_GetDouble(interp, argv[7], &a2) != TCL_OK) {
				opserr << "WARNING invalid a2\n";
				printCommand(argc, argv);
				return 0;	
			}
			if (Tcl_GetDouble(interp, argv[8], &a3) != TCL_OK) {
				opserr << "WARNING invalid a3\n";
				printCommand(argc, argv);
				return 0;	
			}
			if (Tcl_GetDouble(interp, argv[9], &a4) != TCL_OK) {
				opserr << "WARNING invalid a4\n";
				printCommand(argc, argv);
				return 0;	
			}
			theMaterial = new FedeasSteel1Material(tag, fy, E, b, a1, a2, a3, a4);
		}
		else
			theMaterial = new FedeasSteel1Material(tag, fy, E, b);
	}

	else if (strcmp(argv[1],"Steel2") == 0 || strcmp(argv[1],"Steel02") == 0) {
		if (argc < 6) {
			opserr << "WARNING invalid number of arguments\n";
			printCommand(argc,argv);
			opserr << "Want: uniaxialMaterial Steel02 tag? fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
			return 0;
		}    

		double fy, E, b;
		double R0, cR1, cR2;
		double a1, a2, a3, a4;

		if (Tcl_GetDouble(interp, argv[3], &fy) != TCL_OK) {
			opserr << "WARNING invalid fy\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[4], &E) != TCL_OK) {
			opserr << "WARNING invalid E\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (Tcl_GetDouble(interp, argv[5], &b) != TCL_OK) {
			opserr << "WARNING invalid b\n";
			printCommand(argc, argv);
			return 0;	
		}
		if (argc > 8) {
			if (Tcl_GetDouble(interp, argv[6], &R0) != TCL_OK) {
				opserr << "WARNING invalid R0\n";
				printCommand(argc, argv);
				return 0;	
			}
			if (Tcl_GetDouble(interp, argv[7], &cR1) != TCL_OK) {
				opserr << "WARNING invalid cR1\n";
				printCommand(argc, argv);
				return 0;	
			}
			if (Tcl_GetDouble(interp, argv[8], &cR2) != TCL_OK) {
				opserr << "WARNING invalid cR2\n";
				printCommand(argc, argv);
				return 0;	
			}
			if (argc > 12) {
				if (Tcl_GetDouble(interp, argv[9], &a1) != TCL_OK) {
					opserr << "WARNING invalid a1\n";
					printCommand(argc, argv);
					return 0;	
				}
				if (Tcl_GetDouble(interp, argv[10], &a2) != TCL_OK) {
					opserr << "WARNING invalid a2\n";
					printCommand(argc, argv);
					return 0;	
				}
				if (Tcl_GetDouble(interp, argv[11], &a3) != TCL_OK) {
					opserr << "WARNING invalid a3\n";
					printCommand(argc, argv);
					return 0;	
				}
				if (Tcl_GetDouble(interp, argv[12], &a4) != TCL_OK) {
					opserr << "WARNING invalid a4\n";
					printCommand(argc, argv);
					return 0;	
				}
				theMaterial = new FedeasSteel2Material(tag, fy, E, b, R0, cR1, cR2, a1, a2, a3, a4);
			}
			else
				theMaterial = new FedeasSteel2Material(tag, fy, E, b, R0, cR1, cR2);
		}
		else
			theMaterial = new FedeasSteel2Material(tag, fy, E, b);

	}

	return theMaterial;
}
