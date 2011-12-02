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
** With a lot of additions by                                         **
**   Boris Jeremic    (jeremic@ucdavis.edu)                           **
**   Zaohui Yang      (zhyang@ucdavis.edu)                            **
**   Zhao Cheng       (zcheng@ucdavis.edu)                            **
**                                                                    **
**                                                                    **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.40 $
// $Date: 2006-08-11 21:22:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/TclModelBuilderNDMaterialCommand.cpp,v $


// Description: This file contains the function invoked when the user invokes
// the nDMaterial command in the interpreter.
//
// What: "@(#) TclModelBuilderNDMaterialCommand.C, revA"

#include <TclModelBuilder.h>

#include <ElasticIsotropicMaterial.h>
#include <ElasticIsotropic3D.h>
#include <PressureDependentElastic3D.h>
#include <ElasticCrossAnisotropic.h>
#include <J2Plasticity.h>

#include <MultiaxialCyclicPlasticity.h> //Gang Wang

#include <PlaneStressMaterial.h>
#include <PlateFiberMaterial.h>
#include <BeamFiberMaterial.h>

#include <PressureIndependMultiYield.h>
#include <PressureDependMultiYield.h>
#include <PressureDependMultiYield02.h>
#include <FluidSolidPorousMaterial.h>

#include <string.h>

#include <Template3Dep.h>
#include <NewTemplate3Dep.h>
#include <FiniteDeformationElastic3D.h>
#include <FiniteDeformationEP3D.h>

Template3Dep *
TclModelBuilder_addTemplate3Dep(ClientData clientData, Tcl_Interp *interp,  int argc,
				TCL_Char **argv, TclModelBuilder *theTclBuilder, int eleArgStart);

NewTemplate3Dep *
TclModelBuilder_addNewTemplate3Dep(ClientData clientData, Tcl_Interp *interp,  int argc,
				TCL_Char **argv, TclModelBuilder *theTclBuilder, int eleArgStart);				
				
FiniteDeformationElastic3D *
TclModelBuilder_addFiniteDeformationElastic3D(ClientData clientData, Tcl_Interp *interp,  int argc,
				TCL_Char **argv, TclModelBuilder *theTclBuilder, int eleArgStart);

FiniteDeformationEP3D *
TclModelBuilder_addFiniteDeformationEP3D(ClientData clientData, Tcl_Interp *interp,  int argc,
				TCL_Char **argv, TclModelBuilder *theTclBuilder, int eleArgStart);

NDMaterial *
TclModelBuilder_addFeapMaterial(ClientData clientData, Tcl_Interp *interp,
				int argc, TCL_Char **argv,
				TclModelBuilder *theTclBuilder);


static void printCommand(int argc, TCL_Char **argv)
{
    opserr << "Input command: ";
    for (int i=0; i<argc; i++)
	opserr << argv[i] << " ";
    opserr << endln;
}

int
TclModelBuilderNDMaterialCommand (ClientData clientData, Tcl_Interp *interp, int argc,
				  TCL_Char **argv, TclModelBuilder *theTclBuilder)
{
    // Make sure there is a minimum number of arguments
    if (argc < 3) {
	opserr << "WARNING insufficient number of ND material arguments\n";
	opserr << "Want: nDMaterial type? tag? <specific material args>" << endln;
	return TCL_ERROR;
    }

    // Pointer to an ND material that will be added to the model builder
    NDMaterial *theMaterial = 0;

    // Check argv[1] for ND material type

    //Jul. 07, 2001 Boris Jeremic & ZHaohui Yang jeremic|zhyang@ucdavis.edu
    // Pressure dependent elastic material
    if (strcmp(argv[1],"PressureDependentElastic3D") == 0) {
	if (argc < 6) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: nDMaterial PressureDependentElastic3D tag? E? v? rho?" << endln;
	    return TCL_ERROR;
	}

	int tag = 0;
	double E = 0.0;
	double v = 0.0;
	double rho = 0.0;
	double expp = 0.0;
	double prp = 0.0;
	double pop = 0.0;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid PressureDependentElastic3D tag" << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
	    opserr << "WARNING invalid E\n";
	    opserr << "nDMaterial PressureDependentElastic3D: E" << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[4], &v) != TCL_OK) {
	    opserr << "WARNING invalid v\n";
	    opserr << "nDMaterial PressureDependentElastic3D: v" << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[5], &rho) != TCL_OK) {
	    opserr << "WARNING invalid v\n";
	    opserr << "nDMaterial PressureDependentElastic3D: rho" << tag << endln;
	    return TCL_ERROR;
	}

//////////////////////////////////////////////////////////////////////////////////
	if( argc == 6 )
	{
     	   theMaterial = new PressureDependentElastic3D (tag, E, v, rho);
	   //opserr << "nDMaterial PressureDependentElastic3D: expp =" << expp << endln;
	}
//////////////////////////////////////////////////////////////////////////////////
	else if( argc == 7 )
	{
	   //get the exponent of the pressure sensitive elastic material)
	   if (Tcl_GetDouble(interp, argv[6], &expp) != TCL_OK) {
	       opserr << "WARNING invalid v\n";
	       opserr << "nDMaterial PressureDependentElastic3D: " << tag << endln;
	       return TCL_ERROR;
	   }
     	   theMaterial = new PressureDependentElastic3D (tag, E, v, rho, expp);
	   //opserr << "nDMaterial PressureDependentElastic3D: expp =" << expp << endln;
	}
//////////////////////////////////////////////////////////////////////////////////
	else if (argc == 8 )
	{
	   //get the exponent pressure of the pressure sensitive elastic material)
	   if (Tcl_GetDouble(interp, argv[6], &expp) != TCL_OK) {
	       opserr << "WARNING invalid v\n";
	       opserr << "nDMaterial PressureDependentElastic3D: expp" << tag << endln;
	       return TCL_ERROR;
	   }
	   //get the reference pressure of the pressure sensitive elastic material)
	   if (Tcl_GetDouble(interp, argv[7], &prp) != TCL_OK) {
	       opserr << "WARNING invalid v\n";
	       opserr << "nDMaterial PressureDependentElastic3D: prp " << tag << endln;
	       return TCL_ERROR;
	   }
	   //opserr << "nDMaterial ElasticIsotropic3D: prp =" << prp << endln;
     	   theMaterial = new PressureDependentElastic3D (tag, E, v, rho, expp, prp);
	}
//////////////////////////////////////////////////////////////////////////////////
	else if (argc >= 9 )
	{
	   //get the exponent of the pressure sensitive elastic material)
	   if (Tcl_GetDouble(interp, argv[6], &expp) != TCL_OK) {
	       opserr << "WARNING invalid v\n";
	       opserr << "nDMaterial PressureDependentElastic3D: expp" << tag << endln;
	       return TCL_ERROR;
	   }
	   //get the reference pressure of the pressure sensitive elastic material)
	   if (Tcl_GetDouble(interp, argv[7], &prp) != TCL_OK) {
	       opserr << "WARNING invalid v\n";
	       opserr << "nDMaterial PressureDependentElastic3D: prp" << tag << endln;
	       return TCL_ERROR;
	   }
	   //get the cutoff pressure po of the pressure sensitive elastic material)
	   if (Tcl_GetDouble(interp, argv[8], &pop) != TCL_OK) {
	       opserr << "WARNING invalid v\n";
	       opserr << "nDMaterial PressureDependentElastic3D: pop" << tag << endln;
	       return TCL_ERROR;
	   }
	   //opserr << "nDMaterial PressureDependentElastic3D: pop =" << pop << endln;
     	   theMaterial = new PressureDependentElastic3D (tag, E, v, rho, expp, prp, pop);
	}

    }
    //Jul. 07, 2001 Boris Jeremic & ZHaohui Yang jeremic|zhyang@ucdavis.edu
    // Linear Elastic Material (non-pressure dependent)
    else if ( strcmp(argv[1],"ElasticIsotropic3D") == 0 )
    {
	if (argc < 6) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: nDMaterial ElasticIsotropic3D tag? E? v? <rho?>" << endln;
	    return TCL_ERROR;
	}

	int tag = 0;
	double E = 0.0;
	double v = 0.0;
	double rho = 0.0;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid ElasticIsotropic3D tag" << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
	    opserr << "WARNING invalid E\n";
	    opserr << "nDMaterial ElasticIsotropic3D: " << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[4], &v) != TCL_OK) {
	    opserr << "WARNING invalid v\n";
	    opserr << "nDMaterial ElasticIsotropic3D: " << tag << endln;
	    return TCL_ERROR;
	}


 if (argc > 5 && Tcl_GetDouble(interp, argv[5], &rho) != TCL_OK)
   {
     opserr << "WARNING invalid rho\n";
     opserr << "nDMaterial ElasticIsotropic: " << tag << endln;
     return TCL_ERROR;
   }

	theMaterial = new ElasticIsotropic3D (tag, E, v, rho);

    }

    else if (strcmp(argv[1],"ElasticIsotropic") == 0) {
	if (argc < 5) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: nDMaterial ElasticIsotropic tag? E? v? <rho?>" << endln;
	    return TCL_ERROR;
	}

	int tag;
	double E, v;
	double rho = 0.0;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid ElasticIsotropic tag" << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
	    opserr << "WARNING invalid E\n";
	    opserr << "nDMaterial ElasticIsotropic: " << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[4], &v) != TCL_OK) {
	    opserr << "WARNING invalid v\n";
	    opserr << "nDMaterial ElasticIsotropic: " << tag << endln;
	    return TCL_ERROR;
	}

 	if (argc > 5 && Tcl_GetDouble(interp, argv[5], &rho) != TCL_OK)
	   {
	     opserr << "WARNING invalid rho\n";
	     opserr << "nDMaterial ElasticIsotropic: " << tag << endln;
	     return TCL_ERROR;
	   }

	theMaterial = new ElasticIsotropicMaterial (tag, E, v, rho);
    }

//March 20, 2003 Zhaohui Yang, Yi Bian, Boris Jeremic Anisotropic Elastic Material Model
    else if (strcmp(argv[1],"ElasticCrossAnisotropic") == 0) {
 //cout << "argc" << argc;
 if (argc < 8) {
     opserr << "WARNING insufficient arguments\n";
     printCommand(argc,argv);
     opserr << "Want: nDMaterial ElasticCrossAnisotropic tag? Ehh? Ehv? nuhv? nuvv? Ghv? <rho?>" << endln;
     return TCL_ERROR;
 }

 int tag;
 double Eh, Ev, nuhv, nuhh, Ghv;
 double rho = 0.0;

 if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
     opserr << "WARNING invalid ElasticCrossAnisotropic tag" << endln;
     return TCL_ERROR;
 }

 if (Tcl_GetDouble(interp, argv[3], &Eh) != TCL_OK) {
     opserr << "WARNING invalid Eh\n";
     opserr << "nDMaterial ElasticCrossAnisotropic: " << tag << endln;
     return TCL_ERROR;
 }

 if (Tcl_GetDouble(interp, argv[4], &Ev) != TCL_OK) {
     opserr << "WARNING invalid Ev\n";
     opserr << "nDMaterial ElasticCrossAnisotropic: " << tag << endln;
     return TCL_ERROR;
 }

 if (Tcl_GetDouble(interp, argv[5], &nuhv) != TCL_OK) {
     opserr << "WARNING invalid nuhv\n";
     opserr << "nDMaterial ElasticCrossAnisotropic: " << tag << endln;
     return TCL_ERROR;
 }

 if (Tcl_GetDouble(interp, argv[6], &nuhh) != TCL_OK) {
     opserr << "WARNING invalid nuhh\n";
     opserr << "nDMaterial ElasticCrossAnisotropic: " << tag << endln;
     return TCL_ERROR;
 }

 if (Tcl_GetDouble(interp, argv[7], &Ghv) != TCL_OK) {
     opserr << "WARNING invalid Ghv\n";
     opserr << "nDMaterial ElasticCrossAnisotropic: " << tag << endln;
     return TCL_ERROR;
 }

        if (argc > 8 && Tcl_GetDouble(interp, argv[8], &rho) != TCL_OK) {
            opserr << "WARNING invalid rho\n";
            opserr << "nDMaterial ElasticCrossAnisotropic: " << tag << endln;
            return TCL_ERROR;
        }

 //cout << "Ev" << Ev << " Eh " << Eh << " nuhv " << nuhv << " nuhh "<< nuhh << " Ghv " << Ghv << "rho " <<rho << "\n";
 theMaterial = new ElasticCrossAnisotropic(tag, Eh, Ev, nuhv, nuhh, Ghv, rho);
    }


    // Check argv[1] for J2PlaneStrain material type
    else if ((strcmp(argv[1],"J2Plasticity") == 0)  ||
	     (strcmp(argv[1],"J2") == 0)) {
	if (argc < 9) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: nDMaterial J2Plasticity tag? K? G? sig0? sigInf? delta? H? <eta?>" << endln;
	    return TCL_ERROR;
	}

	int tag;
	double K, G, sig0, sigInf, delta, H;
	double eta = 0.0;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid J2Plasticity tag" << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[3], &K) != TCL_OK) {
	    opserr << "WARNING invalid K\n";
	    opserr << "nDMaterial J2Plasticity: " << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[4], &G) != TCL_OK) {
	    opserr << "WARNING invalid G\n";
	    opserr << "nDMaterial J2Plasticity: " << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[5], &sig0) != TCL_OK) {
	    opserr << "WARNING invalid sig0\n";
	    opserr << "nDMaterial J2Plasticity: " << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[6], &sigInf) != TCL_OK) {
	    opserr << "WARNING invalid sigInf\n";
	    opserr << "nDMaterial J2Plasticity: " << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[7], &delta) != TCL_OK) {
	    opserr << "WARNING invalid delta\n";
	    opserr << "nDMaterial J2Plasticity: " << tag << endln;
	    return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[8], &H) != TCL_OK) {
	    opserr << "WARNING invalid H\n";
	    opserr << "nDMaterial J2Plasticity: " << tag << endln;
	    return TCL_ERROR;
	}
	if (argc > 9 && Tcl_GetDouble(interp, argv[9], &eta) != TCL_OK) {
	    opserr << "WARNING invalid eta\n";
	    opserr << "nDMaterial J2Plasticity: " << tag << endln;
	    return TCL_ERROR;
	}

	theMaterial = new J2Plasticity (tag, 0, K, G, sig0, sigInf,
					delta, H, eta);
    }


    //
    //  MultiAxialCyclicPlasticity Model   by Gang Wang
    //
    //  nDMaterial MultiaxialCyclicPlasticity $tag, $rho, $K, $G,
    //      $Su , $Ho , $h, $m, $beta, $KCoeff
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // Check argv[1] for MultiaxialCyclicPlasticity material type
    else if ((strcmp(argv[1],"MultiaxialCyclicPlasticity") == 0)  ||
	     (strcmp(argv[1],"MCP") == 0)) {
      if (argc < 12) {
	opserr << "WARNING insufficient arguments\n";
	printCommand(argc,argv);
	opserr << "Want: nDMaterial MultiaxialCyclicPlasticity tag? rho? K? G? Su? Ho? h? m? beta? KCoeff? <eta?>" << endln;
	return TCL_ERROR;
      }

      int tag;
      double K, G, rho, Su, Ho, h, m, beta, Kcoeff;
      double eta = 0.0;

      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING invalid MultiaxialCyclicPlasticity tag" << endln;
	return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[3], &rho) != TCL_OK) {
	opserr << "WARNING invalid rho\n";
	opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << endln;
	return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[4], &K) != TCL_OK) {
	opserr << "WARNING invalid K\n";
	opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << endln;
	return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[5], &G) != TCL_OK) {
	opserr << "WARNING invalid G\n";
	opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << endln;
	return TCL_ERROR;
      }


      if (Tcl_GetDouble(interp, argv[6], &Su) != TCL_OK) {
	opserr << "WARNING invalid alpha1\n";
	opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << endln;
	return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[7], &Ho) != TCL_OK) {
	opserr << "WARNING invalid Ho\n";
	opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << endln;
	return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[8], &h) != TCL_OK) {
	opserr << "WARNING invalid h\n";
	opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << endln;
	return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[9], &m) != TCL_OK) {
	opserr << "WARNING invalid m\n";
	opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << endln;
	return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[10], &beta) != TCL_OK) {
	opserr << "WARNING invalid beta\n";
	opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[11], &Kcoeff) != TCL_OK) {
	opserr << "WARNING invalid Kcoeff\n";
	opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << endln;
	return TCL_ERROR;
      }


      if (argc > 12 && Tcl_GetDouble(interp, argv[12], &eta) != TCL_OK) {
	opserr << "WARNING invalid eta\n";
	opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << endln;
	return TCL_ERROR;
      }

      theMaterial = new MultiaxialCyclicPlasticity (tag, 0, rho, K, G, Su, Ho, h,m,
						    beta, Kcoeff, eta);
    }


    // Pressure Independend Multi-yield, by ZHY
    else if (strcmp(argv[1],"PressureIndependMultiYield") == 0) {
	const int numParam = 6;
	const int totParam = 10;
	int tag;  double param[totParam];
	param[6] = 0.0;
	param[7] = 100.;
	param[8] = 0.0;
	param[9] = 20;

	char * arg[] = {"nd", "rho", "refShearModul", "refBulkModul",
			"cohesi", "peakShearStra",
			"frictionAng (=0)", "refPress (=100)", "pressDependCoe (=0.0)",
	    "numberOfYieldSurf (=20)"};
	if (argc < (3+numParam)) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: nDMaterial PressureIndependMultiYield tag? " << arg[0];
	    opserr << "? "<< "\n";
	    opserr << arg[1] << "? "<< arg[2] << "? "<< arg[3] << "? "<< "\n";
	    opserr << arg[4] << "? "<< arg[5] << "? "<< arg[6] << "? "<< "\n";
	    opserr << arg[7] << "? "<< arg[8] << "? "<< arg[9] << "? "<<endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid PressureIndependMultiYield tag" << endln;
	    return TCL_ERROR;
	}

	for (int i=3; (i<argc && i<13); i++)
	    if (Tcl_GetDouble(interp, argv[i], &param[i-3]) != TCL_OK) {
		    opserr << "WARNING invalid " << arg[i-3] << "\n";
		    opserr << "nDMaterial PressureIndependMultiYield: " << tag << endln;
		    return TCL_ERROR;
	    }

	static double * gredu = 0;
	// user defined yield surfaces
	if (param[9] < 0 && param[9] > -40) {
     param[9] = -int(param[9]);
     gredu = new double[int(2*param[9])];
		 for (int i=0; i<2*param[9]; i++)
	      if (Tcl_GetDouble(interp, argv[i+13], &gredu[i]) != TCL_OK) {
		      opserr << "WARNING invalid " << arg[i-3] << "\n";
		      opserr << "nDMaterial PressureIndependMultiYield: " << tag << endln;
		      return TCL_ERROR;
				}
  }

	PressureIndependMultiYield * temp =
	    new PressureIndependMultiYield (tag, param[0], param[1], param[2],
					    param[3], param[4], param[5], param[6],
					    param[7], param[8], param[9], gredu);
	theMaterial = temp;

	   if (gredu != 0) {
		   delete [] gredu;
		   gredu = 0;
	   }
    }

    // Pressure Dependend Multi-yield, by ZHY
    else if (strcmp(argv[1],"PressureDependMultiYield") == 0) {
	const int numParam = 15;
	const int totParam = 24;
	int tag;
	double param[totParam];
 	param[15] = 20;
 	param[16] = 0.6;
	param[17] = 0.9;
	param[18] = 0.02;
	param[19] = 0.7;
	param[20] = 101.;
	param[21] = .3;
	param[22] = 0.;
	param[23] = 1.;

	char * arg[] = {"nd", "rho", "refShearModul",
		  "refBulkModul", "frictionAng",
			"peakShearStra", "refPress", "pressDependCoe",
			"phaseTransformAngle", "contractionParam1",
			"dilationParam1", "dilationParam2",
			"liquefactionParam1", "liquefactionParam2",
			"liquefactionParam4", "numberOfYieldSurf (=20)",
			"e (=0.6)", "volLimit1 (=0.9)", "volLimit2 (=0.02)",
			"volLimit3 (=0.7)", "Atmospheric pressure (=101)", "cohesi (=.5)",
	        "Hv (=0)", "Pv (=1.)" };
	if (argc < (3+numParam)) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: nDMaterial PressureDependMultiYield tag? "<< arg[0];
	    opserr << "? "<< "\n";
	    opserr << arg[1] << "? "<< arg[2] << "? "<< arg[3] << "? "<< "\n";
	    opserr << arg[4] << "? "<< arg[5] << "? "<< arg[6] << "? "<< "\n";
	    opserr << arg[7] << "? "<< arg[8] << "? "<< arg[9] << "? "<< "\n";
	    opserr << arg[10] << "? "<< arg[11] << "? "<< arg[12] << "? "<< "\n";
	    opserr << arg[13] << "? "<< arg[14] << "? "<< arg[15] << "? "<< "\n";
	    opserr << arg[16] << "? "<< arg[17] << "? "<< arg[18] << "? "<< "\n";
	    opserr << arg[19] << "? "<< arg[20] << "? "<< arg[21] << "? "<< endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid PressureDependMultiYield tag" << endln;
	    return TCL_ERROR;
	}

	for (int i=3; (i<argc && i<19); i++)
	  if (Tcl_GetDouble(interp, argv[i], &param[i-3]) != TCL_OK) {
		    opserr << "WARNING invalid " << arg[i-3] << "\n";
		    opserr << "nDMaterial PressureDependMultiYield: " << tag << endln;
		    return TCL_ERROR;
	  }

	static double * gredu = 0;
	// user defined yield surfaces
	if (param[15] < 0 && param[15] > -40) {
     param[15] = -int(param[15]);
     gredu = new double[int(2*param[15])];

		 for (int i=0; i<2*param[15]; i++)
	      if (Tcl_GetDouble(interp, argv[i+19], &gredu[i]) != TCL_OK) {
		      opserr << "WARNING invalid " << arg[i-3] << "\n";
		      opserr << "nDMaterial PressureIndependMultiYield: " << tag << endln;
		      return TCL_ERROR;
		  }
	}

	if (gredu != 0) {
	  for (int i=19+int(2*param[15]); i<argc; i++)
	    if (Tcl_GetDouble(interp, argv[i], &param[i-3-int(2*param[15])]) != TCL_OK) {
		      opserr << "WARNING invalid " << arg[i-3-int(2*param[15])] << "\n";
		      opserr << "nDMaterial PressureDependMultiYield: " << tag << endln;
		      return TCL_ERROR;
			}
	} else {
	  for (int i=19; i<argc; i++)
	    if (Tcl_GetDouble(interp, argv[i], &param[i-3]) != TCL_OK) {
		      opserr << "WARNING invalid " << arg[i-3-int(2*param[15])] << "\n";
		      opserr << "nDMaterial PressureDependMultiYield: " << tag << endln;
		      return TCL_ERROR;
		}
	}

	PressureDependMultiYield * temp =
	    new PressureDependMultiYield (tag, param[0], param[1], param[2],
					  param[3], param[4], param[5],
					  param[6], param[7], param[8],
					  param[9], param[10], param[11],
					  param[12], param[13], param[14],
					  param[15], gredu, param[16], param[17],
					  param[18], param[19], param[20], param[21], param[22], param[23]);

	   theMaterial = temp;
	   if (gredu != 0) {
		   delete [] gredu;
		   gredu = 0;
	   }
	}

    // Pressure Dependend Multi-yield, by ZHY
    else if (strcmp(argv[1],"PressureDependMultiYield02") == 0) {
	const int numParam = 13;
	const int totParam = 26;
	int tag;
	double param[totParam];
	param[numParam] = 20;
 	param[numParam+1] = 5.0;
 	param[numParam+2] = 3.;
 	param[numParam+3] = 1.;
	param[numParam+4] = 0.;
 	param[numParam+5] = 0.6;
	param[numParam+6] = 0.9;
	param[numParam+7] = 0.02;
	param[numParam+8] = 0.7;
	param[numParam+9] = 101.;
	param[numParam+10] = 0.1;
	param[numParam+11] = 0.;
	param[numParam+12] = 1.;

	char * arg[] = {"nd", "rho", "refShearModul",
		  "refBulkModul", "frictionAng",
			"peakShearStra", "refPress", "pressDependCoe",
			"phaseTransformAngle", "contractionParam1",
			"contractionParam3","dilationParam1","dilationParam3",
			"numberOfYieldSurf (=20)",
			"contractionParam2=5.0", "dilationParam2=3.0",
			"liquefactionParam1=1.0", "liquefactionParam2=0.0",
			"e (=0.6)", "volLimit1 (=0.9)", "volLimit2 (=0.02)",
			"volLimit3 (=0.7)", "Atmospheric pressure (=101)", "cohesi (=.1)",
	        "Hv (=0)", "Pv (=1.)" };
	if (argc < (3+numParam)) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: nDMaterial PressureDependMultiYield02 tag? "<< arg[0];
	    opserr << "? "<< "\n";
	    opserr << arg[1] << "? "<< arg[2] << "? "<< arg[3] << "? "<< "\n";
	    opserr << arg[4] << "? "<< arg[5] << "? "<< arg[6] << "? "<< "\n";
	    opserr << arg[7] << "? "<< arg[8] << "? "<< arg[9] << "? "<< "\n";
	    opserr << arg[10] << "? "<< arg[11] << "? "<< arg[12] << "? "<< "\n";
	    opserr << arg[13] << "? "<< arg[14] << "? "<< arg[15] << "? "<< "\n";
	    opserr << arg[16] << "? "<< arg[17] << "? "<< arg[18] << "? "<< "\n";
	    opserr << arg[19] << "? "<< arg[20] << "? "<< arg[21] << "? "<< "\n";
		opserr << arg[22] << "? "<< arg[23] << "? " << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid PressureDependMultiYield02 tag" << endln;
	    return TCL_ERROR;
	}

	int in = 17;
	for (int i=3; (i<argc && i<in); i++)
	  if (Tcl_GetDouble(interp, argv[i], &param[i-3]) != TCL_OK) {
		    opserr << "WARNING invalid " << arg[i-3] << "\n";
		    opserr << "nDMaterial PressureDependMultiYield02: " << tag << endln;
		    return TCL_ERROR;
	  }

	static double * gredu = 0;

	// user defined yield surfaces
	if (param[numParam] < 0 && param[numParam] > -100) {
     param[numParam] = -int(param[numParam]);
     gredu = new double[int(2*param[numParam])];

		 for (int i=0; i<2*param[numParam]; i++)
	      if (Tcl_GetDouble(interp, argv[i+in], &gredu[i]) != TCL_OK) {
		      opserr << "WARNING invalid " << arg[i-3] << "\n";
		      opserr << "nDMaterial PressureIndependMultiYield: " << tag << endln;
		      return TCL_ERROR;
		  }
	}

	if (gredu != 0) {
	  for (int i=in+int(2*param[numParam]); i<argc; i++)
	    if (Tcl_GetDouble(interp, argv[i], &param[i-3-int(2*param[numParam])]) != TCL_OK) {
		      opserr << "WARNING invalid " << arg[i-3-int(2*param[numParam])] << "\n";
		      opserr << "nDMaterial PressureDependMultiYield02: " << tag << endln;
		      return TCL_ERROR;
			}
	} else {
	  for (int i=in; i<argc; i++)
	    if (Tcl_GetDouble(interp, argv[i], &param[i-3]) != TCL_OK) {
		      opserr << "WARNING invalid " << arg[i-3-int(2*param[numParam])] << "\n";
		      opserr << "nDMaterial PressureDependMultiYield02: " << tag << endln;
		      return TCL_ERROR;
		}
	}


	PressureDependMultiYield02 * temp =
	    new PressureDependMultiYield02 (tag, param[0], param[1], param[2],
					  param[3], param[4], param[5],
					  param[6], param[7], param[8],
					  param[9], param[10], param[11],
					  param[12], param[13], gredu, param[14],
					  param[15], param[16], param[17],
					  param[18], param[19], param[20], param[21],
					  param[22], param[23], param[24], param[25]);

	   theMaterial = temp;
	   if (gredu != 0) {
		   delete [] gredu;
		   gredu = 0;
	   }
  }

    // Fluid Solid Porous, by ZHY
    else if (strcmp(argv[1],"FluidSolidPorous") == 0) {

	int tag;  double param[4];
	char * arg[] = {"nd", "soilMatTag", "combinedBulkModul", "Atmospheric pressure"};
	if (argc < 6) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: nDMaterial FluidSolidPorous tag? "<< arg[0];
	    opserr << "? "<< "\n";
	    opserr << arg[1] << "? "<< arg[2] << "? "<< endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid FluidSolidPorous tag" << endln;
	    return TCL_ERROR;
	}

	for (int i=3; i<6; i++)
	  if (Tcl_GetDouble(interp, argv[i], &param[i-3] ) != TCL_OK) {
	      opserr << "WARNING invalid " << arg[i-3] << "\n";
	      opserr << "nDMaterial FluidSolidPorous: " << tag << endln;
	      return TCL_ERROR;
	  }

	NDMaterial *soil = theTclBuilder->getNDMaterial(param[1]);
	if (soil == 0) {
	      opserr << "WARNING FluidSolidPorous: couldn't get soil material ";
	      opserr << "tagged: " << param[1] << "\n";
	      return TCL_ERROR;
	}

	param[3] = 101.;
	if (argc == 7) {
	  if (Tcl_GetDouble(interp, argv[6], &param[3] ) != TCL_OK) {
	      opserr << "WARNING invalid " << arg[3] << "\n";
	      opserr << "nDMaterial FluidSolidPorous: " << tag << endln;
	      return TCL_ERROR;
	  }
	}

	theMaterial = new FluidSolidPorousMaterial (tag, param[0], *soil,
						    param[2],param[3]);
    }

    else if (strcmp(argv[1],"Template3Dep") == 0) {
      theMaterial = TclModelBuilder_addTemplate3Dep(clientData, interp, argc, argv,
						    theTclBuilder, 2);
    }

    else if (strcmp(argv[1],"NewTemplate3Dep") == 0) {
      theMaterial = TclModelBuilder_addNewTemplate3Dep(clientData, interp, argc, argv,
						    theTclBuilder, 2);
    }
        
    else if (strcmp(argv[1],"FiniteDeformationElastic3D") == 0 ||
             strcmp(argv[1],"FDElastic3D" ) == 0) {
      theMaterial = TclModelBuilder_addFiniteDeformationElastic3D(clientData, interp, argc, argv,
						    theTclBuilder, 1);
    }

    else if (strcmp(argv[1],"FiniteDeformationEP3D") == 0 ||
             strcmp(argv[1],"FDEP3D" ) == 0) {
      theMaterial = TclModelBuilder_addFiniteDeformationEP3D(clientData, interp, argc, argv,
						    theTclBuilder, 2);
    }


     else if (strcmp(argv[1],"PlaneStressMaterial") == 0 ||
 	     strcmp(argv[1],"PlaneStress") == 0) {
 	if (argc < 4) {
 	    opserr << "WARNING insufficient arguments\n";
 	    printCommand(argc,argv);
 	    opserr << "Want: nDMaterial PlaneStress tag? matTag?" << endln;
 	    return TCL_ERROR;
 	}

 	int tag, matTag;

 	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
 	    opserr << "WARNING invalid nDMaterial PlaneStress tag" << endln;
 	    return TCL_ERROR;
 	}

 	if (Tcl_GetInt (interp, argv[3], &matTag) != TCL_OK) {
 	    opserr << "WARNING invalid matTag" << endln;
 	    opserr << "PlaneStress: " << matTag << endln;
 	    return TCL_ERROR;
 	}

 	NDMaterial *threeDMaterial = theTclBuilder->getNDMaterial(matTag);
 	if (threeDMaterial == 0) {
 	    opserr << "WARNING nD material does not exist\n";
 	    opserr << "nD material: " << matTag;
 	    opserr << "\nPlaneStress nDMaterial: " << tag << endln;
 	    return TCL_ERROR;
 	}

 	theMaterial = new PlaneStressMaterial( tag, *threeDMaterial );
     }


     else if (strcmp(argv[1],"PlateFiberMaterial") == 0 ||
 	     strcmp(argv[1],"PlateFiber") == 0) {
 	if (argc < 4) {
 	    opserr << "WARNING insufficient arguments\n";
 	    printCommand(argc,argv);
 	    opserr << "Want: nDMaterial PlateFiber tag? matTag?" << endln;
 	    return TCL_ERROR;
 	}

 	int tag, matTag;

 	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
 	    opserr << "WARNING invalid nDMaterial PlateFiber tag" << endln;
 	    return TCL_ERROR;
 	}

 	if (Tcl_GetInt (interp, argv[3], &matTag) != TCL_OK) {
 	    opserr << "WARNING invalid matTag" << endln;
 	    opserr << "PlateFiber: " << matTag << endln;
 	    return TCL_ERROR;
 	}

 	NDMaterial *threeDMaterial = theTclBuilder->getNDMaterial(matTag);
 	if (threeDMaterial == 0) {
 	    opserr << "WARNING nD material does not exist\n";
 	    opserr << "nD material: " << matTag;
 	    opserr << "\nPlateFiber nDMaterial: " << tag << endln;
 	    return TCL_ERROR;
 	}

 	theMaterial = new PlateFiberMaterial( tag, *threeDMaterial );
     }

     else if (strcmp(argv[1],"BeamFiberMaterial") == 0 ||
 	     strcmp(argv[1],"BeamFiber") == 0) {
 	if (argc < 4) {
 	    opserr << "WARNING insufficient arguments\n";
 	    printCommand(argc,argv);
 	    opserr << "Want: nDMaterial BeamFiber tag? matTag?" << endln;
 	    return TCL_ERROR;
 	}

 	int tag, matTag;

 	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
 	    opserr << "WARNING invalid nDMaterial BeamFiber tag" << endln;
 	    return TCL_ERROR;
 	}

 	if (Tcl_GetInt (interp, argv[3], &matTag) != TCL_OK) {
 	    opserr << "WARNING invalid matTag" << endln;
 	    opserr << "BeamFiber: " << matTag << endln;
 	    return TCL_ERROR;
 	}

 	NDMaterial *threeDMaterial = theTclBuilder->getNDMaterial(matTag);
 	if (threeDMaterial == 0) {
 	    opserr << "WARNING nD material does not exist\n";
 	    opserr << "nD material: " << matTag;
 	    opserr << "\nBeamFiber nDMaterial: " << tag << endln;
 	    return TCL_ERROR;
 	}

 	theMaterial = new BeamFiberMaterial( tag, *threeDMaterial );
     }

    else if (strcmp(argv[1],"Bidirectional") == 0) {
      opserr << "nDMaterial Bidirectional is now a section model, please "
	   << "change to \'section Bidirectional\'" << endln;
      return TCL_ERROR;
    }

    else {
      theMaterial = TclModelBuilder_addFeapMaterial(clientData,
						    interp,
						    argc,
						    argv,
						    theTclBuilder);
    }

    if (theMaterial == 0) {
	opserr << "WARNING count not create nDMaterial: " << argv[1];
	return TCL_ERROR;
    }

    // Now add the material to the modelBuilder
    if (theTclBuilder->addNDMaterial(*theMaterial) < 0) {
	opserr << "WARNING could not add material to the domain\n";
	opserr << *theMaterial << endln;
	delete theMaterial; // invoke the material objects destructor, otherwise mem leak
	return TCL_ERROR;
    }

    return TCL_OK;
}


