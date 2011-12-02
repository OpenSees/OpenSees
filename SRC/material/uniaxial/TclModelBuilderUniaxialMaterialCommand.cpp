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
                                                                        
// $Revision: 1.18 $
// $Date: 2003-03-06 18:34:15 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/TclModelBuilderUniaxialMaterialCommand.cpp,v $
                                                                        
                                                                        
// Written: fmk, MHS 
// Created: 07/99
//
// Description: This file contains the function invoked when the user invokes
// the uniaxialMaterial command in the interpreter. 
//
// What: "@(#) TclModelBuilderUniaxialMaterialCommand.C, revA"

#include <TclModelBuilder.h>

#include <ElasticMaterial.h>	// fmk
#include <ElasticPPMaterial.h>	// fmk
#include <ParallelMaterial.h>	// fmk
#include <HardeningMaterial.h>	// MHS
#include <Steel01.h>			// MHS
#include <Concrete01.h>			// MHS
#include <HystereticMaterial.h>	// MHS
#include <EPPGapMaterial.h>		// Mackie
#include <ViscousMaterial.h>	// Sasani
#include <PathIndependentMaterial.h>	// MHS
#include <MinMaxMaterial.h>	// MHS
#include <SeriesMaterial.h>		// MHS
#include <ENTMaterial.h>		// MHS
#include <CableMaterial.h>	// CC
#include <BoucWenMaterial.h>	// Terje

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
				  TCL_Char **argv, TclModelBuilder *theTclBuilder);

UniaxialMaterial *
TclModelBuilder_addDrainMaterial(ClientData clientData, Tcl_Interp *interp, int argc, 
				 TCL_Char **argv, TclModelBuilder *theTclBuilder);

UniaxialMaterial *
TclModelBuilder_addSnapMaterial(ClientData clientData, Tcl_Interp *interp, int argc, 
				TCL_Char **argv, TclModelBuilder *theTclBuilder);

UniaxialMaterial *
TclModelBuilder_addPyTzQzMaterial(ClientData clientData, Tcl_Interp *interp, int argc, 
				  TCL_Char **argv, TclModelBuilder *theTclBuilder);

int
TclModelBuilderUniaxialMaterialCommand (ClientData clientData, Tcl_Interp *interp, int argc,
					TCL_Char **argv, TclModelBuilder *theTclBuilder)
{
    // Make sure there is a minimum number of arguments
    if (argc < 3) {
	opserr << "WARNING insufficient number of uniaxial material arguments\n";
	opserr << "Want: uniaxialMaterial type? tag? <specific material args>" << endln;
	return TCL_ERROR;
    }
    
    // Pointer to a uniaxial material that will be added to the model builder
    UniaxialMaterial *theMaterial = 0;

    // Check argv[2] for uniaxial material type
    if (strcmp(argv[1],"Elastic") == 0) {
	if (argc < 4 || argc > 5) {
	    opserr << "WARNING invalid number of arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: uniaxialMaterial Elastic tag? E? <eta?>" << endln;
	    return TCL_ERROR;
	}    

	int tag;
	double E;
        double eta = 0.0;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid uniaxialMaterial Elastic tag" << endln;
	    return TCL_ERROR;		
	}

	if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
	    opserr << "WARNING invalid E\n";
	    opserr << "uniaxiaMaterial Elastic: " << tag << endln;
	    return TCL_ERROR;	
	}

	if (argc == 5) {
	    if (Tcl_GetDouble(interp,argv[4], &eta) != TCL_OK) {
               opserr << "WARNING invalid eta\n";
               opserr << "uniaxialMaterial Elastic: " << tag << endln;
               return TCL_ERROR;
            }
	}

	// Parsing was successful, allocate the material
	theMaterial = new ElasticMaterial(tag, E, eta);       
        //dum << tag " << tag << " E " << E << " eta " << eta <<endln;
    }

	else if (strcmp(argv[1],"ENT") == 0) {
	if (argc < 4) {
	    opserr << "WARNING invalid number of arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: uniaxialMaterial ENT tag? E?" << endln;
	    return TCL_ERROR;
	}    

	int tag;
	double E;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid uniaxialMaterial ENT tag" << endln;
	    return TCL_ERROR;		
	}

	if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
	    opserr << "WARNING invalid E\n";
	    opserr << "uniaxiaMaterial ENT: " << tag << endln;
	    return TCL_ERROR;	
	}

	// Parsing was successful, allocate the material
	theMaterial = new ENTMaterial(tag, E);       
     
	}

    else if (strcmp(argv[1],"ElasticPP") == 0) {
	if (argc < 5) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: uniaxialMaterial ElasticPP tag? E? epsy?" << endln;
	    return TCL_ERROR;
	}

	int tag;
	double E, ep;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid uniaxialMaterial ElasticPP tag" << endln;
	    return TCL_ERROR;		
	}

	if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
	    opserr << "WARNING invalid E\n";
	    opserr << "uniaxialMaterial ElasticPP: " << tag << endln;
	    return TCL_ERROR;	
	}

	if (Tcl_GetDouble(interp, argv[4], &ep) != TCL_OK) {
	    opserr << "WARNING invalid epsy\n";
	    opserr << "uniaxialMaterial ElasticPP: " << tag << endln;
	    return TCL_ERROR;
	}
	
	// read in additional parameters eyn and ezero
	if (argc > 6) {
	    double eyn, ez;
	    if (Tcl_GetDouble(interp, argv[5], &eyn) != TCL_OK) {
		opserr << "WARNING invalid eyn\n";
		opserr << "uniaxialMaterial ElasticPP: " << tag << endln;
		return TCL_ERROR;
	    }	    
	    if (Tcl_GetDouble(interp, argv[6], &ez) != TCL_OK) {
		opserr << "WARNING invalid ez\n";
		opserr << "uniaxialMaterial ElasticPP: " << tag << endln;
		return TCL_ERROR;
	    }	    	    
	    theMaterial = new ElasticPPMaterial(tag, E, ep, eyn, ez);
	} else {
	    // Parsing was successful, allocate the material
	    theMaterial = new ElasticPPMaterial(tag, E, ep);       
	}
    }

	else if (strcmp(argv[1],"ElasticPPGap") == 0) {
	if (argc < 6) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: uniaxialMaterial ElasticPPGap tag? E? fy? gap? <damage>" << endln;
	    return TCL_ERROR;
	}

	int tag;
	int damage = 0;
	double E, fy, gap;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid uniaxialMaterial ElasticPPGap tag" << endln;
	    return TCL_ERROR;		
	}

	if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
	    opserr << "WARNING invalid E\n";
	    opserr << "uniaxialMaterial ElasticPPGap: " << tag << endln;
	    return TCL_ERROR;	
	}

	if (Tcl_GetDouble(interp, argv[4], &fy) != TCL_OK) {
	    opserr << "WARNING invalid fy\n";
	    opserr << "uniaxialMaterial ElasticPPGap: " << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[5], &gap) != TCL_OK) {
	    opserr << "WARNING invalid gap\n";
	    opserr << "uniaxialMaterial ElasticPPGap: " << tag << endln;
	    return TCL_ERROR;
	}

	if (argc > 6 && strcmp(argv[6],"damage") == 0)
	  damage = 1;

	// Parsing was successful, allocate the material
	theMaterial = new EPPGapMaterial(tag, E, fy, gap, damage);       
    }

    else if (strcmp(argv[1],"Hardening") == 0) {
	if (argc < 7) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: uniaxialMaterial Hardening tag? E? sigmaY? H_iso? H_kin? <eta?>" << endln;
	    return TCL_ERROR;
	}

	int tag;
	double E, sigmaY, Hiso, Hkin;
	double eta = 0.0;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid uniaxialMaterial Hardening tag" << endln;
	    return TCL_ERROR;		
	}

	if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
	    opserr << "WARNING invalid E\n";
	    opserr << "uniaxialMaterial Hardening: " << tag << endln;
	    return TCL_ERROR;	
	}

	if (Tcl_GetDouble(interp, argv[4], &sigmaY) != TCL_OK) {
	    opserr << "WARNING invalid sigmaY\n";
	    opserr << "uniaxialMaterial Hardening: " << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[5], &Hiso) != TCL_OK) {
	    opserr << "WARNING invalid H_iso\n";
	    opserr << "uniaxialMaterial Hardening: " << tag << endln;
	    return TCL_ERROR;	
	}

	if (Tcl_GetDouble(interp, argv[6], &Hkin) != TCL_OK) {
	    opserr << "WARNING invalid H_kin\n";
	    opserr << "uniaxialMaterial Hardening: " << tag << endln;
	    return TCL_ERROR;	
	}

	if (argc > 7 && Tcl_GetDouble(interp, argv[7], &eta) != TCL_OK) {
	    opserr << "WARNING invalid eta\n";
	    opserr << "uniaxialMaterial Hardening: " << tag << endln;
	    return TCL_ERROR;	
	}

	// Parsing was successful, allocate the material
	theMaterial = new HardeningMaterial (tag, E, sigmaY, Hiso, Hkin, eta);
    }

    else if (strcmp(argv[1],"BoucWen") == 0) {
      if (argc < 12) {
	opserr << "WARNING insufficient arguments\n";
	printCommand(argc,argv);
	opserr << "Want: uniaxialMaterial BoucWen tag? alpha? ko? n? gamma?" << endln 
	       << " beta? Ao? deltaA? deltaNu? deltaEta?" << endln;
	return TCL_ERROR;
      }
      
      int tag;
      double alpha, ko, n, gamma, beta, Ao, deltaA, deltaNu, deltaEta;
      
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING invalid uniaxialMaterial BoucWen tag" << endln;
	return TCL_ERROR;		
      }
      if (Tcl_GetDouble(interp, argv[3], &alpha) != TCL_OK) {
	opserr << "WARNING invalid alpha\n";
	opserr << "uniaxialMaterial BoucWen: " << tag << endln;
	return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[4], &ko) != TCL_OK) {
	opserr << "WARNING invalid ko\n";
	opserr << "uniaxialMaterial BoucWen: " << tag << endln;
	return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[5], &n) != TCL_OK) {
	opserr << "WARNING invalid n\n";
	opserr << "uniaxialMaterial BoucWen: " << tag << endln;
	return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[6], &gamma) != TCL_OK) {
	opserr << "WARNING invalid gamma\n";
	opserr << "uniaxialMaterial BoucWen: " << tag << endln;
	return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[7], &beta) != TCL_OK) {
	opserr << "WARNING invalid beta\n";
	opserr << "uniaxialMaterial BoucWen: " << tag << endln;
	return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[8], &Ao) != TCL_OK) {
	opserr << "WARNING invalid Ao\n";
	opserr << "uniaxialMaterial BoucWen: " << tag << endln;
	return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[9], &deltaA) != TCL_OK) {
	opserr << "WARNING invalid deltaA\n";
	opserr << "uniaxialMaterial BoucWen: " << tag << endln;
	return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[10], &deltaNu) != TCL_OK) {
	opserr << "WARNING invalid deltaNu\n";
	opserr << "uniaxialMaterial BoucWen: " << tag << endln;
	return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[11], &deltaEta) != TCL_OK) {
	opserr << "WARNING invalid deltaEta\n";
	opserr << "uniaxialMaterial BoucWen: " << tag << endln;
	return TCL_ERROR;	
      }
      
      // Check if the user has given a tolerance for the Newton scheme
      double tolerance = 1.0e-8;
      if (argc > 12) {
	if (Tcl_GetDouble(interp, argv[12], &tolerance) != TCL_OK) {
	  opserr << "WARNING invalid tolerance\n";
	  opserr << "uniaxialMaterial BoucWen: " << tolerance << endln;
	  return TCL_ERROR;	
	}
      }
      
      // Check if the user has given a maxNumIter for the Newton scheme
      int maxNumIter = 20;
      if (argc > 13) {
	if (Tcl_GetInt(interp, argv[13], &maxNumIter) != TCL_OK) {
	  opserr << "WARNING invalid maxNumIter\n";
	  opserr << "uniaxialMaterial BoucWen: " << maxNumIter << endln;
	  return TCL_ERROR;	
	}
      }
      
      // Parsing was successful, allocate the material
      theMaterial = new BoucWenMaterial(tag, alpha, ko, n, gamma, beta, 
					Ao, deltaA, deltaNu, deltaEta,tolerance,maxNumIter);       
    }
    
    else if (strcmp(argv[1],"Parallel") == 0) {
	if (argc < 4) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: uniaxialMaterial Parallel tag? tag1? tag2? ...";
	    opserr << " <-min min?> <-max max?>" << endln;
	    return TCL_ERROR;
	}
 
	int tag;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid uniaxialMaterial Parallel tag" << endln;
	    return TCL_ERROR;		
	}

	int numMaterials = argc-3;
	
	if (numMaterials == 0) {
	    opserr << "WARNING no component material(s) provided\n";
	    opserr << "uniaxialMaterial Parallel: " << tag << endln;
	    return TCL_ERROR;
	}
    
	// Create an array to hold pointers to component materials
	UniaxialMaterial **theMats = new UniaxialMaterial *[numMaterials];
	
	// For each material get the tag and ensure it exists in model already
	for (int i=0; i<numMaterials; i++) {
	    int tagI;
	    if (Tcl_GetInt(interp, argv[i+3], &tagI) != TCL_OK) {
		opserr << "WARNING invalid component tag\n";
		opserr << "uniaxialMaterial Parallel: " << tag << endln;
		return TCL_ERROR;
	    }
	    
	    UniaxialMaterial *theMat = theTclBuilder->getUniaxialMaterial(tagI);
	    
	    if (theMat == 0) {
		opserr << "WARNING component material does not exist\n";
		opserr << "Component material: " << argv[i+3]; 
		opserr << "\nuniaxialMaterial Parallel: " << tag << endln;
		delete [] theMats;
		return TCL_ERROR;
	    }
	    else
		theMats[i] = theMat;
	}	
	
	// Parsing was successful, allocate the material
	theMaterial = new ParallelMaterial(tag, numMaterials, theMats);
	
	// Deallocate the temporary pointers
	delete [] theMats;
    }

    else if (strcmp(argv[1],"Series") == 0) {
	if (argc < 4) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: uniaxialMaterial Series tag? tag1? tag2? ..." << endln;
	    return TCL_ERROR;
	}
 
	int tag;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid uniaxialMaterial Series tag" << endln;
	    return TCL_ERROR;		
	}

	int numMaterials = argc - 3;
	UniaxialMaterial **theMats = new UniaxialMaterial *[numMaterials];

	// For each material get the tag and ensure it exists in model already
	for (int i = 3; i < argc; i++) {
	    int tagI;
	    if (Tcl_GetInt(interp, argv[i], &tagI) != TCL_OK) {
			opserr << "WARNING invalid component tag\n";
			opserr << "uniaxialMaterial Series: " << tag << endln;
			delete [] theMats;
			return TCL_ERROR;
	    }
	    
	    UniaxialMaterial *theMat = theTclBuilder->getUniaxialMaterial(tagI);
	    
	    if (theMat == 0) {
			opserr << "WARNING component material does not exist\n";
			opserr << "Component material: " << tagI;
			opserr << "\nuniaxialMaterial Series: " << tag << endln;
			delete [] theMats;
			return TCL_ERROR;
	    }
	    else
			theMats[i-3] = theMat;
	}	
	
	// Parsing was successful, allocate the material
	theMaterial = new SeriesMaterial(tag, numMaterials, theMats);
	
	// Deallocate the temporary pointers
	delete [] theMats;
    }

    else if (strcmp(argv[1],"Steel01") == 0) {
	// Check that there is the minimum number of arguments
	if (argc < 6) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: uniaxialMaterial Steel01 tag? fy? E0? b?";
	    opserr << " <a1? a2? a3? a4?>" << endln;	
	    return TCL_ERROR;
	}

	int tag;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid uniaxialMaterial Steel01 tag" << endln;
	    return TCL_ERROR;
	}

	// Read required Steel01 material parameters
	double fy, E, b;
	
	if (Tcl_GetDouble(interp, argv[3], &fy) != TCL_OK) {
	    opserr << "WARNING invalid fy\n";
	    opserr << "uniaxialMaterial Steel01: " << tag << endln;
	    return TCL_ERROR;
	}
	
	if (Tcl_GetDouble(interp, argv[4], &E) != TCL_OK) {
	    opserr << "WARNING invalid E0\n";
	    opserr << "uniaxialMaterial Steel01: " << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[5], &b) != TCL_OK) {
	    opserr << "WARNING invalid b\n";
	    opserr << "uniaxialMaterial Steel01: " << tag << endln;
	    return TCL_ERROR;
	}

	// Read optional Steel01 material parameters
	double a1, a2, a3, a4;

	if (argc > 6) {
	    if (argc < 10) {
		opserr << "WARNING insufficient number of hardening parameters\n";
		opserr << "uniaxialMaterial Steel01: " << tag << endln;
		return TCL_ERROR;
	    }
	    
	    if (Tcl_GetDouble(interp, argv[6], &a1) != TCL_OK) {
		opserr << "WARNING invalid a1\n";
		opserr << "uniaxialMaterial Steel01: " << tag << endln;
		return TCL_ERROR;
	    }
	    
	    if (Tcl_GetDouble(interp, argv[7], &a2) != TCL_OK) {
		opserr << "WARNING invalid a2\n";
		opserr << "uniaxialMaterial Steel01: " << tag << endln;
		return TCL_ERROR;
	    }
	    
	    if (Tcl_GetDouble(interp, argv[8], &a3) != TCL_OK) {
		opserr << "WARNING invalid a3\n";
		opserr << "uniaxialMaterial Steel01: " << tag << endln;
		return TCL_ERROR;
	    }

	    if (Tcl_GetDouble(interp, argv[9], &a4) != TCL_OK) {
		opserr << "WARNING invalid a4\n";
		opserr << "uniaxialMaterial Steel01: " << tag << endln;
		return TCL_ERROR;
	    }

		// Parsing was successful, allocate the material
		theMaterial = new Steel01 (tag, fy, E, b, a1, a2, a3, a4);
	}
	else
		// Parsing was successful, allocate the material
		theMaterial = new Steel01 (tag, fy, E, b);
    }

    else if (strcmp(argv[1],"Concrete01") == 0) {
	if (argc < 7) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: uniaxialMaterial Concrete01 tag? fpc? epsc0? fpcu? epscu?" << endln;
	    return TCL_ERROR;
	}

	int tag;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid uniaxialMaterial Concrete01 tag" << endln;
	    return TCL_ERROR;
	}

	// Read required Concrete01 material parameters
	double fpc, epsc0, fpcu, epscu;

	if (Tcl_GetDouble(interp, argv[3], &fpc) != TCL_OK) {
	    opserr << "WARNING invalid fpc\n";
	    opserr << "Concrete01 material: " << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[4], &epsc0) != TCL_OK) {
	    opserr << "WARNING invalid epsc0\n";
	    opserr << "Concrete01 material: " << tag << endln;
	    return TCL_ERROR;
	}
	
	if (Tcl_GetDouble(interp, argv[5], &fpcu) != TCL_OK) {
	    opserr << "WARNING invalid fpcu\n";
	    opserr << "Concrete01 material: " << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[6], &epscu) != TCL_OK) {
	    opserr << "WARNING invalid epscu\n";
	    opserr << "Concrete01 material: " << tag << endln;
	    return TCL_ERROR;
	}

	// Parsing was successful, allocate the material
	theMaterial = new Concrete01(tag, fpc, epsc0, fpcu, epscu);
    }

	else if (strcmp(argv[1],"Hysteretic") == 0) {
		if (argc != 20 && argc != 19 && argc != 16 && argc != 15) {
			opserr << "WARNING insufficient arguments\n";
			printCommand(argc,argv);
			opserr << "Want: uniaxialMaterial Hysteretic tag? mom1p? rot1p? mom2p? rot2p? <mom3p? rot3p?> "
				<< "\nmom1n? rot1n? mom2n? rot2n? <mom3n? rot3n?> pinchX? pinchY? damfc1? damfc2? <beta?>";
			return TCL_ERROR;
		}

		int tag;
		double mom1p, mom2p, mom3p;
		double rot1p, rot2p, rot3p;
		double mom1n, mom2n, mom3n;
		double rot1n, rot2n, rot3n;
		double pinchX, pinchY;
		double damfc1, damfc2;
		double beta = 0.0;

		int i = 2;

		if (Tcl_GetInt(interp, argv[i++], &tag) != TCL_OK) {
			opserr << "WARNING invalid uniaxialMaterial Hysteretic tag" << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &mom1p) != TCL_OK) {
			opserr << "WARNING invalid mom1p\n";
			opserr << "Hysteretic material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &rot1p) != TCL_OK) {
			opserr << "WARNING invalid rot1p\n";
			opserr << "Hysteretic material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &mom2p) != TCL_OK) {
			opserr << "WARNING invalid mom2p\n";
			opserr << "Hysteretic material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &rot2p) != TCL_OK) {
			opserr << "WARNING invalid rot2p\n";
			opserr << "Hysteretic material: " << tag << endln;
			return TCL_ERROR;
		}

		if (argc > 16) {
			if (Tcl_GetDouble(interp, argv[i++], &mom3p) != TCL_OK) {
				opserr << "WARNING invalid mom3p\n";
				opserr << "Hysteretic material: " << tag << endln;
				return TCL_ERROR;
			}
			
			if (Tcl_GetDouble(interp, argv[i++], &rot3p) != TCL_OK) {
				opserr << "WARNING invalid rot3p\n";
				opserr << "Hysteretic material: " << tag << endln;
				return TCL_ERROR;
			}
		}

		if (Tcl_GetDouble(interp, argv[i++], &mom1n) != TCL_OK) {
			opserr << "WARNING invalid mom1n\n";
			opserr << "Hysteretic material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &rot1n) != TCL_OK) {
			opserr << "WARNING invalid rot1n\n";
			opserr << "Hysteretic material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &mom2n) != TCL_OK) {
			opserr << "WARNING invalid mom2n\n";
			opserr << "Hysteretic material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &rot2n) != TCL_OK) {
			opserr << "WARNING invalid rot2n\n";
			opserr << "Hysteretic material: " << tag << endln;
			return TCL_ERROR;
		}

		if (argc > 16) {
			if (Tcl_GetDouble(interp, argv[i++], &mom3n) != TCL_OK) {
				opserr << "WARNING invalid mom3n\n";
				opserr << "Hysteretic material: " << tag << endln;
				return TCL_ERROR;
			}
	
			if (Tcl_GetDouble(interp, argv[i++], &rot3n) != TCL_OK) {
				opserr << "WARNING invalid rot3n\n";
				opserr << "Hysteretic material: " << tag << endln;
				return TCL_ERROR;
			}
		}

		if (Tcl_GetDouble(interp, argv[i++], &pinchX) != TCL_OK) {
			opserr << "WARNING invalid pinchX\n";
			opserr << "Hysteretic material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &pinchY) != TCL_OK) {
			opserr << "WARNING invalid pinchY\n";
			opserr << "Hysteretic material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &damfc1) != TCL_OK) {
			opserr << "WARNING invalid damfc1\n";
			opserr << "Hysteretic material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &damfc2) != TCL_OK) {
			opserr << "WARNING invalid damfc2\n";
			opserr << "Hysteretic material: " << tag << endln;
			return TCL_ERROR;
		}

		if (argc == 20 || argc == 16) {
			if (Tcl_GetDouble(interp, argv[i++], &beta) != TCL_OK) {
				opserr << "WARNING invalid beta\n";
				opserr << "Hysteretic material: " << tag << endln;
				return TCL_ERROR;
			}
		}

		// Parsing was successful, allocate the material

		if (argc > 16)		
       			theMaterial = new HystereticMaterial (tag, 
				mom1p, rot1p, mom2p, rot2p, mom3p, rot3p,
				mom1n, rot1n, mom2n, rot2n, mom3n, rot3n,
				pinchX, pinchY, damfc1, damfc2, beta);

		else
			theMaterial = new HystereticMaterial (tag,
				mom1p, rot1p, mom2p, rot2p,
				mom1n, rot1n, mom2n, rot2n,
				pinchX, pinchY, damfc1, damfc2, beta);
	}

	else if (strcmp(argv[1],"Viscous") == 0) {
		if (argc < 5)
		{
			opserr << "WARNING insufficient arguments\n";
			printCommand(argc,argv);
			opserr << "Want: uniaxialMaterial Viscous tag? C? alpha?" << endln;
			return TCL_ERROR;
		}

		int tag;
		double C, a;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
		{
			opserr << "WARNING invalid tag\n";
			opserr << "Viscous material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[3], &C) != TCL_OK)
		{
			opserr << "WARNING invalid C\n";
			opserr << "Viscous material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[4], &a) != TCL_OK)
		{
			opserr << "WARNING invalid alpha\n";
			opserr << "Viscous material: " << tag << endln;
			return TCL_ERROR;
		}

		theMaterial = new ViscousMaterial (tag, C, a);
	}

	else if (strcmp(argv[1],"PathIndependent") == 0) {
		if (argc < 4)
		{
			opserr << "WARNING insufficient arguments\n";
			printCommand(argc,argv);
			opserr << "Want: uniaxialMaterial PathIndependent tag? matTag?" << endln;
			return TCL_ERROR;
		}

		int tag, matTag;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
		{
			opserr << "WARNING invalid tag\n";
			opserr << "PathIndependent material: " << tag << endln;
			return TCL_ERROR;
		}
		
		if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK)
		{
			opserr << "WARNING invalid matTag\n";
			opserr << "PathIndependent material: " << tag << endln;
			return TCL_ERROR;
		}

		UniaxialMaterial *material = theTclBuilder->getUniaxialMaterial(matTag);
		
		if (material == 0) {
		    opserr << "WARNING material does not exist\n";
		    opserr << "material: " << matTag; 
		    opserr << "\nuniaxialMaterial PathIndependent: " << tag << endln;
		    return TCL_ERROR;
		}

		theMaterial = new PathIndependentMaterial (tag, *material);
	}

    else if (strcmp(argv[1],"MinMax") == 0) {
      if (argc < 4) {
	opserr << "WARNING insufficient arguments\n";
	printCommand(argc,argv);
	opserr << "Want: uniaxialMaterial MinMax tag? matTag?";
	opserr << " <-min min?> <-max max?>" << endln;
	return TCL_ERROR;
      }
      
      int tag, matTag;
      
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING invalid uniaxialMaterial MinMax tag" << endln;
	return TCL_ERROR;		
      }

      if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
	opserr << "WARNING invalid component tag\n";
	opserr << "uniaxialMaterial MinMax: " << tag << endln;
	return TCL_ERROR;
      }

      // Search for min and max strains
      double epsmin = NEG_INF_STRAIN;
      double epsmax = POS_INF_STRAIN;
	
      for (int j = 4; j < argc; j++) {
	if (strcmp(argv[j],"-min") == 0) {
	  if ((j+1) >= argc || Tcl_GetDouble (interp, argv[j+1], &epsmin) != TCL_OK) {
	    opserr << "WARNING invalid min\n";
	    opserr << "uniaxialMaterial MinMax: " << tag << endln;
	    return TCL_ERROR;
	  }
	  j++;
	}
	if (strcmp(argv[j],"-max") == 0) {
	  if ((j+1) >= argc || Tcl_GetDouble (interp, argv[j+1], &epsmax) != TCL_OK) {
	    opserr << "WARNING invalid max\n";
	    opserr << "uniaxialMaterial MinMax: " << tag << endln;
	    return TCL_ERROR;
	  }
	  j++;
	}
      }
	
      UniaxialMaterial *theMat = theTclBuilder->getUniaxialMaterial(matTag);
	    
      if (theMat == 0) {
	opserr << "WARNING component material does not exist\n";
	opserr << "Component material: " << matTag; 
	opserr << "\nuniaxialMaterial MinMax: " << tag << endln;
	return TCL_ERROR;
      }
	
      // Parsing was successful, allocate the material
      theMaterial = new MinMaxMaterial(tag, *theMat, epsmin, epsmax);
      
    }
	else if (strcmp(argv[1],"Cable") == 0) 
	{
		if (argc != 7) {
			opserr << "WARNING invalid number of arguments\n";
			printCommand(argc,argv);
			opserr << "Want: uniaxialMaterial Cable tag? Prestress? E? effUnitWeight? L_element?>" << endln;
			return TCL_ERROR;
		}    

		int tag;
		double Ps, E, unitWeight, L_element;
    
		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
			opserr << "WARNING invalid uniaxialMaterial Cable tag" << endln;
			return TCL_ERROR;		
		}

		if (Tcl_GetDouble(interp, argv[3], &Ps) != TCL_OK) {
			opserr << "WARNING invalid Prestress\n";
			opserr << "uniaxiaMaterial Cable: " << tag << endln;
			return TCL_ERROR;	
		}
        
		if (Tcl_GetDouble(interp, argv[4], &E) != TCL_OK) {
			opserr << "WARNING invalid E\n";
			opserr << "uniaxiaMaterial Cable: " << tag << endln;
			return TCL_ERROR;	
		}
        
		if (Tcl_GetDouble(interp, argv[5], &unitWeight) != TCL_OK) {
			opserr << "WARNING invalid unit weight\n";
			opserr << "uniaxiaMaterial Cable: " << tag << endln;
			return TCL_ERROR;	
		}
        
		if (Tcl_GetDouble(interp, argv[6], &L_element) != TCL_OK) {
			opserr << "WARNING invalid uniaxialMaterial Cable Length " << endln;
			return TCL_ERROR;		
		}
 

		theMaterial = new CableMaterial(tag, Ps, E, unitWeight, L_element);       
    }

	else {
		// Fedeas
		theMaterial = TclModelBuilder_addFedeasMaterial(clientData, interp, argc, argv, theTclBuilder);
		
		// Drain
		if (theMaterial == 0)
			theMaterial = TclModelBuilder_addDrainMaterial(clientData, interp, argc, argv, theTclBuilder);

		// SNAP
		if (theMaterial == 0)
			theMaterial = TclModelBuilder_addSnapMaterial(clientData, interp, argc, argv, theTclBuilder);
		
		// Py, Tz, Qz models
		if (theMaterial == 0)
			theMaterial = TclModelBuilder_addPyTzQzMaterial(clientData, interp, argc, argv, theTclBuilder);
	}

    if (theMaterial == 0) {
		opserr << "WARNING could not create uniaxialMaterial " << argv[1] << endln;
		return TCL_ERROR;
    }

    // Now add the material to the modelBuilder
    if (theTclBuilder->addUniaxialMaterial(*theMaterial) < 0) {
	opserr << "WARNING could not add uniaxialMaterial to the domain\n";
	opserr << *theMaterial << endln;
	delete theMaterial; // invoke the material objects destructor, otherwise mem leak
	return TCL_ERROR;
    }

    
    return TCL_OK;
}
