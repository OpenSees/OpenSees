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
                                                                        
// $Revision: 1.46 $
// $Date: 2007-07-27 19:12:42 $
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
#include <Elastic2Material.h>	// ZHY
#include <ElasticPPMaterial.h>	// fmk
#include <ParallelMaterial.h>	// fmk
#include <HardeningMaterial.h>	// MHS
#include <Steel01.h>			// MHS
#include <Steel02.h>                    // FF 
#include <Steel03.h>			// KM
//#include <Concrete.h>			// FF
#include <Concrete01.h>			// MHS
#include <Concrete01WithSITC.h>		// Won Lee
#include <ECC01.h>                      // Won Lee
#include <Concrete02.h>			// MHS
#include <Concrete04.h>
#include <Concrete05.h>
#include <Concrete06.h>			// LMS
#include <Concrete07.h>			// JDW
#include <HystereticMaterial.h>	// MHS
#include <EPPGapMaterial.h>		// Mackie
#include <ViscousMaterial.h>	// Sasani
#include <PathIndependentMaterial.h>	// MHS
#include <MinMaxMaterial.h>	// MHS
#include <FatigueMaterial.h>	// Patxi
#include <SeriesMaterial.h>		// MHS
#include <ENTMaterial.h>		// MHS
#include <CableMaterial.h>	// CC
#include <BoucWenMaterial.h>	// Terje
#include <Pinching4Material.h>   // NM
#include <ShearPanelMaterial.h>  // NM
#include <BarSlipMaterial.h>     // NM
#include <Bond_SP01.h>	// JZ

#include <Domain.h>			    // RWB bringing in Domain for PyLiq
#include <PySimple1.h>			// RWB
#include <PyLiq1.h>			    // RWB
#include <TzSimple1.h>			// RWB
#include <TzLiq1.h>			    // RWB

#include <Vector.h>
#include <string.h>


#ifdef _LIMITSTATEMATERIAL
extern UniaxialMaterial *
Tcl_AddLimitStateMaterial(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
#endif


#include <packages.h>



typedef struct uniaxialPackageCommand {
  char *funcName;
  int (*funcPtr)(ClientData clientData, Tcl_Interp *interp,  int argc, 
		 TCL_Char **argv, TclModelBuilder *); 
  struct uniaxialPackageCommand *next;
} UniaxialPackageCommand;

static UniaxialPackageCommand *theUniaxialPackageCommands = NULL;

static void printCommand(int argc, TCL_Char **argv)
{
    opserr << "Input command: ";
    for (int i=0; i<argc; i++)
	opserr << argv[i] << " ";
    opserr << endln;
} 


// external functions

int
TclCommand_ReinforcingSteel(ClientData clientData, Tcl_Interp *interp, int argc, 
			    TCL_Char **argv, TclModelBuilder *theTclBuilder);


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
				  TCL_Char **argv, TclModelBuilder *theTclBuilder, Domain *theDomain);

/*
int
TclCommand_KinematicHardening(ClientData clientData, Tcl_Interp *interp, int argc, 
			      TCL_Char **argv, TclModelBuilder *theTclBuilder);
*/

int
TclModelBuilderUniaxialMaterialCommand (ClientData clientData, Tcl_Interp *interp, int argc,
					TCL_Char **argv, TclModelBuilder *theTclBuilder, Domain *theDomain)
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

    } else if (strcmp(argv[1],"Elastic2") == 0) {
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
	theMaterial = new Elastic2Material(tag, E, eta);       
      
    } else if (strcmp(argv[1], "ReinforcingSteel") == 0) {
      return TclCommand_ReinforcingSteel(clientData,interp,argc,argv,theTclBuilder);
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
        opserr << "Want: uniaxialMaterial ElasticPPGap tag? E? fy? gap? <eta?> <damage>" << endln;
        return TCL_ERROR;
      }
      
      int tag;
      int damage = 0;
      double eta = 0.0;
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
      else if (argc > 7 && strcmp(argv[7],"damage") == 0) {
	damage = 1;
	if (Tcl_GetDouble(interp, argv[6], &eta) != TCL_OK) {
	  opserr << "WARNING invalid eta\n";
	  opserr << "uniaxialMaterial ElasticPPGap: " << tag << endln;
	  return TCL_ERROR;
	}
      }
      
      // Parsing was successful, allocate the material
      theMaterial = new EPPGapMaterial(tag, E, fy, gap, eta, damage); 
      
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

    else if (strcmp(argv[1],"Steel2") == 0 || strcmp(argv[1],"Steel02") == 0) {
      if (argc < 6) {
	opserr << "WARNING invalid number of arguments\n";
	printCommand(argc,argv);
	opserr << "Want: uniaxialMaterial Steel02 tag? fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
	return 0;
      }    
      
      double fy, E, b;
      double R0, cR1, cR2;
      double a1, a2, a3, a4, sigini;

      int tag;
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING invalid tag\n";
	printCommand(argc, argv);
	return 0;	
      }
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

	  if (argc > 13) {
	    if (Tcl_GetDouble(interp, argv[13], &sigini) != TCL_OK) {
	      opserr << "WARNING invalid a1\n";
	      printCommand(argc, argv);
	      return 0;	
	    }

	    theMaterial = new Steel02(tag, fy, E, b, R0, cR1, cR2, a1, a2, a3, a4, sigini);

	  } else
	  theMaterial = new Steel02(tag, fy, E, b, R0, cR1, cR2, a1, a2, a3, a4);
	}
	else
	  theMaterial = new Steel02(tag, fy, E, b, R0, cR1, cR2);
      }
      else
	theMaterial = new Steel02(tag, fy, E, b);
      
    }


    else if (strcmp(argv[1],"Steel03") == 0) {
      // Check that there is the minimum number of arguments
      if (argc < 9) {
        opserr << "WARNING insufficient arguments\n";
        printCommand(argc,argv);
        opserr << "Want: uniaxialMaterial Steel03 tag? fy? E0? b? r? cR1 cR2?";
        opserr << " <a1? a2? a3? a4?>" << endln;    
        return TCL_ERROR;
      }
      
      int tag;
      
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
        opserr << "WARNING invalid uniaxialMaterial Steel03 tag" << endln;
        return TCL_ERROR;
      }
      
      
      // Read required Steel01 material parameters
      double fy, E, b, r, cR1, cR2;
      
      if (Tcl_GetDouble(interp, argv[3], &fy) != TCL_OK) {
        opserr << "WARNING invalid fy\n";
        opserr << "uniaxialMaterial Steel03: " << tag << endln;
        return TCL_ERROR;
      }
      
      
      if (Tcl_GetDouble(interp, argv[4], &E) != TCL_OK) {
        opserr << "WARNING invalid E0\n";
        opserr << "uniaxialMaterial Steel03: " << tag << endln;
        return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[5], &b) != TCL_OK) {
        opserr << "WARNING invalid b\n";
        opserr << "uniaxialMaterial Steel03: " << tag << endln;
        return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[6], &r) != TCL_OK) {
	opserr << "WARNING invalid r\n";
	opserr << "uniaxialMaterial Steel03: " << tag << endln;
	return TCL_ERROR;
      }
       

      if (Tcl_GetDouble(interp, argv[7], &cR1) != TCL_OK) {
	opserr << "WARNING invalid cR1\n";
	opserr << "uniaxialMaterial Steel03: " << tag << endln;
	return TCL_ERROR;
      }
      
      
      if (Tcl_GetDouble(interp, argv[8], &cR2) != TCL_OK) {
	opserr << "WARNING invalid cR2\n";
	opserr << "uniaxialMaterial Steel03: " << tag << endln;
	return TCL_ERROR;
      }
      
      // Read optional Steel01 material parameters
      double a1, a2, a3, a4;
      if (argc > 9) {
	if (argc < 13) {
	  opserr << "WARNING insufficient number of hardening parameters\n";
	  opserr << "uniaxialMaterial Steel03: " << tag << endln;
	  return TCL_ERROR;
	}
	
	if (Tcl_GetDouble(interp, argv[9], &a1) != TCL_OK) {
	  opserr << "WARNING invalid a1\n";
	  opserr << "uniaxialMaterial Steel03: " << tag << endln;
	  return TCL_ERROR;
	}
	
	if (Tcl_GetDouble(interp, argv[10], &a2) != TCL_OK) {
	  opserr << "WARNING invalid a2\n";
	  opserr << "uniaxialMaterial Steel03: " << tag << endln;
	  return TCL_ERROR;
	}
	
	if (Tcl_GetDouble(interp, argv[11], &a3) != TCL_OK) {
	  opserr << "WARNING invalid a3\n";
	  opserr << "uniaxialMaterial Steel03: " << tag << endln;
	  return TCL_ERROR;
	}
	
	
	if (Tcl_GetDouble(interp, argv[12], &a4) != TCL_OK) {
	  opserr << "WARNING invalid a4\n";
	  opserr << "uniaxialMaterial Steel03: " << tag << endln;
	  return TCL_ERROR;
	}
	
	
        // Parsing was successful, allocate the material
	theMaterial = new Steel03 (tag, fy, E, b, r, cR1, cR2, a1, a2, a3, a4);
      }
      else
	// Parsing was successful, allocate the material
	theMaterial = new Steel03 (tag, fy, E, b, r, cR1, cR2);
      
    }

    /* ****************************************************************************
    else if (strcmp(argv[1],"Concrete") == 0 || strcmp(argv[1],"Concrete") == 0) {
      if (argc < 9) {
	opserr << "WARNING invalid number of arguments\n";
	printCommand(argc,argv);
	opserr << "Want: uniaxialMaterial Concrete tag? fc? epsc0? Ec? Kfc? ft? Ets?" << endln;
	return 0;
      }    
      
      double fc, eps0, Ec, Kfc;
      double ft, Ets;
      int tag;
      
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING invalid uniaxialMaterial Concrete02 tag" << endln;
	return TCL_ERROR;
      }      

      if (Tcl_GetDouble(interp, argv[3], &fc) != TCL_OK) {
	opserr << "WARNING invalid fc\n";
	printCommand(argc, argv);
	return 0;	
      }

      if (Tcl_GetDouble(interp, argv[4], &eps0) != TCL_OK) {
	opserr << "WARNING invalid eps0\n";
	printCommand(argc, argv);
	return 0;	
      }

      if (Tcl_GetDouble(interp, argv[5], &Ec) != TCL_OK) {
	opserr << "WARNING invalid Ec\n";
	printCommand(argc, argv);
	return 0;	
      }

      if (Tcl_GetDouble(interp, argv[6], &Kfc) != TCL_OK) {
	opserr << "WARNING invalid Kfc\n";
	printCommand(argc, argv);
	return 0;	
      }

      if (Tcl_GetDouble(interp, argv[7], &ft) != TCL_OK) {
	opserr << "WARNING invalid ft\n";
	printCommand(argc, argv);
	return 0;	
      }

      if (Tcl_GetDouble(interp, argv[8], &Ets) != TCL_OK) {
	opserr << "WARNING invalid Ets\n";
	printCommand(argc, argv);
	return 0;	
      }
      
      theMaterial = new Concrete(tag, fc, eps0, Ec, Kfc, ft, Ets);
    }
    *************************************************************** */
  
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

    else if (strcmp(argv[1],"Concrete2") == 0 || strcmp(argv[1],"Concrete02") == 0) {
      if (argc < 10) {
	opserr << "WARNING invalid number of arguments\n";	
	printCommand(argc,argv);
	opserr << "Want: uniaxialMaterial Concrete02 tag? fpc? epsc0? fpcu? epscu? rat? ft? Ets?" << endln;
	return 0;
      }
    
      double fpc, epsc0, fpcu, epscu;
      double rat, ft, Ets;
      int tag;

 
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING invalid uniaxialMaterial Concrete02 tag" << endln;
	return TCL_ERROR;
      }      
      
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
      
      theMaterial = new Concrete02(tag, fpc, epsc0, fpcu, epscu, rat, ft, Ets);
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
    
    else if (strcmp(argv[1],"Concrete04") == 0) {
      //        opserr << argc << endln;
      if (argc != 10 && argc != 9 && argc != 7) {
	opserr << "WARNING insufficient arguments\n";
	printCommand(argc,argv);
	opserr << "Want: uniaxialMaterial Concrete04 tag? fpc? epsc0? epscu? Ec0? <ft? etu? <beta?> >" << endln;
	return TCL_ERROR;
      }
      
      int tag;
      
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING invalid uniaxialMaterial Concrete04 tag"  
	       << endln;
	return TCL_ERROR;
      }
      
      // Read required Concrete04 material parameters
      double fpc, epsc0, ft, epscu, Ec0, etu, beta;
      
      if (Tcl_GetDouble(interp, argv[3], &fpc) != TCL_OK) {
	opserr << "WARNING invalid fpc\n";
	opserr << "Concrete04 material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[4], &epsc0) != TCL_OK) {
	opserr << "WARNING invalid epsc0\n";
	opserr << "Concrete04 material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[5], &epscu) != TCL_OK) {
	opserr << "WARNING invalid epscu\n";
	opserr << "Concrete04 material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[6], &Ec0) != TCL_OK) {
	opserr << "WARNING invalid Ec0\n";
	opserr << "Concrete04 material: " << tag << endln;
	return TCL_ERROR;
      }
      if (argc == 9 || argc == 10) {
	if (Tcl_GetDouble(interp, argv[7], &ft) != TCL_OK) {
	  opserr << "WARNING invalid ft\n";
	  opserr << "Concrete04 material: " << tag << endln;
	  return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[8], &etu) != TCL_OK) {
	  opserr << "WARNING invalid etu\n";
	  opserr << "Concrete04 material: " << tag << endln;
	  return TCL_ERROR;
	}
      }
      if (argc == 10) {
	if (Tcl_GetDouble(interp, argv[9], &beta) != TCL_OK) {
	  opserr << "WARNING invalid beta\n";
	  opserr << "Concrete04 material: " << tag << endln;
	  return TCL_ERROR;
	}
      }
      
      
      // Parsing was successful, allocate the material
      if (argc == 10) {
	theMaterial = new Concrete04(tag, fpc, epsc0, epscu, Ec0,  
				     ft, etu, beta);
      }
      else if (argc == 9) {
	theMaterial = new Concrete04(tag, fpc, epsc0, epscu, Ec0,  
				     ft, etu);
      }
      else if (argc == 7) {
	theMaterial = new Concrete04(tag, fpc, epsc0, epscu, Ec0);
      }
    }



    else if (strcmp(argv[1],"Concrete06") == 0) {
	if (argc < 12) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: uniaxialMaterial Concrete06 tag? fc? eo? r? k? alphaC? fcr? ecr? b? alphaT?" << endln;
	    return TCL_ERROR;
	}

	int tag;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid uniaxialMaterial Concrete06 tag" << endln;
	    return TCL_ERROR;
	}

	// Read required Concrete01 material parameters
	double fc, eo, r, k, fcr, ecr, b, alphaC, alphaT;

	if (Tcl_GetDouble(interp, argv[3], &fc) != TCL_OK) {
	    opserr << "WARNING invalid fc\n";
	    opserr << "Concrete06 material: " << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[4], &eo) != TCL_OK) {
	    opserr << "WARNING invalid eo\n";
	    opserr << "Concrete06 material: " << tag << endln;
	    return TCL_ERROR;
	}
	
	if (Tcl_GetDouble(interp, argv[5], &r) != TCL_OK) {
	    opserr << "WARNING invalid r\n";
	    opserr << "Concrete06 material: " << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[6], &k) != TCL_OK) {
	    opserr << "WARNING invalid k\n";
	    opserr << "Concrete06 material: " << tag << endln;
	    return TCL_ERROR;
	}
	
	if (Tcl_GetDouble(interp, argv[7], &alphaC) != TCL_OK) {
	  opserr << "WARNING invalid alphaC\n";
	  opserr << "Concrete06 material: " << tag << endln;
	  return TCL_ERROR;
	}
	
	if (Tcl_GetDouble(interp, argv[8], &fcr) != TCL_OK) {
	  opserr << "WARNING invalid fcr\n";
	  opserr << "Concrete06 material: " << tag << endln;
	  return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[9], &ecr) != TCL_OK) {
	  opserr << "WARNING invalid ecr\n";
	  opserr << "Concrete06 material: " << tag << endln;
	  return TCL_ERROR;
	}
	
	if (Tcl_GetDouble(interp, argv[10], &b) != TCL_OK) {
	  opserr << "WARNING invalid b\n";
	  opserr << "Concrete06 material: " << tag << endln;
	  return TCL_ERROR;
	}
	
	if (Tcl_GetDouble(interp, argv[11], &alphaT) != TCL_OK) {
	  opserr << "WARNING invalid alphaT\n";
	  opserr << "Concrete06 material: " << tag << endln;
	  return TCL_ERROR;
	}
	
	// Parsing was successful, allocate the material
	theMaterial = new Concrete06(tag, fc, eo, r, k, alphaC, fcr, ecr, b, alphaT);
    }

    
    
    else if (strcmp(argv[1], "Concrete07") == 0) {
      // Check to see if there are enough arquements
      if (argc < 11) {
	opserr << "WARNING: Insufficient arguements\n";
	printCommand(argc, argv);
	opserr << "Want: uniaxialMaterial Concrete07 tag? fpc? epsc0? Ec? fpt? epst0? xcrp? xcrn? r?\n";
	return TCL_ERROR;
      }
      
      int tag;

      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING: Invalid uniaxial Concrete07 tag\n";
	return TCL_ERROR;
      }
      
      // Read in the faluves required for the model
      double fpc, epsc0, Ec, fpt, epst0, xcrp, xcrn, r;
      
      if (Tcl_GetDouble(interp, argv[3], &fpc) != TCL_OK) {
	opserr << "WARNING: Invalid peak compression stress\n";
	opserr << "uniaxialMaterial Concrete07: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[4], &epsc0) != TCL_OK) {
	opserr << "WARNING: Invalid peak compression strain\n";
	opserr << "uniaxialMaterial Concrete07: " << tag <<endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[5], &Ec) != TCL_OK) {
	opserr << "WARNING: Invalid Young's Modulus\n";
	opserr << "uniaxialMaterial Concrete07: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[6], &fpt) != TCL_OK) {
	opserr << "WARNING: Invalid peak tension stress\n";
	opserr << "uniaxialMaterial Concrete07: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[7], &epst0) != TCL_OK) {
	opserr << "WARNING: Invalid peak tension strain\n";
	opserr << "uniaxialMaterial Concrete07: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[8], &xcrp) != TCL_OK) {
	opserr << "WARNING: Invalid critical nondimensional strain in tension\n";
	opserr << "uniaxialMaterial Concrete07: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[9], &xcrn) != TCL_OK) {
	opserr << "WARNING: Invalid critical nondimensional strain in compression\n";
	opserr << "uniaxialMaterial Concrete07: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[10], &r) != TCL_OK) {
	opserr << "WARNING: Invalid value for r\n";
	opserr << "uniaxialMaterial Concrete07: " << tag << endln;
      }
      
      //		opserr << "fpc: " << fpc << endln << "epsc0: " << epsc0 << endln << "Ec: " << Ec << endln;
      //		opserr << "fpt: " << fpt << endln << "epst0: " << epst0 << endln << "xcrp: " << xcrp << endln;
      //		opserr << "xcrn: " << xcrn << endln << "r: " << r << endln;
      
      // Parsing was successful, allocate the material
      
      theMaterial = new Concrete07(tag, fpc, epsc0, Ec, fpt, epst0, xcrp, xcrn, r);
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

    else if (strcmp(argv[1],"Fatigue") == 0) {
      if (argc < 4) {
	opserr << "WARNING insufficient arguments\n";
	printCommand(argc,argv);
	opserr << "Want: uniaxialMaterial Fatigue tag? matTag?";
	opserr << " <-D_max dmax?> <-e0 e0?> <-m m?>" << endln;
	opserr << " <-min min?> <-max max?>" << endln;
	return TCL_ERROR;
      }
      
      int tag, matTag;
      
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING invalid uniaxialMaterial Fatigue tag" << endln;
	return TCL_ERROR;		
      }
      
      if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
	opserr << "WARNING invalid component tag\n";
	opserr << "uniaxialMaterial Fatigue: " << tag << endln;
	return TCL_ERROR;
      }
      
      double Dmax      =  1.0;
      double E0        =  0.191;
      double m         = -0.458;
      double epsmin    = NEG_INF_STRAIN;
      double epsmax    = POS_INF_STRAIN;
      
      for (int j = 4; j < argc; j++) {
	if (strcmp(argv[j],"-Dmax") == 0) {
	  if ((j+1 >= argc) || 
	      (Tcl_GetDouble (interp, argv[j+1], &Dmax) != TCL_OK)) {
	    opserr << "WARNING invalid -Dmax";
	    opserr << "uniaxialMaterial Fatigue: " << tag << endln;
	    return TCL_ERROR;
	  }
	} else if (strcmp(argv[j],"-E0") == 0) {
	  if ((j+1 >= argc) || 
	      (Tcl_GetDouble (interp, argv[j+1], &E0) != TCL_OK)) {
	    opserr << "WARNING invalid -E0";	 
	    opserr << "uniaxialMaterial Fatigue: " << tag << endln;
	    return TCL_ERROR;
	  }
	} else if (strcmp(argv[j],"-m") == 0) {
	  if ((j+1 >= argc) || 
	      (Tcl_GetDouble (interp, argv[j+1], &m) != TCL_OK)) {
	    opserr << "WARNING invalid -m"; 
	    opserr << "uniaxialMaterial Fatigue: " << tag << endln;
	    return TCL_ERROR;
	  }
	} else if (strcmp(argv[j],"-min") == 0) {
	  if ((j+1 >= argc) || 
	      (Tcl_GetDouble (interp, argv[j+1], &epsmin) != TCL_OK)) {
	    opserr << "WARNING invalid -min ";
	    opserr << "uniaxialMaterial Fatigue: " << tag << endln;
	    return TCL_ERROR;
	  }
	} else if (strcmp(argv[j],"-max") == 0) {
	  if ((j+1 >= argc) || 
	      (Tcl_GetDouble (interp, argv[j+1], &epsmax) != TCL_OK)) {
	    opserr << "WARNING invalid -max";
	    opserr << "uniaxialMaterial Fatigue: " << tag << endln;
	    return TCL_ERROR;
	  }
	}
	j++;
      }
      
      
      UniaxialMaterial *theMat = theTclBuilder->getUniaxialMaterial(matTag);
      
      if (theMat == 0) {
	opserr << "WARNING component material does not exist\n";
	opserr << "Component material: " << matTag; 
	opserr << "\nuniaxialMaterial Fatigue: " << tag << endln;
	return TCL_ERROR;
      }
      
      // Parsing was successful, allocate the material
      theMaterial = new FatigueMaterial(tag, *theMat, Dmax, E0, 
					m, epsmin, epsmax);
      
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

	else if (strcmp(argv[1],"Pinching4") == 0) {
		if (argc != 42 && argc != 31 ) {
			opserr << "WARNING insufficient arguments\n";
			printCommand(argc,argv);
			opserr << "Want: uniaxialMaterial Pinching4 tag? stress1p? strain1p? stress2p? strain2p? stress3p? strain3p? stress4p? strain4p? "
				<< "\n<stress1n? strain1n? stress2n? strain2n? stress3n? strain3n? stress4n? strain4n?> rDispP? rForceP? uForceP? "
				<< "\n<rDispN? rForceN? uForceN?> gammaK1? gammaK2? gammaK3? gammaK4? gammaKLimit? gammaD1? gammaD2? gammaD3? gammaD4? "
				<< "\ngammaDLimit? gammaF1? gammaF2? gammaF3? gammaF4? gammaFLimit? gammaE? CycleOrEnergyDamage? ";
			return TCL_ERROR;
		}

		int tag, tDmg;
		double stress1p, stress2p, stress3p, stress4p;
		double strain1p, strain2p, strain3p, strain4p;
		double stress1n, stress2n, stress3n, stress4n;
		double strain1n, strain2n, strain3n, strain4n;
		double rDispP, rForceP, uForceP, rDispN, rForceN, uForceN;
		double gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit;
		double gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit;
		double gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit;
		double gammaE;

		int i = 2;

		if (Tcl_GetInt(interp, argv[i++], &tag) != TCL_OK) {
			opserr << "WARNING invalid uniaxialMaterial Pinching4 tag" << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &stress1p) != TCL_OK) {
			opserr << "WARNING invalid stress1p\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &strain1p) != TCL_OK) {
			opserr << "WARNING invalid strain1p\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &stress2p) != TCL_OK) {
			opserr << "WARNING invalid stress2p\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &strain2p) != TCL_OK) {
			opserr << "WARNING invalid strain2p\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &stress3p) != TCL_OK) {
			opserr << "WARNING invalid stress3p\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &strain3p) != TCL_OK) {
			opserr << "WARNING invalid strain3p\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &stress4p) != TCL_OK) {
			opserr << "WARNING invalid stress4p\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &strain4p) != TCL_OK) {
			opserr << "WARNING invalid strain4p\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

		if (argc == 42) {
			if (Tcl_GetDouble(interp, argv[i++], &stress1n) != TCL_OK) {
				opserr << "WARNING invalid stress1n\n";
				opserr << "Pinching4 material: " << tag << endln;
				return TCL_ERROR;
			}

			if (Tcl_GetDouble(interp, argv[i++], &strain1n) != TCL_OK) {
				opserr << "WARNING invalid strain1n\n";
				opserr << "Pinching4 material: " << tag << endln;
				return TCL_ERROR;
			}

			if (Tcl_GetDouble(interp, argv[i++], &stress2n) != TCL_OK) {
				opserr << "WARNING invalid stress2n\n";
				opserr << "Pinching4 material: " << tag << endln;
				return TCL_ERROR;
			}

			if (Tcl_GetDouble(interp, argv[i++], &strain2n) != TCL_OK) {
				opserr << "WARNING invalid strain2n\n";
				opserr << "Pinching4 material: " << tag << endln;
				return TCL_ERROR;
			}

			if (Tcl_GetDouble(interp, argv[i++], &stress3n) != TCL_OK) {
				opserr << "WARNING invalid stress3n\n";
				opserr << "Pinching4 material: " << tag << endln;
				return TCL_ERROR;
			}

			if (Tcl_GetDouble(interp, argv[i++], &strain3n) != TCL_OK) {
				opserr << "WARNING invalid strain3n\n";
				opserr << "Pinching4 material: " << tag << endln;
				return TCL_ERROR;
			}

			if (Tcl_GetDouble(interp, argv[i++], &stress4n) != TCL_OK) {
				opserr << "WARNING invalid stress4n\n";
				opserr << "Pinching4 material: " << tag << endln;
				return TCL_ERROR;
			}

			if (Tcl_GetDouble(interp, argv[i++], &strain4n) != TCL_OK) {
				opserr << "WARNING invalid strain4n\n";
				opserr << "Pinching4 material: " << tag << endln;
				return TCL_ERROR;
			}
		
		}
		

		if (Tcl_GetDouble(interp, argv[i++], &rDispP) != TCL_OK) {
			opserr << "WARNING invalid rDispP\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &rForceP) != TCL_OK) {
			opserr << "WARNING invalid rForceP\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &uForceP) != TCL_OK) {
			opserr << "WARNING invalid uForceP\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

		if (argc == 42) {
			if (Tcl_GetDouble(interp, argv[i++], &rDispN) != TCL_OK) {
				opserr << "WARNING invalid rDispN\n";
				opserr << "Pinching4 material: " << tag << endln;
				return TCL_ERROR;
			}

			if (Tcl_GetDouble(interp, argv[i++], &rForceN) != TCL_OK) {
				opserr << "WARNING invalid rForceN\n";
				opserr << "Pinching4 material: " << tag << endln;
				return TCL_ERROR;
			}

			if (Tcl_GetDouble(interp, argv[i++], &uForceN) != TCL_OK) {
				opserr << "WARNING invalid uForceN\n";
				opserr << "Pinching4 material: " << tag << endln;
				return TCL_ERROR;
			}
		}

		if (Tcl_GetDouble(interp, argv[i++], &gammaK1) != TCL_OK) {
			opserr << "WARNING invalid gammaK1\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[i++], &gammaK2) != TCL_OK) {
			opserr << "WARNING invalid gammaK2\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[i++], &gammaK3) != TCL_OK) {
			opserr << "WARNING invalid gammaK3\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[i++], &gammaK4) != TCL_OK) {
			opserr << "WARNING invalid gammaK4\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[i++], &gammaKLimit) != TCL_OK) {
			opserr << "WARNING invalid gammaKLimit\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[i++], &gammaD1) != TCL_OK) {
			opserr << "WARNING invalid gammaD1\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}										   
		if (Tcl_GetDouble(interp, argv[i++], &gammaD2) != TCL_OK) {
			opserr << "WARNING invalid gammaD2\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[i++], &gammaD3) != TCL_OK) {
			opserr << "WARNING invalid gammaD3\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[i++], &gammaD4) != TCL_OK) {
			opserr << "WARNING invalid gammaD4\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[i++], &gammaDLimit) != TCL_OK) {
			opserr << "WARNING invalid gammaDLimit\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[i++], &gammaF1) != TCL_OK) {
			opserr << "WARNING invalid gammaF1\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[i++], &gammaF2) != TCL_OK) {
			opserr << "WARNING invalid gammaF2\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[i++], &gammaF3) != TCL_OK) {
			opserr << "WARNING invalid gammaF3\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[i++], &gammaF4) != TCL_OK) {
			opserr << "WARNING invalid gammaF4\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[i++], &gammaFLimit) != TCL_OK) {
			opserr << "WARNING invalid gammaFLimit\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &gammaE) != TCL_OK) {
			opserr << "WARNING invalid gammaE\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

		int y; 
		y = i;

		if ((strcmp(argv[y],"cycle") == 0) || (strcmp(argv[y],"Cycle") == 0) || (strcmp(argv[y],"DamageCycle") == 0) || (strcmp(argv[y],"damageCycle") == 0))
		{ tDmg = 1; }
		else if ((strcmp(argv[y],"energy") == 0) || (strcmp(argv[y],"Energy") == 0) || (strcmp(argv[y],"DamageEnergy") == 0) || (strcmp(argv[y],"damageEnergy") == 0))
		{ tDmg = 0; }
		else
		{
			opserr << "WARNING invalid type of damage calculation specified\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

	// allocate the pinching material
		if (argc == 42) {
		theMaterial = new Pinching4Material (tag,
			stress1p, strain1p, stress2p, strain2p, stress3p, strain3p, stress4p, strain4p,
			stress1n, strain1n, stress2n, strain2n, stress3n, strain3n, stress4n, strain4n,
			rDispP, rForceP, uForceP, rDispN, rForceN, uForceN, 
			gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit,
			gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit,
			gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit, gammaE, tDmg);
		}
		if (argc == 31) {
		theMaterial = new Pinching4Material (tag,
			stress1p, strain1p, stress2p, strain2p, stress3p, strain3p, stress4p, strain4p,
			rDispP, rForceP, uForceP,  
			gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit,
			gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit,
			gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit, gammaE, tDmg);		
		}
   }
   
  else if (strcmp(argv[1],"BarSlip") == 0)
   {
     if (argc != 17 && argc != 15)
       {
	 opserr << "WARNING insufficient arguments\n";
	 printCommand(argc,argv);
	 opserr << "Want: uniaxialMaterial BarSlip tag? fc? fy? Es? fu? Eh? db? ld? nb? width? depth? bsflag? type? <damage? unit?>"  << endln;
	 return TCL_ERROR;
       }
     
     int tag, nb, bsf, typ, dmg, unt;
     double fc, fy, Es, fu, Eh, ld, width, depth, db;
     
     int argStart = 2;
     
     if (Tcl_GetInt(interp, argv[argStart++], &tag) != TCL_OK)
       {
	 opserr << "WARNING invalid tag\n";
	 opserr << "BarSlip: " << tag << endln;
	 return TCL_ERROR;
       }
     if (Tcl_GetDouble(interp, argv[argStart++], &fc) != TCL_OK)
       {
	 opserr << "WARNING invalid fc\n";
	 opserr << "BarSlip: " << tag << endln;
	 return TCL_ERROR;
       }
     if (Tcl_GetDouble(interp, argv[argStart++], &fy) != TCL_OK)
       {
	 opserr << "WARNING invalid fy\n";
	 opserr << "BarSlip: " << tag << endln;
	 return TCL_ERROR;
       }
     if (Tcl_GetDouble(interp, argv[argStart++], &Es) != TCL_OK)
       {
	 opserr << "WARNING invalid Es\n";
	 opserr << "BarSlip: " << tag << endln;
	 return TCL_ERROR;
       }
     if (Tcl_GetDouble(interp, argv[argStart++], &fu) != TCL_OK)
       {
	 opserr << "WARNING invalid fu\n";
	 opserr << "BarSlip: " << tag << endln;
	 return TCL_ERROR;
       }
     if (Tcl_GetDouble(interp, argv[argStart++], &Eh) != TCL_OK)
       {
	 opserr << "WARNING invalid Eh\n";
	 opserr << "BarSlip: " << tag << endln;
	 return TCL_ERROR;
       }
     if (Tcl_GetDouble(interp, argv[argStart++], &db) != TCL_OK)
       {
	 opserr << "WARNING invalid db\n";
	 opserr << "BarSlip: " << tag << endln;
	 return TCL_ERROR;
       }
     if (Tcl_GetDouble(interp, argv[argStart++], &ld) != TCL_OK)
       {
	 opserr << "WARNING invalid ld\n";
	 opserr << "BarSlip: " << tag << endln;
	 return TCL_ERROR;
       }
     if (Tcl_GetInt(interp, argv[argStart++], &nb) != TCL_OK)
       {
	 opserr << "WARNING invalid nbars\n";
	 opserr << "BarSlip: " << tag << endln;
	 return TCL_ERROR;
       }
     if (Tcl_GetDouble(interp, argv[argStart++], &width) != TCL_OK)
       {
	 opserr << "WARNING invalid width\n";
	 opserr << "BarSlip: " << tag << endln;
	 return TCL_ERROR;
       }
     if (Tcl_GetDouble(interp, argv[argStart++], &depth) != TCL_OK)
       {
	 opserr << "WARNING invalid depth\n";
	 opserr << "BarSlip: " << tag << endln;
	 return TCL_ERROR;
       }
     
     int y;
     y = argStart;
     
     
     if ((strcmp(argv[y],"strong") == 0) || (strcmp(argv[y],"Strong") == 0) || (strcmp(argv[y],"weak") == 0) || (strcmp(argv[y],"Weak") == 0))
       {
	 if ((strcmp(argv[y],"strong") == 0) || (strcmp(argv[y],"Strong") == 0))
	   {
	     bsf = 0;
	   }
	 
	 if ((strcmp(argv[y],"weak") == 0) || (strcmp(argv[y],"Weak") == 0))
	   {
	     bsf = 1;
	   }
       }
     else
       {
	 opserr << "WARNING invalid bond strength specified\n";
	 opserr << "BarSlip: " << tag << endln;
	 return TCL_ERROR;
       }
     y ++;
     
     if ((strcmp(argv[y],"beamtop") == 0) || (strcmp(argv[y],"beamTop") == 0) || 
	 (strcmp(argv[y],"beambot") == 0) || (strcmp(argv[y],"beamBot") == 0) || (strcmp(argv[y],"beambottom") == 0) || (strcmp(argv[y],"beamBottom") == 0) ||
	 (strcmp(argv[y],"beam") == 0) || (strcmp(argv[y],"Beam") == 0) || (strcmp(argv[y],"Column") == 0) || (strcmp(argv[y],"column") == 0))
       {
	 if ((strcmp(argv[y],"beamtop") == 0) || (strcmp(argv[y],"beamTop") == 0) || (strcmp(argv[y],"beam") == 0) || (strcmp(argv[y],"Beam") == 0))
	   {
	     typ = 0;
	   }
	 
	 if ((strcmp(argv[y],"beambot") == 0) || (strcmp(argv[y],"beamBot") == 0) || (strcmp(argv[y],"beambottom") == 0) || (strcmp(argv[y],"beamBottom") == 0))
	   {
	     typ = 1;
	   }
	 
	 if ((strcmp(argv[y],"column") == 0) || (strcmp(argv[y],"Column") == 0))
	   {
	     typ = 2;
	   }
       }
     else
       {
	 opserr << "WARNING invalid location of bar specified\n";
	 opserr << "BarSlip: " << tag << endln;
	 return TCL_ERROR;
       }
     if (argc == 17) {
       y ++;
       
       if ((strcmp(argv[y],"damage1") == 0) || (strcmp(argv[y],"Damage1") == 0) || (strcmp(argv[y],"damage2") == 0) || (strcmp(argv[y],"Damage2") == 0) || 
	   (strcmp(argv[y],"nodamage") == 0) || (strcmp(argv[y],"Nodamage") == 0) || (strcmp(argv[y],"NoDamage") == 0) || (strcmp(argv[y],"noDamage") == 0))
	 {
	   if ((strcmp(argv[y],"damage1") == 0) || (strcmp(argv[y],"Damage1") == 0))
	     {
	       dmg = 1;
	     }
	   else if ((strcmp(argv[y],"damage2") == 0) || (strcmp(argv[y],"Damage2") == 0))
	     {
	       dmg = 2;
	     }
	   else if ((strcmp(argv[y],"nodamage") == 0) || (strcmp(argv[y],"Nodamage") == 0) || (strcmp(argv[y],"NoDamage") == 0) || (strcmp(argv[y],"noDamage") == 0))
	     {
	       dmg = 0;
	     }
	   
	 }
       else
	 {
	   opserr << "WARNING invalid damage specified\n";
	   opserr << "BarSlip: " << tag << endln;
	   return TCL_ERROR;
	 }
       
       y ++;
       
       if ((strcmp(argv[y],"mpa") == 0) || (strcmp(argv[y],"MPa") == 0) || (strcmp(argv[y],"mPa") == 0) || (strcmp(argv[y],"Mpa") == 0) ||
	   (strcmp(argv[y],"psi") == 0) || (strcmp(argv[y],"Psi") == 0) || (strcmp(argv[y],"PSI") == 0) || (strcmp(argv[y],"Pa") == 0) ||
	   (strcmp(argv[y],"pa") == 0) ||  (strcmp(argv[y],"psf") == 0) || (strcmp(argv[y],"Psf") == 0) || (strcmp(argv[y],"PSF") == 0) ||
	   (strcmp(argv[y],"ksi") == 0) || (strcmp(argv[y],"Ksi") == 0) || (strcmp(argv[y],"KSI") == 0) || (strcmp(argv[y],"ksf") == 0) ||
	   (strcmp(argv[y],"Ksf") == 0) || (strcmp(argv[y],"KSF") == 0))
	 {
	   if ((strcmp(argv[y],"mpa") == 0) || (strcmp(argv[y],"MPa") == 0) || (strcmp(argv[y],"mPa") == 0) || (strcmp(argv[y],"Mpa") == 0))
	     {
	       unt = 1;
	     }
	   else if ((strcmp(argv[y],"psi") == 0) || (strcmp(argv[y],"Psi") == 0) || (strcmp(argv[y],"PSI") == 0))
	     {
	       unt = 2;
	     }
	   else if ((strcmp(argv[y],"Pa") == 0) || (strcmp(argv[y],"pa") == 0))
	     {
	       unt = 3;
	     }
	   else if ((strcmp(argv[y],"psf") == 0) || (strcmp(argv[y],"Psf") == 0) || (strcmp(argv[y],"PSF") == 0))
	     {
	       unt = 4;
	     }
	   else if ((strcmp(argv[y],"ksi") == 0) || (strcmp(argv[y],"Ksi") == 0) || (strcmp(argv[y],"KSI") == 0))
	     {
	       unt = 5;
	     }
	   else if ((strcmp(argv[y],"ksf") == 0) || (strcmp(argv[y],"Ksf") == 0) || (strcmp(argv[y],"KSF") == 0))
	     {
	       unt = 6;
	     }
	 }
       else
	 {
	   opserr << "WARNING invalid unit specified\n";
	   opserr << "BarSlip: " << tag << endln;
	   return TCL_ERROR;
	 }
     }
     
     // allocate the material
     if (argc == 15 ) {
       theMaterial = new BarSlipMaterial (tag, fc, fy, Es, fu, Eh, db, ld, nb, width, depth, bsf, typ);
     }
     
     if (argc == 17) {
       theMaterial = new BarSlipMaterial (tag, fc, fy, Es, fu, Eh, db, ld, nb, width, depth, bsf, typ, dmg, unt);
     }
     
   } 
    
    
    
  else if (strcmp(argv[1],"ShearPanel") == 0) {
    if (argc != 42 && argc != 31 ) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc,argv);
      opserr << "Want: uniaxialMaterial ShearPanel tag? stress1p? strain1p? stress2p? strain2p? stress3p? strain3p? stress4p? strain4p? "
	     << "\n<stress1n? strain1n? stress2n? strain2n? stress3n? strain3n? stress4n? strain4n?> rDispP? rForceP? uForceP? "
	     << "\n<rDispN? rForceN? uForceN?> gammaK1? gammaK2? gammaK3? gammaK4? gammaKLimit? gammaD1? gammaD2? gammaD3? gammaD4? "
	     << "\ngammaDLimit? gammaF1? gammaF2? gammaF3? gammaF4? gammaFLimit? gammaE? YieldStress? ";
      return TCL_ERROR;
    }
    
    int tag;
    double stress1p, stress2p, stress3p, stress4p;
    double strain1p, strain2p, strain3p, strain4p;
    double stress1n, stress2n, stress3n, stress4n;
    double strain1n, strain2n, strain3n, strain4n;
    double rDispP, rForceP, uForceP, rDispN, rForceN, uForceN;
    double gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit;
    double gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit;
    double gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit;
    double gammaE, yStr;
    
    int i = 2;
    
    if (Tcl_GetInt(interp, argv[i++], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial ShearPanel tag" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &stress1p) != TCL_OK) {
      opserr << "WARNING invalid stress1p\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &strain1p) != TCL_OK) {
      opserr << "WARNING invalid strain1p\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &stress2p) != TCL_OK) {
      opserr << "WARNING invalid stress2p\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &strain2p) != TCL_OK) {
      opserr << "WARNING invalid strain2p\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &stress3p) != TCL_OK) {
      opserr << "WARNING invalid stress3p\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &strain3p) != TCL_OK) {
      opserr << "WARNING invalid strain3p\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &stress4p) != TCL_OK) {
      opserr << "WARNING invalid stress4p\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &strain4p) != TCL_OK) {
      opserr << "WARNING invalid strain4p\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (argc == 42) {
      if (Tcl_GetDouble(interp, argv[i++], &stress1n) != TCL_OK) {
	opserr << "WARNING invalid stress1n\n";
	opserr << "ShearPanel material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[i++], &strain1n) != TCL_OK) {
	opserr << "WARNING invalid strain1n\n";
	opserr << "ShearPanel material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[i++], &stress2n) != TCL_OK) {
	opserr << "WARNING invalid stress2n\n";
	opserr << "ShearPanel material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[i++], &strain2n) != TCL_OK) {
	opserr << "WARNING invalid strain2n\n";
	opserr << "ShearPanel material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[i++], &stress3n) != TCL_OK) {
	opserr << "WARNING invalid stress3n\n";
	opserr << "ShearPanel material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[i++], &strain3n) != TCL_OK) {
	opserr << "WARNING invalid strain3n\n";
	opserr << "ShearPanel material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[i++], &stress4n) != TCL_OK) {
	opserr << "WARNING invalid stress4n\n";
	opserr << "ShearPanel material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[i++], &strain4n) != TCL_OK) {
	opserr << "WARNING invalid strain4n\n";
	opserr << "ShearPanel material: " << tag << endln;
	return TCL_ERROR;
      }
      
    }
    
    
    if (Tcl_GetDouble(interp, argv[i++], &rDispP) != TCL_OK) {
      opserr << "WARNING invalid rDispP\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &rForceP) != TCL_OK) {
      opserr << "WARNING invalid rForceP\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &uForceP) != TCL_OK) {
      opserr << "WARNING invalid uForceP\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (argc == 42) {
      if (Tcl_GetDouble(interp, argv[i++], &rDispN) != TCL_OK) {
	opserr << "WARNING invalid rDispN\n";
	opserr << "ShearPanel material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[i++], &rForceN) != TCL_OK) {
	opserr << "WARNING invalid rForceN\n";
	opserr << "ShearPanel material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[i++], &uForceN) != TCL_OK) {
	opserr << "WARNING invalid uForceN\n";
	opserr << "ShearPanel material: " << tag << endln;
	return TCL_ERROR;
      }
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &gammaK1) != TCL_OK) {
      opserr << "WARNING invalid gammaK1\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaK2) != TCL_OK) {
      opserr << "WARNING invalid gammaK2\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaK3) != TCL_OK) {
      opserr << "WARNING invalid gammaK3\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaK4) != TCL_OK) {
      opserr << "WARNING invalid gammaK4\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaKLimit) != TCL_OK) {
      opserr << "WARNING invalid gammaKLimit\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaD1) != TCL_OK) {
      opserr << "WARNING invalid gammaD1\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }										   
    if (Tcl_GetDouble(interp, argv[i++], &gammaD2) != TCL_OK) {
      opserr << "WARNING invalid gammaD2\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaD3) != TCL_OK) {
      opserr << "WARNING invalid gammaD3\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaD4) != TCL_OK) {
      opserr << "WARNING invalid gammaD4\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaDLimit) != TCL_OK) {
      opserr << "WARNING invalid gammaDLimit\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaF1) != TCL_OK) {
      opserr << "WARNING invalid gammaF1\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaF2) != TCL_OK) {
      opserr << "WARNING invalid gammaF2\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaF3) != TCL_OK) {
      opserr << "WARNING invalid gammaF3\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaF4) != TCL_OK) {
      opserr << "WARNING invalid gammaF4\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaFLimit) != TCL_OK) {
      opserr << "WARNING invalid gammaFLimit\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &gammaE) != TCL_OK) {
      opserr << "WARNING invalid gammaE\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &yStr) != TCL_OK) {
      opserr << "WARNING invalid yield stress\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    // allocate the pinching material
    if (argc == 42) {
      theMaterial = new ShearPanelMaterial (tag,
					    stress1p, strain1p, stress2p, strain2p, stress3p, strain3p, stress4p, strain4p,
					    stress1n, strain1n, stress2n, strain2n, stress3n, strain3n, stress4n, strain4n,
					    rDispP, rForceP, uForceP, rDispN, rForceN, uForceN, 
					    gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit,
					    gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit,
					    gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit, gammaE, yStr);
    }
    if (argc == 31) {
      theMaterial = new ShearPanelMaterial (tag,
					    stress1p, strain1p, stress2p, strain2p, stress3p, strain3p, stress4p, strain4p,
					    rDispP, rForceP, uForceP,  
					    gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit,
					    gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit,
					    gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit, gammaE, yStr);		
    }
  }
    
  else if ((strcmp(argv[1],"Bond_SP01") == 0) || (strcmp(argv[1],"Bond") == 0)) {  //%strain penetration material
      // Check that there is the minimum number of arguments
      if (argc < 9) {
	opserr << "WARNING insufficient arguments\n";
	printCommand(argc,argv);
	opserr << "Want: uniaxialMaterial Bond_SP01 tag? fy? sy? fu? su? b? R?";
	opserr << " <Cd? db? fc? la?>" << endln;	
	return TCL_ERROR;
      }
      
      int tag;
      
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING invalid uniaxialMaterial Bond_SP01 tag" << endln;
	return TCL_ERROR;
      }
      
      // Read required Bond_SP01 material parameters
      double fy, sy, fu, su, Kz, R;
      
      if (Tcl_GetDouble(interp, argv[3], &fy) != TCL_OK) {
	opserr << "WARNING invalid bar yield strength (ksi): fy\n";
	opserr << "uniaxialMaterial Bond_SP01: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[4], &sy) != TCL_OK) {
	opserr << "WARNING invalid slip (in.) @ bar yield: sy\n";
	opserr << "uniaxialMaterial Bond_SP01: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[5], &fu) != TCL_OK) {
	opserr << "WARNING invalid bar failure strength (1.57fy? ksi) fu\n";
	opserr << "uniaxialMaterial Bond_SP01: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[6], &su) != TCL_OK) {
	opserr << "WARNING invalid slip (in.) @ bar failure: su\n";
	opserr << "uniaxialMaterial Bond_SP01: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[7], &Kz) != TCL_OK) {
	opserr << "WARNING invalid hardening ratio for envelop(<0.25<b(0.3)<0.5?): Cr\n";
	opserr << "uniaxialMaterial Bond_SP01: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[8], &R) != TCL_OK) {
	opserr << "WARNING invalid pinching factor (0.5<R<1.0?): R\n";
	opserr << "uniaxialMaterial Bond_SP01: " << tag << endln;
	return TCL_ERROR;
      }
      
      // Read optional Bond_SP01 material parameters (reserved for non-fully anchored cases
      double Cd, db, fc, la;
      
      if (argc > 9) {
	if (argc < 13) {
	  opserr << "WARNING insufficient number of Bond_SP01 parameters\n";
	  opserr << "uniaxialMaterial Bond_SP01: " << tag << endln;
	  return TCL_ERROR;
	}
	
	if (Tcl_GetDouble(interp, argv[9], &Cd) != TCL_OK) {
	  opserr << "WARNING invalid bond damage factor (0<Cd<1?): Cd\n";
	  opserr << "uniaxialMaterial Bond_SP01: " << tag << endln;
	  return TCL_ERROR;
	}
	
	if (Tcl_GetDouble(interp, argv[10], &db) != TCL_OK) {
	  opserr << "WARNING invalid bar diameter (in.): db\n";
	  opserr << "uniaxialMaterial Bond_SP01: " << tag << endln;
	  return TCL_ERROR;
	}
	
	if (Tcl_GetDouble(interp, argv[11], &fc) != TCL_OK) {
	  opserr << "WARNING invalid concrete strength (ksi): fc\n";
	  opserr << "uniaxialMaterial Bond_SP01: " << tag << endln;
	  return TCL_ERROR;
	}
	
	if (Tcl_GetDouble(interp, argv[12], &la) != TCL_OK) {
	  opserr << "WARNING invalid embedded length (in.): la\n";
	  opserr << "uniaxialMaterial Bond_SP01: " << tag << endln;
	  return TCL_ERROR;
	}
	
	// Parsing was successful, allocate the material
	theMaterial = new Bond_SP01 (tag, fy, sy, fu, su, Kz, R, Cd, db, fc, la);
      }
      else
	// Parsing was successful, allocate the material
	theMaterial = new Bond_SP01 (tag, fy, sy, fu, su, Kz, R);
      
    }		//%strain penetration material


    else if (strcmp(argv[1],"Concrete01WithSITC") == 0) {
      if (argc < 7) {
	  opserr << "WARNING insufficient arguments\n";
	  printCommand(argc,argv);
	  opserr << "Want: uniaxialMaterial Concrete01 tag? fpc? epsc0? fpcu? epscu? <endStrainSITC?>" << endln;
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

	if (argc == 7)
	  // Parsing was successful, allocate the material
	  theMaterial = new Concrete01WithSITC(tag, fpc, epsc0, fpcu, epscu);
	else {
	  double endStrainSITC;
	  if (Tcl_GetDouble(interp, argv[7], &endStrainSITC) != TCL_OK) {
	    opserr << "WARNING invalid epscu\n";
	    opserr << "Concrete01 material: " << tag << endln;
	    return TCL_ERROR;
	  }
	  theMaterial = new Concrete01WithSITC(tag, fpc, epsc0, fpcu, epscu, endStrainSITC);
	}
    }


    else if (strcmp(argv[1],"ECC01") == 0) {
      if (argc < 16) {
	opserr << "WARNING insufficient arguments\n";
	printCommand(argc,argv);
	opserr << "Want: uniaxialMaterial ECC01 TAG? SIGT0? EPST0? SIGT1? EPST1? EPST2? SIGC0? EPSC0? EPSC1? ";
	opserr << "ALPHAT1? ALPHAT2? ALPHAC? ALPHACU? BETAT? BETAC\n";
	return TCL_ERROR;
      }
      
      int tag;
      double SIGT0, EPST0, SIGT1, EPST1, EPST2, SIGC0, EPSC0, EPSC1, ALPHAT1, ALPHAT2, ALPHAC, ALPHACU, BETAT, BETAC;
      
      
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING invalid uniaxialMaterial ECC01 tag" << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[3], &SIGT0) != TCL_OK) {
	opserr << "WARNING invalid SIGTO\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[4], &EPST0) != TCL_OK) {
	opserr << "WARNING invalid EPSTO\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[5], &SIGT1) != TCL_OK) {
	opserr << "WARNING invalid SIGT1\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[6], &EPST1) != TCL_OK) {
	opserr << "WARNING invalid EPST1\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[7], &EPST2) != TCL_OK) {
	opserr << "WARNING invalid epscu\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[8], &SIGC0) != TCL_OK) {
	opserr << "WARNING invalid epscu\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[9], &EPSC0) != TCL_OK) {
	opserr << "WARNING invalid epscu\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[10], &EPSC1) != TCL_OK) {
	opserr << "WARNING invalid epscu\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[11], &ALPHAT1) != TCL_OK) {
	opserr << "WARNING invalid epscu\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[12], &ALPHAT2) != TCL_OK) {
	opserr << "WARNING invalid epscu\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[13], &ALPHAC) != TCL_OK) {
	opserr << "WARNING invalid epscu\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[14], &ALPHACU) != TCL_OK) {
	opserr << "WARNING invalid epscu\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[15], &BETAT) != TCL_OK) {
	opserr << "WARNING invalid epscu\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[16], &BETAC) != TCL_OK) {
	opserr << "WARNING invalid epscu\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      
      theMaterial = new ECC01(tag, SIGT0, EPST0, SIGT1, EPST1, EPST2, SIGC0, EPSC0, EPSC1, 
			      ALPHAT1, ALPHAT2, ALPHAC, ALPHACU, BETAT, BETAC);
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
	theMaterial = TclModelBuilder_addPyTzQzMaterial(clientData, interp, argc, argv, theTclBuilder, theDomain);
      
#ifdef _LIMITSTATEMATERIAL
      // LimitState
      if (theMaterial == 0)
	theMaterial = Tcl_AddLimitStateMaterial(clientData, interp, argc, argv);
#endif
    }
    
    if (theMaterial == 0) {
      
      //
      // maybe element in a package already loaded
      //  loop through linked list of loaded functions comparing names & if find call it
      //
      
      UniaxialPackageCommand *matCommands = theUniaxialPackageCommands;
      bool found = false;
      while (matCommands != NULL && found == false) {
	if (strcmp(argv[1], matCommands->funcName) == 0) {
	  
	  int result = (*(matCommands->funcPtr))(clientData, interp, argc, argv, theTclBuilder);
	  return result;
	} else
	  matCommands = matCommands->next;
      }
      
      
      //
      // maybe element command exists in a dll in library path
      //  so try loading package of same name as material name containg
      //  a c function "TclCommand_MatName"
      //
      
      void *libHandle;
      int (*funcPtr)(ClientData clientData, Tcl_Interp *interp,  int argc, 
		     TCL_Char **argv, TclModelBuilder *);       
      int matNameLength = strlen(argv[1]);
      char *tclFuncName = new char[matNameLength+12];
      strcpy(tclFuncName, "TclCommand_");
      strcpy(&tclFuncName[11], argv[1]);    
      int res = getLibraryFunction(argv[1], tclFuncName, &libHandle, (void **)&funcPtr);
      
      delete [] tclFuncName;
      
      if (res == 0) {
	
	//
	// add loaded function to list of functions
	//
	char *matName = new char[matNameLength+1];
	strcpy(matName, argv[1]);
	UniaxialPackageCommand *theMatCommand = new UniaxialPackageCommand;
	theMatCommand->funcPtr = funcPtr;
	theMatCommand->funcName = matName;	
	theMatCommand->next = theUniaxialPackageCommands;
	theUniaxialPackageCommands = theMatCommand;
	
	int result = (*funcPtr)(clientData, interp,
				argc, argv,
				theTclBuilder);	
	
	return result;
      }
    }
    
    //
    // if still here the element command does not exist
    //
    
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
