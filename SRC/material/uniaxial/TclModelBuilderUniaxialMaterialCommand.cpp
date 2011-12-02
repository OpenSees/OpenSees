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
                                                                        
// $Revision: 1.2 $
// $Date: 2000-12-13 05:53:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/TclModelBuilderUniaxialMaterialCommand.cpp,v $
                                                                        
                                                                        
// File: ~/material/uniaxial/TclModelBuilderUniaxialMaterialCommand.C
// 
// Written: fmk, MHS 
// Created: 07/99
// Revision: A
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
#include <SeriesMaterial.h>		// MHS

//#include <BiLinear.h>			// Rohit
//#include <Clough1.h>			// Rohit
//#include <Clough2.h>			// Rohit
//#include <Pinch1.h>				// Rohit

#include <string.h>

static void printCommand(int argc, char **argv)
{
    cerr << "Input command: ";
    for (int i=0; i<argc; i++)
	cerr << argv[i] << " ";
    cerr << endl;
} 

int
TclModelBuilderUniaxialMaterialCommand (ClientData clienData, Tcl_Interp *interp, int argc,
					char **argv, TclModelBuilder *theTclBuilder)
{
    // Make sure there is a minimum number of arguments
    if (argc < 3) {
	cerr << "WARNING insufficient number of uniaxial material arguments\n";
	cerr << "Want: uniaxialMaterial type? tag? <specific material args>" << endl;
	return TCL_ERROR;
    }
    
    // Pointer to a uniaxial material that will be added to the model builder
    UniaxialMaterial *theMaterial = 0;

    // Check argv[2] for uniaxial material type
    if (strcmp(argv[1],"Elastic") == 0) {
	if (argc < 4 || argc > 5) {
	    cerr << "WARNING invalid number of arguments\n";
	    printCommand(argc,argv);
	    cerr << "Want: uniaxialMaterial Elastic tag? E? <eta?>" << endl;
	    return TCL_ERROR;
	}    

	int tag;
	double E;
        double eta = 0.0;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    cerr << "WARNING invalid uniaxialMaterial Elastic tag" << endl;
	    return TCL_ERROR;		
	}

	if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
	    cerr << "WARNING invalid E\n";
	    cerr << "uniaxiaMaterial Elastic: " << tag << endl;
	    return TCL_ERROR;	
	}

	if (argc == 5) {
	    if (Tcl_GetDouble(interp,argv[4], &eta) != TCL_OK) {
               cerr << "WARNING invalid eta\n";
               cerr << "uniaxialMaterial Elastic: " << tag << endl;
               return TCL_ERROR;
            }
	}

	// Parsing was successful, allocate the material
	theMaterial = new ElasticMaterial(tag, E, eta);       
        //dum << tag " << tag << " E " << E << " eta " << eta <<endl;
    }

    else if (strcmp(argv[1],"ElasticPP") == 0) {
	if (argc < 5) {
	    cerr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    cerr << "Want: uniaxialMaterial ElasticPP tag? E? epsy?" << endl;
	    return TCL_ERROR;
	}

	int tag;
	double E, ep;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    cerr << "WARNING invalid uniaxialMaterial ElasticPP tag" << endl;
	    return TCL_ERROR;		
	}

	if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
	    cerr << "WARNING invalid E\n";
	    cerr << "uniaxialMaterial ElasticPP: " << tag << endl;
	    return TCL_ERROR;	
	}

	if (Tcl_GetDouble(interp, argv[4], &ep) != TCL_OK) {
	    cerr << "WARNING invalid epsy\n";
	    cerr << "uniaxialMaterial ElasticPP: " << tag << endl;
	    return TCL_ERROR;
	}
	
	// read in additional parameters eyn and ezero
	if (argc > 6) {
	    double eyn, ez;
	    if (Tcl_GetDouble(interp, argv[5], &eyn) != TCL_OK) {
		cerr << "WARNING invalid eyn\n";
		cerr << "uniaxialMaterial ElasticPP: " << tag << endl;
		return TCL_ERROR;
	    }	    
	    if (Tcl_GetDouble(interp, argv[6], &ez) != TCL_OK) {
		cerr << "WARNING invalid ez\n";
		cerr << "uniaxialMaterial ElasticPP: " << tag << endl;
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
	    cerr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    cerr << "Want: uniaxialMaterial ElasticPPGap tag? E? fy? gap?" << endl;
	    return TCL_ERROR;
	}

	int tag;
	double E, fy, gap;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    cerr << "WARNING invalid uniaxialMaterial ElasticPPGap tag" << endl;
	    return TCL_ERROR;		
	}

	if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
	    cerr << "WARNING invalid E\n";
	    cerr << "uniaxialMaterial ElasticPPGap: " << tag << endl;
	    return TCL_ERROR;	
	}

	if (Tcl_GetDouble(interp, argv[4], &fy) != TCL_OK) {
	    cerr << "WARNING invalid fy\n";
	    cerr << "uniaxialMaterial ElasticPPGap: " << tag << endl;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[5], &gap) != TCL_OK) {
	    cerr << "WARNING invalid gap\n";
	    cerr << "uniaxialMaterial ElasticPPGap: " << tag << endl;
	    return TCL_ERROR;
	}

	// Parsing was successful, allocate the material
	theMaterial = new EPPGapMaterial(tag, E, fy, gap);       
    }

    else if (strcmp(argv[1],"Hardening") == 0) {
	if (argc < 7) {
	    cerr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    cerr << "Want: uniaxialMaterial Hardening tag? E? sigmaY? H_iso? H_kin?" << endl;
	    return TCL_ERROR;
	}

	int tag;
	double E, sigmaY, Hiso, Hkin;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    cerr << "WARNING invalid uniaxialMaterial Hardening tag" << endl;
	    return TCL_ERROR;		
	}

	if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
	    cerr << "WARNING invalid E\n";
	    cerr << "uniaxialMaterial Hardening: " << tag << endl;
	    return TCL_ERROR;	
	}

	if (Tcl_GetDouble(interp, argv[4], &sigmaY) != TCL_OK) {
	    cerr << "WARNING invalid sigmaY\n";
	    cerr << "uniaxialMaterial Hardening: " << tag << endl;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[5], &Hiso) != TCL_OK) {
	    cerr << "WARNING invalid H_iso\n";
	    cerr << "uniaxialMaterial Hardening: " << tag << endl;
	    return TCL_ERROR;	
	}

	if (Tcl_GetDouble(interp, argv[6], &Hkin) != TCL_OK) {
	    cerr << "WARNING invalid H_kin\n";
	    cerr << "uniaxialMaterial Hardening: " << tag << endl;
	    return TCL_ERROR;	
	}
	// Parsing was successful, allocate the material
	theMaterial = new HardeningMaterial (tag, E, sigmaY, Hiso, Hkin);       
    }

    else if (strcmp(argv[1],"Parallel") == 0) {
	if (argc < 4) {
	    cerr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    cerr << "Want: uniaxialMaterial Parallel tag? tag1? tag2? ...";
	    cerr << " <-min min?> <-max max?>" << endl;
	    return TCL_ERROR;
	}
 
	int tag;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    cerr << "WARNING invalid uniaxialMaterial Parallel tag" << endl;
	    return TCL_ERROR;		
	}

	// Search for min and max strains
	double epsmin = NEG_INF_STRAIN;
	double epsmax = POS_INF_STRAIN;
	int numMaterials = argc-3;
	
	for (int j = 3; j < argc; j++) {
	    if (strcmp(argv[j],"-min") == 0) {
		if ((j+1) >= argc || Tcl_GetDouble (interp, argv[j+1], &epsmin) != TCL_OK) {
		    cerr << "WARNING invalid min\n";
		    cerr << "uniaxialMaterial Parallel: " << tag << endl;
		    return TCL_ERROR;
		}
		j++;
		numMaterials -= 2;
	    }
	    if (strcmp(argv[j],"-max") == 0) {
		if ((j+1) >= argc || Tcl_GetDouble (interp, argv[j+1], &epsmax) != TCL_OK) {
		    cerr << "WARNING invalid max\n";
		    cerr << "uniaxialMaterial Parallel: " << tag << endl;
		    return TCL_ERROR;
		}
		j++;
		numMaterials -= 2;
	    }
	}
	
	if (numMaterials == 0) {
	    cerr << "WARNING no component material(s) provided\n";
	    cerr << "uniaxialMaterial Parallel: " << tag << endl;
	    return TCL_ERROR;
	}
    
	// Create an array to hold pointers to component materials
	UniaxialMaterial **theMats = new UniaxialMaterial *[numMaterials];
	
	// For each material get the tag and ensure it exists in model already
	for (int i=0; i<numMaterials; i++) {
	    int tagI;
	    if (Tcl_GetInt(interp, argv[i+3], &tagI) != TCL_OK) {
		cerr << "WARNING invalid component tag\n";
		cerr << "uniaxialMaterial Parallel: " << tag << endl;
		return TCL_ERROR;
	    }
	    
	    UniaxialMaterial *theMat = theTclBuilder->getUniaxialMaterial(tagI);
	    
	    if (theMat == 0) {
		cerr << "WARNING component material does not exist\n";
		cerr << "Component material: " << argv[i+3]; 
		cerr << "\nuniaxialMaterial Parallel: " << tag << endl;
		delete [] theMats;
		return TCL_ERROR;
	    }
	    else
		theMats[i] = theMat;
	}	
	
	// Parsing was successful, allocate the material
	theMaterial = new ParallelMaterial(tag, numMaterials, theMats, epsmin, epsmax);
	
	// Deallocate the temporary pointers
	delete [] theMats;
    }

    else if (strcmp(argv[1],"Series") == 0) {
	if (argc < 4) {
	    cerr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    cerr << "Want: uniaxialMaterial Series tag? tag1? tag2? ..." << endl;
	    return TCL_ERROR;
	}
 
	int tag;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    cerr << "WARNING invalid uniaxialMaterial Series tag" << endl;
	    return TCL_ERROR;		
	}

	int numMaterials = argc - 3;
	UniaxialMaterial **theMats = new UniaxialMaterial *[numMaterials];

	// For each material get the tag and ensure it exists in model already
	for (int i = 3; i < argc; i++) {
	    int tagI;
	    if (Tcl_GetInt(interp, argv[i], &tagI) != TCL_OK) {
			cerr << "WARNING invalid component tag\n";
			cerr << "uniaxialMaterial Series: " << tag << endl;
			delete [] theMats;
			return TCL_ERROR;
	    }
	    
	    UniaxialMaterial *theMat = theTclBuilder->getUniaxialMaterial(tagI);
	    
	    if (theMat == 0) {
			cerr << "WARNING component material does not exist\n";
			cerr << "Component material: " << tagI;
			cerr << "\nuniaxialMaterial Series: " << tag << endl;
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
	    cerr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    cerr << "Want: uniaxialMaterial Steel01 tag? fy? E0? b?";
	    cerr << " <a1? a2? a3? a4?> <-min min?> <-max max?>" << endl;	
	    return TCL_ERROR;
	}

	int tag;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    cerr << "WARNING invalid uniaxialMaterial Steel01 tag" << endl;
	    return TCL_ERROR;
	}

	// Search for min and max strains
	double epsmin = NEG_INF_STRAIN;
	double epsmax = POS_INF_STRAIN;
	int numMatParam = argc - 3;
	
	for (int j = 3; j < argc; j++) {
	    if (strcmp(argv[j],"-min") == 0) {
		if ((j+1) >= argc || Tcl_GetDouble (interp, argv[j+1], &epsmin) != TCL_OK) {
		    cerr << "WARNING invalid min\n";
		    cerr << "uniaxialMaterial Steel01: " << tag << endl;
		    return TCL_ERROR;
		}
		j++;
		numMatParam -= 2;
	    }
	    if (strcmp(argv[j],"-max") == 0) {
		if ((j+1) >= argc || Tcl_GetDouble (interp, argv[j+1], &epsmax) != TCL_OK) {
		    cerr << "WARNING invalid max\n";
		    cerr << "uniaxialMaterial Steel01: " << tag << endl;
		    return TCL_ERROR;
		}
		j++;
		numMatParam -= 2;
	    }
	}	
	
	// Read required Steel01 material parameters
	double fy, E, b;
	
	if (Tcl_GetDouble(interp, argv[3], &fy) != TCL_OK) {
	    cerr << "WARNING invalid fy\n";
	    cerr << "uniaxialMaterial Steel01: " << tag << endl;
	    return TCL_ERROR;
	}
	
	if (Tcl_GetDouble(interp, argv[4], &E) != TCL_OK) {
	    cerr << "WARNING invalid E0\n";
	    cerr << "uniaxialMaterial Steel01: " << tag << endl;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[5], &b) != TCL_OK) {
	    cerr << "WARNING invalid b\n";
	    cerr << "uniaxialMaterial Steel01: " << tag << endl;
	    return TCL_ERROR;
	}

	// Read optional Steel01 material parameters
	double a1 = STEEL_01_DEFAULT_A1;
	double a2 = STEEL_01_DEFAULT_A2;
	double a3 = STEEL_01_DEFAULT_A3;
	double a4 = STEEL_01_DEFAULT_A4;

	if (numMatParam > 3) {
	    if (numMatParam < 7) {
		cerr << "WARNING insufficient number of hardening parameters\n";
		cerr << "uniaxialMaterial Steel01: " << tag << endl;
		return TCL_ERROR;
	    }
	    
	    if (Tcl_GetDouble(interp, argv[6], &a1) != TCL_OK) {
		cerr << "WARNING invalid a1\n";
		cerr << "uniaxialMaterial Steel01: " << tag << endl;
		return TCL_ERROR;
	    }
	    
	    if (Tcl_GetDouble(interp, argv[7], &a2) != TCL_OK) {
		cerr << "WARNING invalid a2\n";
		cerr << "uniaxialMaterial Steel01: " << tag << endl;
		return TCL_ERROR;
	    }
	    
	    if (Tcl_GetDouble(interp, argv[8], &a3) != TCL_OK) {
		cerr << "WARNING invalid a3\n";
		cerr << "uniaxialMaterial Steel01: " << tag << endl;
		return TCL_ERROR;
	    }

	    if (Tcl_GetDouble(interp, argv[9], &a4) != TCL_OK) {
		cerr << "WARNING invalid a4\n";
		cerr << "uniaxialMaterial Steel01: " << tag << endl;
		return TCL_ERROR;
	    }
	}

	// Parsing was successful, allocate the material
	theMaterial = new Steel01 (tag, fy, E, b, a1, a2, a3, a4, epsmin, epsmax);
    }

    else if (strcmp(argv[1],"Concrete01") == 0) {
	if (argc < 7) {
	    cerr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    cerr << "Want: uniaxialMaterial Concrete01 tag? fpc? epsc0? fpcu? epscu?";
	    cerr << " <-min min?> <-max max?>" << endl;
	    return TCL_ERROR;
	}

	int tag;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    cerr << "WARNING invalid uniaxialMaterial Concrete01 tag" << endl;
	    return TCL_ERROR;
	}

	// Search for min and max strain
	double epsmin = NEG_INF_STRAIN;
	double epsmax = POS_INF_STRAIN;
	
	for (int j = 3; j < argc; j++) {
	    if (strcmp(argv[j],"-min") == 0) {
		if ((j+1) >= argc || Tcl_GetDouble(interp, argv[j+1], &epsmin) != TCL_OK) {
		    cerr << "WARNING invalid min\n";
		    cerr << "Concrete01 material: " << tag << endl;
		    return TCL_ERROR;
		}
		j++;
	    }

	    if (strcmp(argv[j],"-max") == 0) {
		if ((j+1) >= argc || Tcl_GetDouble(interp, argv[j+1], &epsmax) != TCL_OK) {
		    cerr << "WARNING invalid max\n";
		    cerr << "Concrete01 material: " << tag << endl;
		    return TCL_ERROR;
		}
		j++;
	    }
	}

	// Read required Concrete01 material parameters
	double fpc, epsc0, fpcu, epscu;

	if (Tcl_GetDouble(interp, argv[3], &fpc) != TCL_OK) {
	    cerr << "WARNING invalid fpc\n";
	    cerr << "Concrete01 material: " << tag << endl;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[4], &epsc0) != TCL_OK) {
	    cerr << "WARNING invalid epsc0\n";
	    cerr << "Concrete01 material: " << tag << endl;
	    return TCL_ERROR;
	}
	
	if (Tcl_GetDouble(interp, argv[5], &fpcu) != TCL_OK) {
	    cerr << "WARNING invalid fpcu\n";
	    cerr << "Concrete01 material: " << tag << endl;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[6], &epscu) != TCL_OK) {
	    cerr << "WARNING invalid epscu\n";
	    cerr << "Concrete01 material: " << tag << endl;
	    return TCL_ERROR;
	}

	// Parsing was successful, allocate the material
	theMaterial = new Concrete01(tag, fpc, epsc0, fpcu, epscu, epsmin, epsmax);
    }

	else if (strcmp(argv[1],"Hysteretic") == 0) {
		if (argc != 20 && argc != 19 && argc != 16 && argc != 15) {
			cerr << "WARNING insufficient arguments\n";
			printCommand(argc,argv);
			cerr << "Want: uniaxialMaterial Hysteretic tag? mom1p? rot1p? mom2p? rot2p? <mom3p? rot3p?> "
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
			cerr << "WARNING invalid uniaxialMaterial Hysteretic tag" << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &mom1p) != TCL_OK) {
			cerr << "WARNING invalid mom1p\n";
			cerr << "Hysteretic material: " << tag << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &rot1p) != TCL_OK) {
			cerr << "WARNING invalid rot1p\n";
			cerr << "Hysteretic material: " << tag << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &mom2p) != TCL_OK) {
			cerr << "WARNING invalid mom2p\n";
			cerr << "Hysteretic material: " << tag << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &rot2p) != TCL_OK) {
			cerr << "WARNING invalid rot2p\n";
			cerr << "Hysteretic material: " << tag << endl;
			return TCL_ERROR;
		}

		if (argc > 16) {
			if (Tcl_GetDouble(interp, argv[i++], &mom3p) != TCL_OK) {
				cerr << "WARNING invalid mom3p\n";
				cerr << "Hysteretic material: " << tag << endl;
				return TCL_ERROR;
			}
			
			if (Tcl_GetDouble(interp, argv[i++], &rot3p) != TCL_OK) {
				cerr << "WARNING invalid rot3p\n";
				cerr << "Hysteretic material: " << tag << endl;
				return TCL_ERROR;
			}
		}

		if (Tcl_GetDouble(interp, argv[i++], &mom1n) != TCL_OK) {
			cerr << "WARNING invalid mom1n\n";
			cerr << "Hysteretic material: " << tag << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &rot1n) != TCL_OK) {
			cerr << "WARNING invalid rot1n\n";
			cerr << "Hysteretic material: " << tag << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &mom2n) != TCL_OK) {
			cerr << "WARNING invalid mom2n\n";
			cerr << "Hysteretic material: " << tag << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &rot2n) != TCL_OK) {
			cerr << "WARNING invalid rot2n\n";
			cerr << "Hysteretic material: " << tag << endl;
			return TCL_ERROR;
		}

		if (argc > 16) {
			if (Tcl_GetDouble(interp, argv[i++], &mom3n) != TCL_OK) {
				cerr << "WARNING invalid mom3n\n";
				cerr << "Hysteretic material: " << tag << endl;
				return TCL_ERROR;
			}
	
			if (Tcl_GetDouble(interp, argv[i++], &rot3n) != TCL_OK) {
				cerr << "WARNING invalid rot3n\n";
				cerr << "Hysteretic material: " << tag << endl;
				return TCL_ERROR;
			}
		}

		if (Tcl_GetDouble(interp, argv[i++], &pinchX) != TCL_OK) {
			cerr << "WARNING invalid pinchX\n";
			cerr << "Hysteretic material: " << tag << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &pinchY) != TCL_OK) {
			cerr << "WARNING invalid pinchY\n";
			cerr << "Hysteretic material: " << tag << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &damfc1) != TCL_OK) {
			cerr << "WARNING invalid damfc1\n";
			cerr << "Hysteretic material: " << tag << endl;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &damfc2) != TCL_OK) {
			cerr << "WARNING invalid damfc2\n";
			cerr << "Hysteretic material: " << tag << endl;
			return TCL_ERROR;
		}

		if (argc == 20 || argc == 16) {
			if (Tcl_GetDouble(interp, argv[i++], &beta) != TCL_OK) {
				cerr << "WARNING invalid beta\n";
				cerr << "Hysteretic material: " << tag << endl;
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
			cerr << "WARNING insufficient arguments\n";
			printCommand(argc,argv);
			cerr << "Want: uniaxialMaterial Viscous tag? C? alpha?" << endl;
			return TCL_ERROR;
		}

		int tag;
		double C, a;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
		{
			cerr << "WARNING invalid tag\n";
			cerr << "Viscous material: " << tag << endl;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[3], &C) != TCL_OK)
		{
			cerr << "WARNING invalid C\n";
			cerr << "Viscous material: " << tag << endl;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[4], &a) != TCL_OK)
		{
			cerr << "WARNING invalid alpha\n";
			cerr << "Viscous material: " << tag << endl;
			return TCL_ERROR;
		}

		theMaterial = new ViscousMaterial (tag, C, a);
	}

	else if (strcmp(argv[1],"PathIndependent") == 0) {
		if (argc < 4)
		{
			cerr << "WARNING insufficient arguments\n";
			printCommand(argc,argv);
			cerr << "Want: uniaxialMaterial PathIndependent tag? matTag?" << endl;
			return TCL_ERROR;
		}

		int tag, matTag;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
		{
			cerr << "WARNING invalid tag\n";
			cerr << "PathIndependent material: " << tag << endl;
			return TCL_ERROR;
		}
		
		if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK)
		{
			cerr << "WARNING invalid matTag\n";
			cerr << "PathIndependent material: " << tag << endl;
			return TCL_ERROR;
		}

		UniaxialMaterial *material = theTclBuilder->getUniaxialMaterial(matTag);
		
		if (material == 0) {
		    cerr << "WARNING material does not exist\n";
		    cerr << "material: " << matTag; 
		    cerr << "\nuniaxialMaterial PathIndependent: " << tag << endl;
		    return TCL_ERROR;
		}

		theMaterial = new PathIndependentMaterial (tag, *material);
	}

	/*
	else if (strcmp(argv[1],"BiLinear") == 0) {
		if (argc < 19)
		{
			cerr << "WARNING insufficient arguments\n";
			printCommand(argc,argv);
			cerr << "Want: uniaxialMaterial BiLinear tag? ..." << endl;
			return TCL_ERROR;
		}
		
		int tag;
		Vector input(16);
		double temp;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
			cerr << "WARNING invalid tag\n";
			cerr << "BiLinear material: " << tag << endl;
			return TCL_ERROR;
		}

		for (int i = 3, j = 0; j < 16; i++, j++) {
			if (Tcl_GetDouble(interp, argv[i], &temp) != TCL_OK) {
				cerr << "WARNING invalid input, data " << i << '\n';
				cerr << "BiLinear material: " << tag << endl;
				return TCL_ERROR;
			}
			input(j) = temp;
		}

		theMaterial = new BiLinear(tag, input);
	}

	else if (strcmp(argv[1],"Clough1") == 0) {
		if (argc < 19)
		{
			cerr << "WARNING insufficient arguments\n";
			printCommand(argc,argv);
			cerr << "Want: uniaxialMaterial Clough1 tag? ..." << endl;
			return TCL_ERROR;
		}
		
		int tag;
		Vector input(16);
		double temp;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
			cerr << "WARNING invalid tag\n";
			cerr << "Clough1 material: " << tag << endl;
			return TCL_ERROR;
		}

		for (int i = 3, j = 0; j < 16; i++, j++) {
			if (Tcl_GetDouble(interp, argv[i], &temp) != TCL_OK) {
				cerr << "WARNING invalid input, data " << i << '\n';
				cerr << "Clough1 material: " << tag << endl;
				return TCL_ERROR;
			}
			input(j) = temp;
		}

		theMaterial = new Clough1(tag, input);
	}

	else if (strcmp(argv[1],"Clough2") == 0) {
		if (argc < 19)
		{
			cerr << "WARNING insufficient arguments\n";
			printCommand(argc,argv);
			cerr << "Want: uniaxialMaterial Clough2 tag? ..." << endl;
			return TCL_ERROR;
		}
		
		int tag;
		Vector input(16);
		double temp;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
			cerr << "WARNING invalid tag\n";
			cerr << "Clough2 material: " << tag << endl;
			return TCL_ERROR;
		}

		for (int i = 3, j = 0; j < 16; i++, j++) {
			if (Tcl_GetDouble(interp, argv[i], &temp) != TCL_OK) {
				cerr << "WARNING invalid input, data " << i << '\n';
				cerr << "Clough2 material: " << tag << endl;
				return TCL_ERROR;
			}
			input(j) = temp;
		}

		theMaterial = new Clough2(tag, input);
	}

	else if (strcmp(argv[1],"Pinch1") == 0) {
		if (argc < 22)
		{
			cerr << "WARNING insufficient arguments\n";
			printCommand(argc,argv);
			cerr << "Want: uniaxialMaterial Pinch1 tag? ..." << endl;
			return TCL_ERROR;
		}
		
		int tag;
		Vector input(19);
		double temp;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
			cerr << "WARNING invalid tag\n";
			cerr << "Pinch1 material: " << tag << endl;
			return TCL_ERROR;
		}

		for (int i = 3, j = 0; j < 19; i++, j++) {
			if (Tcl_GetDouble(interp, argv[i], &temp) != TCL_OK) {
				cerr << "WARNING invalid input, data " << i << '\n';
				cerr << "Pinch1 material: " << tag << endl;
				return TCL_ERROR;
			}
			input(j) = temp;
		}

		theMaterial = new Pinch1(tag, input);
	} */

    else  {
	cerr << "WARNING unknown type of uniaxialMaterial: " << argv[1];
	cerr << "\nValid types: Elastic, ElasticPP, Parallel,\n";
	cerr << "             Hardening, Steel01, Concrete01,\n";
	cerr << "             Hysteretic, Viscous, ElasticPPGap,\n";
	cerr << "             Series, PathIndependent" << endl;
	return TCL_ERROR;
    }

    // Ensure we have created the Material, out of memory if got here and no material
    if (theMaterial == 0) {
	cerr << "WARNING ran out of memory creating uniaxialMaterial\n";
	cerr << argv[1] << endl;
	return TCL_ERROR;
    }

    // Now add the material to the modelBuilder
    if (theTclBuilder->addUniaxialMaterial(*theMaterial) < 0) {
	cerr << "WARNING could not add uniaxialMaterial to the domain\n";
	cerr << *theMaterial << endl;
	delete theMaterial; // invoke the material objects destructor, otherwise mem leak
	return TCL_ERROR;
    }

    
    return TCL_OK;
}

