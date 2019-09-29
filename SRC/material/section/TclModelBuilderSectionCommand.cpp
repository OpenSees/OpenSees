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
                                                                        
// $Revision: 1.31 $
// $Date: 2009-12-17 20:10:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/TclModelBuilderSectionCommand.cpp,v $


// Written: rms, MHS 
// Created: 07/99
//
// Description: This file contains the function invoked when the user invokes
// the section command in the interpreter.
//
// What: "@(#) TclModelBuilderMaterialCommands.C, revA"

#include <TclModelBuilder.h>

#include <tcl.h>
#include <elementAPI.h>

#include <ElasticMaterial.h>

#include <ElasticSection2d.h>
#include <ElasticSection3d.h>
#include <ElasticShearSection2d.h>
#include <ElasticShearSection3d.h>
#include <ElasticWarpingShearSection2d.h>
#include <ElasticTubeSection3d.h>
//#include <GenericSection1d.h>
//#include <GenericSectionNd.h>
#include <SectionAggregator.h>
#include <ParallelSection.h>
//#include <FiberSection.h>
#include <FiberSection2d.h>
#include <NDFiberSection2d.h>
#include <NDFiberSection3d.h>
#include <NDFiberSectionWarping2d.h>
#include <FiberSection2dInt.h>
#include <FiberSection3d.h>
//#include <FiberSectionGJ.h>
#include <FiberSectionRepr.h>

#include <LayeredShellFiberSection.h> // Yuli Huang & Xinzheng Lu 

#include <ElasticPlateSection.h>
#include <ElasticMembranePlateSection.h>
#include <MembranePlateFiberSection.h>

#include <QuadPatch.h>
#include <CircPatch.h>
#include <QuadCell.h>
#include <StraightReinfLayer.h>
#include <CircReinfLayer.h>
#include <ReinfBar.h>

#include <UniaxialFiber2d.h>
#include <UniaxialFiber3d.h>
#include <NDFiber2d.h>
#include <NDFiber3d.h>

#include <Bidirectional.h>
#include <Elliptical2.h>
#include <Isolator2spring.h>

//#include <WSection2d.h>
#include <WideFlangeSectionIntegration.h>
#include <RCSectionIntegration.h>
#include <RCTBeamSectionIntegration.h>
#include <RCCircularSectionIntegration.h>
#include <RCTunnelSectionIntegration.h>
//#include <RCTBeamSectionIntegrationUniMat.h>
#include <TubeSectionIntegration.h>

//#include <McftSection2dfiber.h>

#include <string.h>
#include <fstream>
using std::ifstream;

#include <iostream>
using std::ios;

#include <packages.h>

extern int OPS_ResetInputNoBuilder(ClientData clientData, 
				   Tcl_Interp *interp,  
				   int cArg, 
				   int mArg, 
				   TCL_Char **argv, 
				   Domain *domain);

int
TclCommand_addFiberSection (ClientData clientData, Tcl_Interp *interp, int argc,
			    TCL_Char **argv, TclModelBuilder *theBuilder);

int
TclCommand_addFiberIntSection (ClientData clientData, Tcl_Interp *interp, int argc,
			       TCL_Char **argv, TclModelBuilder *theBuilder);


//--- Adding Thermo-mechanical Sections:[BEGIN]   by UoE OpenSees Group ---//  
#include <FiberSection2dThermal.h>
#include <FiberSection3dThermal.h> //Added by L.Jiang [SIF] 2017
#include <FiberSectionGJThermal.h> //Added by Liming, [SIF] 2017
#include <MembranePlateFiberSectionThermal.h> //Added by Liming, [SIF] 2017
#include <LayeredShellFiberSectionThermal.h> //Added by Liming, [SIF] 2017

int TclCommand_addFiberSectionThermal (ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv, TclModelBuilder *theBuilder);
int buildSectionThermal(Tcl_Interp *interp, TclModelBuilder *theTclModelBuilder,int secTag, bool isTorsion, double GJ);
//--- Adding Thermo-mechanical Sections: [END]   by UoE OpenSees Group ---//  



int
TclCommand_addUCFiberSection (ClientData clientData, Tcl_Interp *interp, int argc,
				   TCL_Char **argv, TclModelBuilder *theBuilder);


SectionForceDeformation *
TclModelBuilderYS_SectionCommand(ClientData clientData, Tcl_Interp *interp, int argc, 
				 TCL_Char **argv, TclModelBuilder *theTclBuilder);

int
TclModelBuilderSectionCommand (ClientData clientData, Tcl_Interp *interp, int argc,
			       TCL_Char **argv, TclModelBuilder *theTclBuilder)
{
  // Make sure there is a minimum number of arguments
    if (argc < 3) {
	opserr << "WARNING insufficient number of section arguments\n";
	opserr << "Want: section type? tag? <specific material args>" << endln;
	return TCL_ERROR;
    }

    //OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, theDomain);	  

    // Pointer to a section that will be added to the model builder
    SectionForceDeformation *theSection = 0;

    int NDM = theTclBuilder->getNDM();  
    
    // Check argv[1] for section type
    if (strcmp(argv[1],"Elastic") == 0) {
      if (argc < 5) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: section Elastic tag? E? A? Iz? <Iy? G? J?>" << endln;
	return TCL_ERROR;
      }
	
	int tag;
	double E, A, Iz, Iy, G, J, alphaY, alphaZ;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid section Elastic tag" << endln;
	    return TCL_ERROR;		
	}

	if (Tcl_GetDouble (interp, argv[3], &E) != TCL_OK) {
	    opserr << "WARNING invalid E" << endln;
	    opserr << "Elastic section: " << tag << endln;	    
	    return TCL_ERROR;
	}	

	if (Tcl_GetDouble (interp, argv[4], &A) != TCL_OK) {
	    opserr << "WARNING invalid A" << endln;
	    opserr << "Elastic section: " << tag << endln;	    
	    return TCL_ERROR;
	}	
	
	if (Tcl_GetDouble (interp, argv[5], &Iz) != TCL_OK) {
	    opserr << "WARNING invalid Iz" << endln;
	    opserr << "Elastic section: " << tag << endln;	    	    
	    return TCL_ERROR;
	}	
	
	if (NDM == 2) {
	  if (argc > 7) {
	    if (Tcl_GetDouble (interp, argv[6], &G) != TCL_OK) {
	      opserr << "WARNING invalid G" << endln;
	      opserr << "Elastic section: " << tag << endln;	    	    
	      return TCL_ERROR;
	    }

	    if (Tcl_GetDouble (interp, argv[7], &alphaY) != TCL_OK) {
	      opserr << "WARNING invalid alpha" << endln;
	      opserr << "Elastic section: " << tag << endln;	    	    
	      return TCL_ERROR;
	    }

	    theSection = new ElasticShearSection2d(tag, E, A, Iz, G, alphaY);
	  }
	  else 
	    theSection = new ElasticSection2d(tag, E, A, Iz);
	} else {
	  // 3D
	  if (argc < 8) {
	    opserr << "WARNING insufficient arguments\n";
	    opserr << "Want: section Elastic tag? E? A? Iz? Iy? G? J?" << endln;
	    return TCL_ERROR;
	  }

	  if (Tcl_GetDouble (interp, argv[6], &Iy) != TCL_OK) {
	    opserr << "WARNING invalid Iy" << endln;
	    opserr << "Elastic section: " << tag << endln;
	    return TCL_ERROR;
	  }
	       
	  if (Tcl_GetDouble (interp, argv[7], &G) != TCL_OK) {
	    opserr << "WARNING invalid G" << endln;
	    opserr << "Elastic section: " << tag << endln;	    
	    return TCL_ERROR;
	  }

	  if (Tcl_GetDouble (interp, argv[8], &J) != TCL_OK) {
	    opserr << "WARNING invalid J" << endln;
	    opserr << "Elastic section: " << tag << endln;	    
	    return TCL_ERROR;
	  }

	  if (argc > 9) {
	    if (Tcl_GetDouble (interp, argv[9], &alphaY) != TCL_OK) {
	      opserr << "WARNING invalid alphaY" << endln;
	      opserr << "Elastic section: " << tag << endln;	    
	      return TCL_ERROR;
	    }

	    if (Tcl_GetDouble (interp, argv[10], &alphaZ) != TCL_OK) {
	      opserr << "WARNING invalid alphaZ" << endln;
	      opserr << "Elastic section: " << tag << endln;	    
	      return TCL_ERROR;
	    }

	    theSection = new ElasticShearSection3d(tag, E, A, Iz, Iy,
                                               G, J, alphaY, alphaZ);
	  }
	  else 
	    theSection = new ElasticSection3d(tag, E, A, Iz, Iy, G, J);
	}
    }	

    else if (strcmp(argv[1],"ElasticWarpingShear") == 0) {
      if (argc < 11) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: section ElasticWarpingShear tag? E? A? Iz? G? alpha? J? B? C?>" << endln;
	return TCL_ERROR;
      }
	
      int tag;
      double E, A, Iz, G, alpha, J, B, C;
      
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING invalid section ElasticWarpingShearSection2d tag" << endln;
	return TCL_ERROR;		
      }
      
      if (Tcl_GetDouble (interp, argv[3], &E) != TCL_OK) {
	opserr << "WARNING invalid E" << endln;
	opserr << "ElasticWarpingShearSection2d section: " << tag << endln;	    
	return TCL_ERROR;
      }	
      
      if (Tcl_GetDouble (interp, argv[4], &A) != TCL_OK) {
	opserr << "WARNING invalid A" << endln;
	opserr << "ElasticWarpingShearSection2d section: " << tag << endln;	    
	return TCL_ERROR;
      }	
      
      if (Tcl_GetDouble (interp, argv[5], &Iz) != TCL_OK) {
	opserr << "WARNING invalid Iz" << endln;
	opserr << "ElasticWarpingShearSection2d section: " << tag << endln;	    	    
	return TCL_ERROR;
      }	
      
      if (Tcl_GetDouble (interp, argv[6], &G) != TCL_OK) {
	opserr << "WARNING invalid G" << endln;
	opserr << "ElasticWarpingShearSection2d section: " << tag << endln;	    	    
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble (interp, argv[7], &alpha) != TCL_OK) {
	opserr << "WARNING invalid alpha" << endln;
	opserr << "ElasticWarpingShearSection2d section: " << tag << endln;	    	    
	return TCL_ERROR;
      }
      if (Tcl_GetDouble (interp, argv[8], &J) != TCL_OK) {
	opserr << "WARNING invalid J" << endln;
	opserr << "ElasticWarpingShearSection2d section: " << tag << endln;	    	    
	return TCL_ERROR;
      }
      if (Tcl_GetDouble (interp, argv[9], &B) != TCL_OK) {
	opserr << "WARNING invalid B" << endln;
	opserr << "ElasticWarpingShearSection2d section: " << tag << endln;	    	    
	return TCL_ERROR;
      }
      if (Tcl_GetDouble (interp, argv[10], &C) != TCL_OK) {
	opserr << "WARNING invalid C" << endln;
	opserr << "ElasticWarpingShearSection2d section: " << tag << endln;	    	    
	return TCL_ERROR;
      }
      
      theSection = new ElasticWarpingShearSection2d(tag, E, A, Iz, G, alpha, J, B, C);
    }
    	
    // Check argv[1] for section type
    else if (strcmp(argv[1],"ElasticTube") == 0) {
      if (argc < 7) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: section ElasticTube tag? E? d? tw? G?" << endln;
	return TCL_ERROR;
      }
	
	int tag;
	double E, d, tw, G;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid section Elastic tag" << endln;
	    return TCL_ERROR;		
	}

	if (Tcl_GetDouble (interp, argv[3], &E) != TCL_OK) {
	    opserr << "WARNING invalid E" << endln;
	    opserr << "ElasticTube section: " << tag << endln;	    
	    return TCL_ERROR;
	}	

	if (Tcl_GetDouble (interp, argv[4], &d) != TCL_OK) {
	    opserr << "WARNING invalid d" << endln;
	    opserr << "ElasticTube section: " << tag << endln;	    
	    return TCL_ERROR;
	}	
	
	if (Tcl_GetDouble (interp, argv[5], &tw) != TCL_OK) {
	    opserr << "WARNING invalid tw" << endln;
	    opserr << "ElasticTube section: " << tag << endln;	    	    
	    return TCL_ERROR;
	}	

	if (Tcl_GetDouble (interp, argv[6], &G) != TCL_OK) {
	    opserr << "WARNING invalid G" << endln;
	    opserr << "ElasticTube section: " << tag << endln;	    	    
	    return TCL_ERROR;
	}

	theSection = new ElasticTubeSection3d(tag, E, d, tw, G);	
    }

    else if (strcmp(argv[1],"Generic1D") == 0 ||
	     strcmp(argv[1],"Generic1d") == 0 ||
	     strcmp(argv[1],"Uniaxial") == 0) {
	if (argc < 5) {
	    opserr << "WARNING insufficient arguments\n";
	    opserr << "Want: section Uniaxial tag? 1DTag? code?" << endln;
	    return TCL_ERROR;
	}

	int tag, uniTag, code;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid section Uniaxial tag" << endln;
	    return TCL_ERROR;		
	}

	if (Tcl_GetInt(interp, argv[3], &uniTag) != TCL_OK) {
	    opserr << "WARNING invalid 1DTag" << endln;
	    opserr << "Uniaxial section: " << tag << endln;	    
	    return TCL_ERROR;		
	}

	if (strcmp(argv[4],"Mz") == 0)
	    code = SECTION_RESPONSE_MZ;
	else if (strcmp(argv[4],"P") == 0)
	    code = SECTION_RESPONSE_P;
	else if (strcmp(argv[4],"Vy") == 0)
	    code = SECTION_RESPONSE_VY;
	else if (strcmp(argv[4],"My") == 0)
	    code = SECTION_RESPONSE_MY;
	else if (strcmp(argv[4],"Vz") == 0)
	    code = SECTION_RESPONSE_VZ;
	else if (strcmp(argv[4],"T") == 0)
	    code = SECTION_RESPONSE_T;
	else {
	    opserr << "WARNING invalid code" << endln;
	    opserr << "Uniaxial section: " << tag << endln;
	    return TCL_ERROR;		
	}
		
	// Retrieve the uniaxial material from the model builder
	UniaxialMaterial *theMat = OPS_getUniaxialMaterial(uniTag);
	
	if (theMat == 0) {
	    opserr << "WARNING uniaxial material does not exist\n";
	    opserr << "uniaxial material: " << uniTag; 
	    opserr << "\nUniaxial section: " << tag << endln;
	    return TCL_ERROR;
	}
	
	// Parsing was successful, allocate the section
	//theSection = new GenericSection1d (tag, *theMat, code);

	UniaxialMaterial *theMats[1];
	theMats[0] = theMat;
	ID codeID(1);
	codeID(0) = code;
	theSection = new SectionAggregator(tag, 1, theMats, codeID);
    }

    else if (strcmp(argv[1],"GenericND") == 0 || strcmp(argv[1],"GenericNd") == 0) {
      opserr << "section GenericND is no longer available" << endln;
      return TCL_ERROR;

      /*
	if (argc < 5) {
	    opserr << "WARNING insufficient arguments\n";
	    opserr << "Want: section GenericNd tag? NDTag? code?" << endln;
	    return TCL_ERROR;
	}
	
	int tag, NDTag;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid section GenericNd tag" << endln;
	    return TCL_ERROR;		
	}
	
	if (Tcl_GetInt(interp, argv[3], &NDTag) != TCL_OK) {
	    opserr << "WARNING invalid NDTag" << endln;
	    opserr << "GenericNd section: " << tag << endln;	    
	    return TCL_ERROR;		
	}
	
	ID code(argc-4);
	
	int i,j;
	
	// Read in the code
	for (i = 4, j = 0; i < argc; i++, j++) {
	    if (strcmp(argv[i],"Mz") == 0)
		code(j) = SECTION_RESPONSE_MZ;
	    else if (strcmp(argv[i],"P") == 0)
		code(j) = SECTION_RESPONSE_P;
	    else if (strcmp(argv[i],"Vy") == 0)
		code(j) = SECTION_RESPONSE_VY;
	    else if (strcmp(argv[i],"My") == 0)
		code(j) = SECTION_RESPONSE_MY;
	    else if (strcmp(argv[i],"Vz") == 0)
		code(j) = SECTION_RESPONSE_VZ;
	    else if (strcmp(argv[i],"T") == 0)
		code(j) = SECTION_RESPONSE_T;
	    else {
		opserr << "WARNING invalid GenericND code" << endln;
		opserr << "\nGenericND section: " << tag << endln;
		return TCL_ERROR;		
	    }
	}
	
	// Retrieve the uniaxial material from the model builder
	NDMaterial *theMat = theTclBuilder->getNDMaterial(NDTag);
	
	if (theMat == 0) {
	    opserr << "WARNING nD material does not exist\n";
	    opserr << "nD material: " << NDTag; 
	    opserr << "\nGenericNd section: " << tag << endln;
	    return TCL_ERROR;
	}
	
	// Parsing was successful, allocate the section
	theSection = new GenericSectionNd (tag, *theMat, code);
	*/
    }	

    else if (strcmp(argv[1],"WFSection2d") == 0 || strcmp(argv[1],"WSection2d") == 0) {
	if (argc < 10) {
	    opserr << "WARNING insufficient arguments\n";
	    opserr << "Want: section WFSection2d tag? matTag? d? tw? bf? tf? nfdw? nftf? <-nd shape?>" << endln;
	    return TCL_ERROR;
	}
	
	int tag, matTag;
	double d, tw, bf, tf;
	int nfdw, nftf;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid section WFSection2d tag" << endln;
	    return TCL_ERROR;		
	}

	if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
	    opserr << "WARNING invalid section WFSection2d matTag" << endln;
	    return TCL_ERROR;		
	}

	if (Tcl_GetDouble (interp, argv[4], &d) != TCL_OK) {
	    opserr << "WARNING invalid d" << endln;
	    opserr << "WFSection2d section: " << tag << endln;	    
	    return TCL_ERROR;
	}	

	if (Tcl_GetDouble (interp, argv[5], &tw) != TCL_OK) {
	    opserr << "WARNING invalid tw" << endln;
	    opserr << "WFSection2d section: " << tag << endln;	    
	    return TCL_ERROR;
	}	
	
	if (Tcl_GetDouble (interp, argv[6], &bf) != TCL_OK) {
	    opserr << "WARNING invalid bf" << endln;
	    opserr << "WFSection2d section: " << tag << endln;	    	    
	    return TCL_ERROR;
	}	

	if (Tcl_GetDouble (interp, argv[7], &tf) != TCL_OK) {
	    opserr << "WARNING invalid tf" << endln;
	    opserr << "WFSection2d section: " << tag << endln;	    	    
	    return TCL_ERROR;
	}	

	if (Tcl_GetInt (interp, argv[8], &nfdw) != TCL_OK) {
	    opserr << "WARNING invalid nfdw" << endln;
	    opserr << "WFSection2d section: " << tag << endln;	    	    
	    return TCL_ERROR;
	}	

	if (Tcl_GetInt (interp, argv[9], &nftf) != TCL_OK) {
	    opserr << "WARNING invalid nftf" << endln;
	    opserr << "WFSection2d section: " << tag << endln;	    	    
	    return TCL_ERROR;
	}	

	WideFlangeSectionIntegration wfsect(d, tw, bf, tf, nfdw, nftf);

	int numFibers = wfsect.getNumFibers();

	if (argc > 10) {

	  double shape = 1.0;
	  if (argc > 11) {
	    if (Tcl_GetDouble(interp, argv[11], &shape) != TCL_OK) {
	      opserr << "WARNING invalid shape" << endln;
	      opserr << "WFSection2d section: " << tag << endln;	    	    
	      return TCL_ERROR;
	    }
	  }

	  NDMaterial *theSteel = OPS_getNDMaterial(matTag);
	
	  if (theSteel == 0) {
	    opserr << "WARNING ND material does not exist\n";
	    opserr << "material: " << matTag; 
	    opserr << "\nWFSection2d section: " << tag << endln;
	    return TCL_ERROR;
	  }
		  
	  NDMaterial **theMats = new NDMaterial *[numFibers];
	  
	  wfsect.arrangeFibers(theMats, theSteel);

	  // Parsing was successful, allocate the section
	  theSection = 0;
	  if (strcmp(argv[10],"-nd") == 0)
	    theSection = new NDFiberSection2d(tag, numFibers, theMats, wfsect, shape);
	  if (strcmp(argv[10],"-ndWarping") == 0)
	    theSection = new NDFiberSectionWarping2d(tag, numFibers, theMats, wfsect, shape);

	  delete [] theMats;	  
	}
	else {
	  UniaxialMaterial *theSteel = OPS_getUniaxialMaterial(matTag);
	
	  if (theSteel == 0) {
	    opserr << "WARNING uniaxial material does not exist\n";
	    opserr << "material: " << matTag; 
	    opserr << "\nWFSection2d section: " << tag << endln;
	    return TCL_ERROR;
	  }
	  
	  UniaxialMaterial **theMats = new UniaxialMaterial *[numFibers];
	  
	  wfsect.arrangeFibers(theMats, theSteel);
	  
	  // Parsing was successful, allocate the section
	  theSection = new FiberSection2d(tag, numFibers, theMats, wfsect);

	  delete [] theMats;
	}
    }

    else if (strcmp(argv[1],"Tube") == 0) {
	if (argc < 8) {
	    opserr << "WARNING insufficient arguments\n";
	    opserr << "Want: section Tube tag? matTag? D? t? nfw? nfr?" << endln;
	    return TCL_ERROR;
	}
	
	int tag, matTag;
	double D, t;
	int nfw, nfr;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid section Tube tag" << endln;
	    return TCL_ERROR;		
	}

	if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
	    opserr << "WARNING invalid section Tube matTag" << endln;
	    return TCL_ERROR;		
	}

	if (Tcl_GetDouble (interp, argv[4], &D) != TCL_OK) {
	    opserr << "WARNING invalid D" << endln;
	    opserr << "Tube section: " << tag << endln;	    
	    return TCL_ERROR;
	}	

	if (Tcl_GetDouble (interp, argv[5], &t) != TCL_OK) {
	    opserr << "WARNING invalid t" << endln;
	    opserr << "Tube section: " << tag << endln;	    
	    return TCL_ERROR;
	}	
	
	if (Tcl_GetInt (interp, argv[6], &nfw) != TCL_OK) {
	    opserr << "WARNING invalid nfw" << endln;
	    opserr << "Tube section: " << tag << endln;	    	    
	    return TCL_ERROR;
	}	

	if (Tcl_GetInt (interp, argv[7], &nfr) != TCL_OK) {
	    opserr << "WARNING invalid nfr" << endln;
	    opserr << "Tube  section: " << tag << endln;	    	    
	    return TCL_ERROR;
	}	

	TubeSectionIntegration tubesect(D, t, nfw, nfr);

	int numFibers = tubesect.getNumFibers();

	if (argc > 8) {

	  double shape = 1.0;
	  if (argc > 9) {
	    if (Tcl_GetDouble(interp, argv[9], &shape) != TCL_OK) {
	      opserr << "WARNING invalid shape" << endln;
	      opserr << "WFSection2d section: " << tag << endln;	    	    
	      return TCL_ERROR;
	    }
	  }

	  NDMaterial *theSteel = OPS_getNDMaterial(matTag);
	
	  if (theSteel == 0) {
	    opserr << "WARNING ND material does not exist\n";
	    opserr << "material: " << matTag; 
	    opserr << "\nTube section: " << tag << endln;
	    return TCL_ERROR;
	  }
		  
	  NDMaterial **theMats = new NDMaterial *[numFibers];
	  
	  tubesect.arrangeFibers(theMats, theSteel);

	  // Parsing was successful, allocate the section
	  theSection = 0;
	  if (strcmp(argv[8],"-nd") == 0)
	    theSection = new NDFiberSection3d(tag, numFibers, theMats, tubesect, shape);
	  if (strcmp(argv[8],"-ndWarping") == 0)
	    theSection = new NDFiberSectionWarping2d(tag, numFibers, theMats, tubesect, shape);

	  delete [] theMats;	  
	}
	else {
	  UniaxialMaterial *theSteel = OPS_getUniaxialMaterial(matTag);
	
	  if (theSteel == 0) {
	    opserr << "WARNING uniaxial material does not exist\n";
	    opserr << "material: " << matTag; 
	    opserr << "\nTube section: " << tag << endln;
	    return TCL_ERROR;
	  }
	  
	  UniaxialMaterial **theMats = new UniaxialMaterial *[numFibers];
	  
	  tubesect.arrangeFibers(theMats, theSteel);
	  
	  // Parsing was successful, allocate the section
	  theSection = new FiberSection2d(tag, numFibers, theMats, tubesect);

	  delete [] theMats;
	}
    }    

    else if (strcmp(argv[1],"RCSection2d") == 0) {
	if (argc < 15) {
	    opserr << "WARNING insufficient arguments\n";
	    opserr << "Want: section RCSection2d tag? coreTag? coverTag? steelTag? d? b? cover? Atop? Abottom? Aside? nfcore? nfcover? nfs?" << endln;
	    return TCL_ERROR;
	}
	
	int tag, coreTag, coverTag, steelTag;
	double d, b, cover, Atop, Abottom, Aside;
	int nfcore, nfcover, nfs;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid section RCSection2d tag" << endln;
	    return TCL_ERROR;		
	}

	if (Tcl_GetInt(interp, argv[3], &coreTag) != TCL_OK) {
	    opserr << "WARNING invalid section RCSection2d coreTag" << endln;
	    return TCL_ERROR;		
	}

	if (Tcl_GetInt(interp, argv[4], &coverTag) != TCL_OK) {
	    opserr << "WARNING invalid section RCSection2d coverTag" << endln;
	    return TCL_ERROR;		
	}

	if (Tcl_GetInt(interp, argv[5], &steelTag) != TCL_OK) {
	    opserr << "WARNING invalid section RCSection2d steelTag" << endln;
	    return TCL_ERROR;		
	}

	if (Tcl_GetDouble (interp, argv[6], &d) != TCL_OK) {
	    opserr << "WARNING invalid d" << endln;
	    opserr << "RCSection2d section: " << tag << endln;	    
	    return TCL_ERROR;
	}	

	if (Tcl_GetDouble (interp, argv[7], &b) != TCL_OK) {
	    opserr << "WARNING invalid b" << endln;
	    opserr << "RCSection2d section: " << tag << endln;	    
	    return TCL_ERROR;
	}	
	
	if (Tcl_GetDouble (interp, argv[8], &cover) != TCL_OK) {
	    opserr << "WARNING invalid cover" << endln;
	    opserr << "RCSection2d section: " << tag << endln;	    	    
	    return TCL_ERROR;
	}	

	if (Tcl_GetDouble (interp, argv[9], &Atop) != TCL_OK) {
	    opserr << "WARNING invalid Atop" << endln;
	    opserr << "RCSection2d section: " << tag << endln;	    	    
	    return TCL_ERROR;
	}	

	if (Tcl_GetDouble (interp, argv[10], &Abottom) != TCL_OK) {
	    opserr << "WARNING invalid Abottom" << endln;
	    opserr << "RCSection2d section: " << tag << endln;	    	    
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble (interp, argv[11], &Aside) != TCL_OK) {
	    opserr << "WARNING invalid Aside" << endln;
	    opserr << "RCSection2d section: " << tag << endln;	    	    
	    return TCL_ERROR;
	}	

	if (Tcl_GetInt (interp, argv[12], &nfcore) != TCL_OK) {
	    opserr << "WARNING invalid nfcore" << endln;
	    opserr << "RCSection2d section: " << tag << endln;	    	    
	    return TCL_ERROR;
	}	

	if (Tcl_GetInt (interp, argv[13], &nfcover) != TCL_OK) {
	    opserr << "WARNING invalid nfcover" << endln;
	    opserr << "RCSection2d section: " << tag << endln;	    	    
	    return TCL_ERROR;
	}	

	if (Tcl_GetInt (interp, argv[14], &nfs) != TCL_OK) {
	    opserr << "WARNING invalid nfs" << endln;
	    opserr << "RCSection2d section: " << tag << endln;	    	    
	    return TCL_ERROR;
	}	

	UniaxialMaterial *theCore = OPS_getUniaxialMaterial(coreTag);
	
	if (theCore == 0) {
	    opserr << "WARNING uniaxial material does not exist\n";
	    opserr << "material: " << coreTag; 
	    opserr << "\nRCSection2d section: " << tag << endln;
	    return TCL_ERROR;
	}
	
	UniaxialMaterial *theCover = OPS_getUniaxialMaterial(coverTag);
	
	if (theCover == 0) {
	    opserr << "WARNING uniaxial material does not exist\4n";
	    opserr << "material: " << coverTag; 
	    opserr << "\nRCSection2d section: " << tag << endln;
	    return TCL_ERROR;
	}
	
	UniaxialMaterial *theSteel = OPS_getUniaxialMaterial(steelTag);

	if (theSteel == 0) {
	    opserr << "WARNING uniaxial material does not exist\n";
	    opserr << "material: " << steelTag; 
	    opserr << "\nRCSection2d section: " << tag << endln;
	    return TCL_ERROR;
	}
	
	RCSectionIntegration rcsect(d, b, Atop, Abottom, Aside, cover, nfcore, nfcover, nfs);

	int numFibers = rcsect.getNumFibers();

	UniaxialMaterial **theMats = new UniaxialMaterial *[numFibers];

	rcsect.arrangeFibers(theMats, theCore, theCover, theSteel);

	// Parsing was successful, allocate the section
	theSection = new FiberSection2d(tag, numFibers, theMats, rcsect);

	delete [] theMats;
    }

    else if (strcmp(argv[1],"RCCircularSection") == 0) {
        if (argc < 13) {
            opserr << "WARNING insufficient arguments\n";
            opserr << "Want: section RCCircularSection tag? coreTag? coverTag? steelTag? d? cover? As? NringsCore? NringsCover Nwedges? Nsteel?" << endln;
            return TCL_ERROR;
        }
        
        int tag, coreTag, coverTag, steelTag;
        double d, cover, As;
        int ncore, ncover, nwedge, nsteel;

        if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
            opserr << "WARNING invalid section RCCircularSection tag" << endln;
            return TCL_ERROR;           
        }

        if (Tcl_GetInt(interp, argv[3], &coreTag) != TCL_OK) {
            opserr << "WARNING invalid section RCCircularSection coreTag" << endln;
            return TCL_ERROR;           
        }

        if (Tcl_GetInt(interp, argv[4], &coverTag) != TCL_OK) {
            opserr << "WARNING invalid section RCCircularSection coverTag" << endln;
            return TCL_ERROR;           
        }

        if (Tcl_GetInt(interp, argv[5], &steelTag) != TCL_OK) {
            opserr << "WARNING invalid section RCCircularSection steelTag" << endln;
            return TCL_ERROR;           
        }

        if (Tcl_GetDouble (interp, argv[6], &d) != TCL_OK) {
            opserr << "WARNING invalid d" << endln;
            opserr << "RCCircularSection section: " << tag << endln;        
            return TCL_ERROR;
        }       

	       if (Tcl_GetDouble (interp, argv[7], &cover) != TCL_OK) {
            opserr << "WARNING invalid cover" << endln;
            opserr << "RCCircularSection section: " << tag << endln;                
            return TCL_ERROR;
        }       

        if (Tcl_GetDouble (interp, argv[8], &As) != TCL_OK) {
            opserr << "WARNING invalid As" << endln;
            opserr << "RCCircularSection section: " << tag << endln;                
            return TCL_ERROR;
        }       

        if (Tcl_GetInt (interp, argv[9], &ncore) != TCL_OK) {
            opserr << "WARNING invalid Ncore" << endln;
            opserr << "RCCircularSection section: " << tag << endln;                
            return TCL_ERROR;
        }

        if (Tcl_GetInt (interp, argv[10], &ncover) != TCL_OK) {
            opserr << "WARNING invalid ncover" << endln;
            opserr << "RCCircularSection section: " << tag << endln;                
            return TCL_ERROR;
        }       

        if (Tcl_GetInt (interp, argv[11], &nwedge) != TCL_OK) {
            opserr << "WARNING invalid nwedge" << endln;
            opserr << "RCCircularSection section: " << tag << endln;                
            return TCL_ERROR;
        }       

        if (Tcl_GetInt (interp, argv[12], &nsteel) != TCL_OK) {
            opserr << "WARNING invalid nsteel" << endln;
            opserr << "RCCircularSection section: " << tag << endln;                
            return TCL_ERROR;
        }       

        UniaxialMaterial *theCore = OPS_getUniaxialMaterial(coreTag);
        
        if (theCore == 0) {
            opserr << "WARNING uniaxial material does not exist\n";
            opserr << "material: " << coreTag; 
            opserr << "\nRCCircularSection section: " << tag << endln;
            return TCL_ERROR;
        }
        
        UniaxialMaterial *theCover = OPS_getUniaxialMaterial(coverTag);
        
        if (theCover == 0) {
            opserr << "WARNING uniaxial material does not exist\4n";
            opserr << "material: " << coverTag; 
            opserr << "\nRCCircularSection section: " << tag << endln;
            return TCL_ERROR;
        }
        
        UniaxialMaterial *theSteel = OPS_getUniaxialMaterial(steelTag);

        if (theSteel == 0) {
            opserr << "WARNING uniaxial material does not exist\n";
            opserr << "material: " << steelTag; 
            opserr << "\nRCCircularSection section: " << tag << endln;
            return TCL_ERROR;
        }
        
        RCCircularSectionIntegration rcsect(d, As, cover, ncore, ncover, nwedge, nsteel);

        int numFibers = rcsect.getNumFibers();

        UniaxialMaterial **theMats = new UniaxialMaterial *[numFibers];

        rcsect.arrangeFibers(theMats, theCore, theCover, theSteel);

        // Parsing was successful, allocate the section
        theSection = new FiberSection3d(tag, numFibers, theMats, rcsect);

        delete [] theMats;
    }

    else if (strcmp(argv[1],"RCTunnelSection") == 0) {
        if (argc < 15) {
            opserr << "WARNING insufficient arguments\n";
            opserr << "Want: section RCTunnelSection tag? concreteTag? steelTag? d? h? coverinner? coverouter? Asinner? Asouter? Nrings? Nwedges? Nbarsinner? Nbarsouter?" << endln;
            return TCL_ERROR;
        }
        
        int tag, concreteTag, steelTag;
        double d, h, coverinner, coverouter, Asinner, Asouter;
        int nring, nwedge, nbarinner, nbarouter;

        if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
            opserr << "WARNING invalid section RCTunnelSection tag" << endln;
            return TCL_ERROR;           
        }

        if (Tcl_GetInt(interp, argv[3], &concreteTag) != TCL_OK) {
            opserr << "WARNING invalid section RCTunnelSection concreteTag" << endln;
            return TCL_ERROR;           
        }

        if (Tcl_GetInt(interp, argv[4], &steelTag) != TCL_OK) {
            opserr << "WARNING invalid section RCTunnelSection steelTag" << endln;
            return TCL_ERROR;           
        }

        if (Tcl_GetDouble (interp, argv[5], &d) != TCL_OK) {
            opserr << "WARNING invalid d" << endln;
            opserr << "RCTunnelSection section: " << tag << endln;        
            return TCL_ERROR;
        }       

        if (Tcl_GetDouble (interp, argv[6], &h) != TCL_OK) {
            opserr << "WARNING invalid h" << endln;
            opserr << "RCTunnelSection section: " << tag << endln;        
            return TCL_ERROR;
        }       

	if (Tcl_GetDouble (interp, argv[7], &coverinner) != TCL_OK) {
            opserr << "WARNING invalid coverinner" << endln;
            opserr << "RCTunnelSection section: " << tag << endln;                
            return TCL_ERROR;
        }       

        if (Tcl_GetDouble (interp, argv[8], &coverouter) != TCL_OK) {
            opserr << "WARNING invalid coverouter" << endln;
            opserr << "RCTunnelSection section: " << tag << endln;                
            return TCL_ERROR;
        }       

	if (Tcl_GetDouble (interp, argv[9], &Asinner) != TCL_OK) {
            opserr << "WARNING invalid Asinner" << endln;
            opserr << "RCTunnelSection section: " << tag << endln;                
            return TCL_ERROR;
        }       

        if (Tcl_GetDouble (interp, argv[10], &Asouter) != TCL_OK) {
            opserr << "WARNING invalid Asouter" << endln;
            opserr << "RCTunnelSection section: " << tag << endln;                
            return TCL_ERROR;
        }       

        if (Tcl_GetInt (interp, argv[11], &nring) != TCL_OK) {
            opserr << "WARNING invalid Nring" << endln;
            opserr << "RCTunnelSection section: " << tag << endln;                
            return TCL_ERROR;
        }

        if (Tcl_GetInt (interp, argv[12], &nwedge) != TCL_OK) {
            opserr << "WARNING invalid nwedge" << endln;
            opserr << "RCTunnelSection section: " << tag << endln;                
            return TCL_ERROR;
        }       

        if (Tcl_GetInt (interp, argv[13], &nbarinner) != TCL_OK) {
            opserr << "WARNING invalid nbarinner" << endln;
            opserr << "RCTunnelSection section: " << tag << endln;                
            return TCL_ERROR;
        }       

        if (Tcl_GetInt (interp, argv[14], &nbarouter) != TCL_OK) {
            opserr << "WARNING invalid nbarouter" << endln;
            opserr << "RCTunnelSection section: " << tag << endln;                
            return TCL_ERROR;
        }       

        UniaxialMaterial *theConcrete = OPS_getUniaxialMaterial(concreteTag);
        
        if (theConcrete == 0) {
            opserr << "WARNING uniaxial material does not exist\n";
            opserr << "material: " << concreteTag; 
            opserr << "\nRCTunnelSection section: " << tag << endln;
            return TCL_ERROR;
        }
        
        UniaxialMaterial *theSteel = OPS_getUniaxialMaterial(steelTag);

        if (theSteel == 0) {
            opserr << "WARNING uniaxial material does not exist\n";
            opserr << "material: " << steelTag; 
            opserr << "\nRCTunnelSection section: " << tag << endln;
            return TCL_ERROR;
        }
        
        RCTunnelSectionIntegration rcsect(d, h, Asinner, Asouter, coverinner, coverouter,
					  nring, nwedge, nbarinner, nbarouter);

        int numFibers = rcsect.getNumFibers();

        UniaxialMaterial **theMats = new UniaxialMaterial *[numFibers];

        rcsect.arrangeFibers(theMats, theConcrete, theSteel);

        // Parsing was successful, allocate the section
        theSection = new FiberSection3d(tag, numFibers, theMats, rcsect);

        delete [] theMats;
    }
	
    else if (strcmp(argv[1],"RCTBeamSection2d") == 0 || strcmp(argv[1],"RCTBeamSectionUniMat2d") == 0) {
      if (argc < 20) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: section RCTBeamSection2d tag? coreTag? coverTag? steelTag? d? bw? beff? hf? Atop? Abottom? flcov? wcov? Nflcover? Nwcover? Nflcore? Nwcore? NsteelTop?  NsteelBottom?" << endln;
	return TCL_ERROR;
      }
      
      int tag, coreTag, coverTag, steelTag;
      double d, bw, beff, hf, Atop, Abottom, flcov, wcov;
      int  Nflcover, Nwcover, Nflcore, Nwcore, NsteelTop, NsteelBottom;
      
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING invalid section RCTBeamSection2d tag" << endln;
	return TCL_ERROR;		
      }
      
      if (Tcl_GetInt(interp, argv[3], &coreTag) != TCL_OK) {
	opserr << "WARNING invalid section RCTBeamSection2d coreTag" << endln;
	return TCL_ERROR;		
      }
      
      if (Tcl_GetInt(interp, argv[4], &coverTag) != TCL_OK) {
	opserr << "WARNING invalid section RCTBeamSection2d coverTag" << endln;
	return TCL_ERROR;		
      }
      
      if (Tcl_GetInt(interp, argv[5], &steelTag) != TCL_OK) {
	opserr << "WARNING invalid section RCTBeamSection2d steelTag" << endln;
	return TCL_ERROR;		
      }
      
      if (Tcl_GetDouble (interp, argv[6], &d) != TCL_OK) {
	opserr << "WARNING invalid d" << endln;
	opserr << "RCTBeamSection2d section: " << tag << endln;	    
	return TCL_ERROR;
      }	
      
      if (Tcl_GetDouble (interp, argv[7], &bw) != TCL_OK) {
	opserr << "WARNING invalid bw" << endln;
	opserr << "RCTBeamSection2d section: " << tag << endln;	    
	return TCL_ERROR;
      }	
      
      if (Tcl_GetDouble (interp, argv[8], &beff) != TCL_OK) {
	opserr << "WARNING invalid beff" << endln;
	opserr << "RCTBeamSection2d section: " << tag << endln;	    	    
	return TCL_ERROR;
      }	
      
      if (Tcl_GetDouble (interp, argv[9], &hf) != TCL_OK) {
	opserr << "WARNING invalid hf" << endln;
	opserr << "RCTBeamSection2d section: " << tag << endln;	    	    
	return TCL_ERROR;
      }	
      
      if (Tcl_GetDouble (interp, argv[10], &Atop) != TCL_OK) {
	opserr << "WARNING invalid Atop" << endln;
	opserr << "RCTBeamSection2d section: " << tag << endln;	    	    
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble (interp, argv[11], &Abottom) != TCL_OK) {
	opserr << "WARNING invalid Abottom" << endln;
	opserr << "RCTBeamSection2d section: " << tag << endln;	    	    
	return TCL_ERROR;
      }	
      
      if (Tcl_GetDouble (interp, argv[12], &flcov) != TCL_OK) {
	opserr << "WARNING invalid flcover" << endln;
	opserr << "RCTBeamSection2d section: " << tag << endln;	    	    
	return TCL_ERROR;
      }	
      
      if (Tcl_GetDouble (interp, argv[13], &wcov) != TCL_OK) {
	opserr << "WARNING invalid wcover" << endln;
	opserr << "RCTBeamSection2d section: " << tag << endln;	    	    
	return TCL_ERROR;
      }	
      
      if (Tcl_GetInt (interp, argv[14], &Nflcover) != TCL_OK) {
	opserr << "WARNING invalid Nflcover" << endln;
	opserr << "RCTBeamSection2d section: " << tag << endln;	    	    
	return TCL_ERROR;
      }	
      
      if (Tcl_GetInt (interp, argv[15], &Nwcover) != TCL_OK) {
	opserr << "WARNING invalid Nwcover" << endln;
	opserr << "RCTBeamSection2d section: " << tag << endln;	    	    
	return TCL_ERROR;
      }	
      
      if (Tcl_GetInt (interp, argv[16], &Nflcore) != TCL_OK) {
	opserr << "WARNING invalid Nflcore" << endln;
	opserr << "RCTBeamSection2d section: " << tag << endln;	    	    
	return TCL_ERROR;
      }	
      
      if (Tcl_GetInt (interp, argv[17], &Nwcore) != TCL_OK) {
	opserr << "WARNING invalid Nwcore" << endln;
	opserr << "RCTBeamSection2d section: " << tag << endln;	    	    
	return TCL_ERROR;
      }
      
      if (Tcl_GetInt (interp, argv[18], &NsteelTop) != TCL_OK) {
	opserr << "WARNING invalid NsteelTop" << endln;
	opserr << "RCTBeamSection2d section: " << tag << endln;	    	    
	return TCL_ERROR;
      }	

      if (Tcl_GetInt (interp, argv[19], &NsteelBottom) != TCL_OK) {
	opserr << "WARNING invalid NsteelBottom" << endln;
	opserr << "RCTBeamSection2d section: " << tag << endln;	    	    
	return TCL_ERROR;
      }	
            
      UniaxialMaterial *theSteel = OPS_getUniaxialMaterial(steelTag);
      if (theSteel == 0) {
	opserr << "WARNING uniaxial material does not exist\n";
	opserr << "material: " << steelTag; 
	opserr << "\nRCTBeamSection2d section: " << tag << endln;
	return TCL_ERROR;
      }

      RCTBeamSectionIntegration
	rctbeamsect(d, bw, beff, hf, Atop, Abottom, flcov, wcov,
		    Nflcover, Nwcover, Nflcore, Nwcore,
		    NsteelTop, NsteelBottom);

      if (strcmp(argv[1],"RCTBeamSectionUniMat2d") == 0) {
	UniaxialMaterial *theCore = OPS_getUniaxialMaterial(coreTag);
	if (theCore == 0) {
	  opserr << "WARNING uniaxial material does not exist\n";
	  opserr << "material: " << coreTag; 
	  opserr << "\nRCTBeamSection2d section: " << tag << endln;
	  return TCL_ERROR;
	}

	UniaxialMaterial *theCover = OPS_getUniaxialMaterial(coverTag);
	if (theCover == 0) {
	  opserr << "WARNING uniaxial material does not exist\n";
	  opserr << "material: " << coreTag; 
	  opserr << "\nRCTBeamSection2d section: " << tag << endln;
	  return TCL_ERROR;
	}

	int numFibers = rctbeamsect.getNumFibers();

	UniaxialMaterial **theUniMat = new UniaxialMaterial *[numFibers];

	rctbeamsect.arrangeFibers(theUniMat, theCore, theCover, theSteel);

	theSection = new FiberSection2d(tag, numFibers, theUniMat, rctbeamsect);

	delete [] theUniMat;
      } 
      else {
	NDMaterial *theCore = OPS_getNDMaterial(coreTag);
	if (theCore == 0) {
	  opserr << "WARNING uniaxial material does not exist\n";
	  opserr << "material: " << coreTag; 
	  opserr << "\nRCTBeamSection2d section: " << tag << endln;
	  return TCL_ERROR;
	}

	NDMaterial *theCover = OPS_getNDMaterial(coverTag);      
	if (theCover == 0) {
	  opserr << "WARNING uniaxial material does not exist\4n";
	  opserr << "material: " << coverTag; 
	  opserr << "\nRCTBeamSection2d section: " << tag << endln;
	  return TCL_ERROR;
	}

	int numCFibers = rctbeamsect.getNumFibers(concrete);
	int numSFibers = rctbeamsect.getNumFibers(steel);

	NDMaterial **theNDMat = new NDMaterial *[numCFibers];
	UniaxialMaterial **theUniMat = new UniaxialMaterial *[numSFibers];

	rctbeamsect.arrangeFibers(theUniMat, theNDMat, theCore, theCover, theSteel);

	//theSection = new McftSection2dfiber(tag, theNDMat, theUniMat, rctbeamsect);



	RCTBeamSectionIntegration
	  steel(d, bw, beff, hf, Atop, Abottom, flcov, wcov,
		0, 0, 0, 0,
		NsteelTop, NsteelBottom);
	steel.arrangeFibers(theUniMat, theNDMat, 0, 0, theSteel);
	FiberSection2d steelSec(0, numSFibers, theUniMat, steel);
	
	RCTBeamSectionIntegration
	  concrete(d, bw, beff, hf, Atop, Abottom, flcov, wcov,
		   Nflcover, Nwcover, Nflcore, Nwcore,
		   0, 0);
	concrete.arrangeFibers(theUniMat, theNDMat, theCore, theCover, 0);
	NDFiberSection2d concSec(0, numCFibers, theNDMat, concrete);
	
	

	SectionForceDeformation *theSections[2];
	theSections[1] = &steelSec;
	theSections[0] = &concSec;
	theSection = new ParallelSection(tag, 2, theSections);
	


	delete [] theNDMat;
	delete [] theUniMat;
      }
    }

    else if (strcmp(argv[1],"Parallel") == 0) {
	if (argc < 4) {
	    opserr << "WARNING insufficient arguments\n";
	    opserr << "Want: section Parallel tag? tag1? tag2? ..." << endln;
	    return TCL_ERROR;
	}
 
	int tag;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid section Parallel tag" << endln;
	    return TCL_ERROR;		
	}

	int numMaterials = argc-3;
	
	if (numMaterials == 0) {
	    opserr << "WARNING no component section(s) provided\n";
	    opserr << "section Parallel: " << tag << endln;
	    return TCL_ERROR;
	}
    
	// Create an array to hold pointers to component materials
	SectionForceDeformation **theMats = new SectionForceDeformation *[numMaterials];
	
	// For each material get the tag and ensure it exists in model already
	for (int i = 0; i < numMaterials; i++) {
	    int tagI;
	    if (Tcl_GetInt(interp, argv[i+3], &tagI) != TCL_OK) {
		opserr << "WARNING invalid component tag\n";
		opserr << "section Parallel: " << tag << endln;
		return TCL_ERROR;
	    }
	    
	    SectionForceDeformation *theMat = theTclBuilder->getSection(tagI);
	    
	    if (theMat == 0) {
		opserr << "WARNING component section does not exist\n";
		opserr << "Component section: " << argv[i+3]; 
		opserr << "\tsection Parallel: " << tag << endln;
		delete [] theMats;
		return TCL_ERROR;
	    }
	    else
		theMats[i] = theMat;
	}	
	
	// Parsing was successful, allocate the material
	theSection = new ParallelSection(tag, numMaterials, theMats);
	
	// Deallocate the temporary pointers
	delete [] theMats;
    }

    else if (strcmp(argv[1],"AddDeformation") == 0 || strcmp(argv[1],"Aggregator") == 0) {
	if (argc < 5) {
	    opserr << "WARNING insufficient arguments\n";
	    opserr << "Want: section Aggregator tag? uniTag1? code1? ... <-section secTag?>" << endln;
	    return TCL_ERROR;
	}
	    
	int tag;
	int secTag;
	SectionForceDeformation *theSec = 0;
	    
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid Aggregator tag" << endln;
	    return TCL_ERROR;		
	}

	int nArgs = argc-3;
	
	for (int ii = 5; ii < argc; ii++) {
	    if (strcmp(argv[ii],"-section") == 0 && ++ii < argc) {
		if (Tcl_GetInt(interp, argv[ii], &secTag) != TCL_OK) {
		    opserr << "WARNING invalid Aggregator tag" << endln;
		    return TCL_ERROR;		
		}
		
		theSec = theTclBuilder->getSection(secTag);
		
		if (theSec == 0) {
		    opserr << "WARNING section does not exist\n";
		    opserr << "section: " << secTag; 
		    opserr << "\nsection Aggregator: " << tag << endln;
		    return TCL_ERROR;
		}
		
		nArgs -= 2;
	    }
	}
	
	int nMats = nArgs / 2;
	
	if (nArgs%2 != 0) {
	    opserr << "WARNING improper number of arguments for Aggregator" << endln;
	    return TCL_ERROR;
	}
	
	UniaxialMaterial **theMats = 0;
	ID codes(nMats);
	
	theMats = new UniaxialMaterial *[nMats];
	
	if (theMats == 0) {
	    opserr << "TclModelBuilderSection (Aggregator) -- unable to create uniaxial array" << endln;
	    return TCL_ERROR;
	}	
	
	int tagI;
	int i, j;
	
	for (i = 3, j = 0; j < nMats; i++, j++) {
	    if (Tcl_GetInt(interp, argv[i], &tagI) != TCL_OK) {
		opserr << "WARNING invalid Aggregator matTag" << endln;
		return TCL_ERROR;		
	    }
	    
	    theMats[j] = OPS_getUniaxialMaterial(tagI);
	    
	    if (theMats[j] == 0) {
		opserr << "WARNING uniaxial material does not exist\n";
		opserr << "uniaxial material: " << tagI; 
		opserr << "\nsection Aggregator: " << tag << endln;
		return TCL_ERROR;
	    }
	    
	    i++;
	    
	    if (strcmp(argv[i],"Mz") == 0)
		codes(j) = SECTION_RESPONSE_MZ;
	    else if (strcmp(argv[i],"P") == 0)
		codes(j) = SECTION_RESPONSE_P;
	    else if (strcmp(argv[i],"Vy") == 0)
		codes(j) = SECTION_RESPONSE_VY;
	    else if (strcmp(argv[i],"My") == 0)
		codes(j) = SECTION_RESPONSE_MY;
	    else if (strcmp(argv[i],"Vz") == 0)
		codes(j) = SECTION_RESPONSE_VZ;
	    else if (strcmp(argv[i],"T") == 0)
		codes(j) = SECTION_RESPONSE_T;
	    else {
		opserr << "WARNING invalid code" << endln;
		opserr << "\nsection Aggregator: " << tag << endln;
		return TCL_ERROR;		
	    }
	}
	
	if (theSec)
	    theSection = new SectionAggregator (tag, *theSec, nMats, theMats, codes);
	else
	    theSection = new SectionAggregator (tag, nMats, theMats, codes);
	
	delete [] theMats;
    }		
    
    else if (strcmp(argv[1],"Fiber") == 0 || 
	     strcmp(argv[1],"fiberSec") == 0 ||
	     strcmp(argv[1],"NDFiberWarping") == 0 ||
	     strcmp(argv[1],"NDFiber") == 0)
	return TclCommand_addFiberSection (clientData, interp, argc, argv,
					   theTclBuilder);

    //--- Adding Thermo-mechanical Sections:[BEGIN]   by UoE OpenSees Group ---//  
    else if (strcmp(argv[1],"FiberThermal") == 0 || strcmp(argv[1],"fiberSecThermal") == 0)
      return TclCommand_addFiberSectionThermal (clientData, interp, argc, argv,
						theTclBuilder);

   else if (strcmp(argv[1],"FiberInt") == 0)
	return TclCommand_addFiberIntSection (clientData, interp, argc, argv,
						theTclBuilder);

    else if (strcmp(argv[1],"UCFiber") == 0)
	return TclCommand_addUCFiberSection (clientData, interp, argc, argv,
						  theTclBuilder);

    else if (strcmp(argv[1],"ElasticPlateSection") == 0) {
	if (argc < 5) {
	    opserr << "WARNING insufficient arguments\n";
	    opserr << "Want: section ElasticPlateSection tag? E? nu? h? " << endln;
	    return TCL_ERROR;
	}
	
	int tag;
	double E, nu, h;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid section ElasticPlateSection tag" << endln;
	    return TCL_ERROR;		
	}

	if (Tcl_GetDouble (interp, argv[3], &E) != TCL_OK) {
	    opserr << "WARNING invalid E" << endln;
	    opserr << "ElasticPlateSection section: " << tag << endln;	    
	    return TCL_ERROR;
	}	

	if (Tcl_GetDouble (interp, argv[4], &nu) != TCL_OK) {
	    opserr << "WARNING invalid nu" << endln;
	    opserr << "ElasticPlateSection section: " << tag << endln;	    
	    return TCL_ERROR;
	}	
	
	if (Tcl_GetDouble (interp, argv[5], &h) != TCL_OK) {
	    opserr << "WARNING invalid h" << endln;
	    opserr << "ElasticPlateSection section: " << tag << endln;	    	    
	    return TCL_ERROR;
	}	

	theSection = new ElasticPlateSection (tag, E, nu, h);
    }	

    else if (strcmp(argv[1],"ElasticMembranePlateSection") == 0) {
	if (argc < 5) {
	    opserr << "WARNING insufficient arguments\n";
	    opserr << "Want: section ElasticMembranePlateSection tag? E? nu? h? <rho?>" << endln;
	    return TCL_ERROR;
	}
	
	int tag;
	double E, nu, h;
	double rho = 0.0;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid section ElasticMembranePlateSection tag" << endln;
	    return TCL_ERROR;		
	}

	if (Tcl_GetDouble (interp, argv[3], &E) != TCL_OK) {
	    opserr << "WARNING invalid E" << endln;
	    opserr << "ElasticMembranePlateSection section: " << tag << endln;	    
	    return TCL_ERROR;
	}	

	if (Tcl_GetDouble (interp, argv[4], &nu) != TCL_OK) {
	    opserr << "WARNING invalid nu" << endln;
	    opserr << "ElasticMembranePlateSection section: " << tag << endln;	    
	    return TCL_ERROR;
	}	
	
	if (Tcl_GetDouble (interp, argv[5], &h) != TCL_OK) {
	    opserr << "WARNING invalid h" << endln;
	    opserr << "ElasticMembranePlateSection section: " << tag << endln;	    	    
	    return TCL_ERROR;
	}	

	if (argc > 6 && Tcl_GetDouble (interp, argv[6], &rho) != TCL_OK) {
	    opserr << "WARNING invalid rho" << endln;
	    opserr << "ElasticMembranePlateSection section: " << tag << endln;	    	    
	    return TCL_ERROR;
	}

	theSection = new ElasticMembranePlateSection (tag, E, nu, h, rho);
    }	

    else if (strcmp(argv[1],"PlateFiber") == 0) {
	if (argc < 5) {
	    opserr << "WARNING insufficient arguments\n";
	    opserr << "Want: section PlateFiber tag? matTag? h? " << endln;
	    return TCL_ERROR;
	}
	
	int tag, matTag;
	double  h;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid section PlateFiber tag" << endln;
	    return TCL_ERROR;		
	}

	if (Tcl_GetInt (interp, argv[3], &matTag) != TCL_OK) {
	    opserr << "WARNING invalid matTag" << endln;
	    opserr << "PlateFiber section: " << matTag << endln;	    	    
	    return TCL_ERROR;
	}	

	if (Tcl_GetDouble (interp, argv[4], &h) != TCL_OK) {
	    opserr << "WARNING invalid h" << endln;
	    opserr << "PlateFiber section: " << tag << endln;	    	    
	    return TCL_ERROR;
	}	

	NDMaterial *theMaterial = OPS_getNDMaterial(matTag);
	if (theMaterial == 0) {
	    opserr << "WARNING nD material does not exist\n";
	    opserr << "nD material: " << matTag; 
	    opserr << "\nPlateFiber section: " << tag << endln;
	    return TCL_ERROR;
	}

	theSection = new MembranePlateFiberSection( tag, h, *theMaterial );
    }	


    //start Yuli Huang & Xinzheng Lu LayeredShellFiberSection
    else if (strcmp(argv[1],"LayeredShell") == 0) {
      if (argc < 6) {
	opserr << "WARNING insufficient arguments" << endln;
	opserr << "Want: section LayeredShell tag? nLayers? matTag1? h1? ... matTagn? hn? " << endln;
	return TCL_ERROR;
	}
      
      int tag, nLayers, matTag;
      double h, *thickness;
      NDMaterial **theMats;
      
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING invalid section LayeredShell tag" << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetInt (interp, argv[3], &nLayers) != TCL_OK) {
	opserr << "WARNING invalid nLayers" << endln;
	opserr << "LayeredShell section: " << tag << endln;	    	    
	return TCL_ERROR;
      }
	
      if (nLayers < 3) {
	opserr << "ERROR number of layers must be larger than 2" << endln;
	opserr << "LayeredShell section: " << tag << endln;	    	    
	return TCL_ERROR;
      }
      
      theMats   = new NDMaterial*[nLayers];
      thickness = new double[nLayers];
      
      for (int iLayer = 0; iLayer < nLayers; iLayer++) {
	if (Tcl_GetInt (interp, argv[4+2*iLayer], &matTag) != TCL_OK) {
	  opserr << "WARNING invalid matTag" << endln;
	  opserr << "LayeredShell section: " << tag << endln;
	  return TCL_ERROR;
	}
	
	theMats[iLayer] = OPS_getNDMaterial(matTag);
	if (theMats[iLayer] == 0) {
	  opserr << "WARNING nD material does not exist" << endln;;
	  opserr << "nD material: " << matTag; 
	  opserr << "LayeredShell section: " << tag << endln;
	  return TCL_ERROR;
	}
	
	if (Tcl_GetDouble (interp, argv[5+2*iLayer], &h) != TCL_OK) {
	  opserr << "WARNING invalid h" << endln;
	  opserr << "LayeredShell section: " << tag << endln;	    	    
	  return TCL_ERROR;
	}
	
	if (h < 0) {
	  opserr << "WARNING invalid h" << endln;
	  opserr << "PlateFiber section: " << tag << endln;	    	    
	  return TCL_ERROR;
	}
	
	thickness[iLayer] = h;
      }
      
      theSection = new LayeredShellFiberSection(tag, nLayers, thickness, theMats);
      if (thickness != 0) delete thickness;
      if (theMats != 0) delete [] theMats;
    }
    //end Yuli Huang & Xinzheng Lu LayeredShellFiberSection

	//-----Thermo-mechanical shell sections added by L.Jiang [SIF] 
	else if (strcmp(argv[1], "PlateFiberThermal") == 0) {
		if (argc < 5) {
			opserr << "WARNING insufficient arguments\n";
			opserr << "Want: section PlateFiberThermal tag? matTag? h? " << endln;
			return TCL_ERROR;
		}

		int tag, matTag;
		double  h;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
			opserr << "WARNING invalid section PlateFiberThermal tag" << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
			opserr << "WARNING invalid matTag" << endln;
			opserr << "PlateFiberThermal section: " << matTag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[4], &h) != TCL_OK) {
			opserr << "WARNING invalid h" << endln;
			opserr << "PlateFiberThermal section: " << tag << endln;
			return TCL_ERROR;
		}

		NDMaterial *theMaterial = OPS_getNDMaterial(matTag);
		if (theMaterial == 0) {
			opserr << "WARNING nD material does not exist\n";
			opserr << "nD material: " << matTag;
			opserr << "\nPlateFiberThermal section: " << tag << endln;
			return TCL_ERROR;
		}

		theSection = new MembranePlateFiberSectionThermal(tag, h, *theMaterial);
	}

	
	// LayeredShellFiberSectionThermal based on the LayeredShellFiberSectionThermal by Yuli Huang & Xinzheng Lu
	else if (strcmp(argv[1], "LayeredShellThermal") == 0) {
		if (argc < 6) {
			opserr << "WARNING insufficient arguments" << endln;
			opserr << "Want: section LayeredShellThermal tag? nLayers? matTag1? h1? ... matTagn? hn? " << endln;
			return TCL_ERROR;
		}

		int tag, nLayers, matTag;
		double h, *thickness;
		NDMaterial **theMats;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
			opserr << "WARNING invalid section LayeredShellThermal tag" << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[3], &nLayers) != TCL_OK) {
			opserr << "WARNING invalid nLayers" << endln;
			opserr << "LayeredShellThermal section: " << tag << endln;
			return TCL_ERROR;
		}

		if (nLayers < 3) {
			opserr << "ERROR number of layers must be larger than 2" << endln;
			opserr << "LayeredShellThermal section: " << tag << endln;
			return TCL_ERROR;
		}

		theMats = new NDMaterial*[nLayers];
		thickness = new double[nLayers];

		for (int iLayer = 0; iLayer < nLayers; iLayer++) {
			if (Tcl_GetInt(interp, argv[4 + 2 * iLayer], &matTag) != TCL_OK) {
				opserr << "WARNING invalid matTag" << endln;
				opserr << "LayeredShellThermal section: " << tag << endln;
				return TCL_ERROR;
			}

			theMats[iLayer] = OPS_getNDMaterial(matTag);
			if (theMats[iLayer] == 0) {
				opserr << "WARNING nD material does not exist" << endln;;
				opserr << "nD material: " << matTag;
				opserr << "LayeredShellThermal section: " << tag << endln;
				return TCL_ERROR;
			}

			if (Tcl_GetDouble(interp, argv[5 + 2 * iLayer], &h) != TCL_OK) {
				opserr << "WARNING invalid h" << endln;
				opserr << "LayeredShellThermal section: " << tag << endln;
				return TCL_ERROR;
			}

			if (h < 0) {
				opserr << "WARNING invalid h" << endln;
				opserr << "LayeredShellThermal section: " << tag << endln;
				return TCL_ERROR;
			}

			thickness[iLayer] = h;
		}

		theSection = new LayeredShellFiberSectionThermal(tag, nLayers, thickness, theMats);
		if (thickness != 0) delete thickness;
		if (theMats != 0) delete[] theMats;
	}
	//end L.Jiang [SIF] added based on LayeredShellFiberSectionThermal section created by Yuli Huang & Xinzheng Lu ----
    
    else if (strcmp(argv[1],"Bidirectional") == 0) {
	if (argc < 7) {
	    opserr << "WARNING insufficient arguments\n";
	    opserr << "Want: section Bidirectional tag? E? sigY? Hiso? Hkin?" << endln;
	    return TCL_ERROR;
	}    

	int tag;
	double E, sigY, Hi, Hk;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid Bidirectional tag" << endln;
	    return TCL_ERROR;		
	}

	if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
	    opserr << "WARNING invalid E\n";
	    opserr << "section Bidirectional: " << tag << endln;
	    return TCL_ERROR;	
	}

	if (Tcl_GetDouble(interp, argv[4], &sigY) != TCL_OK) {
	    opserr << "WARNING invalid sigY\n";
	    opserr << "section Bidirectional: " << tag << endln;
	    return TCL_ERROR;	
	}

	if (Tcl_GetDouble(interp, argv[5], &Hi) != TCL_OK) {
	    opserr << "WARNING invalid Hiso\n";
	    opserr << "section Bidirectional: " << tag << endln;
	    return TCL_ERROR;	
	}

	if (Tcl_GetDouble(interp, argv[6], &Hk) != TCL_OK) {
	    opserr << "WARNING invalid Hkin\n";
	    opserr << "section Bidirectional: " << tag << endln;
	    return TCL_ERROR;	
	}

	if (argc > 8) {
	  int code1, code2;
	  if (strcmp(argv[7],"Mz") == 0)
	    code1 = SECTION_RESPONSE_MZ;
	  else if (strcmp(argv[7],"P") == 0)
	    code1 = SECTION_RESPONSE_P;
	  else if (strcmp(argv[7],"Vy") == 0)
	    code1 = SECTION_RESPONSE_VY;
	  else if (strcmp(argv[7],"My") == 0)
	    code1 = SECTION_RESPONSE_MY;
	  else if (strcmp(argv[7],"Vz") == 0)
	    code1 = SECTION_RESPONSE_VZ;
	  else if (strcmp(argv[7],"T") == 0)
	    code1 = SECTION_RESPONSE_T;
	  else {
	    opserr << "WARNING invalid code 1 " << argv[7] << endln;
	    opserr << "section Bidirectional: " << tag << endln;
	    return TCL_ERROR;		
	  }

	  if (strcmp(argv[8],"Mz") == 0)
	    code2 = SECTION_RESPONSE_MZ;
	  else if (strcmp(argv[8],"P") == 0)
	    code2 = SECTION_RESPONSE_P;
	  else if (strcmp(argv[8],"Vy") == 0)
	    code2 = SECTION_RESPONSE_VY;
	  else if (strcmp(argv[8],"My") == 0)
	    code2 = SECTION_RESPONSE_MY;
	  else if (strcmp(argv[8],"Vz") == 0)
	    code2 = SECTION_RESPONSE_VZ;
	  else if (strcmp(argv[8],"T") == 0)
	    code2 = SECTION_RESPONSE_T;
	  else {
	    opserr << "WARNING invalid code 2 " << argv[8] << endln;
	    opserr << "section Bidirectional: " << tag << endln;
	    return TCL_ERROR;		
	  }
	  theSection = new Bidirectional(tag, E, sigY, Hi, Hk, code1, code2);
	}
	else 
	  theSection = new Bidirectional(tag, E, sigY, Hi, Hk);

	}

        else if (strcmp(argv[1],"Elliptical") == 0 || strcmp(argv[1],"Elliptical2") == 0) {
	if (argc < 10) {
	    opserr << "WARNING insufficient arguments\n";
	    opserr << "Want: section Elliptical tag? E1? E2? sigY1? sigY2? Hiso? Hkin1? Hkin2? <code1? code2?>" << endln;
	    return TCL_ERROR;
	}    

	int tag;
	double E1, E2, sigY1, sigY2, Hi, Hk1, Hk2;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid Elliptical tag" << endln;
	    return TCL_ERROR;		
	}

	if (Tcl_GetDouble(interp, argv[3], &E1) != TCL_OK) {
	    opserr << "WARNING invalid E1\n";
	    opserr << "section Elliptical: " << tag << endln;
	    return TCL_ERROR;	
	}

	if (Tcl_GetDouble(interp, argv[4], &E2) != TCL_OK) {
	    opserr << "WARNING invalid E2\n";
	    opserr << "section Elliptical: " << tag << endln;
	    return TCL_ERROR;	
	}

	if (Tcl_GetDouble(interp, argv[5], &sigY1) != TCL_OK) {
	    opserr << "WARNING invalid sigY1\n";
	    opserr << "section Elliptical: " << tag << endln;
	    return TCL_ERROR;	
	}

	if (Tcl_GetDouble(interp, argv[6], &sigY2) != TCL_OK) {
	    opserr << "WARNING invalid sigY2\n";
	    opserr << "section Elliptical: " << tag << endln;
	    return TCL_ERROR;	
	}

	if (Tcl_GetDouble(interp, argv[7], &Hi) != TCL_OK) {
	    opserr << "WARNING invalid Hiso\n";
	    opserr << "section Elliptical: " << tag << endln;
	    return TCL_ERROR;	
	}

	if (Tcl_GetDouble(interp, argv[8], &Hk1) != TCL_OK) {
	    opserr << "WARNING invalid Hkin1\n";
	    opserr << "section Elliptical: " << tag << endln;
	    return TCL_ERROR;	
	}

	if (Tcl_GetDouble(interp, argv[9], &Hk2) != TCL_OK) {
	    opserr << "WARNING invalid Hkin2\n";
	    opserr << "section Elliptical: " << tag << endln;
	    return TCL_ERROR;	
	}

	if (argc > 11) {
	  int code1, code2;
	  if (strcmp(argv[10],"Mz") == 0)
	    code1 = SECTION_RESPONSE_MZ;
	  else if (strcmp(argv[10],"P") == 0)
	    code1 = SECTION_RESPONSE_P;
	  else if (strcmp(argv[10],"Vy") == 0)
	    code1 = SECTION_RESPONSE_VY;
	  else if (strcmp(argv[10],"My") == 0)
	    code1 = SECTION_RESPONSE_MY;
	  else if (strcmp(argv[10],"Vz") == 0)
	    code1 = SECTION_RESPONSE_VZ;
	  else if (strcmp(argv[10],"T") == 0)
	    code1 = SECTION_RESPONSE_T;
	  else {
	    opserr << "WARNING invalid code 1 " << argv[10] << endln;
	    opserr << "section Elliptical: " << tag << endln;
	    return TCL_ERROR;		
	  }

	  if (strcmp(argv[11],"Mz") == 0)
	    code2 = SECTION_RESPONSE_MZ;
	  else if (strcmp(argv[11],"P") == 0)
	    code2 = SECTION_RESPONSE_P;
	  else if (strcmp(argv[11],"Vy") == 0)
	    code2 = SECTION_RESPONSE_VY;
	  else if (strcmp(argv[11],"My") == 0)
	    code2 = SECTION_RESPONSE_MY;
	  else if (strcmp(argv[11],"Vz") == 0)
	    code2 = SECTION_RESPONSE_VZ;
	  else if (strcmp(argv[11],"T") == 0)
	    code2 = SECTION_RESPONSE_T;
	  else {
	    opserr << "WARNING invalid code 2 " << argv[11] << endln;
	    opserr << "section Elliptical: " << tag << endln;
	    return TCL_ERROR;		
	  }
	  if (strcmp(argv[1],"Elliptical") == 0)
	    //theSection = new Elliptical(tag, E1, E2, sigY1, sigY2, Hi, Hk1, Hk2, code1, code2);
	    theSection = new Elliptical2(tag, E1, E2, sigY1, sigY2, Hi, Hk1, Hk2, code1, code2);
	  else
	    theSection = new Elliptical2(tag, E1, E2, sigY1, sigY2, Hi, Hk1, Hk2, code1, code2);
	}
	else {
	  if (strcmp(argv[1],"Elliptical") == 0)
	    //theSection = new Elliptical(tag, E1, E2, sigY1, sigY2, Hi, Hk1, Hk2);
	    theSection = new Elliptical2(tag, E1, E2, sigY1, sigY2, Hi, Hk1, Hk2);
	  else 
	    theSection = new Elliptical2(tag, E1, E2, sigY1, sigY2, Hi, Hk1, Hk2);
	}
	}

        else if (strcmp(argv[1],"Iso2spring") == 0) {
	  if (argc < 10) {
	    opserr << "WARNING insufficient arguments\n";
	    opserr << "Want: section Iso2spring tag? tol? k1? Fy? k2? kv? hb? Pe? <Po?>" << endln;
	    return TCL_ERROR;
	  }    
	  
	  int tag;
	  double tol, k1, Fy, kb, kvo, hb, Pe;
      double Po = 0.0;

	  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid Iso2spring tag" << endln;
	    return TCL_ERROR;		
	  }

	  if (Tcl_GetDouble(interp, argv[3], &tol) != TCL_OK) {
	    opserr << "WARNING invalid tol\n";
	    opserr << "section Iso2spring: " << tag << endln;
	    return TCL_ERROR;
	  }
	  
	  if (Tcl_GetDouble(interp, argv[4], &k1) != TCL_OK) {
	    opserr << "WARNING invalid k1\n";
	    opserr << "section Iso2spring: " << tag << endln;
	    return TCL_ERROR;	
	  }
	  
	  if (Tcl_GetDouble(interp, argv[5], &Fy) != TCL_OK) {
	    opserr << "WARNING invalid Fy\n";
	    opserr << "section Iso2spring: " << tag << endln;
	    return TCL_ERROR;	
	  }
	  
	  if (Tcl_GetDouble(interp, argv[6], &kb) != TCL_OK) {
	    opserr << "WARNING invalid k2\n";
	    opserr << "section Iso2spring: " << tag << endln;
	    return TCL_ERROR;	
	  }
	  
	  if (Tcl_GetDouble(interp, argv[7], &kvo) != TCL_OK) {
	    opserr << "WARNING invalid kv\n";
	    opserr << "section Iso2spring: " << tag << endln;
	    return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[8], &hb) != TCL_OK) {
	    opserr << "WARNING invalid hb\n";
	    opserr << "section Iso2spring: " << tag << endln;
	    return TCL_ERROR;	
	}
	  
	  if (Tcl_GetDouble(interp, argv[9], &Pe) != TCL_OK) {
	    opserr << "WARNING invalid Pe\n";
	    opserr << "section Iso2spring: " << tag << endln;
	    return TCL_ERROR;	
	  }
	  if (argc > 10) {
	    if (Tcl_GetDouble(interp, argv[10], &Po) != TCL_OK) {
	      opserr << "WARNING invalid Po\n";
	      opserr << "section Iso2spring: " << tag << endln;
	      return TCL_ERROR;	
	    }
	  }	    
	  
	  theSection = new Isolator2spring(tag, tol, k1, Fy, kb, kvo, hb, Pe, Po);
	}
    		
    else {
      theSection = TclModelBuilderYS_SectionCommand(clientData, interp, argc, 
						    argv, theTclBuilder);
    }
    
    // Ensure we have created the Material, out of memory if got here and no section
    if (theSection == 0) {
      opserr << "WARNING could not create section " << argv[1] << endln;
      return TCL_ERROR;
    }
    
    // Now add the material to the modelBuilder
    //    if (theTclBuilder->addSection(*theSection) < 0) {
    if (OPS_addSectionForceDeformation(theSection) != true) {
	opserr << "WARNING could not add section to the domain\n";
	opserr << *theSection << endln;
	delete theSection; // invoke the material objects destructor, otherwise mem leak
	return TCL_ERROR;
    }
    
    return TCL_OK;
}

static int currentSectionTag = 0;
static bool currentSectionIsND = false;
static bool currentSectionIsWarping = false;
    
int
buildSection(Tcl_Interp *interp, TclModelBuilder *theTclModelBuilder,
	     int secTag, bool isTorsion, double GJ);

int
buildSectionInt(Tcl_Interp *interp, TclModelBuilder *theTclModelBuilder,
		int secTag, bool isTorsion, double GJ, 
		int NStrip1, double t1, int NStrip2, double t2, int NStrip3, double t3);

int
TclCommand_addFiberSection (ClientData clientData, Tcl_Interp *interp, int argc,
			    TCL_Char **argv, TclModelBuilder *theTclModelBuilder)
{
    int secTag;
    int maxNumPatches = 30; 
    int maxNumReinfLayers = 30;
    
    if (argc < 4) 
	return TCL_ERROR;
    
    if (Tcl_GetInt(interp, argv[2], &secTag) != TCL_OK) {
      opserr <<  "WARNING bad command - want: \nsection fiberSec secTag { \n\tpatch <patch arguments> \n\tlayer <layer arguments> \n}\n";
	return TCL_ERROR;
    }
    
    currentSectionTag = secTag;
    currentSectionIsND = false;
    currentSectionIsWarping = false;
    if (strcmp(argv[1],"NDFiber") == 0)
      currentSectionIsND = true;
    if (strcmp(argv[1],"NDFiberWarping") == 0) {
      currentSectionIsND = true;
      currentSectionIsWarping = true;
    }

    // create the fiber section representation (with the geometric information) 
      
    SectionRepres *fiberSectionRepr =
	new FiberSectionRepr(secTag, maxNumPatches, maxNumReinfLayers);  

    if (fiberSectionRepr == 0) {
      opserr <<  "WARNING - ran out of memory to create section representation\n";
	return TCL_ERROR;
    }

    if (theTclModelBuilder->addSectionRepres(*fiberSectionRepr) < 0) {
	opserr <<  "WARNING - cannot add section representation\n";
	return TCL_ERROR;
    }	

    int brace = 3; // Start of recursive parse
    double GJ = 1.0;
    bool isTorsion = false;
    if (strcmp(argv[3],"-GJ") == 0) {
      if (Tcl_GetDouble(interp, argv[4], &GJ) != TCL_OK) {
	opserr << "WARNING invalid GJ";
	return TCL_ERROR;
      }
      isTorsion = true;
      brace = 5;
    }

    // parse the information inside the braces (patches and reinforcing layers)
    if (Tcl_Eval(interp, argv[brace]) != TCL_OK) {
	opserr << "WARNING - error reading information in { } \n";
	return TCL_ERROR;
    }

    // build the fiber section (for analysis)
    if (buildSection(interp, theTclModelBuilder, secTag, isTorsion, GJ) != TCL_OK) {
	opserr << "WARNING - error constructing the section\n";
	return TCL_ERROR;
    }
    
//    currentSectionTag = 0;

    return TCL_OK;
}


int
TclCommand_addFiberIntSection (ClientData clientData, Tcl_Interp *interp, int argc,
				 TCL_Char **argv, TclModelBuilder *theTclModelBuilder)
{
    int secTag;
    int maxNumPatches = 30; 
    int maxNumReinfLayers = 30;
    
    if (argc < 4) 
	return TCL_ERROR;
    
    if (Tcl_GetInt(interp, argv[2], &secTag) != TCL_OK) {
      opserr <<  "WARNING bad command - want: \nsection fiberSec secTag { \n\tpatch <patch arguments> \n\tlayer <layer arguments> \n}\n";
	return TCL_ERROR;
    }
    
    currentSectionTag = secTag;
      
    // create the fiber section representation (with the geometric information) 
      
    SectionRepres *fiberSectionRepr =
	new FiberSectionRepr(secTag, maxNumPatches, maxNumReinfLayers);  

    if (fiberSectionRepr == 0) {
      opserr <<  "WARNING - ran out of memory to create section representation\n";
	return TCL_ERROR;
    }

    if (theTclModelBuilder->addSectionRepres(*fiberSectionRepr) < 0) {
	opserr <<  "WARNING - cannot add section representation\n";
	return TCL_ERROR;
    }	

    int brace = 3; // Start of recursive parse
    double GJ = 1.0;
    bool isTorsion = false;
    if (strcmp(argv[3],"-GJ") == 0) {
      if (Tcl_GetDouble(interp, argv[4], &GJ) != TCL_OK) {
	opserr <<  "WARNING invalid GJ";
	return TCL_ERROR;
      }
      isTorsion = true;
      brace = 5;
    }

	int NStrip1, NStrip2, NStrip3;
	double t1, t2, t3;

    if (strcmp(argv[3],"-NStrip") == 0) {
      
	if (Tcl_GetInt(interp, argv[4], &NStrip1) != TCL_OK) {
	opserr <<  "WARNING invalid NStrip1";
	return TCL_ERROR;
      }

	if (Tcl_GetDouble(interp, argv[5], &t1) != TCL_OK) {
	opserr <<  "WARNING invalid t1";
	return TCL_ERROR;
      }

	if (Tcl_GetInt(interp, argv[6], &NStrip2) != TCL_OK) {
	opserr <<  "WARNING invalid NStrip2";
	return TCL_ERROR;
      }

	if (Tcl_GetDouble(interp, argv[7], &t2) != TCL_OK) {
	opserr <<  "WARNING invalid t2";
	return TCL_ERROR;
      }

	if (Tcl_GetInt(interp, argv[8], &NStrip3) != TCL_OK) {
	opserr <<  "WARNING invalid NStrip3";
	return TCL_ERROR;
      }

	if (Tcl_GetDouble(interp, argv[9], &t3) != TCL_OK) {
	opserr <<  "WARNING invalid t3";
	return TCL_ERROR;
      }

      //isTorsion = true;
      brace = 10; //may be 5
    }



    // parse the information inside the braces (patches and reinforcing layers)
    if (Tcl_Eval(interp, argv[brace]) != TCL_OK) {
	opserr << "WARNING - error reading information in { } \n";
	return TCL_ERROR;
    }


    // build the fiber section (for analysis)
    if (buildSectionInt(interp, theTclModelBuilder, secTag, isTorsion, GJ, NStrip1, t1, NStrip2, t2, NStrip3, t3) != TCL_OK) {
	opserr << "WARNING - error constructing the section\n";
	return TCL_ERROR;
    }

//    currentSectionTag = 0;

    return TCL_OK;
}


// add patch to fiber section
int
TclCommand_addPatch(ClientData clientData, Tcl_Interp *interp, int argc, 
			     TCL_Char **argv, TclModelBuilder *theTclModelBuilder)
{
    // check if a section is being processed
    if (currentSectionTag == 0) {
	opserr <<  "WARNING subcommand 'patch' is only valid inside a 'section' command\n";
	return TCL_ERROR;
    }	   
    
    // make sure at least one other argument to contain patch type
    if (argc < 2) {
	opserr <<  "WARNING need to specify a patch type \n";
	return TCL_ERROR;
    }    

    // check argv[1] for type of patch  and create the object
    if (strcmp(argv[1], "quad") == 0 || strcmp(argv[1], "quadr") == 0) {
	int numSubdivIJ, numSubdivJK, matTag, secTag;
	double vertexCoordY, vertexCoordZ;
	static Matrix vertexCoords(4,2);
	int j, argi;

	if (argc < 13) {
	    opserr <<  "WARNING invalid number of parameters: patch quad matTag numSubdivIJ numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL zVertL\n";
	    return TCL_ERROR;
	}
  
	argi = 2;
      
      if (Tcl_GetInt(interp, argv[argi++], &matTag) != TCL_OK)
      {
         opserr <<  "WARNING invalid matTag: patch quad matTag numSubdivIJ numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL zVertL\n";
         return TCL_ERROR;
      }
      //opserr << "\n\tmatTag: " << matTag;

      if (Tcl_GetInt(interp, argv[argi++], &numSubdivIJ) != TCL_OK)
      {
         opserr <<  "WARNING invalid numSubdivIJ: patch quad matTag numSubdivIJ numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL zVertL\n";
         return TCL_ERROR;
      }
      //opserr << "\n\tnumSubdivIJ: " << numSubdivIJ;
 
      if (Tcl_GetInt(interp, argv[argi++], &numSubdivJK) != TCL_OK)
      {
         opserr <<  "WARNING invalid numSubdivJK: patch quad matTag numSubdivIJ numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL zVertL\n";
         return TCL_ERROR;
      }
      //opserr << "\n\tnumSubdivJK: " << numSubdivJK;

      for (j=0; j < 4; j++)
      {
         //opserr << "\n\tVertexCoord: " << j;
         if (Tcl_GetDouble(interp, argv[argi++], &vertexCoordY) != TCL_OK)
         {
            opserr <<  "WARNING invalid Coordinate y: ...yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL zVertL\n";
            return TCL_ERROR;
         }
         //opserr << "\n\t\tvertexCoordY: " << vertexCoordY; 

         if (Tcl_GetDouble(interp, argv[argi++], &vertexCoordZ) != TCL_OK)
         {
            opserr <<  "WARNING invalid Coordinate z: ...yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL zVertL\n";
            return TCL_ERROR;
         }
         //opserr << "\n\t\tvertexCoordZ: " << vertexCoordZ; 

         vertexCoords(j,0) = vertexCoordY;
         vertexCoords(j,1) = vertexCoordZ;
      }
       
      // get section representation
      secTag = currentSectionTag;
      
      SectionRepres *sectionRepres = theTclModelBuilder->getSectionRepres(secTag);
      if (sectionRepres == 0) 
      {
         opserr <<  "WARNING cannot retrieve section\n";
         return TCL_ERROR;
      }    
     
      if (sectionRepres->getType() != SEC_TAG_FiberSection)
      {
         opserr <<  "WARNING section invalid: patch can only be added to fiber sections\n";
         return TCL_ERROR;
      }

      FiberSectionRepr *fiberSectionRepr = (FiberSectionRepr *) sectionRepres;

      // create patch

      QuadPatch *patch = new QuadPatch(matTag, numSubdivIJ, numSubdivJK, vertexCoords);
      if (!patch)
      {
         opserr <<  "WARNING cannot allocate patch\n";
         return TCL_ERROR;
      }

      //opserr << "\n\tpatch: " << *patch;
      
      // add patch to section representation

      int error = fiberSectionRepr->addPatch(*patch);
      delete patch;
      
      if (error)
      {
         opserr <<  "WARNING cannot add patch to section\n";
         return TCL_ERROR;
      }  
  }
    
    
    // check argv[1] for type of patch  and create the object
    else if (strcmp(argv[1], "rect") == 0 || 
	     strcmp(argv[1], "rectangular") == 0) {
	
	int numSubdivIJ, numSubdivJK, matTag, secTag;
	double vertexCoordY, vertexCoordZ;
	static Matrix vertexCoords(4,2);
	int j, argi;

	if (argc < 9) {
	    opserr <<  "WARNING invalid number of parameters: patch quad matTag numSubdivIJ numSubdivJK yVertI zVertI yVertK zVertK\n";
	    return TCL_ERROR;
	}
  
	argi = 2;
      
	if (Tcl_GetInt(interp, argv[argi++], &matTag) != TCL_OK) {
	    opserr <<  "WARNING invalid matTag: patch quad matTag numSubdivIJ numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL zVertL\n";
	    return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[argi++], &numSubdivIJ) != TCL_OK) {
	    opserr <<  "WARNING invalid numSubdivIJ: patch quad matTag numSubdivIJ numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL zVertL\n";
	    return TCL_ERROR;
	}
 
	if (Tcl_GetInt(interp, argv[argi++], &numSubdivJK) != TCL_OK) {
	    opserr <<  "WARNING invalid numSubdivJK: patch quad matTag numSubdivIJ numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL zVertL\n";
	    return TCL_ERROR;
	}

	for (j=0; j < 2; j++) {
	    if (Tcl_GetDouble(interp, argv[argi++], &vertexCoordY) != TCL_OK) {
		opserr <<  "WARNING invalid Coordinate y: ...yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL zVertL\n";
		return TCL_ERROR;
	    }

	    if (Tcl_GetDouble(interp, argv[argi++], &vertexCoordZ) != TCL_OK) {
		opserr <<  "WARNING invalid Coordinate z: ...yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL zVertL\n";
		return TCL_ERROR;
	    }

	    vertexCoords(j*2,0) = vertexCoordY;
	    vertexCoords(j*2,1) = vertexCoordZ;
	}

	vertexCoords(1,0) = vertexCoords(2,0);
	vertexCoords(1,1) = vertexCoords(0,1);	
	vertexCoords(3,0) = vertexCoords(0,0);
	vertexCoords(3,1) = vertexCoords(2,1);		
	    
      // get section representation
      secTag = currentSectionTag;
      
      SectionRepres *sectionRepres = theTclModelBuilder->getSectionRepres(secTag);
      if (sectionRepres == 0) {
         opserr <<  "WARNING cannot retrieve section\n";
         return TCL_ERROR;
      }    
     
      if (sectionRepres->getType() != SEC_TAG_FiberSection) {
         opserr <<  "WARNING section invalid: patch can only be added to fiber sections\n";
         return TCL_ERROR;
      }

      FiberSectionRepr *fiberSectionRepr = (FiberSectionRepr *) sectionRepres;

      // create patch

      QuadPatch *patch = new QuadPatch(matTag, numSubdivIJ, numSubdivJK, vertexCoords);
      if (!patch)
      {
         opserr <<  "WARNING cannot allocate patch\n";
         return TCL_ERROR;
      }

      //opserr << "\n\tpatch: " << *patch;
      
      // add patch to section representation

      int error = fiberSectionRepr->addPatch(*patch);
      delete patch;
      
      if (error)
      {
         opserr <<  "WARNING cannot add patch to section\n";
         return TCL_ERROR;
      }  
  }    
    
    
    
         
    else if (strcmp(argv[1], "circ") == 0) {
	int numSubdivRad, numSubdivCirc, matTag, secTag;
	double yCenter, zCenter;
	static Vector centerPosition(2);
	double intRad, extRad;
	double startAng, endAng;

	int argi;

      argi = 2;
      if (argc < 11)
      {
         opserr <<  "WARNING invalid number of parameters: patch circ matTag numSubdivCirc numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
         return TCL_ERROR;
      }
  
      if (Tcl_GetInt(interp, argv[argi++], &matTag) != TCL_OK)
      {
         opserr <<  "WARNING invalid matTag: patch circ matTag numSubdivCirc numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
         return TCL_ERROR;
      }
      //opserr << "\n\tmatTag: " << matTag;

      if (Tcl_GetInt(interp, argv[argi++], &numSubdivCirc) != TCL_OK)
      {
         opserr <<  "WARNING invalid numSubdivCirc: patch circ matTag numSubdivCirc numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
         return TCL_ERROR;
      }
      //opserr << "\n\tnumSubdivCirc: " << numSubdivCirc;

      if (Tcl_GetInt(interp, argv[argi++], &numSubdivRad) != TCL_OK)
      {
         opserr <<  "WARNING invalid numSubdivRad: patch circ matTag numSubdivCirc numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
         return TCL_ERROR;
      }
      //opserr << "\n\tnumSubdivRad: " << numSubdivRad;

      if (Tcl_GetDouble(interp, argv[argi++], &yCenter) != TCL_OK)
      {
         opserr <<  "WARNING invalid yCenter: patch circ matTag numSubdivCirc numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
         return TCL_ERROR;
      }
      //opserr << "\n\tyCenter: " << yCenter;

      if (Tcl_GetDouble(interp, argv[argi++], &zCenter) != TCL_OK)
      {
         opserr <<  "WARNING invalid zCenter: patch circ matTag numSubdivCirc numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
         return TCL_ERROR;
      }
      //opserr << "\n\tzCenter: " << zCenter;

      if (Tcl_GetDouble(interp, argv[argi++], &intRad) != TCL_OK)
      {
         opserr <<  "WARNING invalid intRad: patch circ matTag numSubdivCirc numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
         return TCL_ERROR;
      }
      //opserr << "\n\tintRad: " << intRad;

      if (Tcl_GetDouble(interp, argv[argi++], &extRad) != TCL_OK)
      {
         opserr <<  "WARNING invalid extRad: patch circ matTag numSubdivCirc numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
         return TCL_ERROR;
      }
      //opserr << "\n\textRad: " << extRad;

      if (Tcl_GetDouble(interp, argv[argi++], &startAng) != TCL_OK)
      {
         opserr <<  "WARNING invalid startAng: patch circ matTag numSubdivCirc numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
         return TCL_ERROR;
      }
      //opserr << "\n\tstartAngle: " << startAng;

      if (Tcl_GetDouble(interp, argv[argi++], &endAng) != TCL_OK)
      {
         opserr <<  "WARNING invalid endAng: patch circ matTag numSubdivCirc numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
         return TCL_ERROR;
      }
      //opserr << "\n\tendAng: " << endAng;


      // get section 
      secTag = currentSectionTag;
      
      SectionRepres *sectionRepres = theTclModelBuilder->getSectionRepres(secTag);
      if (sectionRepres == 0) 
      {
         opserr <<  "WARNING cannot retrieve section\n";
         return TCL_ERROR;
      }    
     
      if (sectionRepres->getType() != SEC_TAG_FiberSection)
      {
         opserr <<  "WARNING section invalid: patch can only be added to fiber sections\n";
         return TCL_ERROR;
      }

      FiberSectionRepr *fiberSectionRepr = (FiberSectionRepr *) sectionRepres;

      centerPosition(0) = yCenter; 
      centerPosition(1) = zCenter; 
   
      // create patch

      CircPatch *patch = new CircPatch(matTag, numSubdivCirc, numSubdivRad,
                                       centerPosition, intRad, extRad, 
                                       startAng, endAng);
      if (!patch)
      {
         opserr <<  "WARNING cannot allocate patch\n";
         return TCL_ERROR;
      }

      //opserr << "\n\tpatch: " << *patch;
      
      // add patch to section

      int error = fiberSectionRepr->addPatch(*patch);
      delete patch;
      
      if (error)
      {
         opserr <<  "WARNING cannot add patch to section\n";
         return TCL_ERROR;
      }
   }

   else
   {
      opserr <<  "WARNING patch type is not available\n";
      return TCL_ERROR;
   }
  
   return TCL_OK;
}



// add patch to fiber section
int
TclCommand_addFiber(ClientData clientData, Tcl_Interp *interp, int argc, 
		    TCL_Char **argv, TclModelBuilder *theTclModelBuilder)
{
    // check if a section is being processed
    if (currentSectionTag == 0) {
	opserr <<  "WARNING subcommand 'fiber' is only valid inside a 'section' command\n";
	return TCL_ERROR;
    }	   
    
    // make sure at least one other argument to contain patch type
    if (argc < 5) {
	opserr <<  "WARNING invalid num args: fiber yLoc zLoc area matTag\n";
	return TCL_ERROR;
    }    

    SectionRepres *sectionRepres = 
	theTclModelBuilder->getSectionRepres(currentSectionTag);
    
    if (sectionRepres == 0) {
	opserr <<  "WARNING cannot retrieve section\n";
	return TCL_ERROR;
    }    
	
    if (sectionRepres->getType() != SEC_TAG_FiberSection) {
	opserr <<  "WARNING section invalid: patch can only be added to fiber sections\n";
	return TCL_ERROR;
    }

    FiberSectionRepr *fiberSectionRepr = (FiberSectionRepr *) sectionRepres;
    int numFibers = fiberSectionRepr->getNumFibers();    
    
    Fiber *theFiber =0;
      
    int matTag;
    double yLoc, zLoc, area;

    
    if (Tcl_GetDouble(interp, argv[1], &yLoc) != TCL_OK) {
      opserr <<  "WARNING invalid yLoc: fiber yLoc zLoc area matTag\n";
      return TCL_ERROR;
    }    
    if (Tcl_GetDouble(interp, argv[2], &zLoc) != TCL_OK) {
      opserr <<  "WARNING invalid zLoc: fiber yLoc zLoc area matTag\n";
      return TCL_ERROR;
    }        
    if (Tcl_GetDouble(interp, argv[3], &area) != TCL_OK) {
      opserr <<  "WARNING invalid area: fiber yLoc zLoc area matTag\n";
      return TCL_ERROR;
    }            
    if (Tcl_GetInt(interp, argv[4], &matTag) != TCL_OK) {
      opserr <<  "WARNING invalid matTag: fiber yLoc zLoc area matTag\n";
      return TCL_ERROR;
    }                

    int NDM = theTclModelBuilder->getNDM();  
        
    // creates 2d section      
    if (NDM == 2) {

      if (currentSectionIsND) {
	NDMaterial *material = OPS_getNDMaterial(matTag);
	if (material == 0) {
	  opserr <<  "WARNING invalid NDMaterial ID for patch\n";
	  return TCL_ERROR;
	}  
	theFiber = new NDFiber2d(numFibers, *material, area, yLoc);
      }
      else {
	UniaxialMaterial *material = OPS_getUniaxialMaterial(matTag);
	if (material == 0) {
	  opserr <<  "WARNING invalid UniaxialMaterial ID for patch\n";
	  return TCL_ERROR;
	}   
	theFiber = new UniaxialFiber2d(numFibers, *material, area, yLoc);
      }
      
      if (theFiber == 0) {
	opserr <<  "WARNING unable to allocate fiber \n";
	return TCL_ERROR;
      }    
    }

    else if (NDM == 3) {
      
      static Vector fiberPosition(2);
      fiberPosition(0) = yLoc;
      fiberPosition(1) = zLoc;
      
      if (currentSectionIsND) {
	NDMaterial *material = OPS_getNDMaterial(matTag);
	if (material == 0) {
	  opserr <<  "WARNING invalid NDMaterial ID for patch\n";
	  return TCL_ERROR;
	}
	theFiber = new NDFiber3d(numFibers, *material, area, yLoc, zLoc);
      }
      else {
	UniaxialMaterial *material = OPS_getUniaxialMaterial(matTag);
	if (material == 0) {
	  opserr <<  "WARNING invalid UniaxialMaterial ID for patch\n";
	  return TCL_ERROR;
	}   
	theFiber = new UniaxialFiber3d(numFibers, *material, area,
				       fiberPosition);
      }

      if (theFiber == 0) {
	opserr <<  "WARNING unable to allocate fiber \n";
	return TCL_ERROR;
      }    
    }

    else {
	opserr <<  "WARNING fiber command for FiberSection only for 2 or 3d \n";
	return TCL_ERROR;
    }    
	
    // add patch to section representation
    int error = fiberSectionRepr->addFiber(*theFiber);

    if (error) {
	opserr <<  "WARNING cannot add patch to section\n";
	return TCL_ERROR;
    }  

    return TCL_OK;
}

// add Hfiber to fiber section
int
TclCommand_addHFiber(ClientData clientData, Tcl_Interp *interp, int argc, 
			 TCL_Char **argv, TclModelBuilder *theTclModelBuilder)
{
    // check if a section is being processed
    if (currentSectionTag == 0) {
	opserr <<  "WARNING subcommand 'Hfiber' is only valid inside a 'section' command\n";
	return TCL_ERROR;
    }	   
    
    // make sure at least one other argument to contain patch type
    if (argc < 5) {
	opserr <<  "WARNING invalid num args: Hfiber yLoc zLoc area matTag\n";
	return TCL_ERROR;
    }    

	SectionRepres *sectionHRepres = 
	theTclModelBuilder->getSectionRepres(currentSectionTag);

	
	if (sectionHRepres == 0) {
	opserr <<  "WARNING cannot retrieve section\n";
	return TCL_ERROR;
    }    
	
    if (sectionHRepres->getType() != SEC_TAG_FiberSection) {
	opserr <<  "WARNING section invalid: patch can only be added to fiber sections\n";
	return TCL_ERROR;
    }

    FiberSectionRepr *fiberSectionHRepr = (FiberSectionRepr *) sectionHRepres;
    int numHFibers = fiberSectionHRepr->getNumHFibers();    
    
    int HNDM = theTclModelBuilder->getNDM();  
    
    Fiber *theHFiber =0;
      
    int matHTag;
    double yHLoc, zHLoc, Harea;

    
    if (Tcl_GetDouble(interp, argv[1], &yHLoc) != TCL_OK) {
         opserr <<  "WARNING invalid yLoc: Hfiber yLoc zLoc area matTag\n";
         return TCL_ERROR;
     }    
    if (Tcl_GetDouble(interp, argv[2], &zHLoc) != TCL_OK) {
         opserr <<  "WARNING invalid zLoc: Hfiber yLoc zLoc area matTag\n";
         return TCL_ERROR;
     }        
    if (Tcl_GetDouble(interp, argv[3], &Harea) != TCL_OK) {
         opserr <<  "WARNING invalid area: Hfiber yLoc zLoc area matTag\n";
         return TCL_ERROR;
     }            
    
    if (Tcl_GetInt(interp, argv[4], &matHTag) != TCL_OK) {
         opserr <<  "WARNING invalid matTag: Hfiber yLoc zLoc area matTag\n";
         return TCL_ERROR;
     }                
    
    UniaxialMaterial *Hmaterial = OPS_getUniaxialMaterial(matHTag);
    
    // creates 2d section      
    if (HNDM == 2) {

	if (Hmaterial == 0) {
	    opserr <<  "WARNING invalid Hmaterial ID for patch\n";
	    return TCL_ERROR;
	}   

	theHFiber = new UniaxialFiber2d(numHFibers, *Hmaterial, Harea, yHLoc);
	if (theHFiber == 0) {
	    opserr <<  "WARNING unable to allocate Hfiber \n";
	    return TCL_ERROR;
	}    
    }

    else if (HNDM == 3) {

      static Vector fiberHPosition(2);
	fiberHPosition(0) = yHLoc;
	fiberHPosition(1) = zHLoc;
	    
	theHFiber = new UniaxialFiber3d(numHFibers, *Hmaterial, Harea, fiberHPosition);
	if (theHFiber == 0) {
	    opserr <<  "WARNING unable to allocate Hfiber \n";
	    return TCL_ERROR;
	}    
    }

    else {
	opserr <<  "WARNING Hfiber command for FiberSection only fo 2 or 3d \n";
	return TCL_ERROR;
    }    
	
    // add patch to section representation
    int error = fiberSectionHRepr->addHFiber(*theHFiber);

    if (error) {
	opserr <<  "WARNING cannot add patch to section\n";
	return TCL_ERROR;
    }  

    return TCL_OK;
}



// add layers of reinforcing bars to fiber section
          
int
TclCommand_addReinfLayer(ClientData clientData, Tcl_Interp *interp, int argc, 
				  TCL_Char **argv, TclModelBuilder *theTclModelBuilder)
{
   //opserr << "\nreading layer:\n";

   // check if a section is being processed
   if (currentSectionTag == 0)
   {
      opserr <<  "WARNING subcommand 'patch' is only valid inside a 'section' command\n";
      return TCL_ERROR;
   }	   

   // make sure at least one other argument to contain layer type
   if (argc < 2) 
   {
      opserr <<  "WARNING need to specify a layer type \n";
      return TCL_ERROR;
   }    

   // check argv[1] for type of layer and create the object
   if (strcmp(argv[1], "straight") == 0) 
   {
      if (argc < 9)
      {
         opserr <<  "WARNING invalid number of parameters: layer straight matTag numReinfBars reinfBarArea yStartPt zStartPt yEndPt zEndPt\n";
         return TCL_ERROR;
      }

      int secTag, matTag, numReinfBars;
      double reinfBarArea;
      double yStartPt, zStartPt, yEndPt, zEndPt;
     
      int argi;

      argi = 2;
      
      if (Tcl_GetInt(interp, argv[argi++], &matTag) != TCL_OK)
      {
         opserr <<  "WARNING invalid matTag: layer straight matTag numReinfBars reinfBarArea  yStartPt zStartPt yEndPt zEndPt\n";
         return TCL_ERROR;
      }
      //opserr << "\n\tmatTag: " << matTag;

      if (Tcl_GetInt(interp, argv[argi++], &numReinfBars) != TCL_OK)
      {
         opserr <<  "WARNING invalid numReinfBars: layer straight matTag numReinfBars reinfBarArea  yStartPt zStartPt yEndPt zEndPt\n";
         return TCL_ERROR;
      }
      //opserr << "\n\tnumReinfBars: " << numReinfBars;

      if (Tcl_GetDouble(interp, argv[argi++], &reinfBarArea) != TCL_OK)
      {
         opserr <<  "WARNING invalid reinfBarArea: layer straight matTag numReinfBars reinfBarArea  yStartPt zStartPt yEndPt zEndPt\n";
         return TCL_ERROR;
      }
      //opserr << "\n\treinfBarArea: " << reinfBarArea;

      if (Tcl_GetDouble(interp, argv[argi++], &yStartPt) != TCL_OK)
      {
         opserr <<  "WARNING invalid yStartPt: layer straight matTag numReinfBars reinfBarArea  yStartPt zStartPt yEndPt zEndPt\n";
         return TCL_ERROR;
      }
      //opserr << "\n\tyStartPt: " << yStartPt;
    
      if (Tcl_GetDouble(interp, argv[argi++], &zStartPt) != TCL_OK)
      {
         opserr <<  "WARNING invalid zStartPt: layer straight matTag numReinfBars reinfBarArea  yStartPt zStartPt yEndPt zEndPt\n";
         return TCL_ERROR;
      }
      //opserr << "\n\tzStartPt: " << zStartPt;
       
      if (Tcl_GetDouble(interp, argv[argi++], &yEndPt) != TCL_OK)
      {
         opserr <<  "WARNING invalid yEndPt: layer straight matTag numReinfBars reinfBarArea  yStartPt zStartPt yEndPt zEndPt\n";
         return TCL_ERROR;
      }
      //opserr << "\n\tyEndPt: " << yEndPt;
    
      if (Tcl_GetDouble(interp, argv[argi++], &zEndPt) != TCL_OK)
      {
         opserr <<  "WARNING invalid zEndPt: layer straight matTag numReinfBars reinfBarArea  yStartPt zStartPt yEndPt zEndPt\n";
         return TCL_ERROR;
      }
      
      //opserr << "\n\tzEndPt: " << zEndPt;
      
      // get section 
      secTag = currentSectionTag;

      SectionRepres *sectionRepres = theTclModelBuilder->getSectionRepres(secTag);
      if (sectionRepres == 0) 
      {
         opserr <<  "WARNING cannot retrieve section\n";
         return TCL_ERROR;
      }    
     
      if (sectionRepres->getType() != SEC_TAG_FiberSection)
      {
         opserr <<  "WARNING section invalid: patch can only be added to fiber sections\n";
         return TCL_ERROR;
      }

      FiberSectionRepr *fiberSectionRepr = (FiberSectionRepr *) sectionRepres;
 
      // create the reinforcing layer

      static Vector startPt(2);
      static Vector endPt(2);

      startPt(0)  = yStartPt;
      startPt(1)  = zStartPt;
      endPt(0) = yEndPt;
      endPt(1) = zEndPt;

      StraightReinfLayer *reinfLayer = new StraightReinfLayer (matTag,
                                                   numReinfBars, reinfBarArea,
                                                   startPt, endPt);
      if (!reinfLayer)
      {
         opserr <<  "WARNING cannot allocate reinfLayer\n";
         return TCL_ERROR;
      }
      //opserr << "\nStraigthReinfLayer: " << *reinfLayer;

      // add reinfLayer to section
      int error = fiberSectionRepr->addReinfLayer(*reinfLayer);
      delete reinfLayer;
      
      if (error)
      {
         opserr <<  "WARNING cannot add reinforcing layer to section\n";
         return TCL_ERROR;
      }
      
   }
   else if (strcmp(argv[1], "circ") == 0) 
   {
      if (argc < 8)
      {
         opserr <<  "WARNING invalid number of parameters: layer circ matTag numReinfBars reinfBarArea yCenter zCenter arcRadius <startAng endAng>\n";
         return TCL_ERROR;
      }

      int secTag, matTag, numReinfBars;
      double reinfBarArea;
      double yCenter, zCenter, radius, startAng, endAng;
     
      int argi;

      argi = 2;
      
      if (Tcl_GetInt(interp, argv[argi++], &matTag) != TCL_OK)
      {
         opserr <<  "WARNING invalid matTag: layer circ matTag numReinfBars reinfBarArea yCenter zCenter radius startAng endAng\n";
         return TCL_ERROR;
      }
      //opserr << "\n\tmatTag: " << matTag;

      if (Tcl_GetInt(interp, argv[argi++], &numReinfBars) != TCL_OK)
      {
         opserr <<  "WARNING invalid numReinfBars: layer circ matTag numReinfBars reinfBarArea yCenter zCenter radius startAng endAng\n";
         return TCL_ERROR;
      }
      //opserr << "\n\tnumReinfBars: " << numReinfBars;

      if (Tcl_GetDouble(interp, argv[argi++], &reinfBarArea) != TCL_OK)
      {
         opserr <<  "WARNING invalid reinfBarArea: layer circ matTag numReinfBars reinfBarArea yCenter zCenter radius startAng endAng\n";
         return TCL_ERROR;
      }
      //opserr << "\n\treinfBarArea: " << reinfBarArea;

      if (Tcl_GetDouble(interp, argv[argi++], &yCenter) != TCL_OK)
      {
         opserr <<  "WARNING invalid yCenter: layer circ matTag numReinfBars reinfBarArea yCenter zCenter radius startAng endAng\n";
         return TCL_ERROR;
      }
      //opserr << "\n\tyCenter: " << yCenter;
    
      if (Tcl_GetDouble(interp, argv[argi++], &zCenter) != TCL_OK)
      {
         opserr <<  "WARNING invalid zCenter: layer circ matTag numReinfBars reinfBarArea yCenter zCenter radius startAng endAng\n";
         return TCL_ERROR;
      }
      //opserr << "\n\tzCenter: " << zCenter;
       
      if (Tcl_GetDouble(interp, argv[argi++], &radius) != TCL_OK)
      {
         opserr <<  "WARNING invalid radius: layer circ matTag numReinfBars reinfBarArea yCenter zCenter radius startAng endAng\n";
         return TCL_ERROR;
      }
      //opserr << "\n\tradius: " << radius;
    
	  bool anglesSpecified = false;

      if (argc > 9) {
		  if (Tcl_GetDouble(interp, argv[argi++], &startAng) != TCL_OK)
		  {
			 opserr <<  "WARNING invalid startAng: layer circ matTag numReinfBars reinfBarArea yCenter zCenter radius startAng endAng\n";
			 return TCL_ERROR;
		  }
		  //opserr << "\n\tstartAng: " << startAng;

		  if (Tcl_GetDouble(interp, argv[argi++], &endAng) != TCL_OK)
		  {
			 opserr <<  "WARNING invalid endAng: layer circ matTag numReinfBars reinfBarArea yCenter zCenter radius startAng endAng\n";
			 return TCL_ERROR;
		  }
		  //opserr << "\n\tendAng: " << endAng;

		  anglesSpecified = true;
	  }

      // get section 
      secTag = currentSectionTag;
      
      SectionRepres *sectionRepres = theTclModelBuilder->getSectionRepres(secTag);
      if (sectionRepres == 0) 
      {
         opserr <<  "WARNING cannot retrieve section\n";
         return TCL_ERROR;
      }    
     
      if (sectionRepres->getType() != SEC_TAG_FiberSection)
      {
         opserr <<  "WARNING section invalid: patch can only be added to fiber sections\n";
         return TCL_ERROR;
      }

      FiberSectionRepr *fiberSectionRepr = (FiberSectionRepr *) sectionRepres;
 
      // create the reinforcing layer

      static Vector center(2);

      center(0) = yCenter; 
      center(1) = zCenter; 

      CircReinfLayer *reinfLayer = 0;
	  if (anglesSpecified)
		  // Construct arc
		  reinfLayer = new CircReinfLayer (matTag, numReinfBars, reinfBarArea,
			center, radius, startAng, endAng);
	  else
		  // Construct circle
		  reinfLayer = new CircReinfLayer (matTag, numReinfBars, reinfBarArea,
			center, radius);

      if (!reinfLayer)
      {
         opserr <<  "WARNING cannot allocate reinfLayer\n";
         return TCL_ERROR;
      }
      //opserr << "\nCircReinfLayer: " << *reinfLayer;

      // add reinfLayer to section
      int error = fiberSectionRepr->addReinfLayer(*reinfLayer);
      delete reinfLayer;
      
      if (error)
      {
         opserr <<  "WARNING cannot add reinforcing layer to section\n";
         return TCL_ERROR;
      }
      
   }
   else
   {
      opserr <<  "WARNING reinforcing layer type is not available\n";
      return TCL_ERROR;
   }
  
   return TCL_OK;
}    



// build the section
int 
buildSection(Tcl_Interp *interp, TclModelBuilder *theTclModelBuilder,
	     int secTag, bool isTorsion, double GJ)
{
   SectionRepres *sectionRepres = theTclModelBuilder->getSectionRepres(secTag);
   if (sectionRepres == 0) 
   {
      opserr <<  "WARNING cannot retrieve section\n";
      return TCL_ERROR;
   }    
     
   if (sectionRepres->getType() == SEC_TAG_FiberSection)
   {
      // build the section
  
      FiberSectionRepr *fiberSectionRepr = (FiberSectionRepr *) sectionRepres;

      int i, j, k;
      int numFibers;
      
      int numPatches;
      Patch **patch;

      int  numReinfLayers;
      ReinfLayer **reinfLayer;

      numPatches     = fiberSectionRepr->getNumPatches();
      patch          = fiberSectionRepr->getPatches();
      numReinfLayers = fiberSectionRepr->getNumReinfLayers();
      reinfLayer     = fiberSectionRepr->getReinfLayers(); 

      int numSectionRepresFibers = fiberSectionRepr->getNumFibers();
      Fiber **sectionRepresFibers = fiberSectionRepr->getFibers();
      
      numFibers = numSectionRepresFibers;
      for (i = 0; i < numPatches; i++)
         numFibers += patch[i]->getNumCells();
      
      for (i = 0; i < numReinfLayers; i++)
         numFibers += reinfLayer[i]->getNumReinfBars();
      
      //opserr << "\nnumFibers: " << numFibers;
      
      static Vector fiberPosition(2);
      int    matTag;
      
      ID     fibersMaterial(numFibers-numSectionRepresFibers);
      Matrix fibersPosition(2,numFibers-numSectionRepresFibers);
      Vector fibersArea(numFibers-numSectionRepresFibers);

      int  numCells;
      Cell **cell;
    
      k = 0;
      for (i = 0; i < numPatches; i++)
      {
         //opserr << "\nPatch :" << i;
      
         numCells   = patch[i]->getNumCells();
         matTag = patch[i]->getMaterialID();

         //opserr << "\nmatTag: " << matTag(k);

         cell = patch[i]->getCells();

         if (cell == 0)
         {
            opserr <<  "WARNING out of run to create fibers\n";
            return TCL_ERROR;
         }    
         
         //opserr << "\n\tnumCells :" << numCells;
      
         for (j = 0; j < numCells; j++)
         {
	    fibersMaterial(k) = matTag;
            fibersArea(k)     = cell[j]->getArea();
            fiberPosition     = cell[j]->getCentroidPosition();

            fibersPosition(0,k) = fiberPosition(0);
	    fibersPosition(1,k) = fiberPosition(1);
	      
            k++;
         }
  
         for (j = 0; j < numCells; j++)
           delete cell[j];
  
         delete [] cell;
      }
         
      ReinfBar *reinfBar;
      int numReinfBars;

      for (i = 0; i < numReinfLayers; i++)
      {
         numReinfBars = reinfLayer[i]->getNumReinfBars();
         reinfBar     = reinfLayer[i]->getReinfBars();
         matTag  = reinfLayer[i]->getMaterialID();
   
         for (j = 0; j < numReinfBars; j++)
         {
	    fibersMaterial(k) = matTag; 
            fibersArea(k)     = reinfBar[j].getArea();
            fiberPosition     = reinfBar[j].getPosition();
     
	    fibersPosition(0,k) = fiberPosition(0);
	    fibersPosition(1,k) = fiberPosition(1);
	
            k++;
         }
         delete [] reinfBar;
      }

      UniaxialMaterial *material;
      NDMaterial *ndmaterial;

      // dimension of the structure (1d, 2d, or 3d)
      int NDM = theTclModelBuilder->getNDM();

      Fiber **fiber = new Fiber *[numFibers];
      if (fiber == 0) {
	  opserr <<  "WARNING unable to allocate fibers \n";
	  return TCL_ERROR;
      }          
      
      // copy the section repres fibers
      for (i=0; i<numSectionRepresFibers; i++)
	  fiber[i] = sectionRepresFibers[i];

      // creates 2d section      


      if (NDM == 2) {
	 k = 0;
	 for (i = numSectionRepresFibers; i < numFibers; i++) {    
	   if (currentSectionIsND) {
	     ndmaterial = OPS_getNDMaterial(fibersMaterial(k));
	     if (ndmaterial == 0) {
               opserr <<  "WARNING invalid NDmaterial ID for patch\n";
               return TCL_ERROR;
	     }
	     fiber[i] = new NDFiber2d(k, *ndmaterial, fibersArea(k), fibersPosition(0,k));
	   }
	   else {
	     material = OPS_getUniaxialMaterial(fibersMaterial(k));
	     if (material == 0) {
               opserr <<  "WARNING invalid UniaxialMaterial ID for patch\n";
               return TCL_ERROR;
	     }   
	     fiber[i] = new UniaxialFiber2d(k, *material, fibersArea(k), fibersPosition(0,k));
	   }
	   if (fiber[i] == 0) {
	     opserr <<  "WARNING unable to allocate fiber \n";
	     return TCL_ERROR;
	   }    
   
	   //opserr << *fiber[k];
	   k++;
	 }
	
	 
	 SectionForceDeformation *section = 0;
	 if (currentSectionIsND) {
           if (currentSectionIsWarping)
  	     section = new NDFiberSectionWarping2d(secTag, numFibers, fiber);
           else 
	     section = new NDFiberSection2d(secTag, numFibers, fiber);
         }
	 else
	   section = new FiberSection2d(secTag, numFibers, fiber);

	 //SectionForceDeformation *section = new FiberSection(secTag, numFibers, fiber);
   
	 // Delete fibers
	 for (i = 0; i < numFibers; i++)
	   delete fiber[i];

         if (section == 0)
         {
            opserr <<  "WARNING - cannot construct section\n";
            return TCL_ERROR;
         }
       
         //if (theTclModelBuilder->addSection (*section) < 0) {
	 if (OPS_addSectionForceDeformation(section) != true) {
            opserr <<  "WARNING - cannot add section\n";
            return TCL_ERROR;
         }

        //opserr << "section: " << *section;
     
      }
      else if (NDM == 3)     
      {

	 static Vector fiberPosition(2);
	 k = 0;
	 for (i = numSectionRepresFibers; i < numFibers; i++) {
	   fiberPosition(0) = fibersPosition(0,k);
	   fiberPosition(1) = fibersPosition(1,k);  
	   if (currentSectionIsND) {
	     ndmaterial = OPS_getNDMaterial(fibersMaterial(k));
	     if (ndmaterial == 0) {
               opserr <<  "WARNING invalid NDmaterial ID for patch\n";
               return TCL_ERROR;
	     }
	     fiber[i] = new NDFiber3d(k, *ndmaterial, fibersArea(k), fiberPosition(0), fiberPosition(1));
	   }
	   else {
	     material = OPS_getUniaxialMaterial(fibersMaterial(k));
	     if (material == 0) {
	       opserr <<  "WARNING invalid UniaxialMaterial ID for patch\n";
	       return TCL_ERROR;
	     }   
	     fiber[i] = new UniaxialFiber3d(k, *material, fibersArea(k), fiberPosition);
	   }
	   if (fiber[k] == 0) {
	     opserr <<  "WARNING unable to allocate fiber \n";
	     return TCL_ERROR;
	   }    
	   k++;
	   //opserr << *fiber[k];
	 }
	
	 //SectionForceDeformation *section = new FiberSection(secTag, numFibers, fiber);
	 SectionForceDeformation *section = 0;
	 if (currentSectionIsND)
	   section = new NDFiberSection3d(secTag, numFibers, fiber);
	 else if (isTorsion) {
           ElasticMaterial theGJ(0, GJ);
           //FiberSection3d theFS(0, numFibers, fiber);
           //section = new SectionAggregator(secTag, theFS, theGJ, SECTION_RESPONSE_T);
           section = new FiberSection3d(secTag, numFibers, fiber, &theGJ);
	 }
	 else
	   section = new FiberSection3d(secTag, numFibers, fiber);
   
	 // Delete fibers
	 for (i = 0; i < numFibers; i++)
	   delete fiber[i];

         if (section == 0)
         {
            opserr <<  "WARNING - cannot construct section\n";
            return TCL_ERROR;
         }
       
         //if (theTclModelBuilder->addSection (*section) < 0) {
	 if (OPS_addSectionForceDeformation(section) != true) {
            opserr <<  "WARNING - cannot add section\n";
            return TCL_ERROR;
         }

        //opserr << "section: " << *section;
       
      }
      else
      {
         opserr << "WARNING NDM = " << NDM << " is imcompatible with available frame elements\n";
         return TCL_ERROR;
      }

      // Delete fiber array
      delete [] fiber;

   }
   else 
   {
      opserr <<  "WARNING section invalid: can only build fiber sections\n";
      return TCL_ERROR;
   }    

   return TCL_OK;
}

// build the section Interaction
int 
buildSectionInt(Tcl_Interp *interp, TclModelBuilder *theTclModelBuilder,
	     int secTag, bool isTorsion, double GJ, int NStrip1, double t1, int NStrip2, double t2, int NStrip3, double t3)
{
   SectionRepres *sectionRepres = theTclModelBuilder->getSectionRepres(secTag);
   if (sectionRepres == 0) 
   {
      opserr <<  "WARNING cannot retrieve section\n";
      return TCL_ERROR;
   }    
     
   if (sectionRepres->getType() == SEC_TAG_FiberSection)
   {
      // build the section
  
      FiberSectionRepr *fiberSectionRepr = (FiberSectionRepr *) sectionRepres;

      int i, j, k;
      int numFibers;
      int numHFibers;
      
      int numPatches;
      Patch **patch;

      int  numReinfLayers;
      ReinfLayer **reinfLayer;

      numPatches     = fiberSectionRepr->getNumPatches();
      patch          = fiberSectionRepr->getPatches();
      numReinfLayers = fiberSectionRepr->getNumReinfLayers();
      reinfLayer     = fiberSectionRepr->getReinfLayers(); 

      int numSectionRepresFibers = fiberSectionRepr->getNumFibers();
      Fiber **sectionRepresFibers = fiberSectionRepr->getFibers();
 
	  int numSectionRepresHFibers = fiberSectionRepr->getNumHFibers();
      Fiber **sectionRepresHFibers = fiberSectionRepr->getHFibers();
 

      numFibers = numSectionRepresFibers;
      for (i = 0; i < numPatches; i++)
         numFibers += patch[i]->getNumCells();
      
      for (i = 0; i < numReinfLayers; i++)
         numFibers += reinfLayer[i]->getNumReinfBars();
      
	  numHFibers = numSectionRepresHFibers;

      static Vector fiberPosition(2);
      int    matTag;
      
      ID     fibersMaterial(numFibers-numSectionRepresFibers);
      Matrix fibersPosition(2,numFibers-numSectionRepresFibers);
      Vector fibersArea(numFibers-numSectionRepresFibers);

      int  numCells;
      Cell **cell;
    
      k = 0;
      for (i = 0; i < numPatches; i++)
      {
         numCells   = patch[i]->getNumCells();
         matTag = patch[i]->getMaterialID();
         cell = patch[i]->getCells();

         if (cell == 0)
         {
            opserr <<  "WARNING out of run to create fibers\n";
            return TCL_ERROR;
         }    
         
         for (j = 0; j < numCells; j++)
         {
			fibersMaterial(k) = matTag;
            fibersArea(k)     = cell[j]->getArea();
            fiberPosition     = cell[j]->getCentroidPosition();
            fibersPosition(0,k) = fiberPosition(0);
			fibersPosition(1,k) = fiberPosition(1);
            k++;
         }
  
         for (j = 0; j < numCells; j++)
           delete cell[j];
  
         delete [] cell;
      }
         
      ReinfBar *reinfBar;
      int numReinfBars;

      for (i = 0; i < numReinfLayers; i++)
      {
         numReinfBars = reinfLayer[i]->getNumReinfBars();
         reinfBar     = reinfLayer[i]->getReinfBars();
         matTag  = reinfLayer[i]->getMaterialID();
   
         for (j = 0; j < numReinfBars; j++)
         {
	    fibersMaterial(k) = matTag; 
            fibersArea(k)     = reinfBar[j].getArea();
            fiberPosition     = reinfBar[j].getPosition();
     
	    fibersPosition(0,k) = fiberPosition(0);
	    fibersPosition(1,k) = fiberPosition(1);
	
            k++;
         }
         delete [] reinfBar;
      }

      UniaxialMaterial *material;

      int NDM = theTclModelBuilder->getNDM();  


      Fiber **fiber = new Fiber *[numFibers];
      if (fiber == 0) {
	  opserr <<  "WARNING unable to allocate fibers \n";
	  return TCL_ERROR;
      }          
      
      // copy the section repres fibers
      for (i=0; i<numSectionRepresFibers; i++)
	  fiber[i] = sectionRepresFibers[i];

      Fiber **Hfiber = new Fiber *[numHFibers];
      if (Hfiber == 0) {
	  opserr <<  "WARNING unable to allocate Hfibers \n";
	  return TCL_ERROR;
      }          
      
      // copy the section repres fibers
      for (i=0; i<numSectionRepresHFibers; i++)
	  Hfiber[i] = sectionRepresHFibers[i];

	  
	  // creates 2d section      

      if (NDM == 2)     
      {
	 k = 0;
	 for (i = numSectionRepresFibers; i < numFibers; i++)
	 {    
            material = OPS_getUniaxialMaterial(fibersMaterial(k));
            if (material == 0)
            {
               opserr <<  "WARNING invalid material ID for patch\n";
               return TCL_ERROR;
            }   
	    
	    fiber[i] = new UniaxialFiber2d(k, *material, fibersArea(k), fibersPosition(0,k));
            if (fiber[i] == 0) 
            {
               opserr <<  "WARNING unable to allocate fiber \n";
               return TCL_ERROR;
            }    
   
	    k++;
	 }

	 SectionForceDeformation *section = new FiberSection2dInt(secTag, numFibers, fiber, numHFibers, Hfiber, NStrip1, t1, NStrip2, t2, NStrip3, t3);

	 // Delete fibers
	 for (i = 0; i < numFibers; i++)
	   delete fiber[i];

	 for (i = 0; i < numHFibers; i++)
	   delete Hfiber[i];

         if (section == 0)
         {
            opserr <<  "WARNING - cannot construct section\n";
            return TCL_ERROR;
         }
       
         //if (theTclModelBuilder->addSection (*section) < 0) {
	 if (OPS_addSectionForceDeformation(section) != true) {
            opserr <<  "WARNING - cannot add section\n";
            return TCL_ERROR;
         }

      }
      else if (NDM == 3)     
      {

	 static Vector fiberPosition(2);
	 k = 0;
	 for (i = numSectionRepresFibers; i < numFibers; i++)
	 {    
            material = OPS_getUniaxialMaterial(fibersMaterial(k));
            if (material == 0)
            {
               opserr <<  "WARNING invalid material ID for patch\n";
               return TCL_ERROR;
            }   
	    
	    fiberPosition(0) = fibersPosition(0,k);
	    fiberPosition(1) = fibersPosition(1,k);
	    
	    fiber[i] = new UniaxialFiber3d(k, *material, fibersArea(k), fiberPosition);
	    if (fibersArea(k) < 0) opserr << "ERROR: " << fiberPosition(0) << " " << fiberPosition(1) << endln;
            if (fiber[k] == 0) 
            {
               opserr <<  "WARNING unable to allocate fiber \n";
               return TCL_ERROR;
            }    
	    k++;
	 }
	
	 SectionForceDeformation *section = 0;
	 if (isTorsion) {
           ElasticMaterial theGJ(0, GJ);
           //FiberSection3d theFS(0, numFibers, fiber);
           //section = new SectionAggregator(secTag, theFS, theGJ, SECTION_RESPONSE_T);
           section = new FiberSection3d(secTag, numFibers, fiber, &theGJ);
	 }
	 else
	   section = new FiberSection3d(secTag, numFibers, fiber);
   
	 // Delete fibers
	 for (i = 0; i < numFibers; i++)
	   delete fiber[i];

         if (section == 0)
         {
            opserr <<  "WARNING - cannot construct section\n";
            return TCL_ERROR;
         }
       
	 // if (theTclModelBuilder->addSection (*section) < 0) {
	 if (OPS_addSectionForceDeformation(section) != true) {
            opserr <<  "WARNING - cannot add section\n";
            return TCL_ERROR;
         }

      }
      else
      {
         opserr << "WARNING NDM = " << NDM << " is imcompatible with available frame elements\n";
         return TCL_ERROR;
      }

      // Delete fiber array
      delete [] fiber;
   //   delete [] Hfiber;

   }
   else 
   {
      opserr <<  "WARNING section invalid: can only build fiber sections\n";
      return TCL_ERROR;
   }    

   return TCL_OK;
}


int
TclCommand_addUCFiberSection (ClientData clientData, Tcl_Interp *interp, int argc,
				   TCL_Char **argv, TclModelBuilder *theTclModelBuilder)
{
    int secTag;
    
    if (argc < 4) 
	return TCL_ERROR;
    
    if (Tcl_GetInt(interp, argv[2], &secTag) != TCL_OK) {
      opserr <<  "could not read section tag\n";
      return TCL_ERROR;
    }

    currentSectionTag = secTag;

    // first create an empty FiberSection
    int NDM = theTclModelBuilder->getNDM();   // dimension of the structure (1d, 2d, or 3d)

    SectionForceDeformation *section = 0;
    FiberSection2d *section2d =0;
    FiberSection3d *section3d =0;

    if (NDM == 2) {
      section2d = new FiberSection2d(secTag, 0, 0);
      section = section2d;
      //SectionForceDeformation *section = new FiberSection(secTag, 0, 0);
    } else if (NDM == 3) {
      ElasticMaterial theGJ(0, 1e10);
      section3d = new FiberSection3d(secTag, 0, 0, &theGJ);
      section = section3d;
    } 

    if (section == 0) {
      return TCL_ERROR;
    }

    //
    // now parse the output file containing the fiber data, 
    // create fibers and add them to the section
    //

    // open the file
    TCL_Char *fileName = argv[3];
    ifstream theFile;
    theFile.open(fileName, ios::in);
    if (!theFile) {
      opserr << "section UCFiber - could not open file named " << fileName;
      return TCL_ERROR;
    } else {
      int foundStart = 0;
      static char garbage[100];

      // parse through until find start of fiber data
      while (foundStart == 0 && theFile >> garbage) 
	if (strcmp(garbage, "#FIBERS") == 0) 
	  foundStart = 1;

      if (foundStart == 0) {
	theFile.close();
	return TCL_ERROR;
      }

      // parse the fiber data until eof, creating a fiber and adding to section as go
      double ycoord, zcoord, area, prestrain;
      int matTag;
      int fiberCount = 0;
      
      while (theFile >> ycoord >> zcoord >> area >> prestrain >> garbage >> matTag) {

	UniaxialMaterial *theMaterial = OPS_getUniaxialMaterial(matTag);
	if (theMaterial == 0) {
	  opserr << "section UCFiber - no material exists with tag << " << matTag << endln;
	  return TCL_ERROR;
	}
	
	Fiber *theFiber = 0;
	if (NDM == 2) {
	  theFiber = new UniaxialFiber2d(fiberCount++, *theMaterial, area, zcoord);
	  if (theFiber != 0) {
	    section2d->addFiber(*theFiber);
	    delete theFiber;
	  }
	} else {
	  static Vector pos(2);
	  pos(0) = ycoord; pos(1) = zcoord;
	  theFiber = new UniaxialFiber3d(fiberCount++, *theMaterial, area, pos);
	  if (theFiber != 0) {
	    section3d->addFiber(*theFiber);
	    delete theFiber;
	  }
	}   
      }

      // close the file
      theFile.close();
    }

    // finally add the section to our modelbuilder
    //if (theTclModelBuilder->addSection (*section) < 0) {
    if (OPS_addSectionForceDeformation(section) != true) {
      opserr <<  "WARNING - cannot add section\n";
      return TCL_ERROR;
    }


    return TCL_OK;
}

////Changes made by L.Jiang [SIF] 2017
///--Adding Tclcommand for FiberSectionThermal:[BEGIN] by UoE OpenSees Group --///  
int TclCommand_addFiberSectionThermal(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv, TclModelBuilder *theTclModelBuilder)
{
	int secTag;
	int maxNumPatches = 30;
	int maxNumReinfLayers = 30;

	if (argc < 4)
		return TCL_ERROR;

	if (Tcl_GetInt(interp, argv[2], &secTag) != TCL_OK) {
		opserr << "WARNING bad command - want: \nsection fiberSec secTag { \n\tpatch <patch arguments> \n\tlayer <layer arguments> \n}\n";
		return TCL_ERROR;
	}
	currentSectionTag = secTag;
	// create the fiber section representation (with the geometric information) 
	SectionRepres *fiberSectionRepr =
		new FiberSectionRepr(secTag, maxNumPatches, maxNumReinfLayers);

	if (fiberSectionRepr == 0) {
		opserr << "WARNING - ran out of memory to create section representation\n";
		return TCL_ERROR;
	}

	if (theTclModelBuilder->addSectionRepres(*fiberSectionRepr) < 0) {
		opserr << "WARNING - cannot add section representation\n";
		return TCL_ERROR;
	}

	int brace = 3; // Start of recursive parse
	double GJ = 1.0;
	bool isTorsion = false;
	if (strcmp(argv[3], "-GJ") == 0) {
		if (Tcl_GetDouble(interp, argv[4], &GJ) != TCL_OK) {
			opserr << "WARNING invalid GJ";
			return TCL_ERROR;
		}
		isTorsion = true;
		brace = 5;
	}
	// parse the information inside the braces (patches and reinforcing layers)
	if (Tcl_Eval(interp, argv[brace]) != TCL_OK) {
		opserr << "WARNING - error reading information in { } \n";
		return TCL_ERROR;
	}
	// build the fiber section (for analysis)
	if (buildSectionThermal(interp, theTclModelBuilder, secTag, isTorsion, GJ) != TCL_OK) {
		opserr << "WARNING - error constructing the section\n";
		return TCL_ERROR;
	}
	//    currentSectionTag = 0;
	return TCL_OK;
}

///--Adding function for building FiberSectionThermal:[BEGIN] by UoE OpenSees Group --///  
int buildSectionThermal(Tcl_Interp *interp, TclModelBuilder *theTclModelBuilder, int secTag, bool isTorsion, double GJ)
{
	SectionRepres *sectionRepres = theTclModelBuilder->getSectionRepres(secTag);
	if (sectionRepres == 0)
	{
		opserr << "WARNING cannot retrieve section\n";
		return TCL_ERROR;
	}

	if (sectionRepres->getType() == SEC_TAG_FiberSection)
	{
		// build the section
		FiberSectionRepr *fiberSectionRepr = (FiberSectionRepr *)sectionRepres;
		int i, j, k;
		int numFibers;
		int numPatches;
		Patch **patch;

		int  numReinfLayers;
		ReinfLayer **reinfLayer;

		numPatches = fiberSectionRepr->getNumPatches();
		patch = fiberSectionRepr->getPatches();
		numReinfLayers = fiberSectionRepr->getNumReinfLayers();
		reinfLayer = fiberSectionRepr->getReinfLayers();

		int numSectionRepresFibers = fiberSectionRepr->getNumFibers();
		Fiber **sectionRepresFibers = fiberSectionRepr->getFibers();

		numFibers = numSectionRepresFibers;
		for (i = 0; i < numPatches; i++)
			numFibers += patch[i]->getNumCells();

		for (i = 0; i < numReinfLayers; i++)
			numFibers += reinfLayer[i]->getNumReinfBars();
		//opserr << "\nnumFibers: " << numFibers;

		static Vector fiberPosition(2);
		int    matTag;

		ID     fibersMaterial(numFibers - numSectionRepresFibers);
		Matrix fibersPosition(2, numFibers - numSectionRepresFibers);
		Vector fibersArea(numFibers - numSectionRepresFibers);

		int  numCells;
		Cell **cell;

		k = 0;
		for (i = 0; i < numPatches; i++)
		{
			//opserr << "\nPatch :" << i;
			numCells = patch[i]->getNumCells();
			matTag = patch[i]->getMaterialID();
			//opserr << "\nmatTag: " << matTag(k);
			cell = patch[i]->getCells();
			if (cell == 0)
			{
				opserr << "WARNING out of run to create fibers\n";
				return TCL_ERROR;
			}
			//opserr << "\n\tnumCells :" << numCells;
			for (j = 0; j < numCells; j++)
			{
				fibersMaterial(k) = matTag;
				fibersArea(k) = cell[j]->getArea();
				fiberPosition = cell[j]->getCentroidPosition();

				fibersPosition(0, k) = fiberPosition(0);
				fibersPosition(1, k) = fiberPosition(1);
				k++;
			}
			for (j = 0; j < numCells; j++)
				delete cell[j];
			delete[] cell;
		}
		ReinfBar *reinfBar;
		int numReinfBars;
		for (i = 0; i < numReinfLayers; i++)
		{
			numReinfBars = reinfLayer[i]->getNumReinfBars();
			reinfBar = reinfLayer[i]->getReinfBars();
			matTag = reinfLayer[i]->getMaterialID();

			for (j = 0; j < numReinfBars; j++)
			{
				fibersMaterial(k) = matTag;
				fibersArea(k) = reinfBar[j].getArea();
				fiberPosition = reinfBar[j].getPosition();
				fibersPosition(0, k) = fiberPosition(0);
				fibersPosition(1, k) = fiberPosition(1);
				k++;
			}
			delete[] reinfBar;
		}
		UniaxialMaterial *material;

		int NDM = theTclModelBuilder->getNDM();   // dimension of the structure (1d, 2d, or 3d)
		Fiber **fiber = new Fiber *[numFibers];
		if (fiber == 0) {
			opserr << "WARNING unable to allocate fibers \n";
			return TCL_ERROR;
		}
		// copy the section repres fibers
		for (i = 0; i<numSectionRepresFibers; i++)
			fiber[i] = sectionRepresFibers[i];
		// creates 2d section      
		if (NDM == 2)
		{
			k = 0;
			for (i = numSectionRepresFibers; i < numFibers; i++)
			{
				material = OPS_getUniaxialMaterial(fibersMaterial(k));
				if (material == 0)
				{
					opserr << "WARNING invalid material ID for patch\n";
					return TCL_ERROR;
				}

				fiber[i] = new UniaxialFiber2d(k, *material, fibersArea(k), fibersPosition(0, k));
				if (!fiber[i])
				{
					opserr << "WARNING unable to allocate fiber \n";
					return TCL_ERROR;
				}
				//opserr << *fiber[k];
				k++;
			}

			SectionForceDeformation *section = new FiberSection2dThermal(secTag, numFibers, fiber);

			// Delete fibers
			for (i = 0; i < numFibers; i++)
				delete fiber[i];

			if (section == 0)
			{
				opserr << "WARNING - cannot construct section\n";
				return TCL_ERROR;
			}

			if (theTclModelBuilder->addSection(*section) < 0)
			{
				opserr << "WARNING - cannot add section\n";
				return TCL_ERROR;
			}
			//opserr << "section: " << *section;
		}
		else if (NDM == 3)
		{
			static Vector fiberPosition(2);
			k = 0;
			for (i = numSectionRepresFibers; i < numFibers; i++)
			{
				material = OPS_getUniaxialMaterial(fibersMaterial(k));
				if (material == 0)
				{
					opserr << "WARNING invalid material ID for patch\n";
					return TCL_ERROR;
				}
				fiberPosition(0) = fibersPosition(0, k);
				fiberPosition(1) = fibersPosition(1, k);

				fiber[i] = new UniaxialFiber3d(k, *material, fibersArea(k), fiberPosition);
				if (fibersArea(k) < 0) opserr << "ERROR: " << fiberPosition(0) << " " << fiberPosition(1) << endln;
				if (!fiber[k])
				{
					opserr << "WARNING unable to allocate fiber \n";
					return TCL_ERROR;
				}
				k++;
				//opserr << *fiber[k];
			}
			//SectionForceDeformation *section = new FiberSection(secTag, numFibers, fiber);

			SectionForceDeformation *section = 0;
			if (isTorsion)
				section = new FiberSectionGJThermal(secTag, numFibers, fiber, GJ);
			else
				section = new FiberSection3dThermal(secTag, numFibers, fiber);

			// Delete fibers
			for (i = 0; i < numFibers; i++)
				delete fiber[i];
			if (section == 0)
			{
				opserr << "WARNING - cannot construct section\n";
				return TCL_ERROR;
			}

			if (theTclModelBuilder->addSection(*section) < 0)
			{
				opserr << "WARNING - cannot add section\n";
				return TCL_ERROR;
			}
			//opserr << "section: " << *section;
		}
		else
		{
			opserr << "WARNING NDM = " << NDM << " is imcompatible with available frame elements\n";
			return TCL_ERROR;
		}
		// Delete fiber array
		delete[] fiber;
	}
	else
	{
		opserr << "WARNING section invalid: can only build fiber sections\n";
		return TCL_ERROR;
	}
	return TCL_OK;
}
///--Adding function for building FiberSectionThermal:[END] by UoE OpenSees Group --///  
//Changes made by L.Jiang [SIF]


