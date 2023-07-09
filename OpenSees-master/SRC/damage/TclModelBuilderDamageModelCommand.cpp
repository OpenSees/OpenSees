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
// $Date: 2004-09-03 15:07:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/damage/TclModelBuilderDamageModelCommand.cpp,v $                                                                        
// Written: Arash Altoontash, Gregory Deierlein
// Created: 08/02
// Revision: AA
//
// Description: This file contains the function invoked when the user invokes
// the DamageModel command in the interpreter. 
//

#include <DamageModel.h>

#include <NormalizedPeak.h>
#include <Kratzig.h>
#include <Mehanny.h>
#include <ParkAng.h>
#include <HystereticEnergy.h>
#include <tcl.h>

#include <Vector.h>
#include <string.h>

static void printCommand(int argc, TCL_Char **argv)
{
  opserr << "Input command: ";
  for (int i=0; i<argc; i++)
    opserr << argv[i] << " ";
  opserr << endln;
} 

int
TclModelBuilderDamageModelCommand (ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
				   
{
  
  // Make sure there is a minimum number of arguments
  if (argc < 3) {
    opserr << "WARNING insufficient number of damage model arguments\n";
    opserr << "Want: damageModel type? tag? <specific material args>" << endln;
    return TCL_ERROR;
  }
  
  // Pointer to a damage model that will be added to the model builder
  DamageModel *theDamage = 0;
  
  
  // Check argv[2] for damage model type
  
  // Mehanny damage model
  if (strcmp(argv[1],"Mehanny") == 0) {
    if ( argc != 8 && argc != 10 && argc != 12 ) {
      opserr << "WARNING invalid number of arguments\n";
      printCommand(argc,argv);
      opserr << "Want: damageModel Mehanny tag? Alpha? Beta? Gamma? ultimatePosDisp? ultimateNegDisp? AbsTol? RelTol?" << endln;
      return TCL_ERROR;
    }    
    
    int tag;
    double alpha,beta,gamma,ultimatePosDisp, ultimateNegDisp, abstol, reltol, posmodifier,negmodifier;
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid tag\n";
      opserr << "Damage model Mehanny " << endln;
      return TCL_ERROR;		
    }
    
    if (Tcl_GetDouble(interp, argv[3], &alpha) != TCL_OK) {
      opserr << "WARNING invalid Alpha\n";
      opserr << "Damage model Mehanny : " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp, argv[4], &beta) != TCL_OK) {
      opserr << "WARNING invalid Beta\n";
      opserr << "Damage Mehanny : " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp, argv[5], &gamma) != TCL_OK) {
      opserr << "WARNING invalid Gamma\n";
      opserr << "Damage Mehanny : " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp, argv[6], &ultimatePosDisp) != TCL_OK) {
      opserr << "WARNING invalid ultimatePosDisp\n";
      opserr << "Damage Mehanny : " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp,argv[7], &ultimateNegDisp) != TCL_OK) {
      opserr << "WARNING invalid ultimateNegDisp\n";
      opserr << "Damage Mehanny : " << tag << endln;
      return TCL_ERROR;
    }
    
    if ( argc > 8 )
      {
	if (Tcl_GetDouble(interp, argv[8], &abstol) != TCL_OK) {
	  opserr << "WARNING invalid AbsTol\n";
	  opserr << "Damage Mehanny : " << tag << endln;
	  return TCL_ERROR;	
	}
	
	if (Tcl_GetDouble(interp,argv[9], &reltol) != TCL_OK) {
	  opserr << "WARNING invalid RelTol\n";
	  opserr << "Damage Mehanny : " << tag << endln;
	  return TCL_ERROR;
	}
	
	
	if ( argc == 12 )
	  {
	    
	    if (Tcl_GetDouble(interp,argv[10], &posmodifier) != TCL_OK) {
	      opserr << "WARNING invalid posmodifier\n";
	      opserr << "Damage Mehanny : " << tag << endln;
	      return TCL_ERROR;
	    }
	    
	    if (Tcl_GetDouble(interp,argv[11], &negmodifier) != TCL_OK) {
	      opserr << "WARNING invalid negmodifier\n";
	      opserr << "Damage Mehanny : " << tag << endln;
	      return TCL_ERROR;
	    }
	  }	else
	    {
	      posmodifier = negmodifier = 1.0;
	    }
      } else
	{
	  abstol = reltol = 0.0;
	  posmodifier = negmodifier = 1.0;		
	}
    
    // Parsing was successful, allocate the damage model
    theDamage = new Mehanny(tag, alpha , beta, gamma , ultimatePosDisp, ultimateNegDisp, abstol, reltol, posmodifier, negmodifier);
  }
  // end of Mehanny damage model
  
  // Kratzig damage model
  else if (strcmp(argv[1],"Kratzig") == 0) {
    if ( argc != 7 && argc != 8) {
      opserr << "WARNING invalid number of arguments\n";
      printCommand(argc,argv);
      opserr << "Want: damageModel Kratzig tag? Alpha? Beta? <Gamma?> ultimatePosDisp? ultimateNegDisp?" << endln;
      return TCL_ERROR;
    }    
    
    int tag;
    double ultimatePosDisp,ultimateNegDisp;
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid tag\n";
      opserr << "Damage model Kratzig " << endln;
      return TCL_ERROR;		
    }
    
    if (Tcl_GetDouble(interp, argv[3], &ultimatePosDisp) != TCL_OK) {
      opserr << "WARNING invalid ultimatePosDisp\n";
      opserr << "Damage Kratzig : " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp,argv[4], &ultimateNegDisp) != TCL_OK) {
      opserr << "WARNING invalid ultimateNegDisp\n";
      opserr << "Damage Kratzig : " << tag << endln;
      return TCL_ERROR;
    }
    
    // Parsing was successful, allocate the damage model
    theDamage = new Kratzig(tag, ultimatePosDisp, ultimateNegDisp);       
  }
  // end of Kratzig damage model
  
  
  // NormalizedPeak damage model
  else if (strcmp(argv[1],"NormalizedPeak") == 0) {
    if (argc != 6 ) {
      opserr << "WARNING invalid number of arguments\n";
      printCommand(argc,argv);
      opserr << "Want: damageModel NormalizedPeak tag? maxValue? minValue? response?" << endln;
      return TCL_ERROR;
    }    
    
    int tag;
    double maxVal, minVal;
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid tag\n";
      opserr << "Damage model NormalizedPeak " << endln;
      return TCL_ERROR;		
    }
    
    if (Tcl_GetDouble(interp, argv[3], &maxVal) != TCL_OK) {
      opserr << "WARNING invalid maxVal\n";
      opserr << "Damage model NormalizedPeak : " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp, argv[4], &minVal) != TCL_OK) {
      opserr << "WARNING invalid Beta\n";
      opserr << "Damage NormalizedPeak : " << tag << endln;
      return TCL_ERROR;	
    }
    
    
    const char *responsetype = argv[5];
    
    // Parsing was successful, allocate the damage model
    theDamage = new NormalizedPeak( tag, maxVal, minVal, responsetype);
  }
  // end of NormalizedPeak damage model
  
  
  // HystereticEnergy damage model
  else if (strcmp(argv[1],"HystereticEnergy") == 0) {
    if (argc != 5 ) {
      opserr << "WARNING invalid number of arguments\n";
      printCommand(argc,argv);
      opserr << "Want: damageModel HystereticEnergy tag?  Etotal? Cpower?" << endln;
      return TCL_ERROR;
    }    
    
    int tag;
    double Etot,Cpow;
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid tag\n";
      opserr << "Damage model HystereticEnergy " << endln;
      return TCL_ERROR;		
    }
    
    if (Tcl_GetDouble(interp, argv[3], &Etot) != TCL_OK) {
      opserr << "WARNING invalid Total energy\n";
      opserr << "Damage model HystereticEnergy : " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp, argv[4], &Cpow) != TCL_OK) {
      opserr << "WARNING invalid Constant exponent\n";
      opserr << "Damage HystereticEnergy : " << tag << endln;
      return TCL_ERROR;	
    }
    
    // Parsing was successful, allocate the damage model
    theDamage = new HystereticEnergy( tag, Etot, Cpow );
  }
  // end of NormalizedPeak damage model
  
  
  // ParkAng damage model
  else if (strcmp(argv[1],"ParkAng") == 0) {
    if (argc != 6 ) {
      opserr << "WARNING invalid number of arguments\n";
      printCommand(argc,argv);
      opserr << "Want: damageModel ParkAng tag?  deltaU?  beta?  sigmaY?" << endln;
      return TCL_ERROR;
    }    
    
    int tag;
    double deltaU, beta, sigmaY;
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid tag\n";
      opserr << "Damage model NormalizedPeak " << endln;
      return TCL_ERROR;		
    }
    
    if (Tcl_GetDouble(interp, argv[3], &deltaU) != TCL_OK) {
      opserr << "WARNING invalid deltaU\n";
      opserr << "Damage model ParkAng : " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp, argv[4], &beta) != TCL_OK) {
      opserr << "WARNING invalid beta\n";
      opserr << "Damage ParkAng : " << tag << endln;
      return TCL_ERROR;	
    }
    
    if (Tcl_GetDouble(interp, argv[5], &sigmaY) != TCL_OK) {
      opserr << "WARNING invalid sigmaY\n";
      opserr << "Damage ParkAng : " << tag << endln;
      return TCL_ERROR;	
    }
    
    
    // Parsing was successful, allocate the damage model
    theDamage = new ParkAng( tag, deltaU,  beta,  sigmaY );
  }
  // end of NormalizedPeak damage model
  
  else {
    opserr << "WARNING the damage model specified does not exist " << argv[1] << endln;
    return TCL_ERROR;
  }
  
  
  // check if the damage model is constructed or not
  if (theDamage == 0) {
    opserr << "WARNING could not create DamageModel, out of memory " << argv[1] << endln;
    return TCL_ERROR;
  }
  
  // Now add the damage model to the modelBuilder
  if (OPS_addDamageModel(theDamage) == false) {
    opserr << "WARNING could not add DamageModel to the domain\n";
    opserr << *theDamage << endln;
    delete theDamage; // invoke the material objects destructor, otherwise mem leak
    return TCL_ERROR;
  }
  
  return TCL_OK;
}
