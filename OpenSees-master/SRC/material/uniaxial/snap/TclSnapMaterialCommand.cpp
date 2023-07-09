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

// $Revision: 1.5 $
// $Date: 2005-01-14 23:45:58 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/snap/TclSnapMaterialCommand.cpp,v $

// Written: Arash Altoontash, Gregory Deierlein,
// Created: Feb 2002
// Modified: Arash June 2004
//
// Description: This file contains the implementation of the
// TclModelBuilder_addSnapMaterial() function. 

#include <Pinching.h>
#include <Clough.h>
#include <CloughHenry.h>
#include <PinchingDamage.h>
#include <CloughDamage.h>
#include <Bilinear.h>
#include <DamageModel.h>

#include <tcl.h>
#include <Vector.h>
#include <string.h>
#include <stdlib.h>

static void printCommand(int argc, TCL_Char **argv)
{
  opserr << "Input command: ";
  for (int i=0; i<argc; i++)
    opserr << argv[i] << " ";
  opserr << endln;
} 

UniaxialMaterial *
TclModelBuilder_addSnapMaterial(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	if (argc < 3) {
		opserr << "WARNING insufficient number of arguments for the Snap material model\n";
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

	if (strcmp(argv[1],"Bilinear") == 0) {
		if (argc < 15) {
			opserr << "WARNING insufficient arguments\n";
			opserr << "Want: uniaxialMaterial Bilinear tag? ..." << endln;
			return 0;
		}
		
		Vector input(12);
		double temp;
		
		for (int i = 3, j = 0; j < 12; i++, j++) {
			if (Tcl_GetDouble(interp, argv[i], &temp) != TCL_OK) {
				opserr << "WARNING invalid input, data " << i << endln;
				printCommand(argc, argv);
				return 0;
			}
			input(j) = temp;
		}
		
		DamageModel *strength;
		if ( (int) input(9) == 0 )
		{
			strength = NULL;
		}
		else
		{
		  strength = OPS_getDamageModel((int) input(9));
		  
		  if (strength == 0) {
		    opserr << "WARNING damage model for strength deterioration not found\n";
		    opserr << "Damage Model: " << input(9);
		    opserr << "\nBinilear material: " << tag << endln;
		    exit (-1);
		  }
		}
		
		DamageModel *stiffness;
		if ( (int) input(10) == 0 )
		{
			stiffness = NULL;
		}
		else
		{
			stiffness = OPS_getDamageModel( (int) input(10) );
			
			if (stiffness == 0) {
				opserr << "WARNING damage model for stiffness deterioration not found\n";
				opserr << "Damage Model: " << input(10);
				opserr << "\nBinilear material: " << tag << endln;
				exit (-1);
			}
		}

		DamageModel *capping;
		if ( (int) input(11) == 0 )
		{
			capping = NULL;
		}
		else
		{
			capping = OPS_getDamageModel( (int) input(11) );
			
			if (capping == 0) {
				opserr << "WARNING damage model for capping deterioration not found\n";
				opserr << "Damage Model: " << input(11);
				opserr << "\nBinilear material: " << tag << endln;
				exit (-1);
			}
		}

		theMaterial = new Bilinear(tag, input, strength, stiffness,capping);
	} 
	
	else if ((strcmp(argv[1],"Clough") == 0) || 
		 (strcmp(argv[1],"clough") == 0) ||
		 (strcmp(argv[1],"CloughHenry") == 0)) {
	    
	    if ( argc < 19 ) {
	      opserr << "WARNING insufficient arguments\n";
	      printCommand(argc,argv);
	      opserr << "Want: uniaxialMaterial Clough tag? 17 args" << endln;
	      return 0;
	    }
	    
	    Vector input(16);
	    double temp;
	    
	    for (int i = 3, j = 0; j < 16; i++, j++) {
	      if (Tcl_GetDouble(interp, argv[i], &temp) != TCL_OK) {
		opserr << "WARNING invalid input, data " << i << endln;
		printCommand(argc, argv);
		return 0;
	      }
	      input(j) = temp;
	    }
	    
	    if ((strcmp(argv[1],"Clough") == 0) || (strcmp(argv[1],"clough") == 0))
	      theMaterial = new Clough(tag, input);
	    else
	      theMaterial = new CloughHenry(tag, input);
	  }
		 
	else if ( strcmp(argv[1],"Clough_Damage") == 0 || strcmp(argv[1],"CloughDamage") == 0 ||
		strcmp(argv[1],"Clough_Damage") == 0 || strcmp(argv[1],"CloughDamage") == 0 ) {
		if ( argc < 15 ) {
			opserr << "WARNING insufficient arguments\n";
			printCommand(argc,argv);
			opserr << "Want: uniaxialMaterial Clough tag? ..." << endln;
			return 0;
		}
		
		Vector input(12);
		double temp;
		
		for (int i = 3, j = 0; j < 12; i++, j++) {
		  if (Tcl_GetDouble(interp, argv[i], &temp) != TCL_OK) {
		    opserr << "WARNING invalid input, data " << i << endln;
		    printCommand(argc, argv);
		    return 0;
		  }
		  input(j) = temp;
		}
		
		DamageModel *strength;
		if ( (int) input(8) == 0 )
		  {
		    strength = NULL;
		  }
		else
		  {
		    strength = OPS_getDamageModel( (int) input(8) );
		    
		    if (strength == 0) {
		      opserr << "WARNING damage model for strength deterioration not found\n";
		      opserr << "Damage Model: " << input(8);
		      opserr << "\nClough material: " << tag << endln;
		      exit (-1);
		    }
		  }
		
		DamageModel *stiffness;
		if ( (int) input(9) == 0 )
		  {
		    stiffness = NULL;
		  }
		else
		  {
		    stiffness = OPS_getDamageModel( (int) input(9) );
		    
		    if (stiffness == 0) {
		      opserr << "WARNING damage model for stiffness deterioration not found\n";
		      opserr << "Damage Model: " << input(9);
		      opserr << "\nClough material: " << tag << endln;
		      exit (-1);
		    }
		  }
		
		DamageModel *accelerated;
		if ( (int) input(10) == 0 )
		  {
		    accelerated = NULL;
		  }
		else
		  {
		    accelerated = OPS_getDamageModel( (int) input(10) );
		    
		    if (accelerated == 0) {
		      opserr << "WARNING damage model for accelerated stiffness deterioration not found\n";
		      opserr << "Damage Model: " << input(10);
		      opserr << "\nClough material: " << tag << endln;
		      exit (-1);
			}
		  }
		
		DamageModel *capping;
		if ( (int) input(11) == 0 )
		  {
		    capping = NULL;
		  }
		else
		  {
		    capping = OPS_getDamageModel( (int) input(11) );
		    
		    if (capping == 0) {
		      opserr << "WARNING damage model for capping deterioration not found\n";
		      opserr << "Damage Model: " << input(11);
		      opserr << "\nClough material: " << tag << endln;
		      exit (-1);
		    }
		  }
		theMaterial = new CloughDamage (tag, input, strength, stiffness, accelerated, capping);
	}
	
	else if ( strcmp(argv[1],"Pinching") == 0 || strcmp(argv[1],"pinching") == 0 ) {
		if ( argc < 22 ) {
			opserr << "WARNING insufficient arguments\n";
			printCommand(argc,argv);
			opserr << "Want: uniaxialMaterial Pinching tag? ..." << endln;
			return 0;
		}
		
		Vector input(19);
		double temp;
		
		for (int i = 3, j = 0; j < 19; i++, j++) {
			if (Tcl_GetDouble(interp, argv[i], &temp) != TCL_OK) {
				opserr << "WARNING invalid input, data " << i << endln;
				printCommand(argc, argv);
				return 0;
			}
			input(j) = temp;
		}
		theMaterial = new Pinching(tag, input);
	}
	
	else if ( strcmp(argv[1],"Pinching_Damage") == 0 || strcmp(argv[1],"pinching_Damage") == 0 ||
		strcmp(argv[1],"PinchingDamage") == 0 || strcmp(argv[1],"pinchingDamage") == 0 ) {
		if ( argc < 18 ) {
			opserr << "WARNING insufficient arguments\n";
			printCommand(argc,argv);
			opserr << "Want: uniaxialMaterial Pinching tag? ..." << endln;
			return 0;
		}
		
		Vector input(15);
		double temp;
		
		for (int i = 3, j = 0; j < 15; i++, j++) {
			if (Tcl_GetDouble(interp, argv[i], &temp) != TCL_OK) {
				opserr << "WARNING invalid input, data " << i << endln;
				printCommand(argc, argv);
				return 0;
			}
			input(j) = temp;
		}
		
		DamageModel *strength;
		if ( (int) input(11) == 0 )
		{
			strength = NULL;
		}
		else
		{
			strength = OPS_getDamageModel( (int) input(11) );
			
			if (strength == 0) {
				opserr << "WARNING damage model for strength deterioration not found\n";
				opserr << "Damage Model: " << input(11);
				opserr << "\nPinching material: " << tag << endln;
				exit (-1);
			}
		}
			
		DamageModel *stiffness;
		if ( (int) input(12) == 0 )
		{
			stiffness = NULL;
		}
		else
		{
			stiffness = OPS_getDamageModel( (int) input(12) );
			if (stiffness == 0) {
				opserr << "WARNING damage model for stiffness deterioration not found\n";
				opserr << "Damage Model: " << input(12);
				opserr << "\nPinching material: " << tag << endln;
				exit (-1);
			}
		}
			
		DamageModel *accelerated;
		if ( (int) input(13) == 0 )
		{
			accelerated = NULL;					
		}
		else
		{
			accelerated = OPS_getDamageModel( (int) input(13) );
			if (accelerated == 0) {
				opserr << "WARNING damage model for accelerated stiffness deterioration not found\n";
				opserr << "Damage Model: " << input(13);
				opserr << "\nPinching material: " << tag << endln;
				exit (-1);
			}
		}

		DamageModel *capping;
		if ( (int) input(14) == 0 )
		{
			capping = NULL;
		}
		else
		{
			capping = OPS_getDamageModel( (int) input(14) );
			
			if (capping == 0) {
				opserr << "WARNING damage model for capping deterioration not found\n";
				opserr << "Damage Model: " << input(14);
				opserr << "\nPinching material: " << tag << endln;
				exit (-1);
			}
		}
		theMaterial = new PinchingDamage(tag, input, strength, stiffness, accelerated, capping);
	}
	
	return theMaterial;
}

