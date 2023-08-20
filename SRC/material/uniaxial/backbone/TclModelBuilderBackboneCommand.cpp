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
                                                                        
// $Revision: 1.3 $
// $Date: 2008-11-24 17:12:12 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/backbone/TclModelBuilderBackboneCommand.cpp,v $

// Written: MHS
// Created: Aug 2000
//
// Description: This file contains the parsing routines for the
// TCL hystereticBackbone command.

#include <OPS_Globals.h>
#include <UniaxialMaterial.h>

#include <tcl.h>
#include <elementAPI.h>
extern "C" int         OPS_ResetInputNoBuilder(ClientData clientData, Tcl_Interp * interp, int cArg, int mArg, TCL_Char * *argv, Domain * domain);

#include <string.h>

#include <HystereticBackbone.h>

extern void *OPS_ArctangentBackbone(void);
extern void *OPS_ManderBackbone(void);
extern void *OPS_TrilinearBackbone(void);
extern void *OPS_BilinearBackbone(void);
extern void *OPS_MultilinearBackbone(void);
extern void* OPS_MaterialBackbone();
extern void* OPS_ReeseStiffClayBelowWS();
extern void* OPS_ReeseStiffClayAboveWS();
extern void* OPS_VuggyLimestone();
extern void* OPS_CementedSoil();
extern void* OPS_WeakRock();
extern void* OPS_LiquefiedSand();
extern void* OPS_RaynorBackbone();
extern void* OPS_ReeseSandBackbone();
extern void* OPS_ReeseSoftClayBackbone();
extern void* OPS_CappedBackbone();
extern void* OPS_LinearCappedBackbone();

#include <packages.h>

static void printCommand(int argc, TCL_Char **argv)
{
  opserr << "Input command: ";
  for (int i=0; i<argc; i++)
    opserr << argv[i] << " ";
  opserr << endln;
} 

int
TclModelBuilderHystereticBackboneCommand(ClientData clientData,
					 Tcl_Interp *interp,
					 int argc, TCL_Char **argv, Domain *theDomain)
{
  if (argc < 3) {
    opserr << "WARNING insufficient number of hystereticBackbone arguments\n";
    opserr << "Want: hystereticBackbone type? tag? <specific hystereticBackbone args>" << endln;
    return TCL_ERROR;
  }
  
  OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, theDomain);	  

  // Pointer to a hysteretic backbone that will be added to the model builder
  HystereticBackbone *theBackbone = 0;
  
  // Check argv[1] for backbone type
  if (strcmp(argv[1],"Bilinear") == 0) {
    void *theBB = OPS_BilinearBackbone();
    if (theBB != 0)
      theBackbone = (HystereticBackbone *)theBB;
    else
      return TCL_ERROR;
  }
  
  else if (strcmp(argv[1],"Trilinear") == 0) {
    void *theBB = OPS_TrilinearBackbone();
    if (theBB != 0)
      theBackbone = (HystereticBackbone *)theBB;
    else
      return TCL_ERROR;
  }
  
  else if (strcmp(argv[1],"Multilinear") == 0) {
    void *theBB = OPS_MultilinearBackbone();
    if (theBB != 0)
      theBackbone = (HystereticBackbone *)theBB;
    else
      return TCL_ERROR;
  }
  
  else if (strcmp(argv[1],"Arctangent") == 0) {
    void *theBB = OPS_ArctangentBackbone();
    if (theBB != 0)
      theBackbone = (HystereticBackbone *)theBB;
    else
      return TCL_ERROR;
  }
  
  else if (strcmp(argv[1],"Material") == 0) {
    void *theBB = OPS_MaterialBackbone();
    if (theBB != 0)
      theBackbone = (HystereticBackbone *)theBB;
    else
      return TCL_ERROR;
  }
  
  else if (strcmp(argv[1],"ReeseStiffClayBelowWS") == 0) {
    void *theBB = OPS_ReeseStiffClayBelowWS();
    if (theBB != 0)
      theBackbone = (HystereticBackbone *)theBB;
    else
      return TCL_ERROR;
  }
  
  else if (strcmp(argv[1],"ReeseStiffClayAboveWS") == 0) {
    void *theBB = OPS_ReeseStiffClayAboveWS();
    if (theBB != 0)
      theBackbone = (HystereticBackbone *)theBB;
    else
      return TCL_ERROR;
  }
  
  else if (strcmp(argv[1],"VuggyLimestone") == 0) {
    void *theBB = OPS_VuggyLimestone();
    if (theBB != 0)
      theBackbone = (HystereticBackbone *)theBB;
    else
      return TCL_ERROR;
  }
  
  else if (strcmp(argv[1],"CementedSoil") == 0) {
    void *theBB = OPS_CementedSoil();
    if (theBB != 0)
      theBackbone = (HystereticBackbone *)theBB;
    else
      return TCL_ERROR;
  }
  
  else if (strcmp(argv[1],"WeakRock") == 0) {
    void *theBB = OPS_WeakRock();
    if (theBB != 0)
      theBackbone = (HystereticBackbone *)theBB;
    else
      return TCL_ERROR;
  }
  
  else if (strcmp(argv[1],"LiquefiedSand") == 0) {
    void *theBB = OPS_LiquefiedSand();
    if (theBB != 0)
      theBackbone = (HystereticBackbone *)theBB;
    else
      return TCL_ERROR;
  }
  
  else if (strcmp(argv[1],"ReeseSoftClay") == 0) {
    void *theBB = OPS_ReeseSoftClayBackbone();
    if (theBB != 0)
      theBackbone = (HystereticBackbone *)theBB;
    else
      return TCL_ERROR;
  }
  
  else if (strcmp(argv[1],"ReeseSand") == 0) {
    void *theBB = OPS_ReeseSandBackbone();
    if (theBB != 0)
      theBackbone = (HystereticBackbone *)theBB;
    else
      return TCL_ERROR;
  }

  else if (strcmp(argv[1],"Mander") == 0) {
    void *theBB = OPS_ManderBackbone();
    if (theBB != 0)
      theBackbone = (HystereticBackbone *)theBB;
    else
      return TCL_ERROR;
  }
  
  else if (strcmp(argv[1],"Raynor") == 0) {
    void *theBB = OPS_RaynorBackbone();
    if (theBB != 0)
      theBackbone = (HystereticBackbone *)theBB;
    else
      return TCL_ERROR;
  }

  else if (strcmp(argv[1],"Capped") == 0) {
    void *theBB = OPS_CappedBackbone();
    if (theBB != 0)
      theBackbone = (HystereticBackbone *)theBB;
    else
      return TCL_ERROR;
  }
  
  else if (strcmp(argv[1],"LinearCapped") == 0) {
    void *theBB = OPS_LinearCappedBackbone();
    if (theBB != 0)
      theBackbone = (HystereticBackbone *)theBB;
    else
      return TCL_ERROR;
  }
  
  else  {
    opserr << "WARNING unknown type of hystereticBackbone: " << argv[1];
    opserr << "\nValid types: Bilinear, Trilinear, Arctangent," << endln;
    opserr << "\tCapped, LinearCapped, Material" << endln;
    return TCL_ERROR;
  }
  
  // Ensure we have created the Material, out of memory if got here and no material
  if (theBackbone == 0) {
    opserr << "WARNING ran out of memory creating hystereticBackbone\n";
    opserr << argv[1] << endln;
    return TCL_ERROR;
  }
  
  // Now add the material to the modelBuilder
  if (OPS_addHystereticBackbone(theBackbone) == false) {
    opserr << "WARNING could not add hystereticBackbone to the domain\n";
    opserr << *theBackbone << endln;
    delete theBackbone; // invoke the material objects destructor, otherwise mem leak
    return TCL_ERROR;
  }
  
  return TCL_OK;
}
