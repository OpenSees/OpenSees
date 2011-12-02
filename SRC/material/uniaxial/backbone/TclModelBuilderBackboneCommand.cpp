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

#include <TclModelBuilder.h>

#include <TrilinearBackbone.h>
//#include <MultilinearBackbone.h>
#include <ArctangentBackbone.h>
//#include <CappedBackbone.h>
//#include <MaterialBackbone.h>
//#include <LinearCappedBackbone.h>
#include <ReeseSoftClayBackbone.h>
#include <ReeseStiffClayBelowWS.h>
#include <ReeseSandBackbone.h>
#include <ManderBackbone.h>
#include <RaynorBackbone.h>

#include <string.h>

static void printCommand(int argc, TCL_Char **argv)
{
  opserr << "Input command: ";
  for (int i=0; i<argc; i++)
    opserr << argv[i] << " ";
  opserr << endln;
} 

int
TclModelBuilderHystereticBackboneCommand(ClientData clienData,
					 Tcl_Interp *interp,
					 int argc, TCL_Char **argv,
					 TclModelBuilder *theTclBuilder)
{
  if (argc < 2) {
    opserr << "WARNING insufficient number of hystereticBackbone arguments\n";
    opserr << "Want: hystereticBackbone type? tag? <specific hystereticBackbone args>" << endln;
    return TCL_ERROR;
  }
  
  // Pointer to a hysteretic backbone that will be added to the model builder
  HystereticBackbone *theBackbone = 0;
  
  // Check argv[1] for backbone type
  if (strcmp(argv[1],"Bilinear") == 0) {
    if (argc < 7) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc,argv);
      opserr << "Want: hystereticBackbone Bilinear tag? e1? s1? e2? s2?" << endln;
      return TCL_ERROR;
    }
    
    int tag;
    double e1, e2, s1, s2;
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Bilinear tag" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[3], &e1) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Bilinear e1" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[4], &s1) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Bilinear s1" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[5], &e2) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Bilinear e2" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[6], &s2) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Bilinear s2" << endln;
      return TCL_ERROR;
    }
    
    // Parsing was successful, allocate the material
    theBackbone = new TrilinearBackbone (tag, e1, s1, e2, s2);
  }
  
  else if (strcmp(argv[1],"Trilinear") == 0) {
    if (argc < 9) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc,argv);
      opserr << "Want: hystereticBackbone Trilinear tag? e1? s1? e2? s2? e3? s3?" << endln;
      return TCL_ERROR;
    }
    
    int tag;
    double e1, e2, e3, s1, s2, s3;
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Trilinear tag" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[3], &e1) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Trilinear e1" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[4], &s1) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Trilinear s1" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[5], &e2) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Trilinear e2" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[6], &s2) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Trilinear s2" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[7], &e3) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Trilinear e3" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[8], &s3) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Trilinear s3" << endln;
      return TCL_ERROR;
    }
    
    // Parsing was successful, allocate the material
    theBackbone = new TrilinearBackbone (tag, e1, s1, e2, s2, e3, s3);
  }
  
  else if (strcmp(argv[1],"Multilinear") == 0) {
    if (argc < 5 || (argc-3)%2 != 0) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc,argv);
      opserr << "Want: hystereticBackbone Multilinear tag? e1? s1? e2? s2? ..." << endln;
      return TCL_ERROR;
    }
    
    int tag;
    int numPoints = (argc-3)/2;
    Vector e(numPoints);
    Vector s(numPoints);
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Trilinear tag" << endln;
      return TCL_ERROR;
    }
    
    double def,force;
    for (int i = 3, j = 0; i < argc; i++, j++) {
      if (Tcl_GetDouble(interp, argv[i++], &def) != TCL_OK) {
	opserr << "WARNING invalid hystereticBackbone Trilinear e" << endln;
	return TCL_ERROR;
      }
      e(j) = def;
      
      if (Tcl_GetDouble(interp, argv[i], &force) != TCL_OK) {
	opserr << "WARNING invalid hystereticBackbone Trilinear s" << endln;
	return TCL_ERROR;
      }
      s(j) = force;
    }
    
    // Parsing was successful, allocate the material
    //theBackbone = new MultilinearBackbone (tag, numPoints, e, s);
    theBackbone = 0;
  }
  
  else if (strcmp(argv[1],"Arctangent") == 0) {
    if (argc < 6) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc,argv);
      opserr << "Want: hystereticBackbone Arctangent tag? K1? gammaY? alpha?" << endln;
      return TCL_ERROR;
    }
    
    int tag;
    double K1, gammaY, alpha;
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Arctangent tag" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[3], &K1) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Arctangent K1" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[4], &gammaY) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Arctangent gammaY" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[5], &alpha) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Arctangent alpha" << endln;
      return TCL_ERROR;
    }
    
    theBackbone = new ArctangentBackbone (tag, K1, gammaY, alpha);
  }
  
  else if (strcmp(argv[1],"ReeseSoftClay") == 0) {
    if (argc < 6) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc,argv);
      opserr << "Want: hystereticBackbone ReeseSoftClay tag? pu? y50? n?" << endln;
      return TCL_ERROR;
    }
    
    int tag;
    double pu, y50, n;
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone ReeseSoftClay tag" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[3], &pu) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone ReeseSoftClay pu" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[4], &y50) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone ReeseSoftClay y50" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &n) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone ReeseSoftClay n" << endln;
      return TCL_ERROR;
    }
    
    theBackbone = new ReeseSoftClayBackbone(tag, pu, y50, n);
  }
  
  else if (strcmp(argv[1],"ReeseSand") == 0) {
    if (argc < 8) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc,argv);
      opserr << "Want: hystereticBackbone ReeseSand tag? kx? ym? pm? yu? pu?" << endln;
      return TCL_ERROR;
    }
    
    int tag;
    double kx, ym, pm, yu, pu;
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone ReeseSand tag" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[3], &kx) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone ReeseSand kx" << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[4], &ym) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone ReeseSand ym" << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[5], &pm) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone ReeseSand pm" << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[6], &yu) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone ReeseSand yu" << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[7], &pu) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone ReeseSand pu" << endln;
      return TCL_ERROR;
    }
    
    theBackbone = new ReeseSandBackbone(tag, kx, ym, pm, yu, pu);
  }

  else if (strcmp(argv[1],"ReeseStiffClayBelowWS") == 0)
  {
    if (argc < 7)
    {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc,argv);
      opserr << "Want: hystereticBackbone ReeseStiffClayBelowWS tag? Esi? y50? As?  Pc?" <<endln;
      return TCL_ERROR;
    }

    int tag;
    double Esi, y50, As, Pc;

    if( Tcl_GetInt(interp, argv[2], &tag ) != TCL_OK)
    {
      opserr << "WARNING invalid hystereticBackbone ReeseStiffClayBelowWS tag" << endln;
      return TCL_ERROR;
    }
    if( Tcl_GetDouble(interp, argv[3], &Esi ) != TCL_OK)
    {
      opserr << "WARNING invalid hystereticBackbone ReeseStiffClayBelowWS Esi" << endln;
      return TCL_ERROR;
    }
    if( Tcl_GetDouble(interp, argv[4], &y50 ) != TCL_OK)
    {
      opserr << "WARNING invalid hystereticBackbone ReeseStiffClayBelowWS y50" << endln;
      return TCL_ERROR;
    }
    if( Tcl_GetDouble(interp, argv[5], &As ) != TCL_OK)
    {
      opserr << "WARNING invalid hystereticBackbone ReeseStiffClayBelowWS As" << endln;
      return TCL_ERROR;
    }
    if( Tcl_GetDouble(interp, argv[6], &Pc ) != TCL_OK)
    {
      opserr << "WARNING invalid hystereticBackbone ReeseStiffClayBelowWS Pc" << endln;
      return TCL_ERROR;
    }

    theBackbone = new ReeseStiffClayBelowWS (tag, Esi, y50, As, Pc);

  }

  else if (strcmp(argv[1],"Mander") == 0) {
    if (argc < 6) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc,argv);
      opserr << "Want: hystereticBackbone Mander tag? fc? epsc? Ec?" << endln;
      return TCL_ERROR;
    }
    
    int tag;
    double fc, epsc, Ec;
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Mander tag" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[3], &fc) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Mander fc" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[4], &epsc) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Mander epsc" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[5], &Ec) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Mander Ec" << endln;
      return TCL_ERROR;
    }
    
    theBackbone = new ManderBackbone (tag, fc, epsc, Ec);
  }
  
  else if (strcmp(argv[1],"Raynor") == 0)
  {
    if (argc < 10)
    {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc,argv);
      opserr << "Want: hystereticBackbone Raynor tag? Es? fy? fsu? Epsilonsh? Epsilonsm? C1? Ey?" <<endln;
      return TCL_ERROR;
    }

    int tag;
    double Es, fy, fsu, Epsilonsh, Epsilonsm, C1, Ey;

    if( Tcl_GetInt(interp, argv[2], &tag ) != TCL_OK)
    {
      opserr << "WARNING invalid hystereticBackbone Raynor tag" << endln;
      return TCL_ERROR;
    }
    if( Tcl_GetDouble(interp, argv[3], &Es ) != TCL_OK)
    {
      opserr << "WARNING invalid hystereticBackbone Raynor Es" << endln;
      return TCL_ERROR;
    }
    if( Tcl_GetDouble(interp, argv[4], &fy ) != TCL_OK)
    {
      opserr << "WARNING invalid hystereticBackbone Raynor fy" << endln;
      return TCL_ERROR;
    }
    if( Tcl_GetDouble(interp, argv[5], &fsu ) != TCL_OK)
    {
      opserr << "WARNING invalid hystereticBackbone Raynor fsu" << endln;
      return TCL_ERROR;
    }
    if( Tcl_GetDouble(interp, argv[6], &Epsilonsh ) != TCL_OK)
    {
      opserr << "WARNING invalid hystereticBackbone Raynor Epsilonsh" << endln;
      return TCL_ERROR;
    }
    if( Tcl_GetDouble(interp, argv[7], &Epsilonsm ) != TCL_OK)
    {
      opserr << "WARNING invalid hystereticBackbone Raynor Epsilonsm" << endln;
      return TCL_ERROR;
    }
    if( Tcl_GetDouble(interp, argv[8], &C1 ) != TCL_OK)
    {
      opserr << "WARNING invalid hystereticBackbone Raynor fy" << endln;
      return TCL_ERROR;
    }
    if( Tcl_GetDouble(interp, argv[9], &Ey ) != TCL_OK)
    {
      opserr << "WARNING invalid hystereticBackbone Raynor fsu" << endln;
      return TCL_ERROR;
    }

    theBackbone = new RaynorBackbone (tag, Es, fy, fsu, Epsilonsh, Epsilonsm, C1, Ey);

  }

  else if (strcmp(argv[1],"Capped") == 0) {
    if (argc < 5) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc,argv);
      opserr << "Want: hystereticBackbone Capped tag? hystereticBackboneTag? capTag?" << endln;
      return TCL_ERROR;
    }
    
    int tag, bTag, cTag;
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Capped tag" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetInt(interp, argv[3], &bTag) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Capped hystereticBackboneTag" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetInt(interp, argv[4], &cTag) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone Capped capTag" << endln;
      return TCL_ERROR;
    }
    
    HystereticBackbone *hystereticBackbone = 0;
    
    //hystereticBackbone = theTclBuilder->getHystereticBackbone(bTag);
    
    if (hystereticBackbone == 0) {
      opserr << "WARNING hystereticBackbone does not exist\n";
      opserr << "hystereticBackbone: " << bTag; 
      opserr << "\nhystereticBackbone Capped: " << tag << endln;
      return TCL_ERROR;
    }
    
    HystereticBackbone *cap = 0;
    
    //cap = theTclBuilder->getHystereticBackbone(cTag);
    
    if (cap == 0) {
      opserr << "WARNING cap does not exist\n";
      opserr << "cap: " << cTag; 
      opserr << "\nhystereticBackbone Capped: " << tag << endln;
      return TCL_ERROR;
    }
    
    //theBackbone = new CappedBackbone (tag, *hystereticBackbone, *cap);
    theBackbone = 0;
  }
  
  else if (strcmp(argv[1],"LinearCapped") == 0) {
    if (argc < 7) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc,argv);
      opserr << "Want: hystereticBackbone LinearCapped tag? backboneTag? eCap? E? sRes?" << endln;
      return TCL_ERROR;
    }
    
    int tag, bTag;
    double e, E, s;
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone LinearCapped tag" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetInt(interp, argv[3], &bTag) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone LinearCapped backboneTag" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[4], &e) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone LinearCapped eCap" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[5], &E) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone LinearCapped E" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[6], &s) != TCL_OK) {
      opserr << "WARNING invalid hystereticBackbone LinearCapped sRes" << endln;
      return TCL_ERROR;
    }
    
    HystereticBackbone *backbone = 0;
    //backbone = theTclBuilder->getHystereticBackbone(bTag);
    
    if (backbone == 0) {
      opserr << "WARNING hystereticBackbone does not exist\n";
      opserr << "hystereticBackbone: " << bTag; 
      opserr << "\nhystereticBackbone LinearCapped: " << tag << endln;
      return TCL_ERROR;
    }
    
    //theBackbone = new LinearCappedBackbone (tag, *backbone, e, E, s);
    theBackbone = 0;
  }
  
  else if (strcmp(argv[1],"Material") == 0) {
    if (argc < 4)
      {
	opserr << "WARNING insufficient arguments\n";
	printCommand(argc,argv);
	opserr << "Want: hystereticBackbone Material tag? matTag?" << endln;
	return TCL_ERROR;
      }
    
    int tag, matTag;
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid tag\n";
      opserr << "hystereticBackbone Material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
      opserr << "WARNING invalid matTag\n";
      opserr << "hystereticBackbone Material: " << tag << endln;
      return TCL_ERROR;
    }
    
    UniaxialMaterial *material = theTclBuilder->getUniaxialMaterial(matTag);
    
    if (material == 0) {
      opserr << "WARNING material does not exist\n";
      opserr << "material: " << matTag; 
      opserr << "\nhystereticBackbone Material: " << tag << endln;
      return TCL_ERROR;
    }
    
    //theBackbone = new MaterialBackbone(tag, *material);
    theBackbone = 0;
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
  if (theTclBuilder->addHystereticBackbone(*theBackbone) < 0) {
    opserr << "WARNING could not add hystereticBackbone to the domain\n";
    opserr << *theBackbone << endln;
    delete theBackbone; // invoke the material objects destructor, otherwise mem leak
    return TCL_ERROR;
  }
  
  return TCL_OK;
}
