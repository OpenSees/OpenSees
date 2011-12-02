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
                                                                        
// $Revision: 1.39 $
// $Date: 2006-01-17 19:32:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/TclElementCommands.cpp,v $
                                                                        
// Written: fmk 
// Created: 07/99
// Revision: A
//
// Description: This file contains the implementation of the TclElementCommands.
// The file contains the routine TclElementCommands which is invoked by the
// TclModelBuilder.
//
// What: "@(#) TclModelBuilder.C, revA"

#include <stdlib.h>
#include <string.h>
#include <OPS_Stream.h>
#include <Domain.h>

#include <ElasticBeam2d.h>
#include <ElasticBeam3d.h>

//Zhaohui Yang (UCD)
#include <EightNodeBrick.h>
#include <TwentyNodeBrick.h>

#include <CrdTransf2d.h>
#include <CrdTransf3d.h>

#include <TclModelBuilder.h>
#include <packages.h>

//
// SOME STATIC POINTERS USED IN THE FUNCTIONS INVOKED BY THE INTERPRETER
//


typedef struct elementPackageCommand {
  char *funcName;
  int (*funcPtr)(ClientData clientData, Tcl_Interp *interp,  int argc, 
		 TCL_Char **argv, Domain*, TclModelBuilder *); 
  struct elementPackageCommand *next;
} ElementPackageCommand;

static ElementPackageCommand *theElementPackageCommands = NULL;

extern void printCommand(int argc, TCL_Char **argv);

// 
// THE PROTOTYPES OF THE FUNCTIONS INVOKED BY THE INTERPRETER
//

extern int
TclModelBuilder_addFeapTruss(ClientData clientData, Tcl_Interp *interp,  int argc, 
			     TCL_Char **argv, Domain*, TclModelBuilder *, int argStart);

extern int
TclModelBuilder_addTruss(ClientData clientData, Tcl_Interp *interp,  int argc, 
			 TCL_Char **argv, Domain*, TclModelBuilder *, int argStart); 

extern int
TclModelBuilder_addElasticBeam(ClientData clientData, Tcl_Interp *interp,  int argc, 
			       TCL_Char **argv, Domain*, TclModelBuilder *, int argStart);

extern int
TclModelBuilder_addBrick(ClientData clientData, Tcl_Interp *interp,
			 int argc, TCL_Char **argv, Domain*, 
			 TclModelBuilder *, int argStart);

extern int
TclModelBuilder_addShellMITC4(ClientData clientData, Tcl_Interp *interp,
			      int argc, TCL_Char **argv, Domain*, 
			      TclModelBuilder *, int argStart);

extern int 
TclModelBuilder_addConstantPressureVolumeQuad(ClientData, Tcl_Interp *, int, TCL_Char **,
					      Domain*, TclModelBuilder *);

extern int 
TclModelBuilder_addJoint2D(ClientData, Tcl_Interp *, int, TCL_Char **,
			   Domain*, TclModelBuilder *);

extern int 
TclModelBuilder_addJoint3D(ClientData, Tcl_Interp *, int, TCL_Char **,
			   Domain*, TclModelBuilder *);

extern int 
TclModelBuilder_addEnhancedQuad(ClientData, Tcl_Interp *, int, TCL_Char **,
				Domain*, TclModelBuilder *);

extern int 
TclModelBuilder_addNineNodeMixedQuad(ClientData, Tcl_Interp *, int, TCL_Char **,
				     Domain*, TclModelBuilder *);


// GLF			       
extern int 
TclModelBuilder_addZeroLength(ClientData, Tcl_Interp *, int, TCL_Char **,
			      Domain*, TclModelBuilder *);

// MHS			       
extern int 
TclModelBuilder_addZeroLengthSection(ClientData, Tcl_Interp *, int, TCL_Char **,
				     Domain*, TclModelBuilder *);

// MHS
extern int 
TclModelBuilder_addZeroLengthND(ClientData, Tcl_Interp *, int, TCL_Char **,
				Domain*, TclModelBuilder *);


// REMO
extern int 
TclModelBuilder_addNLBeamColumn(ClientData, Tcl_Interp *, int, TCL_Char **,
				Domain*, TclModelBuilder *);
			
// MHS
extern int 
TclModelBuilder_addBeamWithHinges(ClientData, Tcl_Interp *, int, TCL_Char **,
				  Domain*, TclModelBuilder *);
extern int 
TclModelBuilder_addFourNodeQuad(ClientData, Tcl_Interp *, int, TCL_Char **,
				Domain*, TclModelBuilder *);
extern int 
TclModelBuilder_addDispBeamColumn(ClientData, Tcl_Interp *, int, TCL_Char **,
				  Domain*, TclModelBuilder *);
extern int 
TclModelBuilder_addForceBeamColumn(ClientData, Tcl_Interp *, int, TCL_Char **,
				   Domain*, TclModelBuilder *);
		   
// NM
extern int
TclModelBuilder_addBeamColumnJoint(ClientData, Tcl_Interp *, int, TCL_Char **,
				   Domain*, TclModelBuilder *, int);

//Boris Jeremic & Zhaohui
extern int TclModelBuilder_addEightNodeBrick(ClientData, 
                                             Tcl_Interp *,  
					     int, 
					     TCL_Char **,
					     Domain*, 
					     TclModelBuilder *, 
					     int);
//Boris Jeremic & Zhaohui
extern int TclModelBuilder_addTwentyNodeBrick(ClientData, 
                                              Tcl_Interp *,  
					      int, 
					      TCL_Char **,
					      Domain*, 
					      TclModelBuilder *, 
					      int);

//Boris Jeremic & Xiaoyan 01/07/2002
extern int TclModelBuilder_addEightNodeBrick_u_p_U(ClientData, 
                                                   Tcl_Interp *,  
						   int, 
						   TCL_Char **,
						   Domain*, 
						   TclModelBuilder *, 
						   int);
//Boris Jeremic & Xiaoyan 01/07/2002
extern int TclModelBuilder_addTwentyNodeBrick_u_p_U(ClientData, 
                                                    Tcl_Interp *,  
						    int, 
						    TCL_Char **,
						    Domain*, 
						    TclModelBuilder *, 
						    int);

//Boris Jeremic & Guanzhou Jie 10/30/2003
extern int TclModelBuilder_addTwentySevenNodeBrick(ClientData, 
                                                   Tcl_Interp *, 
                                                   int, 
                                                   TCL_Char **, 
                                                   Domain*, 
                                                   TclModelBuilder *, 
                                                   int);



//Rohit Kraul
extern int
TclModelBuilder_addElastic2dGNL(ClientData, Tcl_Interp *, int, TCL_Char **,
				Domain *,TclModelBuilder *);
extern int
TclModelBuilder_addElement2dYS(ClientData, Tcl_Interp *, int, TCL_Char **,
			       Domain *,TclModelBuilder *);

// Zhaohui Yang
extern int
TclModelBuilder_addFourNodeQuadUP(ClientData, Tcl_Interp *, int, TCL_Char **,
				Domain*, TclModelBuilder *);

// Zhaohui Yang
extern int
TclModelBuilder_addBrickUP(ClientData, Tcl_Interp *, int, TCL_Char **,
				Domain*, TclModelBuilder *);

// Zhaohui Yang
extern int
TclModelBuilder_addNineFourNodeQuadUP(ClientData, Tcl_Interp *, int, TCL_Char **,
				Domain*, TclModelBuilder *);
// Jinchi Lu
extern int
TclModelBuilder_addTwentyEightNodeBrickUP(ClientData, Tcl_Interp *, int, TCL_Char **,
				Domain*, TclModelBuilder *);
// Jinchi Lu
extern int
TclModelBuilder_addTwentyNodeBrick(ClientData, Tcl_Interp *, int, TCL_Char **,
				Domain*, TclModelBuilder *);

// Boris Jeremic and Zhao Cheng
extern int
TclModelBuilder_addTLFD20nBrick(ClientData, Tcl_Interp *, int, TCL_Char **,
				Domain*, TclModelBuilder *, int);

int
TclModelBuilderElementCommand(ClientData clientData, Tcl_Interp *interp,
			      int argc, TCL_Char **argv, 
			      Domain *theTclDomain, TclModelBuilder *theTclBuilder)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed\n";    
    return TCL_ERROR;
  }

  // check at least two arguments so don't segemnt fault on strcmp  
  if (argc < 2) {
    opserr << "WARNING need to specify an element type\n";
    opserr << "Want: element eleType <specific element args> .. see manual for valid eleTypes & arguments\n";
    return TCL_ERROR;
  }

  if (strcmp(argv[1],"fTruss") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addFeapTruss(clientData, interp, argc, argv,
					      theTclDomain, theTclBuilder, eleArgStart);
    return result;
  } else if (strcmp(argv[1],"truss") == 0 || strcmp(argv[1],"corotTruss") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addTruss(clientData, interp, argc, argv,
					  theTclDomain, theTclBuilder, eleArgStart);
    return result;
  } else if (strcmp(argv[1],"elasticBeamColumn") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addElasticBeam(clientData, interp, argc, argv,
						theTclDomain, theTclBuilder, eleArgStart);    
    return result;
  } else if (strcmp(argv[1],"nonlinearBeamColumn") == 0) {
    int result = TclModelBuilder_addNLBeamColumn(clientData, interp, argc, argv,
						 theTclDomain, theTclBuilder);
    return result;
  } else if (strcmp(argv[1],"dispBeamColumn") == 0) {
    int result = TclModelBuilder_addDispBeamColumn(clientData, interp, argc, argv,
						   theTclDomain, theTclBuilder);
    return result;
  } else if (strcmp(argv[1],"forceBeamColumn") == 0) {
    int result = TclModelBuilder_addForceBeamColumn(clientData, interp, argc, argv,
						    theTclDomain, theTclBuilder);
    return result;
  } else if (strcmp(argv[1],"beamWithHinges") == 0 ||
	     strcmp(argv[1],"beamWithHinges1") == 0 ||
	     strcmp(argv[1],"beamWithHinges2") == 0 ||
	     strcmp(argv[1],"beamWithHinges3") == 0) {
	  int result = TclModelBuilder_addBeamWithHinges(clientData, interp, argc, argv,
						 theTclDomain, theTclBuilder);
    return result;
  } else if ((strcmp(argv[1],"quad") == 0) || (strcmp(argv[1],"stdQuad") == 0)) {
    int result = TclModelBuilder_addFourNodeQuad(clientData, interp, argc, argv,
						 theTclDomain, theTclBuilder);
	  return result;
  } else if (strcmp(argv[1],"enhancedQuad") == 0) {
    int result = TclModelBuilder_addEnhancedQuad(clientData, interp, argc, argv,
						 theTclDomain, theTclBuilder);
    return result;
  } else if ((strcmp(argv[1],"bbarQuad") == 0) || (strcmp(argv[1],"mixedQuad") == 0)) {
    int result = TclModelBuilder_addConstantPressureVolumeQuad(clientData, interp, 
							       argc, argv,
							       theTclDomain, 
							       theTclBuilder);
    return result;
  } else if ((strcmp(argv[1],"nineNodeMixedQuad") == 0) 
	     || (strcmp(argv[1],"nineNodeQuad") == 0)) {
    int result = TclModelBuilder_addNineNodeMixedQuad(clientData, interp, 
						      argc, argv,
						      theTclDomain, 
						      theTclBuilder);
    return result;
  } else if (strcmp(argv[1],"quadUP") == 0) {
          int result = TclModelBuilder_addFourNodeQuadUP(clientData, interp, argc, argv,
						       theTclDomain, theTclBuilder);
	  return result;
  } else if (strcmp(argv[1],"brickUP") == 0) {
          int result = TclModelBuilder_addBrickUP(clientData, interp, argc, argv,
						       theTclDomain, theTclBuilder);
	  return result;
  } else if (strcmp(argv[1],"9_4_QuadUP") == 0) {
          int result = TclModelBuilder_addNineFourNodeQuadUP(clientData, interp, argc, argv,
						       theTclDomain, theTclBuilder);
	  return result;
  } else if (strcmp(argv[1],"20_8_BrickUP") == 0) {
          int result = TclModelBuilder_addTwentyEightNodeBrickUP(clientData, interp, argc, argv,
						       theTclDomain, theTclBuilder);
	  return result;
  } else if (strcmp(argv[1],"20NodeBrick") == 0) {
          int result = TclModelBuilder_addTwentyNodeBrick(clientData, interp, argc, argv,
						       theTclDomain, theTclBuilder);
	  return result;
  } else if ((strcmp(argv[1],"shell") == 0) || (strcmp(argv[1],"shellMITC4") == 0) ||
	     (strcmp(argv[1],"Shell") == 0) || (strcmp(argv[1],"ShellMITC4") == 0)) {


    int eleArgStart = 1;
    int result = TclModelBuilder_addShellMITC4(clientData, interp, 
					       argc, argv,
					       theTclDomain, 
					       theTclBuilder, 
					       eleArgStart);
    return result;
  }

  //Boris Jeremic & Zhaohui
  else if (strcmp(argv[1],"Brick8N") == 0) {

    int eleArgStart = 1;
    int result = TclModelBuilder_addEightNodeBrick(clientData,
						   interp,
						   argc,
						   argv,
						   theTclDomain,
						   theTclBuilder,
						   eleArgStart);
    return result;
  }

  //Boris Jeremic & Zhaohui
  else if (strcmp(argv[1],"Brick20N") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addTwentyNodeBrick(clientData, 
						    interp, 
						    argc, 
						    argv,
						    theTclDomain, 
						    theTclBuilder, 
						    eleArgStart);
    return result;
  } 

  //Boris Jeremic & Guanzhou Jie
  else if (strcmp(argv[1],"Brick27N") == 0) 
    {
      int eleArgStart = 1;
      int result = TclModelBuilder_addTwentySevenNodeBrick(clientData, 
                                                           interp, 
                                                           argc, 
                                                           argv,
                                                           theTclDomain, 
                                                           theTclBuilder, 
                                                           eleArgStart);
                                                           return result;
    }


  // Boris jeremic & Zhao Cheng
  else if (strcmp(argv[1],"TLFD20nBrick") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addTLFD20nBrick(clientData,
                                                   interp,
                                                   argc,
                                                   argv,
                                                   theTclDomain,
                                                   theTclBuilder,
                                                   eleArgStart);
    return result;
  }


  //Boris Jeremic & Zhaohui  
  else if (strcmp(argv[1],"Brick8N_u_p_U") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addEightNodeBrick_u_p_U(clientData, 
							 interp, 
							 argc, 
							 argv,
							 theTclDomain, 
							 theTclBuilder, 
							 eleArgStart);
    return result;
  } 
  //Boris Jeremic & Zhaohui  
  else if (strcmp(argv[1],"Brick20N_u_p_U") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addTwentyNodeBrick_u_p_U(clientData, 
							  interp, 
							  argc, 
							  argv,
							  theTclDomain, 
							  theTclBuilder, 
							  eleArgStart);
    return result;
  } 
  else if (strcmp(argv[1],"stdBrick") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addBrick(clientData, interp, argc, argv,
					  theTclDomain, theTclBuilder, eleArgStart);
    return result;
  } else if (strcmp(argv[1],"bbarBrick") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addBrick(clientData, interp, argc, argv,
					  theTclDomain, theTclBuilder, eleArgStart);
    return result;
  } else if (strcmp(argv[1],"flBrick") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addBrick(clientData, interp, argc, argv,
					  theTclDomain, theTclBuilder, eleArgStart);
    return result;
  } else if (strcmp(argv[1],"zeroLength") == 0) {
    int result = TclModelBuilder_addZeroLength(clientData, interp, argc, argv,
					       theTclDomain, theTclBuilder);
    return result;
  } else if (strcmp(argv[1],"zeroLengthSection") == 0) {
    int result = TclModelBuilder_addZeroLengthSection(clientData, interp, argc, argv,
					       theTclDomain, theTclBuilder);
    return result;
  } else if (strcmp(argv[1],"zeroLengthND") == 0) {
    opserr << "element zeroLengthND is no longer available, please use "
	 << "the zeroLengthSection element instead" << endln;
    return TCL_ERROR;
  } else if ((strcmp(argv[1],"Joint2D") == 0) ||
	     (strcmp(argv[1],"Joint2d") == 0)) {
    int result = TclModelBuilder_addJoint2D(clientData, interp, argc, argv,
					    theTclDomain, theTclBuilder);
    return result;
  } else if ((strcmp(argv[1],"Joint3D") == 0) ||
	     (strcmp(argv[1],"Joint3d") == 0)) {
    int result = TclModelBuilder_addJoint3D(clientData, interp, argc, argv,
					    theTclDomain, theTclBuilder);
    return result;
  } else if ((strcmp(argv[1], "inelastic2dYS01")== 0) ||
	     (strcmp(argv[1], "inelastic2dYS02")== 0) ||
	     (strcmp(argv[1], "inelastic2dYS03")== 0) ||
	     (strcmp(argv[1], "inelastic2dYS04")== 0) ||
	     (strcmp(argv[1], "inelastic2dYS05")== 0)) {
    int result = TclModelBuilder_addElement2dYS (clientData, interp,
						 argc, argv,
						 theTclDomain, theTclBuilder);
    return result;	
  } else if ((strcmp(argv[1],"element2dGNL") == 0) ||
	     (strcmp(argv[1],"elastic2dGNL") == 0)) {
    int result = TclModelBuilder_addElastic2dGNL(clientData, interp, argc, argv,
						 theTclDomain, theTclBuilder);
    return result;
  }

  else if (strcmp(argv[1],"beamColumnJoint") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addBeamColumnJoint(clientData, interp,
						    argc, argv, theTclDomain,
						    theTclBuilder, eleArgStart);
    return result;
  }

  else {

    //
    // maybe element in a package
    //

    // try existing loaded packages

    ElementPackageCommand *eleCommands = theElementPackageCommands;
    bool found = false;
    while (eleCommands != NULL && found == false) {
      if (strcmp(argv[1], eleCommands->funcName) == 0) {
	int result = (*(eleCommands->funcPtr))(clientData, interp, argc, argv, theTclDomain, theTclBuilder);
	return result;
      } else
	eleCommands = eleCommands->next;
    }

    // load new package
    void *libHandle;
    int (*funcPtr)(ClientData clientData, Tcl_Interp *interp,  int argc, 
		   TCL_Char **argv, Domain*, TclModelBuilder *);       
    int eleNameLength = strlen(argv[1]);
    char *tclFuncName = new char[eleNameLength+12];
    strcpy(tclFuncName, "TclCommand_");
    strcpy(&tclFuncName[11], argv[1]);    

    int res = getLibraryFunction(argv[1], tclFuncName, &libHandle, (void **)&funcPtr);

    delete [] tclFuncName;

    if (res == 0) {

      char *eleName = new char[eleNameLength+1];
      strcpy(eleName, argv[1]);
      ElementPackageCommand *theEleCommand = new ElementPackageCommand;
      theEleCommand->funcPtr = funcPtr;
      theEleCommand->funcName = eleName;	
      theEleCommand->next = theElementPackageCommands;
      theElementPackageCommands = theEleCommand;
      
      int result = (*funcPtr)(clientData, interp,
			      argc, argv,
			      theTclDomain, theTclBuilder);	
      return result;
    }

    // element type not recognized
    opserr << "WARNING unknown element type: " <<  argv[1] << " :check the manual\n";
    return TCL_ERROR;
  }    
}
