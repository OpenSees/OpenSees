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

// $Revision: 1.69 $
// $Date: 2010-07-28 22:42:43 $
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
#include <CrdTransf.h>

#include <TclModelBuilder.h>
#include <packages.h>
#include <elementAPI.h>

extern 
#ifdef _WIN32
int __cdecl
#else
int
#endif
httpGET_File(char const *URL, char const *page, unsigned int port, const char *filename);

//
// SOME STATIC POINTERS USED IN THE FUNCTIONS INVOKED BY THE INTERPRETER
//

extern int OPS_ResetInput(ClientData clientData,
			  Tcl_Interp *interp,
			  int cArg,
			  int mArg,
			  TCL_Char **argv,
			  Domain *domain,
			  TclModelBuilder *builder);

typedef struct elementPackageCommand {
  char *funcName;
  void *(*funcPtr)();
  struct elementPackageCommand *next;
} ElementPackageCommand;

static ElementPackageCommand *theElementPackageCommands = NULL;

extern void printCommand(int argc, TCL_Char **argv);

//
// THE PROTOTYPES OF THE FUNCTIONS INVOKED BY THE INTERPRETER
//

extern  void *OPS_NewTrussElement(void);
extern  void *OPS_NewTrussSectionElement(void);
extern  void *OPS_NewCorotTrussElement(void);
extern  void *OPS_NewCorotTrussSectionElement(void);
extern  void *OPS_ElasticTubularJoint(void);
extern Element *OPS_NewZeroLengthContactNTS2D(void);
extern Element *OPS_NewZeroLengthInterface2D(void);
extern "C" void *OPS_PY_Macro2D(void);
extern void *OPS_SimpleContact2D(void);
extern void *OPS_SimpleContact3D(void);
extern void *OPS_BeamContact2D(void);
extern void *OPS_BeamContact2Dp(void);
extern void *OPS_BeamContact3D(void);
extern void *OPS_BeamContact3Dp(void);
extern void *OPS_SurfaceLoad(void);
extern void *OPS_ModElasticBeam2d(void);
extern void *OPS_TPB1D(void);
extern void *OPS_BeamEndContact3D(void);
extern void *OPS_BeamEndContact3Dp(void);
extern void *OPS_TFP_Bearing(void);
extern void *OPS_MultiFP2d(void);
extern void *OPS_CoupledZeroLength(void);
extern void *OPS_FourNodeQuad3d(void);
extern void *OPS_Tri31(void);
extern void *OPS_SSPquad(void);
extern void *OPS_SSPquadUP(void);
extern void *OPS_SSPbrick(void);
extern void *OPS_SSPbrickUP(void);
extern void *OPS_NewShellMITC4(void);
extern void *OPS_NewShellNL(void);
extern void *OPS_Quad4FiberOverlay(void);
extern void *OPS_Brick8FiberOverlay(void);
extern void *OPS_TripleFrictionPendulum(void);


extern int TclModelBuilder_addFeapTruss(ClientData clientData, Tcl_Interp *interp,  int argc,
					TCL_Char **argv, Domain*, TclModelBuilder *, int argStart);

extern int
Tcl_addWrapperElement(eleObj *, ClientData clientData, Tcl_Interp *interp,  int argc,
		      TCL_Char **argv, Domain*, TclModelBuilder *);

extern int
TclModelBuilder_addElasticBeam(ClientData clientData, Tcl_Interp *interp,  int argc,
			       TCL_Char **argv, Domain*, TclModelBuilder *, int argStart);

extern int
TclModelBuilder_addBrick(ClientData clientData, Tcl_Interp *interp,
			 int argc, TCL_Char **argv, Domain*,
			 TclModelBuilder *, int argStart);

extern int
TclModelBuilder_addConstantPressureVolumeQuad(ClientData, Tcl_Interp *, int, TCL_Char **,
					      Domain*, TclModelBuilder *);

extern int
TclModelBuilder_addJoint2D(ClientData, Tcl_Interp *, int, TCL_Char **, Domain*);
			   

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


// add by Gang Wang for Contact Element
extern int
TclModelBuilder_addZeroLengthContact2D(ClientData, Tcl_Interp *, int, TCL_Char **,
				       Domain*, TclModelBuilder *);

// add by Gang Wang for Contact Element
extern int
TclModelBuilder_addZeroLengthContact3D(ClientData, Tcl_Interp *, int, TCL_Char **,
				       Domain*, TclModelBuilder *);

// KRM added for rocking element
extern int
TclModelBuilder_addZeroLengthRocking(ClientData, Tcl_Interp *, int, TCL_Char **,
                                       Domain*, TclModelBuilder *);

// MHS
extern int
TclModelBuilder_addZeroLengthSection(ClientData, Tcl_Interp *, int, TCL_Char **,
				     Domain*, TclModelBuilder *);

// MHS
extern int
TclModelBuilder_addZeroLengthND(ClientData, Tcl_Interp *, int, TCL_Char **,
				Domain*, TclModelBuilder *);


// MHS
extern int
TclModelBuilder_addBeamWithHinges(ClientData, Tcl_Interp *, int, TCL_Char **,
				  Domain*, TclModelBuilder *);
// Minjie Zhu
extern int
TclModelBuilder_addPFEMElement2D(ClientData clientData, Tcl_Interp *interp,  int argc,
                                 TCL_Char **argv, Domain*, TclModelBuilder *);

// Quan
extern int
TclModelBuilder_addFourNodeQuadWithSensitivity(ClientData, Tcl_Interp *, int, TCL_Char **,
					       Domain*, TclModelBuilder *);

extern int
TclModelBuilder_addFourNodeQuad(ClientData, Tcl_Interp *, int, TCL_Char **,
				Domain*, TclModelBuilder *);

extern int
TclModelBuilder_addDispBeamColumnInt(ClientData, Tcl_Interp *, int, TCL_Char **,
				     Domain*, TclModelBuilder *);

extern int
TclModelBuilder_addForceBeamColumn(ClientData, Tcl_Interp *, int, TCL_Char **,
				   Domain*, TclModelBuilder *);

// NM
extern int
TclModelBuilder_addBeamColumnJoint(ClientData, Tcl_Interp *, int, TCL_Char **, Domain*, int);


// Andreas Schellenberg
extern int
TclModelBuilder_addGenericClient(ClientData , Tcl_Interp *,  int, TCL_Char **,
				 Domain*, TclModelBuilder *, int argStart);

extern int
TclModelBuilder_addGenericCopy(ClientData , Tcl_Interp *,  int, TCL_Char **,
				 Domain*, TclModelBuilder *, int argStart);


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

// Zhaohui Yang
extern int
TclModelBuilder_addBBarFourNodeQuadUP(ClientData, Tcl_Interp *, int, TCL_Char **,
				      Domain*, TclModelBuilder *);

// Zhaohui Yang
extern int
TclModelBuilder_addBBarBrickUP(ClientData, Tcl_Interp *, int, TCL_Char **,
			       Domain*, TclModelBuilder *);

// Jinchi Lu
extern int
TclModelBuilder_addTwentyEightNodeBrickUP(ClientData, Tcl_Interp *, int, TCL_Char **,
					  Domain*, TclModelBuilder *);
// Jinchi Lu
extern int
TclModelBuilder_addTwentyNodeBrick(ClientData, Tcl_Interp *, int, TCL_Char **,
				   Domain*, TclModelBuilder *);

// Andreas Schellenberg
extern int
TclModelBuilder_addActuator(ClientData clientData, Tcl_Interp *interp,  int argc,
			    TCL_Char **argv, Domain*, TclModelBuilder *, int argStart);

extern int
TclModelBuilder_addActuatorCorot(ClientData clientData, Tcl_Interp *interp,  int argc,
				 TCL_Char **argv, Domain*, TclModelBuilder *, int argStart);

extern int
TclModelBuilder_addAdapter(ClientData clientData, Tcl_Interp *interp,  int argc,
			   TCL_Char **argv, Domain*, TclModelBuilder *, int argStart);

extern int
TclModelBuilder_addFlatSliderBearing(ClientData clientData, Tcl_Interp *interp,  int argc,
				     TCL_Char **argv, Domain*, TclModelBuilder *, int argStart);

extern int
  TclModelBuilder_addSingleFPBearing(ClientData clientData, Tcl_Interp *interp,  int argc,
  TCL_Char **argv, Domain*, TclModelBuilder *, int argStart);
/*  
  extern int
  TclModelBuilder_addDoubleFPBearing(ClientData clientData, Tcl_Interp *interp,  int argc,
  TCL_Char **argv, Domain*, TclModelBuilder *, int argStart);
*/
extern int
TclModelBuilder_addElastomericBearing(ClientData clientData, Tcl_Interp *interp,  int argc,
				      TCL_Char **argv, Domain*, TclModelBuilder *, int argStart);

extern int
TclModelBuilder_addTwoNodeLink(ClientData clientData, Tcl_Interp *interp,  int argc,
			       TCL_Char **argv, Domain*, TclModelBuilder *, int argStart);


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


  OPS_ResetInput(clientData, interp, 2, argc, argv, theTclDomain, theTclBuilder);

  // check at least two arguments so don't segemnt fault on strcmp
  if (argc < 2) {
    opserr << "WARNING need to specify an element type\n";
    opserr << "Want: element eleType <specific element args> .. see manual for valid eleTypes & arguments\n";
    return TCL_ERROR;
  }

  Element *theElement = 0;

  if ((strcmp(argv[1],"truss") == 0) || (strcmp(argv[1],"Truss") == 0)) {
    
    void *theEle = OPS_NewTrussElement();
    // for backward compatability
	if (theEle == 0) {
      theEle = OPS_NewTrussSectionElement(); 
	}

    if (theEle != 0) 
      theElement = (Element *)theEle;

    else  {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;      
    }

  } else if ((strcmp(argv[1],"trussSection") == 0) || (strcmp(argv[1],"TrussSection") == 0)) {

    void *theEle = OPS_NewTrussSectionElement(); 
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else  {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;      
    }
  }

  else if ((strcmp(argv[1],"corotTruss") == 0) || (strcmp(argv[1],"CorotTruss") == 0)) {
    
    void *theEle = OPS_NewCorotTrussElement();
    
    // for backward compatability
    if (theEle == 0)
      theEle = OPS_NewCorotTrussSectionElement(); 
    
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;      
    }

  } else if ((strcmp(argv[1],"corotTrussSection") == 0) || (strcmp(argv[1],"CorotTrussSection") == 0)) {

    void *theEle = OPS_NewCorotTrussSectionElement(); 
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else  {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;      
    }

  } else if (strcmp(argv[1],"zeroLengthContactNTS2D") == 0) {
    Element *theEle = OPS_NewZeroLengthContactNTS2D();
    if (theEle != 0) 
      theElement = theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if (strcmp(argv[1],"zeroLengthInterface2D") == 0) {
    Element *theEle = OPS_NewZeroLengthInterface2D();
    if (theEle != 0) 
      theElement = theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if ((strcmp(argv[1],"ModElasticBeam2d") == 0) || (strcmp(argv[1],"modElasticBeam2d")) == 0) {
    Element *theEle = (Element *)OPS_ModElasticBeam2d();
    if (theEle != 0) 
      theElement = theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if ((strcmp(argv[1],"pyMacro2D") == 0) || (strcmp(argv[1],"PY_Macro2D") == 0)) {
    
    void *theEle = OPS_PY_Macro2D();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if ((strcmp(argv[1],"SimpleContact2d") == 0) || (strcmp(argv[1],"SimpleContact2D") == 0)) {
    
    void *theEle = OPS_SimpleContact2D();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if ((strcmp(argv[1],"SimpleContact3d") == 0) || (strcmp(argv[1],"SimpleContact3D") == 0)) {
    
    void *theEle = OPS_SimpleContact3D();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if ((strcmp(argv[1],"BeamContact3d") == 0) || (strcmp(argv[1],"BeamContact3D") == 0)) {
    
    void *theEle = OPS_BeamContact3D();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if ((strcmp(argv[1],"BeamContact3dp") == 0) || (strcmp(argv[1],"BeamContact3Dp") == 0)) {
    
    void *theEle = OPS_BeamContact3Dp();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if ((strcmp(argv[1],"TFPbearing") == 0) || (strcmp(argv[1],"TFP") == 0)) {
    
    void *theEle = OPS_TFP_Bearing();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }


  } else if (strcmp(argv[1],"TripleFrictionPendulum") == 0) {
    
    void *theEle = OPS_TripleFrictionPendulum();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }


  } else if ((strcmp(argv[1],"MultiFP2d") == 0) || (strcmp(argv[1],"MultiFPB2d") == 0)){
    
    void *theEle = OPS_MultiFP2d();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if ((strcmp(argv[1],"shell") == 0) || (strcmp(argv[1],"shellMITC4") == 0) ||
	     (strcmp(argv[1],"Shell") == 0) || (strcmp(argv[1],"ShellMITC4") == 0)) {
    
    void *theEle = OPS_NewShellMITC4();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }    

  } else if ((strcmp(argv[1],"shellNL") == 0) || (strcmp(argv[1],"ShellNL") == 0)) {
    
    void *theEle = OPS_NewShellNL();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }    


  } else if ((strcmp(argv[1],"CoupledZeroLength") == 0) || (strcmp(argv[1],"ZeroLengthCoupled") == 0)) {
    
    void *theEle = OPS_CoupledZeroLength();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if ((strcmp(argv[1],"BeamContact2d") == 0) || (strcmp(argv[1],"BeamContact2D") == 0)) {
    
    void *theEle = OPS_BeamContact2D();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if ((strcmp(argv[1],"BeamContact2dp") == 0) || (strcmp(argv[1],"BeamContact2Dp") == 0)) {
    
    void *theEle = OPS_BeamContact2Dp();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if ((strcmp(argv[1],"BeamEndContact3d") == 0) || (strcmp(argv[1],"BeamEndContact3D") == 0)) {
    
    void *theEle = OPS_BeamEndContact3D();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if ((strcmp(argv[1],"BeamEndContact3dp") == 0) || (strcmp(argv[1],"BeamEndContact3Dp") == 0)) {
    
    void *theEle = OPS_BeamEndContact3Dp();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if ((strcmp(argv[1],"Tri31") == 0) || (strcmp(argv[1],"tri31") == 0)) {
    
    void *theEle = OPS_Tri31();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if ((strcmp(argv[1],"SSPquad") == 0) || (strcmp(argv[1],"SSPQuad") == 0)) {
    
    void *theEle = OPS_SSPquad();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if ((strcmp(argv[1],"SSPquadUP") == 0) || (strcmp(argv[1],"SSPQuadUP") == 0)) {
    
    void *theEle = OPS_SSPquadUP();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if ((strcmp(argv[1],"SSPbrick") == 0) || (strcmp(argv[1],"SSPBrick") == 0)) {
    
    void *theEle = OPS_SSPbrick();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if ((strcmp(argv[1],"SSPbrickUP") == 0) || (strcmp(argv[1],"SSPBrickUP") == 0)) {
    
    void *theEle = OPS_SSPbrickUP();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if ((strcmp(argv[1],"SurfaceLoad") == 0)) {
    
    void *theEle = OPS_SurfaceLoad();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if ((strcmp(argv[1],"TPB1D") == 0)) {
    
    void *theEle = OPS_TPB1D();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }
  
  } else if ((strcmp(argv[1],"elasticTubularJoint") == 0) || (strcmp(argv[1],"ElasticTubularJoint") == 0)) {
    
    void *theEle = OPS_ElasticTubularJoint();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if ((strcmp(argv[1],"quad3d") == 0) || (strcmp(argv[1],"Quad3d") == 0)) {

    void *theEle = OPS_FourNodeQuad3d();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }	

  } else if ((strcmp(argv[1],"Quad4FiberOverlay") == 0)) { 	//////////////////////// mmc
    
    void *theEle = OPS_Quad4FiberOverlay();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "tclelementcommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }																								
  } else if ((strcmp(argv[1],"Brick8FiberOverlay") == 0)) { //////////////////////// mmc
    
    void *theEle = OPS_Brick8FiberOverlay();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "tclelementcommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }																									
  }

  // if one of the above worked
  if (theElement != 0) {
    if (theTclDomain->addElement(theElement) == false) {
      opserr << "WARNING could not add element of with tag: " << theElement->getTag() << " and of type: " 
	     << theElement->getClassType() << " to the Domain\n";
      delete theElement;
      return TCL_ERROR;
    } else
      return TCL_OK;
  }


  if (strcmp(argv[1],"fTruss") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addFeapTruss(clientData, interp, argc, argv,
					      theTclDomain, theTclBuilder, eleArgStart);
    return result;

  } else if (strcmp(argv[1],"elasticBeamColumn") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addElasticBeam(clientData, interp, argc, argv,
						theTclDomain, theTclBuilder, eleArgStart);
    return result;
  }
  /*
else if (strcmp(argv[1],"nonlinearBeamColumn") == 0) {
    int result = TclModelBuilder_addNLBeamColumn(clientData, interp, argc, argv,
						 theTclDomain, theTclBuilder);
    return result;
  } else if (strcmp(argv[1],"dispBeamColumn") == 0) {
    int result = TclModelBuilder_addDispBeamColumn(clientData, interp, argc, argv,
						   theTclDomain, theTclBuilder);
    return result;

  }
  */
  else if (strcmp(argv[1],"dispBeamColumnInt") == 0) {
    int result = TclModelBuilder_addDispBeamColumnInt(clientData, interp, argc, argv,
						   theTclDomain, theTclBuilder);
    return result;

  } else if (strcmp(argv[1],"forceBeamColumn") == 0 || 
	     strcmp(argv[1],"dispBeamColumn") == 0  || 
	     strcmp(argv[1],"dispBeamColumnThermal") == 0  || 
	     strcmp(argv[1],"elasticForceBeamColumn") == 0 || 
	     strcmp(argv[1],"nonlinearBeamColumn") == 0 || 
	     strcmp(argv[1],"dispBeamColumnWithSensitivity") == 0) {

    int result = TclModelBuilder_addForceBeamColumn(clientData, interp, argc, argv,
						    theTclDomain, theTclBuilder);
    return result;
  } else if (strstr(argv[1],"beamWithHinges") != 0) {
    int result = TclModelBuilder_addBeamWithHinges(clientData, interp, argc, argv,
						   theTclDomain, theTclBuilder);
    return result;
  } else if ((strcmp(argv[1],"quad") == 0) || (strcmp(argv[1],"stdQuad") == 0)) {
    int result = TclModelBuilder_addFourNodeQuad(clientData, interp, argc, argv,
						 theTclDomain, theTclBuilder);
    return result;

  } else if (strcmp(argv[1],"quadWithSensitivity") == 0) {
    int result = TclModelBuilder_addFourNodeQuadWithSensitivity(clientData, interp, argc, argv,
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
  } else if (strcmp(argv[1],"bbarQuadUP") == 0) {
    int result = TclModelBuilder_addBBarFourNodeQuadUP(clientData, interp, argc, argv,
						       theTclDomain, theTclBuilder);
    return result;
  } else if (strcmp(argv[1],"bbarBrickUP") == 0) {
    int result = TclModelBuilder_addBBarBrickUP(clientData, interp, argc, argv,
						theTclDomain, theTclBuilder);
    return result;
  }

  // Andreas Schellenberg
  else if (strcmp(argv[1],"genericClient") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addGenericClient(clientData, interp, argc, argv,
						theTclDomain, theTclBuilder, eleArgStart);
    return result;
  }

  else if (strcmp(argv[1],"genericCopy") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addGenericCopy(clientData, interp, argc, argv,
						theTclDomain, theTclBuilder, eleArgStart);
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
  } else if (strcmp(argv[1],"bbarBrickWithSensitivity") == 0) {
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
  } else if (strcmp(argv[1],"zeroLengthRocking") == 0) {
      int result = TclModelBuilder_addZeroLengthRocking(clientData, interp, argc, argv,
                            theTclDomain, theTclBuilder);
    return result;
  } else if (strcmp(argv[1],"zeroLengthContact2D") == 0) {
    int result = TclModelBuilder_addZeroLengthContact2D(clientData, interp, argc, argv,
							theTclDomain, theTclBuilder);
    return result;
  } else if (strcmp(argv[1],"zeroLengthContact3D") == 0) {
    int result = TclModelBuilder_addZeroLengthContact3D(clientData, interp, argc, argv,
							theTclDomain, theTclBuilder);
    return result;
  } else if (strcmp(argv[1],"zeroLengthND") == 0) {
    int result = TclModelBuilder_addZeroLengthND(clientData, interp, argc, argv,
						 theTclDomain, theTclBuilder);
    return result;
  } else if ((strcmp(argv[1],"Joint2D") == 0) ||
	     (strcmp(argv[1],"Joint2d") == 0)) {
    int result = TclModelBuilder_addJoint2D(clientData, interp, argc, argv,theTclDomain);
					    
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
						    argc, argv, theTclDomain,eleArgStart);

    return result;
  }
  
  // Andreas Schellenberg
  else if (strcmp(argv[1],"actuator") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addActuator(clientData, interp, argc, argv,
					     theTclDomain, theTclBuilder, eleArgStart);
    return result;
  }
  
  else if (strcmp(argv[1],"corotActuator") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addActuatorCorot(clientData, interp, argc, argv,
						  theTclDomain, theTclBuilder, eleArgStart);
    return result;
  }
  
  else if (strcmp(argv[1],"adapter") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addAdapter(clientData, interp, argc, argv,
					    theTclDomain, theTclBuilder, eleArgStart);
    return result;
  }

  else if (strcmp(argv[1],"flatSliderBearing") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addFlatSliderBearing(clientData, interp, argc, argv,
						      theTclDomain, theTclBuilder, eleArgStart);
    return result;
  }

  else if (strcmp(argv[1],"singleFPBearing") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addSingleFPBearing(clientData, interp, argc, argv,
						theTclDomain, theTclBuilder, eleArgStart);
    return result;
  }

  /*
  else if (strcmp(argv[1],"doubleFPBearing") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addDoubleFPBearing(clientData, interp, argc, argv,
						theTclDomain, theTclBuilder, eleArgStart);
    return result;
  }*/

  else if (strcmp(argv[1],"elastomericBearing") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addElastomericBearing(clientData, interp, argc, argv,
						       theTclDomain, theTclBuilder, eleArgStart);
    return result;
  }

  else if (strcmp(argv[1],"twoNodeLink") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addTwoNodeLink(clientData, interp, argc, argv,
						theTclDomain, theTclBuilder, eleArgStart);
    return result;
  }

  else if (strcmp(argv[1],"PFEMElement2D") == 0) {
      int result = TclModelBuilder_addPFEMElement2D(clientData, interp, argc, argv,
                                                    theTclDomain, theTclBuilder);
    return result;
  }

  else {

    //
    // maybe element already loaded as c++ class from a package
    //

    // try existing loaded packages

    ElementPackageCommand *eleCommands = theElementPackageCommands;
    bool found = false;
    int result = TCL_ERROR;
    while (eleCommands != NULL && found == false) {
      if (strcmp(argv[1], eleCommands->funcName) == 0) {

	OPS_ResetInput(clientData, interp, 2, argc, argv, theTclDomain, theTclBuilder);
	void *theRes = (*(eleCommands->funcPtr))();
	if (theRes != 0) {
	  Element *theEle = (Element *)theRes;
	  result = theTclDomain->addElement(theEle);

	  if (result >= 0)
	    return TCL_OK;
	  else
	    return TCL_ERROR;
	}
	return TCL_ERROR;;
      } else
	eleCommands = eleCommands->next;
    }

   //
    // maybe element in a routine, check existing ones or try loading new ones
    //

    char *eleType = new char[strlen(argv[1])+1];
    strcpy(eleType, argv[1]);
    eleObj *eleObject = OPS_GetElementType(eleType, strlen(eleType));

    delete [] eleType;
    
    if (eleObject != 0) {

      int result = Tcl_addWrapperElement(eleObject, clientData, interp,
					 argc, argv,
					 theTclDomain, theTclBuilder);
      
      if (result != 0)
	delete eleObject;
      else
	return result;
    }

    //
    // try loading new dynamic library containg a c+= class
    //
    
    void *libHandle;
    void *(*funcPtr)();
    int eleNameLength = strlen(argv[1]);
    char *tclFuncName = new char[eleNameLength+5];
    strcpy(tclFuncName, "OPS_");

    strcpy(&tclFuncName[4], argv[1]);
    
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
      
      OPS_ResetInput(clientData, interp, 2, argc, argv, theTclDomain, theTclBuilder);
      
      void *theRes = (*funcPtr)();
      
      if (theRes != 0) {
	Element *theEle = (Element *)theRes;
	result = theTclDomain->addElement(theEle);
	if (result >= 0)
	  return TCL_OK;
	else
	  return TCL_ERROR;
      } else {
	return TCL_ERROR;
      }
    }

    
  }

  // If we get here, the element type is unknown
  opserr << "ERROR -- element of type " << argv[1] << " not known" << endln;
  return TCL_ERROR;
}
