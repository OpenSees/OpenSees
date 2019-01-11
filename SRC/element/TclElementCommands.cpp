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

#include <Element.h>
#include <stdlib.h>
#include <string.h>
#include <OPS_Stream.h>
#include <Domain.h>

#include <CrdTransf.h>

#include <TclModelBuilder.h>
#include <packages.h>
#include <elementAPI.h>

#include <UniaxialMaterial.h>
#include <MultipleShearSpring.h>
#include <KikuchiBearing.h>
#include <YamamotoBiaxialHDR.h>
#include<WheelRail.h>

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

extern  void *OPS_ComponentElement2d(void);
//extern  void *OPS_ComponentElementDamp2d(void);
extern  void *OPS_TrussElement(void);
extern  void *OPS_TrussSectionElement(void);
extern  void *OPS_CorotTrussElement(void);
extern  void *OPS_CorotTrussSectionElement(void);
extern  void *OPS_ElasticTubularJoint(void);
extern void *OPS_ZeroLengthContactNTS2D(void);
extern void *OPS_ZeroLengthVG_HG(void);
extern void *OPS_ZeroLengthInterface2D(void);
extern "C" void *OPS_PY_Macro2D(void);
extern void *OPS_SimpleContact2D(void);
extern void *OPS_SimpleContact3D(void);
extern void *OPS_BeamContact2D(void);
extern void *OPS_BeamContact2Dp(void);
extern void *OPS_BeamContact3D(void);
extern void *OPS_BeamContact3Dp(void);
extern void *OPS_PileToe3D(void);
extern void *OPS_SurfaceLoad(void);
extern void *OPS_TriSurfaceLoad(void);
extern void *OPS_ModElasticBeam2d(void);
extern void *OPS_ElasticBeam2d(const ID &info);
extern void *OPS_ElasticBeam3d(void);
extern void *OPS_ElasticTimoshenkoBeam2d(void);
extern void *OPS_ElasticTimoshenkoBeam3d(void);
extern void *OPS_TPB1D(void);
extern void *OPS_BeamEndContact3D(void);
extern void *OPS_BeamEndContact3Dp(void);
extern void *OPS_TFP_Bearing(void);
extern void *OPS_TFP_Bearing(void);
extern void *OPS_FPBearingPTV();
extern void *OPS_MultiFP2d(void);
extern void *OPS_CoupledZeroLength(void);
extern void *OPS_FourNodeQuad3d(void);
extern void *OPS_Tri31(const ID &info);
extern void *OPS_SSPquad(void);
extern void *OPS_SSPquadUP(void);
extern void *OPS_SSPbrick(void);
extern void *OPS_SSPbrickUP(void);
extern void *OPS_ShellMITC4(void);
extern void *OPS_ShellMITC9(void);
extern void *OPS_ShellDKGQ(void);     //Added by Lisha Wang, Xinzheng Lu, Linlin Xie, Song Cen & Quan Gu
extern void *OPS_ShellNLDKGQ(void);   //Added by Lisha Wang, Xinzheng Lu, Linlin Xie, Song Cen & Quan Gu
extern void *OPS_ShellDKGT(void);     //Added by Shuhao Zhang and  Xinzheng Lu 
extern void *OPS_ShellNLDKGT(void);   //Added by Shuhao Zhang and  Xinzheng Lu 
extern void *OPS_Quad4FiberOverlay(void);
extern void *OPS_Brick8FiberOverlay(void);
extern void *OPS_QuadBeamEmbedContact(void);
extern void *OPS_TripleFrictionPendulum(void);
extern void *OPS_Truss2(void);
extern void *OPS_CorotTruss2(void);
extern void *OPS_ZeroLengthImpact3D(void);
extern void *OPS_HDR(void);
extern void *OPS_LeadRubberX(void);
extern void *OPS_ElastomericX(void);
extern void *OPS_N4BiaxialTruss(void);
extern void *OPS_AC3D8HexWithSensitivity(void);
extern void *OPS_ASID8QuadWithSensitivity(void);
extern void *OPS_AV3D4QuadWithSensitivity(void);
extern void *OPS_VS3D4WuadWithSensitivity(void);
extern void *OPS_MVLEM(void);
extern void *OPS_SFI_MVLEM(void);
extern void *OPS_SFI_MVLEM(void);
extern void *OPS_AxEqDispBeamColumn2d(void);
extern void *OPS_ElastomericBearingBoucWenMod3d(void);
extern void *OPS_PFEMElement2DBubble(const ID &info);
extern void *OPS_PFEMElement2Dmini(const ID &info);
extern void *OPS_PFEMElement2D();

extern void *OPS_ShellMITC4Thermal(void);//Added by L.Jiang [SIF]
extern void *OPS_ShellNLDKGQThermal(void);//Added by L.Jiang [SIF]

extern  void *OPS_CatenaryCableElement(void);
extern  void *OPS_ShellANDeS(void);
extern  void *OPS_FourNodeTetrahedron(void);
extern  void *OPS_LysmerTriangle(void);

extern int TclModelBuilder_addFeapTruss(ClientData clientData, Tcl_Interp *interp,  int argc,
					TCL_Char **argv, Domain*, TclModelBuilder *, int argStart);

extern int
Tcl_addWrapperElement(eleObj *, ClientData clientData, Tcl_Interp *interp,  int argc,
		      TCL_Char **argv, Domain*, TclModelBuilder *);

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

extern int
TclModelBuilder_addRJWatsonEqsBearing(ClientData clientData, Tcl_Interp *interp,  int argc,
				     TCL_Char **argv, Domain*, TclModelBuilder *, int argStart);
/*  
  extern int
  TclModelBuilder_addDoubleFPBearing(ClientData clientData, Tcl_Interp *interp,  int argc,
  TCL_Char **argv, Domain*, TclModelBuilder *, int argStart);
*/
extern int
TclModelBuilder_addElastomericBearingPlasticity(ClientData clientData, Tcl_Interp *interp,  int argc,
				      TCL_Char **argv, Domain*, TclModelBuilder *, int argStart);

extern int
TclModelBuilder_addElastomericBearingBoucWen(ClientData clientData, Tcl_Interp *interp,  int argc,
				      TCL_Char **argv, Domain*, TclModelBuilder *, int argStart);

extern int
TclModelBuilder_addElastomericBearingUFRP(ClientData clientData, Tcl_Interp *interp,  int argc,
				      TCL_Char **argv, Domain*, TclModelBuilder *, int argStart);

extern int
TclModelBuilder_addTwoNodeLink(ClientData clientData, Tcl_Interp *interp,  int argc,
			       TCL_Char **argv, Domain*, TclModelBuilder *, int argStart);

extern int
TclModelBuilder_addMultipleShearSpring(ClientData clientData, Tcl_Interp *interp,  int argc,
			               TCL_Char **argv, Domain*, TclModelBuilder *);

extern int
TclModelBuilder_addMultipleNormalSpring(ClientData clientData, Tcl_Interp *interp,  int argc,
				        TCL_Char **argv, Domain*, TclModelBuilder *);

extern int
TclModelBuilder_addKikuchiBearing(ClientData clientData, Tcl_Interp *interp,  int argc,
				  TCL_Char **argv, Domain*, TclModelBuilder *);

extern int
TclModelBuilder_addYamamotoBiaxialHDR(ClientData clientData, Tcl_Interp *interp,  int argc,
				  TCL_Char **argv, Domain*, TclModelBuilder *);


// Added by Quan Gu and Yongdou Liu, et al. on 2018/10/31 (Xiamen University)
extern int
TclModelBuilder_addWheelRail(ClientData clientData, Tcl_Interp *interp, int argc,
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
    
    void *theEle = OPS_TrussElement();
    // for backward compatibility
	if (theEle == 0) {
      theEle = OPS_TrussSectionElement(); 
	}

    if (theEle != 0) 
      theElement = (Element *)theEle;

    else  {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;      
    }

  } else if ((strcmp(argv[1],"trussSection") == 0) || (strcmp(argv[1],"TrussSection") == 0)) {

    void *theEle = OPS_TrussSectionElement(); 
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else  {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;      
    }
  }

  else if ((strcmp(argv[1],"corotTruss") == 0) || (strcmp(argv[1],"CorotTruss") == 0)) {
    
    void *theEle = OPS_CorotTrussElement();
    
    // for backward compatibility
    if (theEle == 0)
      theEle = OPS_CorotTrussSectionElement(); 
    
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;      
    }

  } else if ((strcmp(argv[1],"corotTrussSection") == 0) || (strcmp(argv[1],"CorotTrussSection") == 0)) {

    void *theEle = OPS_CorotTrussSectionElement(); 
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else  {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;      
    }

  } else if (strcmp(argv[1],"zeroLengthContactNTS2D") == 0) {
      Element *theEle = (Element*) OPS_ZeroLengthContactNTS2D();
    if (theEle != 0) 
      theElement = theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if (strcmp(argv[1],"zeroLengthInterface2D") == 0) {
      Element *theEle = (Element*) OPS_ZeroLengthInterface2D();
    if (theEle != 0) 
      theElement = theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }


  } else if (strcmp(argv[1],"componentElement2d") == 0) {
    void *theEle = OPS_ComponentElement2d();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

    /*
  } else if (strcmp(argv[1],"componentElementDamp2d") == 0) {
    void *theEle = OPS_ComponentElementDamp2d();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }
    */
  } else if (strcmp(argv[1],"zeroLengthImpact3D") == 0) {
    void *theEle = OPS_ZeroLengthImpact3D();
    if (theEle != 0) 
      theElement = (Element *)theEle;
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

  } else if ((strcmp(argv[1],"elasticBeamColumn") == 0) || (strcmp(argv[1],"elasticBeam")) == 0) {
    Element *theEle = 0;
    ID info;
    if (OPS_GetNDM() == 2)
      theEle = (Element *)OPS_ElasticBeam2d(info);
    else
      theEle = (Element *)OPS_ElasticBeam3d();
    if (theEle != 0) 
      theElement = theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  }

  // Beginning of WheelRail element TCL command
  //Added by Quan Gu and Yongdou Liu, et al. on 2018/10/31
  else if((strcmp(argv[1], "WheelRail") == 0)) {
  // ------------------------------add------------------------------------------
  int eleArgStart = 1;
  int result = TclModelBuilder_addWheelRail(clientData, interp, argc, argv,
	  theTclDomain, theTclBuilder, eleArgStart);
  return result;

}//End of WheelRail element TCL command*/

  else if ((strcmp(argv[1],"ElasticTimoshenkoBeam") == 0) || (strcmp(argv[1],"elasticTimoshenkoBeam")) == 0) {
    Element *theEle = 0;
    if (OPS_GetNDM() == 2)
      theEle = (Element *)OPS_ElasticTimoshenkoBeam2d();
    else
      theEle = (Element *)OPS_ElasticTimoshenkoBeam3d();
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

  } else if ((strcmp(argv[1],"N4BiaxialTruss") == 0)) {
    
    void *theEle = OPS_N4BiaxialTruss();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "tclelementcommand -- unable to create element of type : " << argv[1] << endln;
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
    
  } else if ((strcmp(argv[1],"PileToe3d") == 0) || (strcmp(argv[1],"PileToe3D") == 0)) {
    
    void *theEle = OPS_PileToe3D();
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

  } else if ((strcmp(argv[1],"FPBearingPTV") == 0)) {
    
    void *theEle = OPS_FPBearingPTV();
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

  } else if (strcmp(argv[1],"HDR") == 0) {
    
    void *theEle = OPS_HDR();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if (strcmp(argv[1],"LeadRubberX") == 0) {
    
    void *theEle = OPS_LeadRubberX();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if (strcmp(argv[1],"ElastomericX") == 0) {
    
    void *theEle = OPS_ElastomericX();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if (strcmp(argv[1], "AxEqDispBeamColumn2d") == 0) {
    
    void *theEle = OPS_AxEqDispBeamColumn2d();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if (strcmp(argv[1],"MVLEM") == 0) {
    
    void *theEle = OPS_MVLEM();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }

  } else if (strcmp(argv[1],"SFI_MVLEM") == 0) {
    
    void *theEle = OPS_SFI_MVLEM();
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
    
    void *theEle = OPS_ShellMITC4();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }    

    //Added by L.Jiang [SIF]
  } else if ((strcmp(argv[1], "shellMITC4Thermal") == 0) 
	     || (strcmp(argv[1], "ShellMITC4Thermal") == 0)) {
      
      void *theEle = OPS_ShellMITC4Thermal();
      if (theEle != 0)
	theElement = (Element *)theEle;
      else {
	opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
	return TCL_ERROR;
      }      
    }

    else if ((strcmp(argv[1], "shellNLDKGQThermal") == 0) || (strcmp(argv[1], "ShellNLDKGQThermal") == 0)) {
      
      void *theEle = OPS_ShellNLDKGQThermal();
      if (theEle != 0)
	theElement = (Element *)theEle;
      else {
	opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
	return TCL_ERROR;
      }
      //end of adding thermo-mechanical shell elements by L.Jiang [SIF]  
      
  } else if ((strcmp(argv[1],"shellNL") == 0) || (strcmp(argv[1],"ShellNL") == 0) ||
	     (strcmp(argv[1],"shellMITC9") == 0) || (strcmp(argv[1],"ShellMITC9") == 0)) {
    
    void *theEle = OPS_ShellMITC9();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }    

  } else if ((strcmp(argv[1],"shellDKGQ") == 0) || (strcmp(argv[1],"ShellDKGQ") == 0)) {     //Lisha Wang & Xinzheng Lu
    
    void *theEle = OPS_ShellDKGQ();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }    

  } else if ((strcmp(argv[1],"shellNLDKGQ") == 0) || (strcmp(argv[1],"ShellNLDKGQ") == 0)) {    //Lisha Wang & Xinzheng Lu
    
    void *theEle = OPS_ShellNLDKGQ();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    } 

  } else if ((strcmp(argv[1],"shellDKGT") == 0) || (strcmp(argv[1],"ShellDKGT") == 0)) {  
    
    void *theEle = OPS_ShellDKGT();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }    

  } else if ((strcmp(argv[1],"shellNLDKGT") == 0) || (strcmp(argv[1],"ShellNLDKGT") == 0)) {    
    
    void *theEle = OPS_ShellNLDKGT();
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
    
    ID info;
    void *theEle = OPS_Tri31(info);
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
  } else if ((strcmp(argv[1],"TriSurfaceLoad") == 0)) {
    
    void *theEle = OPS_TriSurfaceLoad();
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
  } else if ((strcmp(argv[1],"QuadBeamEmbedContact") == 0)) { 	//////////////////////// mmc
    
    void *theEle = OPS_QuadBeamEmbedContact();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "tclelementcommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }																								
  } else if ((strcmp(argv[1],"Truss2") == 0)) { 	//////////////////////// mmc
    
    void *theEle = OPS_Truss2();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "tclelementcommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }																								
  } else if ((strcmp(argv[1],"CorotTruss2") == 0)) { 	//////////////////////// mmc
    
    void *theEle = OPS_CorotTruss2();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "tclelementcommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }																								
  }  

  else if (strcmp(argv[1], "AC3D8") == 0) {
      
    void *theEle = OPS_AC3D8HexWithSensitivity();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "tclelementcommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }																								
  }

  else if (strcmp(argv[1], "ASI3D8") == 0) {
    
    void *theEle = OPS_ASID8QuadWithSensitivity();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "tclelementcommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }
    
  } 

  else if (strcmp(argv[1], "AV3D4") == 0) {
    
    void *theEle = OPS_AV3D4QuadWithSensitivity();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "tclelementcommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }
  } 
  

  else if (strcmp(argv[1],"elastomericBearingBoucWenMod") == 0) {
    void *theEle = OPS_ElastomericBearingBoucWenMod3d();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "tclelementcommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }
  }


  else if (strcmp(argv[1], "VS3D4") == 0) {

    void *theEle = OPS_VS3D4WuadWithSensitivity();
    if (theEle != 0) 
      theElement = (Element *)theEle;
    else {
      opserr << "tclelementcommand -- unable to create element of type : " << argv[1] << endln;
      return TCL_ERROR;
    }
  }

  else if (strcmp(argv[1], "PFEMElement2DBuble") == 0) {
    ID info;
      void *theEle = OPS_PFEMElement2DBubble(info);
      if (theEle != 0) {
	  theElement = (Element*)theEle;
      } else {
	  opserr<<"tclelementcommand -- unable to create element of type : "
		<<argv[1]<<endln;
	  return TCL_ERROR;
      }
  }

  else if (strcmp(argv[1], "PFEMElement2DMini") == 0) {
      ID info;
      void *theEle = OPS_PFEMElement2Dmini(info);
      if (theEle != 0) {
	  theElement = (Element*)theEle;
      } else {
	  opserr<<"tclelementcommand -- unable to create element of type : "
		<<argv[1]<<endln;
	  return TCL_ERROR;
      }
  }

  else if (strcmp(argv[1], "PFEMElement2D") == 0) {
      void *theEle = OPS_PFEMElement2D();
      if (theEle != 0) {
	  theElement = (Element*)theEle;
      } else {
	  opserr<<"tclelementcommand -- unable to create element of type : "
		<<argv[1]<<endln;
	  return TCL_ERROR;
      }
  }

  else if (strcmp(argv[1], "CatenaryCable") == 0) {
      void *theEle = OPS_CatenaryCableElement();
      if (theEle != 0) {
    theElement = (Element*)theEle;
      } else {
    opserr<<"tclelementcommand -- unable to create element of type : "
    <<argv[1]<<endln;
    return TCL_ERROR;
      }
  }

  else if (strcmp(argv[1], "ShellANDeS") == 0) {
      void *theEle = OPS_ShellANDeS();
      if (theEle != 0) {
    theElement = (Element*)theEle;
      } else {
    opserr<<"tclelementcommand -- unable to create element of type : "
    <<argv[1]<<endln;
    return TCL_ERROR;
      }
  }


  else if (strcmp(argv[1], "LysmerTriangle") == 0) {
      void *theEle = OPS_LysmerTriangle();
      if (theEle != 0) {
    theElement = (Element*)theEle;
      } else {
    opserr<<"tclelementcommand -- unable to create element of type : "
    <<argv[1]<<endln;
    return TCL_ERROR;
      }
  }

  else if (strcmp(argv[1], "FourNodeTetrahedron") == 0) {
      void *theEle = OPS_FourNodeTetrahedron();
      if (theEle != 0) 
      {
        theElement = (Element*)theEle;
      } 
      else 
      {
        opserr<<"tclelementcommand -- unable to create element of type : " <<argv[1]<<endln;
        return TCL_ERROR;
      }

  }
  
  else if (strcmp(argv[1],"ZeroLengthVG_HG") == 0) {
      Element *theEle = (Element*) OPS_ZeroLengthVG_HG();
    if (theEle != 0) 
      theElement = theEle;
    else {
      opserr << "TclElementCommand -- unable to create element of type : " << argv[1] << endln;
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

  } else if (strcmp(argv[1],"dispBeamColumnInt") == 0) {
    int result = TclModelBuilder_addDispBeamColumnInt(clientData, interp, argc, argv,
						   theTclDomain, theTclBuilder);
    return result;

  } else if (strcmp(argv[1],"forceBeamColumn") == 0 || 
	     strcmp(argv[1],"dispBeamColumn") == 0  || 
	     strcmp(argv[1],"timoshenkoBeamColumn") == 0  || 
	     strcmp(argv[1],"forceBeamColumnCBDI") == 0  || 
	     strcmp(argv[1],"forceBeamColumnCSBDI") == 0  || 
	     strcmp(argv[1],"forceBeamColumnWarping") == 0  || 
	     strcmp(argv[1],"elasticForceBeamColumnWarping") == 0  || 
	     strcmp(argv[1],"dispBeamColumnNL") == 0  || 
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

   else if (strcmp(argv[1],"RJWatsonEqsBearing") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addRJWatsonEqsBearing(clientData, interp, argc, argv,
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

  else if (strcmp(argv[1],"elastomericBearing") == 0 ||
      strcmp(argv[1],"elastomericBearingPlasticity") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addElastomericBearingPlasticity(clientData, interp, argc, argv,
						       theTclDomain, theTclBuilder, eleArgStart);
    return result;
  }

  else if (strcmp(argv[1],"elastomericBearingBoucWen") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addElastomericBearingBoucWen(clientData, interp, argc, argv,
						       theTclDomain, theTclBuilder, eleArgStart);
    return result;
  }


  else if (strcmp(argv[1],"elastomericBearingUFRP") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addElastomericBearingUFRP(clientData, interp, argc, argv,
						       theTclDomain, theTclBuilder, eleArgStart);
    return result;
  }

  else if (strcmp(argv[1],"twoNodeLink") == 0) {
    int eleArgStart = 1;
    int result = TclModelBuilder_addTwoNodeLink(clientData, interp, argc, argv,
						theTclDomain, theTclBuilder, eleArgStart);
    return result;
  }
  
  else if ((strcmp(argv[1],"multipleShearSpring") == 0) ||
	   (strcmp(argv[1],"MSS") == 0)) {
    int result = TclModelBuilder_addMultipleShearSpring(clientData, interp, argc, argv,
							theTclDomain, theTclBuilder);
    return result;
  }

  else if ((strcmp(argv[1],"multipleNormalSpring") == 0) ||
	   (strcmp(argv[1],"MNS") == 0)) {
    int result = TclModelBuilder_addMultipleNormalSpring(clientData, interp, argc, argv,
    theTclDomain, theTclBuilder);
    return result;
  }

  else if (strcmp(argv[1],"KikuchiBearing") == 0) {
    int result = TclModelBuilder_addKikuchiBearing(clientData, interp, argc, argv,
						   theTclDomain, theTclBuilder);
    return result;
  }

  else if (strcmp(argv[1],"YamamotoBiaxialHDR") == 0) {
    int result = TclModelBuilder_addYamamotoBiaxialHDR(clientData, interp, argc, argv,
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
    // try loading new dynamic library containing a c+= class
    //
    
    void *libHandle;
    void *(*funcPtr)();
    int eleNameLength = strlen(argv[1]);
    char *tclFuncName = new char[eleNameLength+5];
    strcpy(tclFuncName, "OPS_");

    strcpy(&tclFuncName[4], argv[1]);

    opserr << "checking library: " << tclFuncName << endln;    
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


int TclModelBuilder_addMultipleShearSpring(ClientData clientData,
    Tcl_Interp *interp, int argc, TCL_Char **argv, Domain *theTclDomain,
    TclModelBuilder *theTclBuilder)
{

  // ensure the destructor has not been called
  if (theTclBuilder == 0)  {
    opserr << "WARNING builder has been destroyed - multipleShearSpring\n";    
    return TCL_ERROR;
  }

  // 3-dim, 6-dof
  int ndm = theTclBuilder->getNDM();
  int ndf = theTclBuilder->getNDF();

  if (ndm != 3 || ndf != 6) {
    opserr << "ndm=" << ndm << ", ndf=" << ndf << endln;
    opserr << "WARNING multipleShearSpring command only works when ndm is 3 and ndf is 6" << endln;
    return TCL_ERROR;
  }

  //arguments (necessary)
  int eleTag;
  int iNode;
  int jNode;
  int nSpring;
  int matTag;

  //material
  UniaxialMaterial *material = 0;
  UniaxialMaterial **theMaterials = 0;
  int recvMat = 0;

  //arguments (optional)
  double limDisp = 0.0;
  Vector oriX(0);
  Vector oriYp(3); oriYp(0) = 0.0; oriYp(1) = 1.0; oriYp(2) = 0.0;
  double mass = 0.0;

  //
  Element *theElement = 0;


  //error flag
  bool ifNoError = true;


  if (argc < 8)  { //element multipleShearSpring eleTag? iNode? jNode? nSpring? -mat matTag?

    opserr << "WARNING insufficient arguments\n";
    ifNoError = false;

  } else {

    //argv[2~5]
    if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK)  {
      opserr << "WARNING invalid multipleShearSpring eleTag\n";
      ifNoError = false;
    }

    if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK)  {
      opserr << "WARNING invalid iNode\n";
      ifNoError = false;
    }

    if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK)  {
      opserr << "WARNING invalid jNode\n";
      ifNoError = false;
    }

    if (Tcl_GetInt(interp, argv[5], &nSpring) != TCL_OK || nSpring <= 0)  {
      opserr << "WARNING invalid nSpring\n";
      ifNoError = false;
    }

    //argv[6~]
    for (int i=6; i<=(argc-1); i++) {
      
      double value;
      
      if (strcmp(argv[i],"-mat")==0 && (i+1)<=(argc-1)) { // -mat matTag?
	
	if (Tcl_GetInt(interp,argv[i+1], &matTag) != TCL_OK) {
	  opserr << "WARNING invalid matTag\n";
	  ifNoError = false;
	}

	material = OPS_getUniaxialMaterial(matTag);
	if (material == 0)  {
	  opserr << "WARNING material model not found\n";
	  opserr << "uniaxialMaterial: " << matTag << endln;
	  opserr << "multipleShearSpring element: " << eleTag << endln;
	  return TCL_ERROR;
	}

	//opserr << "org material " << material->getClassType() << "\n";
	recvMat++ ;
	i += 1;

      } else if (strcmp(argv[i],"-nMat")==0 && (i+nSpring)<=(argc-1)) { // -mat matTag?

	theMaterials = new UniaxialMaterial *[nSpring];
	for (int j=0; j<nSpring; j++) {
	  if (Tcl_GetInt(interp,argv[j+i+1], &matTag) != TCL_OK) {
	    opserr << "WARNING invalid matTag\n";
	    ifNoError = false;
	  }
	  
	  theMaterials[j] = OPS_getUniaxialMaterial(matTag);
	  if (theMaterials[j] == 0)  {
	    opserr << "WARNING material model not found\n";
	    opserr << "uniaxialMaterial: " << matTag << endln;
	    opserr << "multipleShearSpring element: " << eleTag << endln;
	    return TCL_ERROR;
	  }
	}
	//opserr << "org material " << material->getClassType() << "\n";
	recvMat++ ;
	i += nSpring;	

      } else if (strcmp(argv[i],"-orient")==0 && (i+6)<=(argc-1) && Tcl_GetDouble(interp,argv[i+4], &value) == TCL_OK) { // <-orient x1? x2? x3? yp1? yp2? yp3?>

	oriX.resize(3);

	for (int j=1; j<=3; j++) {
	  if (Tcl_GetDouble(interp, argv[i+j], &value) != TCL_OK )  {
	    opserr << "WARNING invalid -orient value\n";
	    ifNoError = false;
	  } else {
	    oriX(j-1) = value;
	  }
	}

	i += 3;

	for (int j=1; j<=3; j++) {
	  if (Tcl_GetDouble(interp, argv[i+j], &value) != TCL_OK )  {
	    opserr << "WARNING invalid -orient value\n";
	    ifNoError = false;
	  } else {
	    oriYp(j-1) = value;
	  }
	}

	i += 3;


      } else if (strcmp(argv[i],"-orient")==0 && (i+3)<=(argc-1)) { // <-orient yp1? yp2? yp3?> の読み込み	  

	for (int j=1; j<=3; j++) {
	  if (Tcl_GetDouble(interp, argv[i+j], &value) != TCL_OK )  {
	    opserr << "WARNING invalid -orient value\n";
	    ifNoError = false;
	  } else {
	    oriYp(j-1) = value;
	  }
	}

	i += 3;

      } else if (strcmp(argv[i],"-mass")==0 && (i+1)<=(argc-1)) { // <-mass m?> の読み込み
	
	if (Tcl_GetDouble(interp, argv[i+1], &mass) != TCL_OK || mass <= 0)  {
	  opserr << "WARNING invalid mass\n";
	  ifNoError = false;
	}

	i += 1;

      } else if (strcmp(argv[i],"-lim")==0 && (i+1)<=(argc-1)) { // <-lim limDisp?> の読み込み
	
	if (Tcl_GetDouble(interp, argv[i+1], &limDisp) != TCL_OK || limDisp < 0)  {
	  opserr << "WARNING invalid limDisp\n";
	  ifNoError = false;
	}

	i += 1;
	
      } else { //invalid option
	
	opserr << "WARNING invalid optional arguments \n";
	ifNoError = false;
	break;

      }

    }
    
  } //end input
  


  //confirm material
  if (recvMat != 1)  {
    opserr << "WARNING wrong number of -mat inputs\n";
    opserr << "got " << recvMat << " inputs, but want 1 input\n";
    ifNoError = false;
  }

  
  //if error detected
  if (!ifNoError) {
    //input:
    printCommand(argc, argv);
    //want:
    opserr << "Want: element multipleShearSpring eleTag? iNode? jNode? nSpring? -mat matTag? <-lim dsp> <-orient <x1? x2? x3?> yp1? yp2? yp3?> <-mass m?>\n";
    return TCL_ERROR;
  }
  

  // now create the multipleShearSpring
  if (theMaterials == 0) {
    theElement = new MultipleShearSpring(eleTag, iNode, jNode, nSpring, material, limDisp, oriYp, oriX, mass);
  } else {
    theElement = new MultipleShearSpring(eleTag, iNode, jNode, theMaterials, nSpring, limDisp, oriYp, oriX, mass);
    delete [] theMaterials;
  }

  if (theElement == 0)  {
    opserr << "WARNING ran out of memory creating element\n";
    opserr << "multipleShearSpring element: " << eleTag << endln;
    return TCL_ERROR;
  }
  
  // then add the multipleShearSpring to the domain
  if (theTclDomain->addElement(theElement) == false)  {
    opserr << "WARNING could not add element to the domain\n";
    opserr << "multipleShearSpring element: " << eleTag << endln;
    delete theElement;
    return TCL_ERROR;
  }       
  
  // if get here we have successfully created the multipleShearSpring and added it to the domain
  return TCL_OK;
}

static bool errDetected(bool ifNoError,char *msg){
  if (ifNoError){
    opserr << "" << endln;
    opserr << "========================================" << endln;
    opserr << " element : input error detected" << endln;
    opserr << "------------------------------" << endln;
  }
  opserr << "  " << msg << endln;
  return false;
};

int TclModelBuilder_addMultipleNormalSpring(ClientData clientData,
    Tcl_Interp *interp, int argc, TCL_Char **argv, Domain *theTclDomain,
    TclModelBuilder *theTclBuilder)
{

  // ensure the destructor has not been called
  if (theTclBuilder == 0)  {
    opserr << "WARNING builder has been destroyed - multipleNormalSpring\n";    
    return TCL_ERROR;
  }

  // 3-dim, 6-dof
  int ndm = theTclBuilder->getNDM();
  int ndf = theTclBuilder->getNDF();

  if (ndm != 3 || ndf != 6) {
    opserr << "ndm=" << ndm << ", ndf=" << ndf << endln;
    opserr << "WARNING multipleNormalSpring command only works when ndm is 3 and ndf is 6" << endln;
    return TCL_ERROR;
  }


  //arguments (necessary)
  int eleTag;
  int iNode;
  int jNode;
  int nDivide;

  //arguments (necessary, input with -???)
  int matTag;
  UniaxialMaterial *material;
  int shape;
  double size;

  //arguments (optional, input with -???)
  double lambda = -1.0;
  Vector oriX(0);
  Vector oriYp(3); oriYp(0) = 0.0; oriYp(1) = 1.0; oriYp(2) = 0.0;
  double mass = 0.0;

  // input comfirmation
  int recvMat    = 0;
  int recvShape  = 0;
  int recvSize   = 0;
  int recvLambda = 0;
  int recvOrient = 0;
  int recvMass   = 0;


  //
  Element *theElement = 0;


  // error flag
  bool ifNoError = true;



  if (argc < 6)  { //element multipleNormalSpring eleTag? iNode? jNode? nDivide?

    ifNoError = errDetected(ifNoError,"insufficient arguments");

  } else {

    //argv[2~5]
    if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK)  {
      ifNoError = errDetected(ifNoError,"invalid eleTag");
    }

    if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK)  {
      ifNoError = errDetected(ifNoError,"invalid iNode");
    }

    if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK)  {
      ifNoError = errDetected(ifNoError,"invalid jNode");
    }

    if (Tcl_GetInt(interp, argv[5], &nDivide) != TCL_OK || nDivide <= 0)  {
      ifNoError = errDetected(ifNoError,"invalid nDivide");
    }

    //argv[6~]
    for (int i=6; i<=(argc-1); i++) {
      
      double value;
      
      if (strcmp(argv[i],"-mat")==0 && (i+1)<=(argc-1)) { // -mat matTag?
	
	if (Tcl_GetInt(interp,argv[i+1], &matTag) != TCL_OK) {
	  ifNoError = errDetected(ifNoError,"invalid matTag");
	}


	material = OPS_getUniaxialMaterial(matTag);
	if (material == 0)  {
	  ifNoError = errDetected(ifNoError,"material model not found");
	}

	recvMat++ ;
	i += 1;


      } else if (strcmp(argv[i],"-shape")==0 && (i+1)<=(argc-1)) { // -shape shape?
	
	if (strcmp(argv[i+1],"round") == 0) {
	  shape = 1; //round shape
	} else if (strcmp(argv[i+1],"square") == 0) {
	  shape = 2; //square
	} else {
	  ifNoError = errDetected(ifNoError,"invalid shape (\"round\" or \"square\" are available)");
	}
	
	recvShape++ ;
	i += 1;
	
      } else if (strcmp(argv[i],"-size")==0 && (i+1)<=(argc-1)) { // -size size?
	
	if (Tcl_GetDouble(interp, argv[i+1], &size) != TCL_OK || size <= 0)  {
	  ifNoError = errDetected(ifNoError,"invalid size");
	}

	recvSize++ ;
	i += 1;

      } else if (strcmp(argv[i],"-lambda")==0 && (i+1)<=(argc-1)) { // <-lambda lambda?>
	
	if (Tcl_GetDouble(interp, argv[i+1], &lambda) != TCL_OK || lambda < 0)  {
	  ifNoError = errDetected(ifNoError,"invalid lambda");
	}

	recvLambda++ ;
	i += 1;

      } else if (strcmp(argv[i],"-orient")==0 && (i+6)<=(argc-1) && Tcl_GetDouble(interp,argv[i+4], &value) == TCL_OK) { // <-orient x1? x2? x3? yp1? yp2? yp3?>

	oriX.resize(3);

	for (int j=1; j<=3; j++) {
	  if (Tcl_GetDouble(interp, argv[i+j], &value) != TCL_OK )  {
	    ifNoError = errDetected(ifNoError,"invalid orient");
	  } else {
	    oriX(j-1) = value;
	  }
	}

	i += 3;

	for (int j=1; j<=3; j++) {
	  if (Tcl_GetDouble(interp, argv[i+j], &value) != TCL_OK )  {
	    ifNoError = errDetected(ifNoError,"invalid orient");
	  } else {
	    oriYp(j-1) = value;
	  }
	}

	recvOrient++ ;
	i += 3;

      } else if (strcmp(argv[i],"-orient")==0 && (i+3)<=(argc-1)) { // <-orient yp1? yp2? yp3?>

	for (int j=1; j<=3; j++) {
	  if (Tcl_GetDouble(interp, argv[i+j], &value) != TCL_OK )  {
	    ifNoError = errDetected(ifNoError,"invalid orient");
	  } else {
	    oriYp(j-1) = value;
	  }
	}

	recvOrient++ ;
	i += 3;

      } else if (strcmp(argv[i],"-mass")==0 && (i+1)<=(argc-1)) { // <-mass m?> の読み込み
	
	if (Tcl_GetDouble(interp, argv[i+1], &mass) != TCL_OK || mass <= 0)  {
	  ifNoError = errDetected(ifNoError,"invalid mass");
	}

	recvMass++ ;
	i += 1;
	
      } else { //invalid option
	
	ifNoError = errDetected(ifNoError,"invalid optional arguments");
	break;

      }

    }
    
  } // end input
  


  // input cofirmation
  // necessary arguments
  if (recvMat != 1)  {
    char buf[100];
    sprintf(buf,"wrong number of -mat inputs (got %d inputs, but want 1 input)",recvMat);
    ifNoError = errDetected(ifNoError,buf);
  }

  if (recvShape != 1)  {
    char buf[100];
    sprintf(buf,"wrong number of -shape inputs (got %d inputs, but want 1 input)",recvShape);
    ifNoError = errDetected(ifNoError,buf);
  }

  if (recvSize != 1)  {
    char buf[100];
    sprintf(buf,"wrong number of -size inputs (got %d inputs, but want 1 input)",recvSize);
    ifNoError = errDetected(ifNoError,buf);
  }

  // optional arguments
  if (recvLambda >= 2)  {
    char buf[100];
    sprintf(buf,"wrong number of -lambda inputs (got %d inputs, but want 1 input)",recvLambda);
    ifNoError = errDetected(ifNoError,buf);
  }

  if (recvOrient >= 2)  {
    char buf[100];
    sprintf(buf,"wrong number of -ori inputs (got %d inputs, but want 1 input)",recvOrient);
    ifNoError = errDetected(ifNoError,buf);
  }

  if (recvMass >= 2)  {
    char buf[100];
    sprintf(buf,"wrong number of -mass inputs (got %d inputs, but want 1 input)",recvMass);
    ifNoError = errDetected(ifNoError,buf);
  }

  
  // if error detected
  if (!ifNoError) {
    opserr << "------------------------------" << endln;
    //input:
    printCommand(argc, argv);
    //want:
    opserr << "Want: element multipleNormalSpring eleTag? iNode? jNode? nDivide? -mat matTag? -shape shape? -size size? <-lambda lambda?> <-orient <x1? x2? x3?> yp1? yp2? yp3?> <-mass m?>\n";
    opserr << "========================================" << endln;
    opserr << "" << endln;
    return TCL_ERROR;
  }
 
  

  // now create the multipleNormalSpring
  //theElement = new MultipleNormalSpring(eleTag, iNode, jNode, nDivide, material, shape, size, lambda, oriYp, oriX, mass);

  if (theElement == 0)  {
    opserr << "WARNING ran out of memory creating element\n";
    opserr << "multipleNormalSpring element: " << eleTag << endln;
    return TCL_ERROR;
  }
  
  // then add the multipleNormalSpring to the domain
  if (theTclDomain->addElement(theElement) == false)  {
    opserr << "WARNING could not add element to the domain\n";
    opserr << "multipleNormalSpring element: " << eleTag << endln;
    delete theElement;
    return TCL_ERROR;
  }       
  
  // if get here we have successfully created the multipleNormalSpring and added it to the domain
  return TCL_OK;
}

int TclModelBuilder_addKikuchiBearing(ClientData clientData,
    Tcl_Interp *interp, int argc, TCL_Char **argv, Domain *theTclDomain,
    TclModelBuilder *theTclBuilder)
{

  // ensure the destructor has not been called
  if (theTclBuilder == 0)  {
    opserr << "WARNING builder has been destroyed - KikuchiBearing\n";    
    return TCL_ERROR;
  }

  //3-dim, 6dof
  int ndm = theTclBuilder->getNDM();
  int ndf = theTclBuilder->getNDF();

  if (ndm != 3 || ndf != 6) {
    opserr << "ndm=" << ndm << ", ndf=" << ndf << endln;
    opserr << "WARNING KikuchiBearing command only works when ndm is 3 and ndf is 6" << endln;
    return TCL_ERROR;
  }

  //arguments (necessary)
  int eleTag;
  int iNode;
  int jNode;

  //arguments (necessary, input with -???)
  int shape;
  double size;
  double totalRubber;
  int nMSS;
  int matMSSTag;
  UniaxialMaterial *matMSS;
  int nMNS;
  int matMNSTag;
  UniaxialMaterial *matMNS;

  //arguments (optional, input with -???)
  double totalHeight = -1.0; //default: Norm(I->J)
  double limDisp = -1.0; //default: INF
  double lambda = -1.0; //default: INF
  Vector oriX(0); //default: local-x Vec(I->J)
  Vector oriYp(3); oriYp(0) = 0.0; oriYp(1) = 1.0; oriYp(2) = 0.0; //default: global-Y
  double mass = 0.0;
  bool ifPDInput = true;
  bool ifTilt = true;
  double adjCi = 0.5;
  double adjCj = 0.5;
  bool ifBalance = false;
  double limFo = -1.0; //default: INF
  double limFi = -1.0; //default: INF
  int nIter = 1;

  // input comfirmation
  int recvShape  = 0;
  int recvSize   = 0;
  int recvHeight = 0;
  int recvNMSS  = 0;
  int recvMatMSS  = 0;
  int recvLimDisp = 0;
  int recvNMNS  = 0;
  int recvMatMNS  = 0;
  int recvLambda = 0;
  int recvOrient = 0;
  int recvMass   = 0;
  int recvIfPD = 0;
  int recvIfTl = 0;
  int recvAdj = 0;
  int recvBal = 0;

  //
  Element *theElement = 0;


  //error flag
  bool ifNoError = true;



  if (argc < 5)  { //element KikuchiBearing eleTag? iNode? jNode?

    ifNoError = errDetected(ifNoError,"insufficient arguments");

  } else {

    //argv[2~4]
    if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK)  {
      ifNoError = errDetected(ifNoError,"invalid eleTag");
    }

    if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK)  {
      ifNoError = errDetected(ifNoError,"invalid iNode");
    }

    if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK)  {
      ifNoError = errDetected(ifNoError,"invalid jNode");
    }

    //argv[5~]
    for (int i=5; i<=(argc-1); i++) {
      
      double value;

      if (strcmp(argv[i],"-shape")==0 && (i+1)<=(argc-1)) { // -shape shape?
	
      	if (strcmp(argv[i+1],"round") == 0) {
	  shape = 1; //round
	} else if (strcmp(argv[i+1],"square") == 0) {
	  shape = 2; //square
	} else {
	  ifNoError = errDetected(ifNoError,"invalid shape (\"round\" or \"square\" are available)");
	}
	
	recvShape++ ;
	i += 1;


      } else if (strcmp(argv[i],"-size")==0 && (i+2)<=(argc-1)) { // -size size? totalRubber?
	
	if (Tcl_GetDouble(interp, argv[i+1], &size) != TCL_OK || size <= 0)  {
	  ifNoError = errDetected(ifNoError,"invalid size");
	}

	if (Tcl_GetDouble(interp, argv[i+2], &totalRubber) != TCL_OK || totalRubber <= 0)  {
	  ifNoError = errDetected(ifNoError,"invalid totalRubber");
	}

	recvSize++ ;
	i += 2;
	
      } else if (strcmp(argv[i],"-totalHeight")==0 && (i+1)<=(argc-1)) { // -totalHeight totalHeight?
	
	if (Tcl_GetDouble(interp, argv[i+1], &totalHeight) != TCL_OK || totalHeight <= 0)  {
	  ifNoError = errDetected(ifNoError,"invalid totalHeight");
	}

	recvHeight++ ;
	i += 1;

      } else if (strcmp(argv[i],"-nMSS")==0 && (i+1)<=(argc-1)) { // -nMSS nMSS?
	
	if (Tcl_GetInt(interp, argv[i+1], &nMSS) != TCL_OK || nMSS <= 0)  {
	  ifNoError = errDetected(ifNoError,"invalid nMSS");
	}

	recvNMSS++ ;
	i += 1;


      } else if (strcmp(argv[i],"-matMSS")==0 && (i+1)<=(argc-1)) { // -matMSS matMSSTag?
	
	if (Tcl_GetInt(interp,argv[i+1], &matMSSTag) != TCL_OK) {
	  ifNoError = errDetected(ifNoError,"invalid matMSSTag");
	}

	matMSS = OPS_getUniaxialMaterial(matMSSTag);
	if (matMSS == 0)  {
	  ifNoError = errDetected(ifNoError,"material for MSS model not found");
	}

	recvMatMSS++ ;
	i += 1;


      } else if (strcmp(argv[i],"-limDisp")==0 && (i+1)<=(argc-1)) { // <-limDisp limDisp?>
	
	if (Tcl_GetDouble(interp, argv[i+1], &limDisp) != TCL_OK || limDisp < 0)  {
	  ifNoError = errDetected(ifNoError,"invalid limDisp");
	}

	recvLimDisp++ ;
	i += 1;


      } else if (strcmp(argv[i],"-nMNS")==0 && (i+1)<=(argc-1)) { // -nMNS nMNS?
	
	if (Tcl_GetInt(interp, argv[i+1], &nMNS) != TCL_OK || nMNS <= 0)  {
	  ifNoError = errDetected(ifNoError,"invalid nMNS");
	}

	recvNMNS++ ;
	i += 1;


      } else if (strcmp(argv[i],"-matMNS")==0 && (i+1)<=(argc-1)) { // -matMNS matMNSTag?
	
	if (Tcl_GetInt(interp,argv[i+1], &matMNSTag) != TCL_OK) {
	  ifNoError = errDetected(ifNoError,"invalid matMNSTag");
	}

	matMNS = OPS_getUniaxialMaterial(matMNSTag);
	if (matMNS == 0)  {
	  ifNoError = errDetected(ifNoError,"material for MNS model not found");
	}

	recvMatMNS++ ;
	i += 1;


      } else if (strcmp(argv[i],"-lambda")==0 && (i+1)<=(argc-1)) { // <-lambda lambda?>
	
	if (Tcl_GetDouble(interp, argv[i+1], &lambda) != TCL_OK || lambda < 0)  {
	  ifNoError = errDetected(ifNoError,"invalid lambda");
	}

	recvLambda++ ;
	i += 1;


      } else if (strcmp(argv[i],"-orient")==0 && (i+6)<=(argc-1) && Tcl_GetDouble(interp,argv[i+4], &value) == TCL_OK) { // <-orient x1? x2? x3? yp1? yp2? yp3?>

	oriX.resize(3);

	for (int j=1; j<=3; j++) {
	  if (Tcl_GetDouble(interp, argv[i+j], &value) != TCL_OK )  {
	    ifNoError = errDetected(ifNoError,"invalid orient");
	  } else {
	    oriX(j-1) = value;
	  }
	}

	i += 3;

	for (int j=1; j<=3; j++) {
	  if (Tcl_GetDouble(interp, argv[i+j], &value) != TCL_OK )  {
	    ifNoError = errDetected(ifNoError,"invalid orient");
	  } else {
	    oriYp(j-1) = value;
	  }
	}

	recvOrient++ ;
	i += 3;

      } else if (strcmp(argv[i],"-orient")==0 && (i+3)<=(argc-1)) { // <-orient yp1? yp2? yp3?>

	for (int j=1; j<=3; j++) {
	  if (Tcl_GetDouble(interp, argv[i+j], &value) != TCL_OK )  {
	    ifNoError = errDetected(ifNoError,"invalid orient");
	  } else {
	    oriYp(j-1) = value;
	  }
	}

	recvOrient++ ;
	i += 3;

      } else if (strcmp(argv[i],"-mass")==0 && (i+1)<=(argc-1)) { // <-mass mass?>
	
	if (Tcl_GetDouble(interp, argv[i+1], &mass) != TCL_OK || mass <= 0)  {
	  ifNoError = errDetected(ifNoError,"invalid mass");
	}

	recvMass++ ;
	i += 1;

      } else if (strcmp(argv[i],"-noPDInput")==0) { // <-noPDInput>
	
	ifPDInput = false;

	recvIfPD++ ;
	i += 0;

      } else if (strcmp(argv[i],"-noTilt")==0) { // <-noTilt>
	
	ifTilt = false;

	recvIfTl++ ;
	i += 0;
	
      } else if (strcmp(argv[i],"-adjustPDOutput")==0 && (i+2)<=(argc-1)) { // -adjustPDOutput ci? cj?
	
	if (Tcl_GetDouble(interp, argv[i+1], &adjCi) != TCL_OK)  {
	  ifNoError = errDetected(ifNoError,"invalid ci");
	}

	if (Tcl_GetDouble(interp, argv[i+2], &adjCj) != TCL_OK)  {
	  ifNoError = errDetected(ifNoError,"invalid cj");
	}

	recvAdj++ ;
	i += 2;

      } else if (strcmp(argv[i],"-doBalance")==0 && (i+3)<=(argc-1)) { // -doBalance limFo? limFi? nIter?
	
	if (Tcl_GetDouble(interp, argv[i+1], &limFo) != TCL_OK || limFo <= 0)  {
	  ifNoError = errDetected(ifNoError,"invalid limFo");
	}

	if (Tcl_GetDouble(interp, argv[i+2], &limFi) != TCL_OK || limFi <= 0)  {
	  ifNoError = errDetected(ifNoError,"invalid limFi");
	}

	if (Tcl_GetInt(interp, argv[i+3], &nIter) != TCL_OK || nIter <= 0)  {
	  ifNoError = errDetected(ifNoError,"invalid nIter");
	}

	ifBalance = true;

	recvBal++ ;
	i += 3;


      } else { //invalid option
	
	ifNoError = errDetected(ifNoError,"invalid optional arguments");
	break;

      }

    }
    
  } //end input
  

  // input cofirmation
  // necessary arguments
  if (recvShape != 1)  {
    char buf[100];
    sprintf(buf,"wrong number of -shape inputs (got %d inputs, but want 1 input)",recvShape);
    ifNoError = errDetected(ifNoError,buf);
  }

  if (recvSize != 1)  {
    char buf[100];
    sprintf(buf,"wrong number of -size inputs (got %d inputs, but want 1 input)",recvSize);
    ifNoError = errDetected(ifNoError,buf);
  }

  if (recvNMSS != 1)  {
    char buf[100];
    sprintf(buf,"wrong number of -NMSS inputs (got %d inputs, but want 1 input)",recvNMSS);
    ifNoError = errDetected(ifNoError,buf);
  }

  if (recvMatMSS != 1)  {
    char buf[100];
    sprintf(buf,"wrong number of -matMSS inputs (got %d inputs, but want 1 input)",recvMatMSS);
    ifNoError = errDetected(ifNoError,buf);
  }

  if (recvNMNS != 1)  {
    char buf[100];
    sprintf(buf,"wrong number of -NMNS inputs (got %d inputs, but want 1 input)",recvNMNS);
    ifNoError = errDetected(ifNoError,buf);
  }

  if (recvMatMNS != 1)  {
    char buf[100];
    sprintf(buf,"wrong number of -matMNS inputs (got %d inputs, but want 1 input)",recvMatMNS);
    ifNoError = errDetected(ifNoError,buf);
  }


  //optional arguments
  if (recvHeight >= 2)  {
    char buf[100];
    sprintf(buf,"wrong number of -totalHeight inputs (got %d inputs, but want 1 input)",recvHeight);
    ifNoError = errDetected(ifNoError,buf);
  }

  if (recvLimDisp >= 2)  {
    char buf[100];
    sprintf(buf,"wrong number of -limDisp inputs (got %d inputs, but want 1 input)",recvLimDisp);
    ifNoError = errDetected(ifNoError,buf);
  }

  if (recvLambda >= 2)  {
    char buf[100];
    sprintf(buf,"wrong number of -lambda inputs (got %d inputs, but want 1 input)",recvLambda);
    ifNoError = errDetected(ifNoError,buf);
  }

  if (recvOrient >= 2)  {
    char buf[100];
    sprintf(buf,"wrong number of -ori inputs (got %d inputs, but want 1 input)",recvOrient);
    ifNoError = errDetected(ifNoError,buf);
  }

  if (recvMass >= 2)  {
    char buf[100];
    sprintf(buf,"wrong number of -mass inputs (got %d inputs, but want 1 input)",recvMass);
    ifNoError = errDetected(ifNoError,buf);
  }

  if (recvIfPD >= 2)  {
    char buf[100];
    sprintf(buf,"wrong number of -noPDInput inputs (got %d inputs, but want 1 input)",recvIfPD);
    ifNoError = errDetected(ifNoError,buf);
  }

  if (recvIfTl >= 2)  {
    char buf[100];
    sprintf(buf,"wrong number of -noTilt inputs (got %d inputs, but want 1 input)",recvIfTl);
    ifNoError = errDetected(ifNoError,buf);
  }

  if (recvAdj >= 2)  {
    char buf[100];
    sprintf(buf,"wrong number of -adjustPDOutput inputs (got %d inputs, but want 1 input)",recvAdj);
    ifNoError = errDetected(ifNoError,buf);
  }

  if (recvBal >= 2)  {
    char buf[100];
    sprintf(buf,"wrong number of -doBalance inputs (got %d inputs, but want 1 input)",recvBal);
    ifNoError = errDetected(ifNoError,buf);
  }

  
  //if error detected
  if (!ifNoError) {
    opserr << "------------------------------" << endln;
    //input:
    printCommand(argc, argv);
    //want:
    opserr << "Want: element KikuchiBearing eleTag? iNode? jNode?\n";
    opserr << "                             -shape shape? -size size? totalRubber? <-totalHeight totalHeight?>\n";
    opserr << "                             -nMSS nMSS? -matMSS matMSSTag? <-lim limDisp?>\n";
    opserr << "                             -nMNS nMNS? -matMNS matMNSTag? <-lambda lambda?>\n";
    opserr << "                             <-orient <x1? x2? x3?> yp1? yp2? yp3?> <-mass m?>\n";
    opserr << "                             <-noPDInput> <-noTilt> <-adjustPDOutput ci? cj?> <-doBalance limFo? limFi? nIter?>\n";
    opserr << "========================================" << endln;
    opserr << "" << endln;
    return TCL_ERROR;
  }

  // now create the KikuchiBearing
  theElement = new KikuchiBearing(eleTag, iNode, jNode, shape, size, totalRubber, totalHeight, nMSS, matMSS, limDisp, nMNS, matMNS, lambda, oriYp, oriX, mass, ifPDInput, ifTilt, adjCi, adjCj, ifBalance, limFo, limFi, nIter);

  if (theElement == 0)  {
    opserr << "WARNING ran out of memory creating element\n";
    opserr << "KikuchiBearing element: " << eleTag << endln;
    return TCL_ERROR;
  }
  
  // then add the KikuchiBearing to the domain
  if (theTclDomain->addElement(theElement) == false)  {
    opserr << "WARNING could not add element to the domain\n";
    opserr << "KikuchiBearing element: " << eleTag << endln;
    delete theElement;
    return TCL_ERROR;
  }       
  
  // if get here we have successfully created the KikuchiBearing and added it to the domain
  return TCL_OK;
}

int TclModelBuilder_addYamamotoBiaxialHDR(ClientData clientData,
				      Tcl_Interp *interp, int argc, TCL_Char **argv, Domain *theTclDomain,
				      TclModelBuilder *theTclBuilder)
{

  // ensure the destructor has not been called
  if (theTclBuilder == 0)  {
    opserr << "WARNING builder has been destroyed - YamamotoBiaxialHDR\n";    
    return TCL_ERROR;
  }


  // 3-dim, 6-dof
  int ndm = theTclBuilder->getNDM();
  int ndf = theTclBuilder->getNDF();

  if (ndm != 3 || ndf != 6) {
    opserr << "ndm=" << ndm << ", ndf=" << ndf << endln;
    opserr << "WARNING YamamotoBiaxialHDR command only works when ndm is 3 and ndf is 6" << endln;
    return TCL_ERROR;
  }


  //arguments (necessary)
  int eleTag;
  int iNode;
  int jNode;

  int Tp = 1;
  double DDo;
  double DDi;
  double Hr;

  //arguments (optional)
  double Cr=1.0;
  double Cs=1.0;
  Vector oriX(0);
  Vector oriYp(3); oriYp(0) = 0.0; oriYp(1) = 1.0; oriYp(2) = 0.0;
  double mass = 0.0;


  //
  Element *theElement = 0;


  //error flag
  bool ifNoError = true;


  if (argc < 9)  { //element YamamotoBiaxialHDR eleTag? iNode? jNode? Tp? DDo? DDi? Hr?
    // argc =            1           2             3      4      5     6   7    8    9
    // argv =       argv[0]      argv[1]      argv[2]  argv[3] ................. argv[8]
    opserr << "WARNING insufficient arguments\n";
    ifNoError = false;

  } else {


    //argv[2~8]
    if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK)  {
      opserr << "WARNING invalid YamamotoBiaxialHDR eleTag\n";
      ifNoError = false;
    }

    // iNode
    if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK)  {
      opserr << "WARNING invalid iNode\n";
      ifNoError = false;
    }

    // jNode
    if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK)  {
      opserr << "WARNING invalid jNode\n";
      ifNoError = false;
    }

    // Tp
    if (strcmp(argv[5],"1") == 0) {
      Tp = 1; // Bridgestone X0.6R (EESD version)
    } else {
      opserr << "WARNING invalid YamamotoBiaxialHDR Tp" << endln;
      ifNoError = false;
    }

    // DDo
    if (Tcl_GetDouble(interp, argv[6], &DDo) != TCL_OK || DDo <= 0.0) {
      opserr << "WARNING invalid YamamotoBiaxialHDR DDo" << endln;
      ifNoError = false;
    }

    // DDi
    if (Tcl_GetDouble(interp, argv[7], &DDi) != TCL_OK || DDi < 0.0) {
      opserr << "WARNING invalid YamamotoBiaxialHDR DDi" << endln;
      ifNoError = false;
    }

    // Hr
    if (Tcl_GetDouble(interp, argv[8], &Hr) != TCL_OK || Hr <= 0.0) {
      opserr << "WARNING invalid YamamotoBiaxialHDR Hr" << endln;
      ifNoError = false;
    }

    // check print--------------------------------------------/
    //  opserr << "   \n";
    //  opserr << "TclModelBuilder_addYamamotoBiaxialHDR()\n";
    //  opserr << "  tp  = " << Tp << endln;
    //  opserr << "  ddo = " << DDo << endln;
    //  opserr << "  ddi = " << DDi << endln;
    //  opserr << "  hr  = " << Hr << endln;
    //------------------------------------------------------

    // argv[9~]
    for (int i=9; i<=(argc-1); i++) {
      double value;

      if (strcmp(argv[i],"-orient")==0 && (i+6)<=(argc-1) && Tcl_GetDouble(interp,argv[i+4], &value) == TCL_OK) { // <-orient x1? x2? x3? yp1? yp2? yp3?>

	oriX.resize(3);

	// x1, x2, x3
	for (int j=1; j<=3; j++) {
	  if (Tcl_GetDouble(interp, argv[i+j], &value) != TCL_OK )  {
	    opserr << "WARNING invalid -orient value\n";
	    ifNoError = false;
	  } else {
	    oriX(j-1) = value;
	  }
	}
	
	i += 3;
	
	// yp1, yp2, yp3
	for (int j=1; j<=3; j++) {
	  if (Tcl_GetDouble(interp, argv[i+j], &value) != TCL_OK )  {
	    opserr << "WARNING invalid -orient value\n";
	    ifNoError = false;
	  } else {
	    oriYp(j-1) = value;
	  }
	}
	
	i += 3;
	
      } else if (strcmp(argv[i],"-orient")==0 && (i+3)<=(argc-1)) { // <-orient yp1? yp2? yp3?>
	
	for (int j=1; j<=3; j++) {
	  if (Tcl_GetDouble(interp, argv[i+j], &value) != TCL_OK )  {
	    opserr << "WARNING invalid -orient value\n";
	    ifNoError = false;
	  } else {
	    oriYp(j-1) = value;
	  }
	}
	
	i += 3;
	
      } else if (strcmp(argv[i],"-mass")==0 && (i+1)<=(argc-1)) { // <-mass m?>
	
	if (Tcl_GetDouble(interp, argv[i+1], &mass) != TCL_OK || mass <= 0)  {
	  opserr << "WARNING invalid mass\n";
	  ifNoError = false;
	}

	i += 1;

      } else if (strcmp(argv[i],"-coRS")==0 && (i+2)<=(argc-1)) { // <-coRS cr? cs?>
	
	if (Tcl_GetDouble(interp, argv[i+1], &Cr) != TCL_OK || Cr <= 0)  {
	  opserr << "WARNING invalid cr\n";
	  ifNoError = false;
	}
	if (Tcl_GetDouble(interp, argv[i+2], &Cs) != TCL_OK || Cs <= 0)  {
	  opserr << "WARNING invalid cs\n";
	  ifNoError = false;
	}

	i += 2;

	
      } else {
	
	opserr << "WARNING invalid optional arguments \n";
	ifNoError = false;
	break;
	
      }
    }

  } //end input

  
  //if error detected
  if (!ifNoError) {
    //input:
    printCommand(argc, argv);
    //want:
    opserr << "Want: element YamamotoBiaxialHDR eleTag? iNode? jNode? Tp? DDo? DDi? Hr?  <-coRS cr? cs?> <-orient <x1? x2? x3?> y1? y2? y3?> <-mass m?>\n";
    return TCL_ERROR;
  }
  

  // now create the YamamotoBiaxialHDR
  theElement = new YamamotoBiaxialHDR(eleTag, iNode, jNode, Tp, DDo, DDi, Hr, Cr, Cs, oriYp, oriX, mass);

  if (theElement == 0)  {
    opserr << "WARNING ran out of memory creating element\n";
    opserr << "YamamotoBiaxialHDR element: " << eleTag << endln;
    return TCL_ERROR;
  }
  
  // then add the YamamotoBiaxialHDR to the domain
  if (theTclDomain->addElement(theElement) == false)  {
    opserr << "WARNING could not add element to the domain\n";
    opserr << "YamamotoBiaxialHDR element: " << eleTag << endln;
    delete theElement;
    return TCL_ERROR;
  }       
  
  // if get here we have successfully created the YamamotoBiaxialHDR and added it to the domain
  return TCL_OK;
  }



  int
	  TclModelBuilder_addWheelRail(ClientData clientData, Tcl_Interp *interp, int argc,
		  TCL_Char **argv, Domain *theTclDomain, TclModelBuilder *theTclBuilder,
		  int eleArgStart)
  {
	  // ensure the destructor has not been called - 
	  if (theTclBuilder == 0) {
		  opserr << "WARNING builder has been destroyed - elasticBeamColumn \n";
		  return TCL_ERROR;
	  }

	  int ndm = theTclBuilder->getNDM();
	  int ndf = theTclBuilder->getNDF();

	  Element *theElement = 0;

	  int pTag, pnLoad;
	  //-------------Beginning of a 2D wheel-rail element(By Quan Gu, Yongdou Liu, et al.) on 2018/10/29
	  if (ndm == 2) {

		  // check plane frame problem has 3 dof per node
		  if (ndf != 3) {
			  opserr << "WARNING invalid ndf: " << ndf;
			  opserr << ", for plane problem need 3 - elasticBeamColumn \n";
			  return TCL_ERROR;
		  }

		  // check the number of arguments
		  if ((argc - eleArgStart) < 8) {
			  opserr << "WARNING bad command - want: elasticBeamColumn beamId iNode jNode A E I <alpha> <d> transTag <-mass m> <-cMass>\n";
			  printCommand(argc, argv);
			  return TCL_ERROR;
		  }

		  // get the id, end nodes, and section properties
			int pNd1,beamId, transTag;
			double A, E, I;

			double pDeltT, pVel, pInitLocation, pRWheel, pI, pE, pA;

			if (Tcl_GetInt(interp, argv[1 + eleArgStart], &pTag) != TCL_OK) {
				opserr << "WARNING invalid elemId: " << argv[1 + eleArgStart];
				opserr << " - WheelRail elemId iNode jNode";
				return TCL_ERROR;
			}

			if (Tcl_GetDouble(interp, argv[2 + eleArgStart], &pDeltT) != TCL_OK) {
				opserr << "WARNING invalid iNode - WheelRail " << pTag << " iNode jNode A E I\n";
				return TCL_ERROR;
			}

			if (Tcl_GetDouble(interp, argv[3 + eleArgStart], &pVel) != TCL_OK) {
				opserr << "WARNING invalid jNode - WheelRail " << pTag << " iNode jNode A E I\n";
				return TCL_ERROR;
			}

			if (Tcl_GetDouble(interp, argv[4 + eleArgStart], &pInitLocation) != TCL_OK) {
				opserr << "WARNING invalid A - WheelRail " << pTag << " iNode jNode A E I\n";
				return TCL_ERROR;
			}

			if (Tcl_GetInt(interp, argv[5 + eleArgStart], &pNd1) != TCL_OK) {
				opserr << "WARNING invalid E - WheelRail " << pTag << " iNode jNode A E I\n";
				return TCL_ERROR;
			}

			if (Tcl_GetDouble(interp, argv[6 + eleArgStart], &pRWheel) != TCL_OK) {
				opserr << "WARNING invalid A - WheelRail " << pTag << " iNode jNode A E I\n";
				return TCL_ERROR;
			}

			if (Tcl_GetDouble(interp, argv[7 + eleArgStart], &pI) != TCL_OK) {
				opserr << "WARNING invalid A - WheelRail " << pTag << " iNode jNode A E I\n";
				return TCL_ERROR;
			}

			if (Tcl_GetDouble(interp, argv[8 + eleArgStart], &pE) != TCL_OK) {
				opserr << "WARNING invalid A - WheelRail " << pTag << " iNode jNode A E I\n";
				return TCL_ERROR;
			}

			if (Tcl_GetDouble(interp, argv[9 + eleArgStart], &pA) != TCL_OK) {
				opserr << "WARNING invalid A - WheelRail " << pTag << " iNode jNode A E I\n";
				return TCL_ERROR;
			}

			if (Tcl_GetInt(interp, argv[10 + eleArgStart], &transTag) != TCL_OK) {
				opserr << "WARNING invalid transTag - elasticBeamColumn " << beamId << " iNode jNode A E I alpha d transTag\n";
				return TCL_ERROR;
			}
			CrdTransf *theTransRWheel = OPS_getCrdTransf(transTag);

			if (Tcl_GetInt(interp, argv[11 + eleArgStart], &pnLoad) != TCL_OK) {
				opserr << "WARNING invalid I - WheelRail " << pTag << " iNode jNode A E I\n";
				return TCL_ERROR;
			}
			//----------------------------------
			Vector * pNodeList = 0;
			Vector * pDeltaYList = 0;
			Vector * pDeltaYLocationList = 0;

			if (strcmp(argv[12 + eleArgStart], "-NodeList") == 0) {
				int pathSize;
				TCL_Char **pathStrings;

				int debug = Tcl_SplitList(interp, argv[13 + eleArgStart], &pathSize, &pathStrings);

				if (Tcl_SplitList(interp, argv[13 + eleArgStart], &pathSize, &pathStrings) != TCL_OK) {
					opserr << "WARNING problem splitting path list " << argv[13 + eleArgStart] << " - ";
					opserr << " NodeList -values {path} ... \n";
					return TCL_OK;
				}
				pNodeList = new Vector(pathSize);
				for (int i = 0; i < pathSize; i++) {
					double value;
					int debug = Tcl_GetDouble(interp, pathStrings[i], &value);
					if (Tcl_GetDouble(interp, pathStrings[i], &value) != TCL_OK) {
						opserr << "WARNING problem reading path data value " <<
							pathStrings[i] << " - ";
						opserr << " -strain {path} ... \n";
						return 0;
					}
					(*pNodeList)(i) = value;
				} //for
			}
			if (strcmp(argv[14 + eleArgStart], "-DeltaYList") == 0) {
				int pathSize;
				TCL_Char **pathStrings;
				if (Tcl_SplitList(interp, argv[15 + eleArgStart], &pathSize, &pathStrings) != TCL_OK) {
					opserr << "WARNING problem splitting path list " << argv[15 + eleArgStart] << " - ";
					opserr << " NodeList -values {path} ... \n";
					return TCL_OK;
				}
				pDeltaYList = new Vector(pathSize);
				for (int i = 0; i < pathSize; i++) {
					double value;
					if (Tcl_GetDouble(interp, pathStrings[i], &value) != TCL_OK) {
						opserr << "WARNING problem reading path data value " <<
							pathStrings[i] << " - ";
						opserr << " -strain {path} ... \n";
						return 0;
					}
					(*pDeltaYList)(i) = value;
				} //for
			}
			if (strcmp(argv[16 + eleArgStart], "-LocationList") == 0) {
				int pathSize;
				TCL_Char **pathStrings;
				if (Tcl_SplitList(interp, argv[17 + eleArgStart], &pathSize, &pathStrings) != TCL_OK) {
					opserr << "WARNING problem splitting path list " << argv[17 + eleArgStart] << " - ";
					opserr << " NodeList -values {path} ... \n";
					return TCL_OK;
				}
				pDeltaYLocationList = new Vector(pathSize);
				for (int i = 0; i < pathSize; i++) {
					double value;
					if (Tcl_GetDouble(interp, pathStrings[i], &value) != TCL_OK) {
						opserr << "WARNING problem reading path data value " <<
							pathStrings[i] << " - ";
						opserr << " -strain {path} ... \n";
						return 0;
					}
					(*pDeltaYLocationList)(i) = value;
				} //for
			}
			  theElement = new WheelRail(pTag, pDeltT, pVel, pInitLocation, pNd1, pRWheel, pI, pE, pA, theTransRWheel, pnLoad,
				  pNodeList, pDeltaYList, pDeltaYLocationList);


			if (theElement == 0) {
				opserr << "WARNING ran out of memory creating beam - WheelRail ";
				opserr << beamId << " iNode jNode A E I\n";
				return TCL_ERROR;
			}

		  } //--------------End of a 2D wheel-rail element(By Quan Gu, Yongdou Liu, et al.) on 2018/10/29 */
		  else if (ndm == 3) {

		  opserr << "Have not developed yet." << endln;
		  return TCL_ERROR;

		  }

		  //add the WheelRail element to the Domain
		  if (theTclDomain->addElement(theElement) == false) {
			  opserr << "WARNING could not add element to the domain\n";
			  opserr << "YamamotoBiaxialHDR element: " << pTag << endln;
			  delete theElement;
			  return TCL_ERROR;
		  }

		  return 0;

  }
