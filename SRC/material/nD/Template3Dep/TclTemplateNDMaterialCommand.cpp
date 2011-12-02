/*
################################################################################
# COPYRIGHT (C):     :-))                                                      #
# PROJECT:           Object Oriented Finite Element Program                    #
# PURPOSE:           General platform for elaso-plastic constitutive model     #
#                    implementation                                            #
# CLASS:             Template3Dep (the base class for all material point)     #
#                                                                              #
# VERSION:                                                                     #
# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
# TARGET OS:         DOS || UNIX || . . .                                      #
# DESIGNER(S):       Boris Jeremic, Zhaohui Yang                               #
# PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                               #
#                                                                              #
#                                                                              #
# DATE:              
# UPDATE HISTORY:    
#                                                                              #
#                                                                              #
#                                                                              #
#                                                                              #
*/

#include <tcl.h>
#include <stdlib.h>
#include <string.h>
#include <iOPS_Stream.h>
#include <ErrorHandler.h>

// includes for Drucker-Prager and von Mises yield surface and potential surface classes
#include <DP_YS.h>
#include <DP_PS.h>
#include <CAM_YS.h>
#include <CAM_PS.h>
#include <VM_YS.h>
#include <VM_PS.h>
#include <MD_YS.h>
#include <MD_PS.h>

// includes for scalar evolution law class
#include <EL_S.h>
#include <EL_LEp.h>
#include <EL_LEeq.h>
#include <EL_NLEeq.h>

// includes for tensorial evolution law class
#include <EL_T.h>
#include <EL_LEij.h>
#include <EL_NLEij.h>

// include for elasto-plastic state class
#include <EPState.h>

// include for Material class
#include <Template3Dep.h>


YieldSurface *
TclYieldSurfaceCommand(ClientData clientData, Tcl_Interp *interp, 
		       int argc, char **argv);

PotentialSurface *
TclPotentialSurfaceCommand(ClientData clientData, Tcl_Interp *interp, 
		       int argc, char **argv); 
		       
EPState *
TclEPStateCommand(ClientData clientData, Tcl_Interp *interp, 
		  int argc, char **argv); 
		       
		       
EvolutionLaw_S *
TclEvolutionLawSCommand(ClientData clientData, Tcl_Interp *interp, 
			int argc, char **argv);	
			
EvolutionLaw_T *
TclEvolutionLawTCommand(ClientData clientData, Tcl_Interp *interp, 
			int argc, char **argv);	

NDMaterial *
TcladdTemplate3Dep(ClientData clientData, Tcl_Interp *interp,  int argc, 
		   char **argv, Domain*theTclDomain, int startArg)
{
  int tag;
  YieldSurface     *YS   =0;        
  PotentialSurface *PS   =0;
  EPState          *EPS  =0;
  EvolutionLaw_S   *ELS1 =0; 
  EvolutionLaw_S   *ELS2 =0; 
  EvolutionLaw_S   *ELS3 =0; 
  EvolutionLaw_S   *ELS4 =0; 
  EvolutionLaw_T   *ELT1 =0;
  EvolutionLaw_T   *ELT2 =0;
  EvolutionLaw_T   *ELT3 =0;
  EvolutionLaw_T   *ELT4 =0;   

  
  if (Tcl_GetInt(interp, argv[startArg+1], &tag) != TCL_OK) {
        g3ErrorHandler->warning("invalid tag: %s for Template3Dep material",
			  argv[startArg+1]);
	return 0;		
  }
  
  int endMarker = startArg+2;
  while (endMarker < argc) {

      int numListArgs =0;
      char **listArgs =0;

      // create the yield surface object on the flag -ys
      if ((strcmp(argv[endMarker],"-ys") == 0) || 
	  (strcmp(argv[endMarker],"-YS") == 0)) {

	  endMarker++;
	  
	  if (Tcl_SplitList(interp, argv[endMarker], 
			    &numListArgs, &listArgs) != TCL_OK) {

	      g3ErrorHandler->warning("problem splitting list: %s %s %d",
				      argv[endMarker],
				      "for Template3Dep material", tag);
	      return 0;
	  }	  
	  
	  YS = TclYieldSurfaceCommand(clientData, interp, numListArgs, listArgs);
	  if (YS == 0) {
	      g3ErrorHandler->warning("problem splitting list: %s %s %d",
				      argv[endMarker],
				      "for Template3Dep material", tag);
	      return 0;
	  }
	  endMarker++;
      }	

      // create the potential surface object on the flag -ys
      if ((strcmp(argv[endMarker],"-ps") == 0) || 
	  (strcmp(argv[endMarker],"-PS") == 0)) {

	  endMarker++;
	  
	  if (Tcl_SplitList(interp, argv[endMarker], 
			    &numListArgs, &listArgs) != TCL_OK) {

	      g3ErrorHandler->warning("problem splitting list: %s %s %d",
				      argv[endMarker],
				      "for Template3Dep material", tag);
	      return 0;
	  }	  
	  
	  PS = TclPotentialSurfaceCommand(clientData, interp, numListArgs, listArgs);
	  if (PS == 0) {
	      g3ErrorHandler->warning("problem splitting list: %s %s %d",
				      argv[endMarker],
				      "for Template3Dep material", tag);
	      return 0;
	  }
	  endMarker++;
      }	


      // create the EPState object on the flag -eps
      if ((strcmp(argv[endMarker],"-eps") == 0) || 
	  (strcmp(argv[endMarker],"-EPS") == 0)) {

	  endMarker++;
	  
	  if (Tcl_SplitList(interp, argv[endMarker], 
			    &numListArgs, &listArgs) != TCL_OK) {

	      g3ErrorHandler->warning("problem splitting list: %s %s %d",
				      argv[endMarker],
				      "for Template3Dep material", tag);
	      return 0;
	  }	  
	  
	  EPS = TclEPStateCommand(clientData, interp, numListArgs, listArgs);
	  if (PS == 0) {
	      g3ErrorHandler->warning("problem splitting list: %s %s %d",
				      argv[endMarker],
				      "for Template3Dep material", tag);
	      return 0;
	  }
	  endMarker++;
      }	
	  
	  
      else {
	  endMarker++;
      }

      if (listArgs != 0) {
	  Tcl_Free((char *) listArgs);
      }
  }

  NDMaterial *theMaterial = new Template3Dep(tag,
					     YS,
					     PS,
					     EPS,
					     ELS1,
					     ELS2,
					     ELS3,
					     ELS4,
					     ELT1,
					     ELT2,
					     ELT3,
					     ELT4);
  
  return theMaterial;
}



YieldSurface *
TclYieldSurfaceCommand(ClientData clientData, Tcl_Interp *interp, 
		       int argc, char **argv)
{
  // parse args and return a DruckerPrager yield surface
  if ((strcmp(argv[0],"DruckerPrager") == 0) || 
      (strcmp(argv[0],"DP") == 0)) {

    return new DPYieldSurface();
  }

  // parse args and return a Cam Clay yield surface
  else if ((strcmp(argv[0],"CamClay") == 0) || 
      (strcmp(argv[0],"Cam") == 0)) {

    double mp = 0.0;
    if (argc == 2) {
      if (Tcl_GetDouble(interp, argv[1], &mp) != TCL_OK) {
        g3ErrorHandler->warning("invalid M: %s for -PS CamClay M",
				argv[1]);
	return 0;		
      }
    }
    
    // create the object & return it
    return new CAMYieldSurface(mp);
  }

  // parse args and return a VonMises yield surface
  else if ((strcmp(argv[0],"VonMises") == 0) || 
      (strcmp(argv[0],"VM") == 0)) {

    return new VMYieldSurface();

  }

  else {
    g3ErrorHandler->warning("unkown Yield Surface type: %s\n",
			    argv[0]);
    return 0;
  }
}



PotentialSurface *
TclPotentialSurfaceCommand(ClientData clientData, Tcl_Interp *interp, 
		       int argc, char **argv)

{  
  // parse args and return a DruckerPrager potential surface
  if ((strcmp(argv[0],"DruckerPrager") == 0) || 
      (strcmp(argv[0],"DP") == 0)) {

    double a2d = 0.0;
    if (argc == 2) {
      if (Tcl_GetDouble(interp, argv[1], &a2d) != TCL_OK) {
        g3ErrorHandler->warning("invalid a2d: %s for -PS DruckerPrage a2d",
				argv[1]);
	return 0;		
      }
    }
    
    // create the object & return it
    return new DPPotentialSurface(a2d);
  }
  
  // parse args and return a CamClay potential surface
  else if ((strcmp(argv[0],"CamClay") == 0) || 
      (strcmp(argv[0],"Cam") == 0)) {

    double mp = 0.0;
    if (argc == 2) {
      if (Tcl_GetDouble(interp, argv[1], &mp) != TCL_OK) {
        g3ErrorHandler->warning("invalid M: %s for -PS CamClay M",
				argv[1]);
	return 0;		
      }
    }
    
    // create the object & return it
    return new CAMPotentialSurface(mp);
  }

  // parse args and return a VonMises potential surface
  else if ((strcmp(argv[0],"VonMises") == 0) || 
      (strcmp(argv[0],"VM") == 0)) {

    return new VMPotentialSurface();

  }

  // unknown type return error
  else {
    g3ErrorHandler->warning("unkown Potential Surface type: %s\n",
			    argv[0]);
    return 0;
  }
}		       


EPState *
TclEPStateCommand(ClientData clientData, Tcl_Interp *interp, 
		  int argc, char **argv)
{
  double 	      Eod;
  double              Ed;
  double              nu;
  double              rho;
  stresstensor        stressp;
  straintensor        strainp;
  straintensor        Estrainp;
  straintensor        Pstrainp;
  straintensor        dEstrainp;
  straintensor        dPstrainp;
  int                 NScalarp;
  double              Scalarp;
  int                 NTensorp;
  stresstensor        Tensorp;
  tensor              Eepp;
  stresstensor        Stress_commitp;
  straintensor        Strain_commitp;	  
  double              Scalar_commitp;
  stresstensor        Tensor_commitp; 
  tensor              Eep_commitp;
  stresstensor        Stress_initp;    
  straintensor        Strain_initp;	 
  double              Scalar_initp;
  stresstensor        Tensor_initp;
  tensor              Eep_initp; 
  bool                Convergedp;

  return 0;
}		       
		       
EvolutionLaw_S *
TclEvolutionLawSCommand(ClientData clientData, Tcl_Interp *interp, 
			int argc, char **argv)
{
  return 0;
}
			
EvolutionLaw_T *
TclEvolutionLawTCommand(ClientData clientData, Tcl_Interp *interp, 
			int argc, char **argv)
{
  return 0;
}











