//===============================================================================
//# COPYRIGHT (C): Woody's license (by BJ):
//                 ``This    source  code is Copyrighted in
//                 U.S.,  for  an  indefinite  period,  and anybody
//                 caught  using it without our permission, will be
//                 mighty good friends of ourn, cause we don't give
//                 a  darn.  Hack it. Compile it. Debug it. Run it.
//                 Yodel  it.  Enjoy it. We wrote it, that's all we
//                 wanted to do.''
//
//# PROJECT:           Object Oriented Finite Element Program
//# PURPOSE:           Finite Deformation Hyper-Elastic classes
//# CLASS:
//#
//# VERSION:           0.6_(1803398874989) (golden section)
//# LANGUAGE:          C++
//# TARGET OS:         all...
//# DESIGN:            Zhao Cheng, Boris Jeremic (jeremic@ucdavis.edu)
//# PROGRAMMER(S):     Zhao Cheng, Boris Jeremic
//#
//#
//# DATE:              July 2004
//# UPDATE HISTORY:
//#
//===============================================================================


#include <stdlib.h>
#include <string.h>

#include <Domain.h>

#include <ErrorHandler.h>
#include <TclModelBuilder.h>

#include <NDMaterial.h>
#include <FiniteDeformationEP3D.h>
#include <FiniteDeformationElastic3D.h>

#include <fdYield.h>
#include <fdYieldVM.h>
#include <fdYieldDP.h>

#include <fdFlow.h>
#include <fdFlowVM.h>
#include <fdFlowDP.h>

#include <fdEvolution_S.h>
#include <fdEvolution_SLS.h>

#include <fdEvolution_T.h>
#include <fdEvolution_TL.h>


// the functions to create the component objects (defined at eof)
fdYield         *EvaluatefdYield(ClientData, Tcl_Interp *, TCL_Char *tclString);
fdFlow          *EvaluatefdFlow(ClientData, Tcl_Interp *, TCL_Char *tclString);
fdEvolution_S   *EvaluatefdEvolution_S(ClientData, Tcl_Interp *, TCL_Char *tclString);
fdEvolution_T   *EvaluatefdEvolution_T(ClientData, Tcl_Interp *, TCL_Char *tclString);


static void cleanup(TCL_Char **argv) {
    Tcl_Free((char *) argv);
}


//**************************************************************************************
//**************************************************************************************
FiniteDeformationEP3D *
TclModelBuilder_addFiniteDeformationEP3D(ClientData clientData, Tcl_Interp *interp,  int argc,
          TCL_Char **argv, TclModelBuilder *theTclBuilder, int eleArgStart)
{
  int tag = 0;
  int tagElasticFD = 0;
  fdYield           *fdY  = 0;
  fdFlow            *fdF  = 0;
  fdEvolution_S  *fdES = 0;
  fdEvolution_T  *fdET = 0;

  NDMaterial	   *matFDElastic = 0;

  int loc = eleArgStart;

  if (Tcl_GetInt(interp, argv[loc++], &tag) != TCL_OK) {
    opserr << "Warning: nDMaterial FDEP3D - invalid tag " << argv[loc] << "\n";
    exit (-1);
  }

  if (Tcl_GetInt(interp, argv[loc++], &tagElasticFD) != TCL_OK) {
    opserr << "Warning: nDMaterial FDEP3D - invalid elastic material tag " << argv[loc] << "\n";
    exit (-1);		
  }

  if (tagElasticFD == tag) {
    opserr << "Error: nDMaterial FDEP3D, elastic matTag is the same with FDEP3D matTag" << argv[loc] << "\n";
    exit (-1);		
  }

  matFDElastic = theTclBuilder->getNDMaterial(tagElasticFD);
  
  if (tagElasticFD == 0) {
 	    opserr << "WARNING: nD FD elastic material does not exist\n";
 	    opserr << "nD FD material: " << tagElasticFD; 
 	    opserr << "\n FDEP3D nDMaterial: " << tag << "\n";
 	    exit (-1);
 	}

  while (loc < argc) {

    if ((strcmp(argv[loc],"-fdYield") == 0) || (strcmp(argv[loc],"-fdY") == 0)) {
      fdY = EvaluatefdYield(clientData, interp, argv[loc+1]);
      if (fdY == 0) {
        opserr << "Warning: nDMaterial FDEP3D - could not create a fdYield from" << argv[loc+1] << "\n";
        exit (-1);
      }
    }

    else if ((strcmp(argv[loc],"-fdFlow") == 0) || (strcmp(argv[loc],"-fdF") == 0)) {
      fdF = EvaluatefdFlow(clientData, interp, argv[loc+1]);
      if (fdF == 0) {
        opserr << "Warning: nDMaterial FDEP3D - could not create a fdFlow from" << argv[loc+1] << "\n";
        exit (-1);
      }
    }

    else if ((strcmp(argv[loc],"-fdEvolution_S") == 0) || (strcmp(argv[loc],"-fdES") == 0)) {
      fdES = EvaluatefdEvolution_S(clientData, interp, argv[loc+1]);
      if (fdES == 0) {
        opserr << "Warning: nDMaterial FDEP3D - could not create a fdES from" << argv[loc+1] << "\n";
        exit (-1);
      }
    }

    else if ((strcmp(argv[loc],"-fdEvolution_T") == 0) || (strcmp(argv[loc],"-fdET") == 0)) {
      fdET = EvaluatefdEvolution_T(clientData, interp, argv[loc+1]);
      if (fdET == 0) {
        opserr << "Warning: nDMaterial FDEP3D - could not create a fdET from" << argv[loc+1] << "\n";
        exit (-1);
      }
    }

    else {
      opserr << "Warning: nDMaterial FDEP3D - don't understand %s\n";
      exit (-1);
    }

    // increment locator by 2 and do next one
    loc += 2;
  }

  FiniteDeformationEP3D *theMaterial = 0;

  if ( (fdY != 0) && (fdF != 0) && (fdES != 0) && (fdET != 0) && (matFDElastic != 0) )
    theMaterial = new FiniteDeformationEP3D(tag, matFDElastic, fdY, fdF, fdES, fdET);
  else if ( (fdY != 0) && (fdF != 0) && (fdES != 0) && (fdET == 0) && (matFDElastic != 0) ) 
    theMaterial = new FiniteDeformationEP3D(tag, matFDElastic, fdY, fdF, fdES);
  else if ( (fdY != 0) && (fdF != 0) && (fdES == 0) && (fdET != 0) && (matFDElastic != 0) ) 
    theMaterial = new FiniteDeformationEP3D(tag, matFDElastic, fdY, fdF, fdET);
  else if ( (fdY != 0) && (fdF != 0) && (fdES == 0) && (fdET == 0) && (matFDElastic != 0) ) 
    theMaterial = new FiniteDeformationEP3D(tag, matFDElastic, fdY, fdF);
  else
    opserr << "Warning: invalid args used to create a FiniteDeformationEP3D material\n";

  return theMaterial;
}


//**************************************************************************************
//**************************************************************************************
// Function - to create a FD Yield Surface
fdYield *EvaluatefdYield(ClientData clientData, Tcl_Interp *interp, TCL_Char *tclString)
{
  int argc;
  TCL_Char **argv;

  // split the list
  if (Tcl_SplitList(interp, tclString, &argc, &argv) != TCL_OK) {
    exit (-1);
  }

  if (argc == 0)
    exit (-1);

  // now parse the list & construct the required object
  fdYield *fdY = 0;

  // 1. von Mises fd Yield Surface
  //
  if ((strcmp(argv[0],"-VM") == 0) || (strcmp(argv[0],"-vM") == 0) || (strcmp(argv[0],"-J2") == 0)) {
    double Y0 = 0.0;

    if (argc == 2) {
      if (Tcl_GetDouble(interp, argv[1], &Y0) != TCL_OK) {
        opserr << "Warning: nDMaterial FDEP3D - invalid Y0 " << argv[1] << "\n";
        exit (-1);
      }
    }   

    fdY = new fdYieldVM(Y0);
  }

  // 2. Druke-Prager fd Yield Surface
  //
  else if ((strcmp(argv[0],"-DP") == 0) || (strcmp(argv[0],"-dp") == 0) ) {
    double FrictionAng_in = 0.0;
    double Cohension_in = 0.0;
    int ConeIndex_in = 0;

    if (argc >= 3) {
      if (Tcl_GetDouble(interp, argv[1], &FrictionAng_in) != TCL_OK) {
        opserr << "Warning: nDMaterial FDEP3D - invalid Friction Angle " << argv[1] << "\n";
        exit (-1);
      }
      if (Tcl_GetDouble(interp, argv[2], &Cohension_in) != TCL_OK) {
        opserr << "Warning: nDMaterial FDEP3D - invalid Conhension " << argv[2] << "\n";
        exit (-1);
      }
    }   

    if (argc == 4) {
      if (Tcl_GetInt(interp, argv[3], &ConeIndex_in) != TCL_OK) {
        opserr << "Warning: nDMaterial FDEP3D - invalid Cone Index " << argv[1] << "\n";
        exit (-1);
      }
    } 

    fdY = new fdYieldDP(FrictionAng_in, Cohension_in, ConeIndex_in);
  }

  else {
    opserr << "Warning: invalid fd yield function: " << argv[0] << "\n";
    exit (-1);
  }

  cleanup(argv);
  return fdY;
}


//**************************************************************************************
//**************************************************************************************
// Function - to create a FD Flow Rule
fdFlow *EvaluatefdFlow(ClientData clientData, Tcl_Interp *interp, TCL_Char *tclString)
{
  int argc;
  TCL_Char **argv;

  // split the list
  if (Tcl_SplitList(interp, tclString, &argc, &argv) != TCL_OK) {
    exit (-1);
  }

  if (argc == 0)
    exit (-1);

  fdFlow *fdF = 0;

  // 1. von Mises fd Yield Surface
  //
  if ((strcmp(argv[0],"-VM") == 0) || (strcmp(argv[0],"-vM") == 0) || (strcmp(argv[0],"-J2") == 0)) {
    double Y0 = 0.0;

    if (argc == 2) {
      if (Tcl_GetDouble(interp, argv[1], &Y0) != TCL_OK) {
        opserr << "Warning: nDMaterial FDEP3D - invalid Y0 " << argv[1] << "\n";
        exit (-1);
      }
    }   

    fdF = new fdFlowVM(Y0);
  }

  // 2. Druke-Prager fd Flow Rule
  //
  else if ((strcmp(argv[0],"-DP") == 0) || (strcmp(argv[0],"-dp") == 0) ) {
    double DilatedAngle_in = 0.0;
    int ConeIndex_in = 0;

    if (argc >= 2) {
      if (Tcl_GetDouble(interp, argv[1], &DilatedAngle_in) != TCL_OK) {
        opserr << "Warning: nDMaterial FDEP3D - invalid Dilated Angle " << argv[1] << "\n";
        exit (-1);
      }
    }   

    if (argc == 3) {
      if (Tcl_GetInt(interp, argv[2], &ConeIndex_in) != TCL_OK) {
        opserr << "Warning: nDMaterial FDEP3D - invalid Cone Index " << argv[2] << "\n";
        exit (-1);
      }
    } 

    fdF = new fdFlowDP(DilatedAngle_in, ConeIndex_in);
  }

  else {
    opserr << "Warning: invalid fd flow rule: " << argv[0] << "\n";
    exit (-1);
  }

  cleanup(argv);
  return fdF;
}


//**************************************************************************************
//**************************************************************************************
// Function - to create an fd EvolutionLaw_S object
fdEvolution_S *EvaluatefdEvolution_S(ClientData clientData, Tcl_Interp *interp, TCL_Char *tclString)
{
  int argc;
  TCL_Char **argv;

  // split the list
  if (Tcl_SplitList(interp, tclString, &argc, &argv) != TCL_OK) {
    exit (-1);
  }

  fdEvolution_S *fdES = 0;

  //1. Linear and Saturation isotropic (scalar) evolution law:
  //
  if ((strcmp(argv[0],"-LS") == 0) || (strcmp(argv[0],"-LinearSaturated") == 0)) {

    double H_linear    = 0.0;
    double q_saturated = 0.0;
    double beta        = 0.0;

    if (argc >= 2) {
      if (Tcl_GetDouble(interp, argv[1], &H_linear) != TCL_OK) {
        opserr << "Warning: nDMaterial FDEP3D - invalid H_linear " << argv[1] << "\n";
        cleanup(argv);
        exit (-1);
      }
    }
    
    if (argc >= 4) {
      if (Tcl_GetDouble(interp, argv[2], &q_saturated) != TCL_OK) {
        opserr << "Warning: nDMaterial FDEP3D - invalid q_saturated " << argv[2] << "\n";
        cleanup(argv);
        exit (-1);
      }
      if (Tcl_GetDouble(interp, argv[3], &beta) != TCL_OK) {
        opserr << "Warning: nDMaterial FDEP3D - invalid beta " << argv[3] << "\n";
        cleanup(argv);
        exit (-1);
      }
    }
  
    fdES = new fdEvolution_SLS(H_linear, q_saturated, beta);
  }
 
  cleanup(argv);
  return fdES;
}

//**************************************************************************************
//**************************************************************************************
// Function - to create an fd EvolutionLaw_T object
fdEvolution_T *EvaluatefdEvolution_T(ClientData clientData, Tcl_Interp *interp, TCL_Char *tclString)
{
  int argc;
  TCL_Char **argv;

  // split the list
  if (Tcl_SplitList(interp, tclString, &argc, &argv) != TCL_OK) {
    exit (-1);
  }

  fdEvolution_T *fdET = 0;

  //1. Linear kinematic (tensor) evolution law:
  //
  if ((strcmp(argv[0],"-Linear") == 0) ) {

    double H_linear    = 0.0;

    if (argc >= 2) {
      if (Tcl_GetDouble(interp, argv[1], &H_linear) != TCL_OK) {
        opserr << "Warning: nDMaterial FDEP3D - invalid H_linear " << argv[1] << "\n";
        cleanup(argv);
        exit (-1);
      }
    }
      
    fdET = new fdEvolution_TL(H_linear);
  }
 
  cleanup(argv);
  return fdET;
}


