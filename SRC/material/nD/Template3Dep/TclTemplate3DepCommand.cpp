//===============================================================================
//# COPYRIGHT (C):     :-))
//# PROJECT:           Object Oriented Finite Element Program
//# PURPOSE:           Tcl based Parser
//# CLASS:
//#
//# VERSION:
//# LANGUAGE:          C++
//# TARGET OS:         DOS || UNIX || . . .
//# PROGRAMMER(S):     Zhao Cheng, Boris Jeremic
//#
//#
//# DATE:              19Feb2003
//# UPDATE HISTORY:    18May2004 ZHao fixed Armstrong-Frederick command
//#		       May2004 Zhao Cheng splitted elastic part 
//#
//#
//#
//===============================================================================

#include <stdlib.h>
#include <string.h>

#include <Domain.h>

#include <ErrorHandler.h>
#include <TclModelBuilder.h>

#include <Template3Dep.h>
#include <straint.h>
#include <stresst.h>

#include <YS.h>
#include <DP_YS.h>
#include <VM_YS.h>
#include <CAM_YS.h>
#include <RMC01_YS.h>
#include <MD_YS.h>

#include <PS.h>
#include <DP_PS.h>
#include <VM_PS.h>
#include <CAM_PS.h>
#include <RMC01_PS.h>
#include <MD_PS.h>

#include <EPState.h>

#include <EL_S.h>
#include <EL_LEeq.h>
#include <EL_NLEp.h>

#include <EL_T.h>
#include <EL_LEij.h>
#include <EL_NLEij.h>
#include <EL_NLEijMD.h>





// the functions to create the component objects (defined at eof)
YieldSurface     *EvaluateYieldSurfaceArgs(ClientData, Tcl_Interp *, TCL_Char *tclString);
PotentialSurface *EvaluatePotentialSurfaceArgs(ClientData, Tcl_Interp *, TCL_Char *tclString);
EPState          *EvaluateEPStateArgs(ClientData, Tcl_Interp *, TCL_Char *tclString);
EvolutionLaw_S   *EvaluateEvolutionLawSArgs(ClientData, Tcl_Interp *, TCL_Char *tclString);
EvolutionLaw_T   *EvaluateEvolutionLawTArgs(ClientData, Tcl_Interp *, TCL_Char *tclString);


// little function to free memory after invoke Tcl_SplitList
//   note Tcl_Split list stores the array of pointers and the strings in
//   one array, which is why Tcl_Free needs only be called on the array.

static void cleanup(TCL_Char **argv) {
    Tcl_Free((char *) argv);
}

Template3Dep *
TclModelBuilder_addTemplate3Dep(ClientData clientData, Tcl_Interp *interp,  int argc,
          TCL_Char **argv, TclModelBuilder *theTclBuilder, int eleArgStart)
{
  // create some empty pointers which we fill in as parse the command line
  int tag =0;
  YieldSurface     *YS =0;
  PotentialSurface *PS =0;
  EPState          *EPS =0;
  EvolutionLaw_S   *ELS1 =0;
  EvolutionLaw_S   *ELS2 =0;
  EvolutionLaw_S   *ELS3 =0;
  EvolutionLaw_S   *ELS4 =0;
  EvolutionLaw_T   *ELT1 =0;
  EvolutionLaw_T   *ELT2 =0;
  EvolutionLaw_T   *ELT3 =0;
  EvolutionLaw_T   *ELT4 =0;
  int ElMattag;

  int loc = eleArgStart;
  if (Tcl_GetInt(interp, argv[loc++], &tag) != TCL_OK) {
    opserr << "nDMaterial Templated3Dep - invalid tag " << argv[loc] << endln;
    return 0;
  }

  if (Tcl_GetInt(interp, argv[loc++], &ElMattag) != TCL_OK) {
    opserr << "nDMaterial Templated3Dep - invalid elastic material tag " << argv[loc] << endln;
    return 0;		
  }

  if (ElMattag == tag) {
    opserr << "nDMaterial Templated3Dep - error: the elastic matTag is the same with Templated3Dep matTag" << argv[loc] << endln;
    return 0;		
  }

  NDMaterial *theElMat = theTclBuilder->getNDMaterial(ElMattag);
  
  if (theElMat == 0) {
 	    opserr << "WARNING nD elastic material does not exist\n";
 	    opserr << "nD material: " << ElMattag; 
 	    opserr << "\nTamplate3Dep nDMaterial: " << tag << endln;
 	    return 0;
 	}

  // parse rest of command, switching on the flags -YS, -PS, -EPS, -ELS, -ELT, ...
  while (loc < argc) {
    if ((strcmp(argv[loc],"-YS") == 0) || (strcmp(argv[loc],"-ys") == 0)) {
      YS = EvaluateYieldSurfaceArgs(clientData, interp, argv[loc+1]);
      if (YS == 0) {
  opserr << "nDMaterial Templated3Dep - could not create a YS from" << argv[loc+1] << endln;
  return 0;
      }
    }

    else if ((strcmp(argv[loc],"-PS") == 0) || (strcmp(argv[loc],"-ps") == 0)) {
      PS = EvaluatePotentialSurfaceArgs(clientData, interp, argv[loc+1]);
      if (PS == 0) {
  opserr << "nDMaterial Templated3Dep - could not create a PS from" << argv[loc+1] << endln;
  return 0;
      }
    }

   else if ((strcmp(argv[loc],"-EPS") == 0) || (strcmp(argv[loc],"-eps") == 0)) {
      EPS = EvaluateEPStateArgs(clientData, interp, argv[loc+1]);
      if (EPS == 0) {
  opserr << "nDMaterial Templated3Dep - could not create an EPS from" << argv[loc+1] << endln;
  return 0;
      }
    }

    else if ((strcmp(argv[loc],"-ELS1") == 0) || (strcmp(argv[loc],"-els1") == 0)) {
      ELS1 = EvaluateEvolutionLawSArgs(clientData, interp, argv[loc+1]);
      if (ELS1 == 0) {
  opserr << "nDMaterial Templated3Dep - could not create an ELS1 from" << argv[loc+1] << endln;
  return 0;
      }
    }

    else if ((strcmp(argv[loc],"-ELS2") == 0) || (strcmp(argv[loc],"-els2") == 0)) {
      ELS2 = EvaluateEvolutionLawSArgs(clientData, interp, argv[loc+1]);
      if (ELS2 == 0) {
  opserr << "nDMaterial Templated3Dep - could not create an ELS2 from" << argv[loc+1] << endln;
  return 0;
      }
    }

    else if ((strcmp(argv[loc],"-ELS3") == 0) || (strcmp(argv[loc],"-els3") == 0)) {
      ELS3 = EvaluateEvolutionLawSArgs(clientData, interp, argv[loc+1]);
      if (ELS3 == 0) {
  opserr << "nDMaterial Templated3Dep - could not create an ELS3 from" << argv[loc+1] << endln;
  return 0;
      }
    }

    else if ((strcmp(argv[loc],"-ELS4") == 0) || (strcmp(argv[loc],"-els4") == 0)) {
      ELS4 = EvaluateEvolutionLawSArgs(clientData, interp, argv[loc+1]);
      if (ELS4 == 0) {
  opserr << "nDMaterial Templated3Dep - could not create an ELS4 from" << argv[loc+1] << endln;
  return 0;
      }
    }

    else if ((strcmp(argv[loc],"-ELT1") == 0) || (strcmp(argv[loc],"-elt1") == 0)) {
      ELT1 = EvaluateEvolutionLawTArgs(clientData, interp, argv[loc+1]);
      if (ELT1 == 0) {
  opserr << "nDMaterial Templated3Dep - could not create an ELT1 from" << argv[loc+1] << endln;
  return 0;
      }
    }

    else if ((strcmp(argv[loc],"-ELT2") == 0) || (strcmp(argv[loc],"-elt2") == 0)) {
      ELT2 = EvaluateEvolutionLawTArgs(clientData, interp, argv[loc+1]);
      if (ELT2 == 0) {
  opserr << "nDMaterial Templated3Dep - could not create an ELT2 from" << argv[loc+1] << endln;
  return 0;
      }
    }

    else if ((strcmp(argv[loc],"-ELT3") == 0) || (strcmp(argv[loc],"-elt3") == 0)) {
      ELT3 = EvaluateEvolutionLawTArgs(clientData, interp, argv[loc+1]);
      if (ELT3 == 0) {
  opserr << "nDMaterial Templated3Dep - could not create an ELT3 from" << argv[loc+1] << endln;
  return 0;
      }
    }

    else if ((strcmp(argv[loc],"-ELT4") == 0) || (strcmp(argv[loc],"-elt4") == 0)) {
      ELT4 = EvaluateEvolutionLawTArgs(clientData, interp, argv[loc+1]);
      if (ELT4 == 0) {
  opserr << "nDMaterial Templated3Dep - could not create an ELT4 from" << argv[loc+1] << endln;
  return 0;
      }
    }

    else {
      opserr << "nDMaterial Templated3Dep - don't understand %s\n";
  return 0;
    }

    // increment locator by 2 and do next one
    loc += 2;
  }

  // now depending on the objects types that are not NULL we use the appropriate
  // constructor to construct or Template3Dep material object
  Template3Dep *theMaterial = 0;
  if ( (YS != 0) && (PS != 0) && (EPS != 0) &&
       (ELS1 != 0) && (ELS2 != 0) && (ELS3 != 0) && (ELS4 != 0) &&
       (ELT1 != 0) && (ELT2 != 0) && (ELT3 != 0) && (ELT4 != 0) )
    theMaterial = new Template3Dep(tag, *theElMat, YS, PS, EPS, ELS1, ELS2, ELS3, ELS4,
           ELT1, ELT2, ELT3, ELT4);

  // constructor 0
  else if ( (YS != 0) && (PS != 0) && (EPS != 0) &&
       (ELS1 == 0) && (ELS2 == 0) && (ELS3 == 0) && (ELS4 == 0) &&
       (ELT1 == 0) && (ELT2 == 0) && (ELT3 == 0) && (ELT4 == 0) )
    theMaterial = new Template3Dep(tag, *theElMat, YS, PS, EPS);

  // constructor 1
  else if ( (YS != 0) && (PS != 0) && (EPS != 0) &&
       (ELS1 != 0) && (ELS2 == 0) && (ELS3 == 0) && (ELS4 == 0) &&
       (ELT1 == 0) && (ELT2 == 0) && (ELT3 == 0) && (ELT4 == 0) )
    theMaterial = new Template3Dep(tag, *theElMat, YS, PS, EPS, ELS1);

  // constructor 2
  else if ( (YS != 0) && (PS != 0) && (EPS != 0) &&
       (ELS1 == 0) && (ELS2 == 0) && (ELS3 == 0) && (ELS4 == 0) &&
       (ELT1 != 0) && (ELT2 == 0) && (ELT3 == 0) && (ELT4 == 0) )
    theMaterial = new Template3Dep(tag, *theElMat, YS, PS, EPS, ELT1);

  // constructor 3
  else if ( (YS != 0) && (PS != 0) && (EPS != 0) &&
       (ELS1 != 0) && (ELS2 == 0) && (ELS3 == 0) && (ELS4 == 0) &&
       (ELT1 != 0) && (ELT2 == 0) && (ELT3 == 0) && (ELT4 == 0) )
    theMaterial = new Template3Dep(tag, *theElMat, YS, PS, EPS, ELS1, ELT1);

  // constructor 4
  else if ( (YS != 0) && (PS != 0) && (EPS != 0) &&
       (ELS1 != 0) && (ELS2 != 0) && (ELS3 == 0) && (ELS4 == 0) &&
       (ELT1 != 0) && (ELT2 == 0) && (ELT3 == 0) && (ELT4 == 0) )
    theMaterial = new Template3Dep(tag, *theElMat, YS, PS, EPS, ELS1, ELS2, ELT1);

  // constructor 5
  else if ( (YS != 0) && (PS != 0) && (EPS != 0) &&
       (ELS1 != 0) && (ELS2 != 0) && (ELS3 == 0) && (ELS4 == 0) &&
       (ELT1 != 0) && (ELT2 != 0) && (ELT3 == 0) && (ELT4 == 0) )
    theMaterial = new Template3Dep(tag, *theElMat, YS, PS, EPS, ELS1, ELS2, ELT1, ELT2);

  else
    opserr << "invalid number of args used to create a Template3Dep material\n";

    return theMaterial;
}


// Function - to create a YieldSurface
YieldSurface *EvaluateYieldSurfaceArgs(ClientData clientData, Tcl_Interp *interp, TCL_Char *tclString)
{
  int argc;
  TCL_Char **argv;

  // split the list
  if (Tcl_SplitList(interp, tclString, &argc, &argv) != TCL_OK) {
    return 0;
  }

  if (argc == 0)
    return 0;

  // now parse the list & construct the required object
  YieldSurface *YS = 0;

  // 1. Drucker-Prager Yield Surface
  if ((strcmp(argv[0],"-DP") == 0) || (strcmp(argv[0],"-dp") == 0)) {
    YS = new DPYieldSurface();
  }

  // 2. von Mises Yield Surface
  else if ( (strcmp(argv[0],"-VM") == 0) || (strcmp(argv[0],"-vM") == 0) || (strcmp(argv[0],"-vm") == 0)) {
    YS = new VMYieldSurface();
  }

  // 3. Cam Clay yield surface
  else if ((strcmp(argv[0],"-CC") == 0) || (strcmp(argv[0],"-cc") == 0)) {
    double mp = 1.2;
    if (argc > 1) {
      if (Tcl_GetDouble(interp, argv[1], &mp) != TCL_OK) {
       opserr << "invalid M: argv[1] for -PS CamClay M\n";
       return 0;
      }
    }
    YS = new CAMYieldSurface(mp);
  }

   // 4. Rounded Mohr-Conlomb 01 Yield Surface
  else if ( (strcmp(argv[0],"-RMC01") == 0) || (strcmp(argv[0],"-rmc01") == 0) || (strcmp(argv[0],"-rmc01") == 0)) {
    YS = new RMC01YieldSurface();
  }

  // 5. Manzari-Dafalias (conic-shaped) Yield Surface
  else if ( (strcmp(argv[0],"-MD") == 0) || (strcmp(argv[0],"-md") == 0) ) {
    YS = new MDYieldSurface();
  }

  else {
    opserr << "invalid yield function: " << argv[0] << endln;
    return 0;
  }

  cleanup(argv);
  return YS;
}

// Function - to create a PotentialSurface object
PotentialSurface *EvaluatePotentialSurfaceArgs(ClientData clientData, Tcl_Interp *interp, TCL_Char *tclString)
{
  int argc;
  TCL_Char **argv;

  // split the list
  if (Tcl_SplitList(interp, tclString, &argc, &argv) != TCL_OK) {
    return 0;
  }

  if (argc == 0)
    return 0;

  // now parse the list & construct the required object
  PotentialSurface *PS = 0;

  //1. Drucker-Prager Potential Surface
  if ((strcmp(argv[0],"-DP") == 0) || (strcmp(argv[0],"-dp") == 0)) {
    double alpha = 0.0;
    if (argc > 1)
      if (Tcl_GetDouble(interp, argv[1], &alpha) != TCL_OK) {
       opserr << "nDMaterial Templated3Dep - invalid alpha " << argv[1] << endln;
       cleanup(argv);
       return 0;
      }
    PS = new DPPotentialSurface(alpha);
  }

  //2. von-Mises Potential Surface
  else if ((strcmp(argv[0],"-VM") == 0) || (strcmp(argv[0],"-vm") == 0)) {
    PS = new VMPotentialSurface();
  }

  // 3. CamClay potential surface
  else if ((strcmp(argv[0],"-CC") == 0) || (strcmp(argv[0],"-cc") == 0)) {
    double mp = 1.2;
    if (argc == 2) {
      if (Tcl_GetDouble(interp, argv[1], &mp) != TCL_OK) {
        opserr << "nDMaterial Templated3Dep - invalid M " << argv[1] << endln;
    return 0;
      }
    }
    PS = new CAMPotentialSurface(mp);
  }

 //4. RMC01 Potential Surface
  else if ((strcmp(argv[0],"-RMC01") == 0) || (strcmp(argv[0],"-rmc01") == 0))
    {
      PS = new RMC01PotentialSurface();
    }

 //5. Manzari-Dafalias Potential Surface
  else if ((strcmp(argv[0],"-MD") == 0) || (strcmp(argv[0],"-md") == 0))
    {
      PS = new MDPotentialSurface();
    }

  else {
    opserr << "invalid potential function: " << argv[0] << endln;
    return 0;
  }

  cleanup(argv);
  return PS;
}

//ZCMay04 // Function - to read in a stress tensor
//ZCMay04 
//ZCMay04 int EvaluateStressTensor(ClientData clientData, Tcl_Interp *interp, TCL_Char *tclString,
//ZCMay04        stresstensor &stress)
//ZCMay04 {
//ZCMay04   int argc;
//ZCMay04   TCL_Char **argv;
//ZCMay04 
//ZCMay04   // split the list
//ZCMay04   if (Tcl_SplitList(interp, tclString, &argc, &argv) != TCL_OK) {
//ZCMay04     return 0;
//ZCMay04   }
//ZCMay04 
//ZCMay04   double *values = new double[9];
//ZCMay04   if (values == 0) {
//ZCMay04     opserr << "nDMaterial Template3DEp -ESP stresstenor - out of memory for size: " << argc << endln;
//ZCMay04     cleanup(argv);
//ZCMay04     return -1;
//ZCMay04   }
//ZCMay04 
//ZCMay04   int i;
//ZCMay04   for (i=0; i<arc; i++) {
//ZCMay04     if (Tcl_GetDouble(interp, argv[i], &values[i]) != TCL_OK) {
//ZCMay04       opserr << "nDMaterial Templated3Dep -ESP stresstensor - invalid value: " << argv[i] << endln;
//ZCMay04       cleanup(argv);
//ZCMay04       //delete [] values;
//ZCMay04       return -1;
//ZCMay04     }
//ZCMay04 
//ZCMay04 
//ZCMay04   }
//ZCMay04 
//ZCMay04   for (i=1; i<9; i++) {
//ZCMay04      values[i] = 0;
//ZCMay04      if ((i == 4)|| (i==8))
//ZCMay04         values[i] = -1.0 * values[0];
//ZCMay04   }
//ZCMay04   values[0] = -1.0 * values[0];
//ZCMay04   stresstensor newStress(values);
//ZCMay04 
//ZCMay04   //newStress.printshort("tcl:");
//ZCMay04 
//ZCMay04   stress = newStress;
//ZCMay04   cleanup(argv);
//ZCMay04   delete [] values;
//ZCMay04   return 0;
//ZCMay04 }

//ZCMay04 // Function - to read in a strain tensor
//ZCMay04 int EvaluateStrainTensor(ClientData clientData, Tcl_Interp *interp, TCL_Char *tclString,
//ZCMay04        straintensor &strain)
//ZCMay04 {
//ZCMay04   int argc;
//ZCMay04   TCL_Char **argv;
//ZCMay04 
//ZCMay04   // split the list
//ZCMay04   if (Tcl_SplitList(interp, tclString, &argc, &argv) != TCL_OK) {
//ZCMay04     return 0;
//ZCMay04   }
//ZCMay04 
//ZCMay04   double *values = new double[argc];
//ZCMay04   if (values == 0) {
//ZCMay04     opserr << "nDMaterial Template3DEp -ESP straintensor - out of memory for size: " << argc << endln;
//ZCMay04     cleanup(argv);
//ZCMay04     return -1;
//ZCMay04   }
//ZCMay04 
//ZCMay04   for (int i=0; i<argc; i++) {
//ZCMay04     if (Tcl_GetDouble(interp, argv[i], &values[i]) != TCL_OK) {
//ZCMay04       opserr << "nDMaterial Templated3Dep -ESP straintensor - invalid value: " << argv[i] << endln;
//ZCMay04       cleanup(argv);
//ZCMay04       //delete [] values;
//ZCMay04       return -1;
//ZCMay04     }
//ZCMay04   }
//ZCMay04 
//ZCMay04   straintensor newStrain(values);
//ZCMay04   strain = newStrain;
//ZCMay04   cleanup(argv);
//ZCMay04   delete [] values;
//ZCMay04   return 0;
//ZCMay04 }

// Function - to create an EPState object
EPState *EvaluateEPStateArgs(ClientData clientData, Tcl_Interp *interp, TCL_Char *tclString)
{
  int argc;
  TCL_Char **argv;

  // split the list
  if (Tcl_SplitList(interp, tclString, &argc, &argv) != TCL_OK) {
    return 0;
  }
  EPState *EPS = 0;

//ZCMay04  double Eod = 0.0;
//ZCMay04  double Ed = 0.0;
//ZCMay04  double nu = 0.0;
//ZCMay04  double rho = 0.0; 
//ZCMay04  
//ZCMay04  if (argc < 5) {
//ZCMay04 
//ZCMay04    cleanup(argv);
//ZCMay04    return 0;
//ZCMay04  }
//ZCMay04 
//ZCMay04  if (Tcl_GetDouble(interp, argv[0], &Eod) != TCL_OK) {
//ZCMay04    opserr << "nDMaterial Templated3Dep -EPS - invalid Eod " << argv[1] << endln;
//ZCMay04    cleanup(argv);
//ZCMay04    return 0;
//ZCMay04  }
//ZCMay04  if (Tcl_GetDouble(interp, argv[1], &Ed) != TCL_OK) {
//ZCMay04    opserr << "nDMaterial Templated3Dep -EPS - invalid Ed " << argv[1] << endln;
//ZCMay04    cleanup(argv);
//ZCMay04    return 0;
//ZCMay04  }
//ZCMay04  if (Tcl_GetDouble(interp, argv[2], &nu) != TCL_OK) {
//ZCMay04    opserr << "nDMaterial Templated3Dep -EPS - invalid nu " << argv[1] << endln;
//ZCMay04    cleanup(argv);
//ZCMay04    return 0;
//ZCMay04  }
//ZCMay04  if (Tcl_GetDouble(interp, argv[3], &rho) != TCL_OK) {
//ZCMay04    opserr << "nDMaterial Templated3Dep -EPS - invalid rho " << argv[1] << endln;
//ZCMay04  }
//ZCMay04 
//ZCMay04  int loc = 4;
  int loc = 0;
  stresstensor stressp(0.0);
  straintensor strainp(0.0);
  straintensor Estrainp(0.0);
  straintensor Pstrainp(0.0);

  int NoS = 0;
  int NoD = 0;
  double *scalars = 0;
  stresstensor *tensors = 0;
  double ep = 0.85;
  double psip = 0.05;
//ZCMay04   int     Elasticflagp = 0;
//ZCMay04   double  Evp   = 0.0;
//ZCMay04   double  nuhvp = 0.0;
//ZCMay04   double  Ghvp = 0.0;  
//ZCMay04   double  eod = 0.80;
//ZCMay04   double  ecd = 0.75;
//ZCMay04   double  lamd=0.025;
//ZCMay04   double  p_refd= 100.0;
//ZCMay04   double  ad=0.5;

  // switch on remaining args;
  //double *values_stp = new double[9];
  double values_stp[9]= {-0.00010, 0.0, 0.0, 0.0, -0.00010, 0.0, 0.0, 0.0, -0.00010};

  while (loc < argc) {

    if ((strcmp(argv[loc],"-stressP") == 0) || (strcmp(argv[loc],"-stressp") == 0)) {
      int loc_i = loc;
      for (int i=0; i<9; i++) {
        if (Tcl_GetDouble(interp, argv[loc+1], &values_stp[i]) != TCL_OK) {
          opserr << "nDMaterial Templated3Dep -EPS - invalid stress tensor: need nine components " << argv[loc_i] << endln;
          return 0;
        }
        loc++;
      }
      stressp.val(1,1) = values_stp[0]; stressp.val(1,2) = values_stp[1]; stressp.val(1,3) = values_stp[2]; 
      stressp.val(2,1) = values_stp[3]; stressp.val(2,2) = values_stp[4]; stressp.val(2,3) = values_stp[5];
      stressp.val(3,1) = values_stp[6]; stressp.val(3,2) = values_stp[7]; stressp.val(3,3) = values_stp[8]; 
      loc++;
    }

    else if ((strcmp(argv[loc],"-strainP") == 0) || (strcmp(argv[loc],"-strainp") == 0)) {
      int loc_i = loc;
      for (int i=0; i<9; i++) {
        if (Tcl_GetDouble(interp, argv[loc+1], &values_stp[i]) != TCL_OK) {
          opserr << "nDMaterial Templated3Dep -EPS - invalid strain tensor: need nine components " << argv[loc_i] << endln;
          return 0;
        }
        loc++;
      }
      strainp.val(1,1) = values_stp[0]; strainp.val(1,2) = values_stp[1]; strainp.val(1,3) = values_stp[2]; 
      strainp.val(2,1) = values_stp[3]; strainp.val(2,2) = values_stp[4]; strainp.val(2,3) = values_stp[5];
      strainp.val(3,1) = values_stp[6]; strainp.val(3,2) = values_stp[7]; strainp.val(3,3) = values_stp[8]; 
      loc++;
    }

    else if ((strcmp(argv[loc],"-EstrainP") == 0) || (strcmp(argv[loc],"-Estrainp") == 0)) {
      int loc_i = loc;
      for (int i=0; i<9; i++) {
        if (Tcl_GetDouble(interp, argv[loc+1], &values_stp[i]) != TCL_OK) {
          opserr << "nDMaterial Templated3Dep -EPS - invalid Estrain tensor: need nine components " << argv[loc_i] << endln;
          return 0;
        }
        loc++;
      }
      Estrainp.val(1,1) = values_stp[0]; Estrainp.val(1,2) = values_stp[1]; Estrainp.val(1,3) = values_stp[2]; 
      Estrainp.val(2,1) = values_stp[3]; Estrainp.val(2,2) = values_stp[4]; Estrainp.val(2,3) = values_stp[5];
      Estrainp.val(3,1) = values_stp[6]; Estrainp.val(3,2) = values_stp[7]; Estrainp.val(3,3) = values_stp[8]; 
      loc++;
    }

    else if ((strcmp(argv[loc],"-PstrainP") == 0) || (strcmp(argv[loc],"-Pstrainp") == 0)) {
      int loc_i = loc;
      for (int i=0; i<9; i++) {
        if (Tcl_GetDouble(interp, argv[loc+1], &values_stp[i]) != TCL_OK) {
          opserr << "nDMaterial Templated3Dep -EPS - invalid Pstrain tensor: need nine components " << argv[loc_i] << endln;
          return 0;
        }
        loc++;
      }
      Pstrainp.val(1,1) = values_stp[0]; Pstrainp.val(1,2) = values_stp[1]; Pstrainp.val(1,3) = values_stp[2]; 
      Pstrainp.val(2,1) = values_stp[3]; Pstrainp.val(2,2) = values_stp[4]; Pstrainp.val(2,3) = values_stp[5];
      Pstrainp.val(3,1) = values_stp[6]; Pstrainp.val(3,2) = values_stp[7]; Pstrainp.val(3,3) = values_stp[8]; 
      loc++;
    }


//ZCMay04     else if ((strcmp(argv[loc],"-strainP") == 0) || (strcmp(argv[loc],"-strainp") == 0)) {
//ZCMay04       if (EvaluateStrainTensor(clientData, interp, argv[loc+1], strainp) < 0) {
//ZCMay04     opserr << "nDMaterial Templated3Dep -EPS - invalid strainp " << argv[loc+1] << endln;
//ZCMay04     cleanup(argv);
//ZCMay04     return 0;
//ZCMay04       }
//ZCMay04     }
//ZCMay04     else if ((strcmp(argv[loc],"-EstrainP") == 0) || (strcmp(argv[loc],"-Estrainp") == 0)) {
//ZCMay04       if (EvaluateStrainTensor(clientData, interp, argv[loc+1], Estrainp) < 0)
//ZCMay04     opserr << "nDMaterial Templated3Dep -EPS - invalid Estrainp " << argv[loc+1] << endln;
//ZCMay04     cleanup(argv);
//ZCMay04     return 0;
//ZCMay04       loc+=2;
//ZCMay04     }
//ZCMay04     else if ((strcmp(argv[loc],"-PstrainP") == 0) || (strcmp(argv[loc],"-Estrainp") == 0)) {
//ZCMay04       if (EvaluateStrainTensor(clientData, interp, argv[loc+1], Pstrainp) < 0) {
//ZCMay04     opserr << "nDMaterial Templated3Dep -EPS - invalid Pstrainp " << argv[loc+1] << endln;
//ZCMay04     cleanup(argv);
//ZCMay04     return 0;
//ZCMay04       }
//ZCMay04       loc+=2;
//ZCMay04     }

    else if ((strcmp(argv[loc],"-NOS") == 0) || (strcmp(argv[loc],"-nos") == 0)) {
      if (Tcl_GetInt(interp, argv[loc+1], &NoS) != TCL_OK) {
  	opserr << "nDMaterial Templated3Dep -EPS - invalid NOS " << argv[loc+1] << endln;
  	return 0;  
      }
      loc+= 2;

      if (NoS > 0) {
  	scalars = new double[NoS];
  	for (int i=0; i<NoS; i++) {
          if (Tcl_GetDouble(interp, argv[loc++], &scalars[i]) != TCL_OK) {
      	    opserr << "nDMaterial Templated3Dep -EPS - invalid scalar: " << argv[loc-1] << endln;
      	    return 0;      
          }
  	} 
      } 
      else 
  	NoS = 0;
    }
    
    else if ((strcmp(argv[loc],"-NOD") == 0) || (strcmp(argv[loc],"-nod") == 0)) {
      if (Tcl_GetInt(interp, argv[loc+1], &NoD) != TCL_OK) {
  	opserr << "nDMaterial Templated3Dep -EPS - invalid NOD: " << argv[loc+1] << endln;
  	return 0;
      }
      loc += 2;
      
      if (NoD > 0) {
	 tensors = new stresstensor[NoD];
  	 if (tensors == 0) {
    	   opserr << "nDMaterial Templated3Dep -EPS - invalid NOD " << argv[loc+1] << endln;
    	   return 0;
  	}
      }
      else 
  	NoD = 0;
    }

//ZCMay04     else if  (strcmp(argv[loc],"-ElasticFlag") == 0 || strcmp(argv[loc],"-elasticflag") == 0) {
//ZCMay04       if (Tcl_GetInt(interp, argv[loc+1], &Elasticflagp) != TCL_OK) {
//ZCMay04       opserr << "nDMaterial Templated3Dep -EPS - invalid e " << argv[loc+1] << endln;
//ZCMay04       cleanup(argv);
//ZCMay04       return 0;
//ZCMay04       }
//ZCMay04       loc+=2;
//ZCMay04     }
//ZCMay04     else if  (strcmp(argv[loc],"-Ev") == 0) {
//ZCMay04       if (Tcl_GetDouble(interp, argv[loc+1], &Evp) != TCL_OK) {
//ZCMay04       opserr << "nDMaterial Templated3Dep -EPS - invalid Ev " << argv[loc+1] << endln;
//ZCMay04       cleanup(argv);
//ZCMay04       return 0;
//ZCMay04       }
//ZCMay04       loc+=2;
//ZCMay04     }
//ZCMay04     else if  (strcmp(argv[loc],"-nuhv") == 0) {
//ZCMay04       if (Tcl_GetDouble(interp, argv[loc+1], &nuhvp) != TCL_OK) {
//ZCMay04       opserr << "nDMaterial Templated3Dep -EPS - invalid nuhv " << argv[loc+1] << endln;
//ZCMay04       cleanup(argv);
//ZCMay04       return 0;
//ZCMay04       }
//ZCMay04       loc+=2;
//ZCMay04     }
//ZCMay04     else if  (strcmp(argv[loc],"-Ghv") == 0) {
//ZCMay04       if (Tcl_GetDouble(interp, argv[loc+1], &Ghvp) != TCL_OK) {
//ZCMay04       opserr << "nDMaterial Templated3Dep -EPS - invalid Ghv " << argv[loc+1] << endln;
//ZCMay04       cleanup(argv);
//ZCMay04       return 0;
//ZCMay04       }
//ZCMay04       loc+=2;
//ZCMay04     }    
//ZCMay04     else if  (strcmp(argv[loc],"-eo") == 0) {
//ZCMay04       if (Tcl_GetDouble(interp, argv[loc+1], &eod) != TCL_OK) {
//ZCMay04       opserr << "nDMaterial Templated3Dep -EPS - invalid e " << argv[loc+1] << endln;
//ZCMay04       cleanup(argv);
//ZCMay04       return 0;
//ZCMay04       }
//ZCMay04       loc+=2;
//ZCMay04     }
//ZCMay04     else if  (strcmp(argv[loc],"-ec") == 0) {
//ZCMay04       if (Tcl_GetDouble(interp, argv[loc+1], &ecd) != TCL_OK) {
//ZCMay04       opserr << "nDMaterial Templated3Dep -EPS - invalid ec " << argv[loc+1] << endln;
//ZCMay04       cleanup(argv);
//ZCMay04       return 0;
//ZCMay04       }
//ZCMay04       loc+=2;
//ZCMay04     }
//ZCMay04     else if  (strcmp(argv[loc],"-lambda") == 0) {
//ZCMay04       if (Tcl_GetDouble(interp, argv[loc+1], &lamd) != TCL_OK) {
//ZCMay04       opserr << "nDMaterial Templated3Dep -EPS - invalid lambda " << argv[loc+1] << endln;
//ZCMay04       cleanup(argv);
//ZCMay04       return 0;
//ZCMay04       }
//ZCMay04       loc+=2;
//ZCMay04     }
//ZCMay04     else if  (strcmp(argv[loc],"-p_ref") == 0) {
//ZCMay04       if (Tcl_GetDouble(interp, argv[loc+1], &p_refd) != TCL_OK) {
//ZCMay04       opserr << "nDMaterial Templated3Dep -EPS - invalid p_ref " << argv[loc+1] << endln;
//ZCMay04       cleanup(argv);
//ZCMay04       return 0;
//ZCMay04       }
//ZCMay04       loc+=2;
//ZCMay04     }
//ZCMay04     else if  (strcmp(argv[loc],"-a") == 0) {
//ZCMay04       if (Tcl_GetDouble(interp, argv[loc+1], &ad) != TCL_OK) {
//ZCMay04       opserr << "nDMaterial Templated3Dep -EPS - invalid a " << argv[loc+1] << endln;
//ZCMay04       cleanup(argv);
//ZCMay04       return 0;
//ZCMay04       }
//ZCMay04       loc+=2;
//ZCMay04     }

   else if  (strcmp(argv[loc],"-e") == 0) {
      if (Tcl_GetDouble(interp, argv[loc+1], &ep) != TCL_OK) {
        opserr << "nDMaterial Templated3Dep -EPS - invalid e " << argv[loc+1] << endln;
        cleanup(argv);
        return 0;
      }
      loc+=2;
    }
    else if  (strcmp(argv[loc],"-psi") == 0) {
      if (Tcl_GetDouble(interp, argv[loc+1], &psip) != TCL_OK) {
        opserr << "nDMaterial Templated3Dep -EPS - invalid psip " << argv[loc+1] << endln;
        cleanup(argv);
        return 0;
      }
      loc+=2;
    }
  } //end of while


  //stresstensor stressp(values_stp);
  //delete [] values_stp;
//  cout << " Tcl check e " << ed << " ec " << ecd << " lambda " << lamd << " p_ref " << p_refd << " a " << ad << "\n";


//   int                 Elasticflagp = 1;
//   double              Evp   = 0.0;
//   double              nuhvp = 0.0;
//   double              Ghvp = 0.0;



//ZCMay04   EPS = new EPState(Eod,            //      PState(double              Eod,
//ZCMay04                     Ed,             //             double              Ed,
//ZCMay04                     nu,             //             double              nu,
//ZCMay04                     rho,            //             double              rho,
//ZCMay04                     stressp,        //             const stresstensor  stressp,
//ZCMay04                     strainp,        //             const straintensor  strainp,
//ZCMay04                     Estrainp,       //             const straintensor  Estrainp,
//ZCMay04                     Pstrainp,       //             const straintensor  Pstrainp,
//ZCMay04                     NoS,            //       int                 NScalarp,
//ZCMay04                     scalars,        //       const double       *Scalarp,
//ZCMay04                     NoD,            //       int                 NTensorp,
//ZCMay04                     tensors,        //       const stresstensor *Tensorp,
//ZCMay04                     Elasticflagp,   //       int                 Elasticflagp = 1,
//ZCMay04                     Evp,            //       double              Evp   = 0.0,
//ZCMay04                     nuhvp,          //       double              nuhvp = 0.0,
//ZCMay04                     Ghvp,           //       double              Ghvp = 0.0,
//ZCMay04                     eod,            //       double              eop = 0.85,
//ZCMay04                     ecd,            //       double              ecp = 0.80,
//ZCMay04                     lamd,           //       double              Lam = 0.025,
//ZCMay04                     p_refd,         //       double              pop = 100.0,
//ZCMay04                     ad);            //       double              ap = 0.5

  EPS = new EPState(stressp,        //             const stresstensor  stressp,
                    strainp,        //             const straintensor  strainp,
                    Estrainp,       //             const straintensor  Estrainp,
                    Pstrainp,       //             const straintensor  Pstrainp,
                    NoS,            //       int                 NScalarp,
                    scalars,        //       const double       *Scalarp,
                    NoD,            //       int                 NTensorp,
                    tensors,        //       const stresstensor *Tensorp,
		    ep,
		    psip);
  if (EPS == 0) {
    opserr << "nDMaterial Templated3Dep -EPS - out of memory\n";
  }

  cleanup(argv);

  //if (NoD > 0)
  //   for (int i=0; i<=NoD; i++) {
  //     delete tensors;
  //   }
  if (NoS > 0)
     delete [] scalars;

  return EPS;
}



// Function - to create an EvolutionLaw_S object
EvolutionLaw_S *EvaluateEvolutionLawSArgs(ClientData clientData, Tcl_Interp *interp, TCL_Char *tclString)
{
  int argc;
  TCL_Char **argv;

  // split the list
  if (Tcl_SplitList(interp, tclString, &argc, &argv) != TCL_OK) {
    return 0;
  }

  //1. Linear scalar evolution law: eq
  EvolutionLaw_S *ELS = 0;

  if ((strcmp(argv[0],"-Leq") == 0) || (strcmp(argv[0],"-leq") == 0)) {
    double alpha = 0.0;
    if (argc > 1)
      if (Tcl_GetDouble(interp, argv[1], &alpha) != TCL_OK) {
  opserr << "nDMaterial Templated3Dep - invalid alpha " << argv[1] << endln;
  cleanup(argv);
  return 0;
      }
    ELS = new EvolutionLaw_L_Eeq(alpha);
  }

  //2. Nonlinear scalar evolution law: p for Cam clay
  //EvolutionLaw_S *ELS = 0;
  if ((strcmp(argv[0],"-NLp") == 0) || (strcmp(argv[0],"-nlp") == 0)) {
    double eod = 0.65, lambdad=0.19, kappad=0.06;
    if (argc ==3 ) {
      if (Tcl_GetDouble(interp, argv[1], &eod) != TCL_OK) {
  opserr << "nDMaterial Templated3Dep - invalid eo " << argv[1] << endln;

  //cleanup(argv);
  return 0;
      }
      if (Tcl_GetDouble(interp, argv[2], &lambdad) != TCL_OK) {
  opserr << "nDMaterial Templated3Dep - invalid lambda " << argv[2] << endln;

  //cleanup(argv);
  return 0;
      }
      if (Tcl_GetDouble(interp, argv[3], &kappad) != TCL_OK) {
  opserr << "nDMaterial Templated3Dep - invalid kappa " << argv[3] << endln;

  cleanup(argv);
  return 0;
      }

    }
    ELS = new EvolutionLaw_NL_Ep(eod, lambdad, kappad);
  }

  cleanup(argv);
  return ELS;
}

// Function - to create an EvolutionLaw_T object
EvolutionLaw_T *EvaluateEvolutionLawTArgs(ClientData clientData, Tcl_Interp *interp, TCL_Char *tclString)
{
  int argc;
  TCL_Char **argv;

  // split the list
  if (Tcl_SplitList(interp, tclString, &argc, &argv) != TCL_OK) {
    return 0;
  }

  EvolutionLaw_T *ELT = 0;
  
  // Linear Tensorial      
  if ((strcmp(argv[0],"-Linear") == 0) || (strcmp(argv[0],"-linear") == 0)) {
    double alpha = 0.0;
    if (argc > 1)
      if (Tcl_GetDouble(interp, argv[1], &alpha) != TCL_OK) {
    opserr << "nDMaterial Templated3Dep - invalid alpha " << argv[1] << endln;
    cleanup(argv);
    return 0;
      }
    ELT = new EvolutionLaw_L_Eij(alpha);
  }
  
  // Armstrong-Freederick Nolinear Tensorial Evolution Law 
  else if ((strcmp(argv[0],"-NLEij") == 0) || (strcmp(argv[0],"-AF") == 0)) {
    double had = 0.0;
    double Crd = 0.0;    
    
      if (Tcl_GetDouble(interp, argv[1], &had) != TCL_OK) {
    opserr << "nDMaterial Templated3Dep - invalid ha " << argv[1] << endln;
    cleanup(argv);
    return 0;
      }
      
      if (Tcl_GetDouble(interp, argv[2], &Crd) != TCL_OK) {
    opserr << "nDMaterial Templated3Dep - invalid Cr " << argv[2] << endln;
    cleanup(argv);
    return 0;
      }      
    
    ELT = new EvolutionLaw_NL_Eij(had, Crd);

  }
  
  //Manzari-Dafalias Two Surface model
  else if ((strcmp(argv[0],"-NLEijMD") == 0) || (strcmp(argv[0],"-MD") == 0 )) {
    double eod = 0.85, ad = 0.5, Lambdad = 0.025, ec_refd = 0.8, p_refd = 160.0;
    double Mcd = 1.14, Med = 1.14, kc_bd = 3.975, kc_dd = 4.200;
    double ke_bd = 2.000, ke_dd = 0.07, hod = 1200, Cmd = 0.00, Aod = 2.64, Fmaxd = 100, Cfd = 100;

    int loc = 1;
    // switch on remaining args;
    while (loc < argc) {

      if  (strcmp(argv[loc],"-eo") == 0) {
        if (Tcl_GetDouble(interp, argv[loc+1], &eod) != TCL_OK) {
        opserr << "nDMaterial Templated3Dep -EPS - invalid eo " << argv[loc+1] << endln;
        cleanup(argv);
        return 0;
        }
        loc+=2;
      }
      else if  (strcmp(argv[loc],"-a") == 0) {
        if (Tcl_GetDouble(interp, argv[loc+1], &ad) != TCL_OK) {
        opserr << "nDMaterial Templated3Dep -EPS - invalid a " << argv[loc+1] << endln;
        cleanup(argv);
        return 0;
        }
        loc+=2;
      }
      else if  (strcmp(argv[loc],"-Mc") == 0) {
        if (Tcl_GetDouble(interp, argv[loc+1], &Mcd) != TCL_OK) {
        opserr << "nDMaterial Templated3Dep -EPS - invalid Mc " << argv[loc+1] << endln;
        cleanup(argv);
        return 0;
        }
        loc+=2;
      }
      else if  (strcmp(argv[loc],"-Me") == 0) {
        if (Tcl_GetDouble(interp, argv[loc+1], &Med) != TCL_OK) {
        opserr << "nDMaterial Templated3Dep -EPS - invalid ec " << argv[loc+1] << endln;
        cleanup(argv);
        return 0;
        }
        loc+=2;
      }
      else if  (strcmp(argv[loc],"-Lambda") == 0) {
        if (Tcl_GetDouble(interp, argv[loc+1], &Lambdad) != TCL_OK) {
        opserr << "nDMaterial Templated3Dep -EPS - invalid lambdad " << argv[loc+1] << endln;
        cleanup(argv);
        return 0;
        }
        loc+=2;
      }
      else if  (strcmp(argv[loc],"-ec_ref") == 0) {
        if (Tcl_GetDouble(interp, argv[loc+1], &ec_refd) != TCL_OK) {
        opserr << "nDMaterial Templated3Dep -EPS - invalid ec_ref " << argv[loc+1] << endln;
        cleanup(argv);
        return 0;
        }
        loc+=2;
      }
      else if  (strcmp(argv[loc],"-p_ref") == 0) {
        if (Tcl_GetDouble(interp, argv[loc+1], &p_refd) != TCL_OK) {
        opserr << "nDMaterial Templated3Dep -EPS - invalid p_ref " << argv[loc+1] << endln;
        cleanup(argv);
        return 0;
        }
        loc+=2;
      }

      else if  (strcmp(argv[loc],"-kc_b") == 0) {
        if (Tcl_GetDouble(interp, argv[loc+1], &kc_bd) != TCL_OK) {
        opserr << "nDMaterial Templated3Dep -EPS - invalid kc_b " << argv[loc+1] << endln;
        cleanup(argv);
        return 0;
        }
        loc+=2;
      }
      else if  (strcmp(argv[loc],"-kc_d") == 0) {
        if (Tcl_GetDouble(interp, argv[loc+1], &kc_dd) != TCL_OK) {
        opserr << "nDMaterial Templated3Dep -EPS - invalid kc_d " << argv[loc+1] << endln;
        cleanup(argv);
        return 0;
        }
        loc+=2;
      }
      else if  (strcmp(argv[loc],"-ke_b") == 0) {
        if (Tcl_GetDouble(interp, argv[loc+1], &ke_bd) != TCL_OK) {
        opserr << "nDMaterial Templated3Dep -EPS - invalid ke_b " << argv[loc+1] << endln;
        cleanup(argv);
        return 0;
        }
        loc+=2;
      }
      else if  (strcmp(argv[loc],"-ke_d") == 0) {
        if (Tcl_GetDouble(interp, argv[loc+1], &ke_dd) != TCL_OK) {
        opserr << "nDMaterial Templated3Dep -EPS - invalid ke_d " << argv[loc+1] << endln;
        cleanup(argv);
        return 0;
        }
        loc+=2;
      }
      else if  (strcmp(argv[loc],"-ho") == 0) {
        if (Tcl_GetDouble(interp, argv[loc+1], &hod) != TCL_OK) {
        opserr << "nDMaterial Templated3Dep -EPS - invalid a " << argv[loc+1] << endln;
        cleanup(argv);
        return 0;
        }
        loc+=2;
      }
      else if  (strcmp(argv[loc],"-Cm") == 0) {
        if (Tcl_GetDouble(interp, argv[loc+1], &Cmd) != TCL_OK) {
        opserr << "nDMaterial Templated3Dep -EPS - invalid a " << argv[loc+1] << endln;
        cleanup(argv);
        return 0;
        }
        loc+=2;
      }
      else if  (strcmp(argv[loc],"-Ao") == 0) {
        if (Tcl_GetDouble(interp, argv[loc+1], &Aod) != TCL_OK) {
        opserr << "nDMaterial Templated3Dep -EPS - invalid a " << argv[loc+1] << endln;
        cleanup(argv);
        return 0;
        }
        loc+=2;
      }
      else if  (strcmp(argv[loc],"-Fmax") == 0) {
        if (Tcl_GetDouble(interp, argv[loc+1], &Fmaxd) != TCL_OK) {
        opserr << "nDMaterial Templated3Dep -EPS - invalid a " << argv[loc+1] << endln;
        cleanup(argv);
        return 0;
        }
        loc+=2;
      }
      else if  (strcmp(argv[loc],"-Cf") == 0) {
        if (Tcl_GetDouble(interp, argv[loc+1], &Cfd) != TCL_OK) {
        opserr << "nDMaterial Templated3Dep -EPS - invalid a " << argv[loc+1] << endln;
        cleanup(argv);
        return 0;
        }
        loc+=2;
      }

    }//end of while

    ELT = new EvolutionLaw_NL_EijMD(eod, ad, Mcd, Med, Lambdad, ec_refd, p_refd, kc_bd, kc_dd, ke_bd, ke_dd, hod, Cmd, Aod, Fmaxd, Cfd);
  }

  cleanup(argv);
  return ELT;
}


