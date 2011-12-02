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
//# DATE:              19AUg2003
//# UPDATE HISTORY:    Sept 2003
//#		       28May2004
//#
//===============================================================================
#include <tcl.h>
#include <OPS_Globals.h>

#include <Domain.h>

#include <ErrorHandler.h>
#include <TclModelBuilder.h>

#include <NeoHookeanCompressible3D.h>
#include <FDdecoupledElastic3D.h>

#include <W.h>
#include <LogWEnergy.h>
#include <MooneyRivlinWEnergy.h>
#include <NeoHookeanWEnergy.h>
#include <OgdenWEnergy.h>
#include <SimoPisterWEnergy.h>
#include <OgdenSimoWEnergy.h>
#include <MooneyRivlinSimoWEnergy.h>

//static void cleanup(TCL_Char **argv)
//{
//    Tcl_Free((char *) argv);
//}

FiniteDeformationElastic3D *
TclModelBuilder_addFiniteDeformationElastic3D(ClientData clientData, Tcl_Interp *interp,  int argc,
          TCL_Char **argv, TclModelBuilder *theTclBuilder, int eleArgStart)
{
  //int argc;
  //TCL_Char **argv;

  int tag = 0;
  
  int loc = eleArgStart;

  FiniteDeformationElastic3D *theMaterial = 0;

  if (Tcl_GetInt(interp, argv[loc+1], &tag) != TCL_OK) {
    opserr << "Warning: nDMaterial FiniteDeformationElastic3D - invalid tag " << argv[loc+1] << "\n";
    exit (-1);
  }
  
  // Neo-Hookean (Compressible)
  if ( (strcmp(argv[loc+2],"NeoHookean3D") == 0) || (strcmp(argv[loc+2],"NeoHookeanCompressible3D") == 0) ) {    

    double rho_in = 0.0;
    double K_in = 0.0;
    double G_in = 0.0;

    if (argc < 7) {
      opserr << "Warning: NeoHookeanCompressible3D -insufficient number of arguments\n";
      exit (-1);
    }

    if (Tcl_GetDouble(interp, argv[loc+3], &K_in) != TCL_OK) {
      opserr << "nDMaterial NeoHookeanCompressible3D - invalid K " << argv[loc+3] << "\n";
      exit (-1);
    }

    if (Tcl_GetDouble(interp, argv[loc+4], &G_in) != TCL_OK) {
      opserr << "Warning: nDMaterial NeoHookeanCompressible3D - invalid G " << argv[loc+4] << "\n";
      exit (-1);
    }

    if (Tcl_GetDouble(interp, argv[loc+5], &rho_in) != TCL_OK) {
      opserr << "Warning: nDMaterial NeoHookeanCompressible3D - invalid rho " << argv[loc+5] << "\n";
      exit (-1);
    }
    
    theMaterial = new NeoHookeanCompressible3D(tag, K_in, G_in, rho_in);
  }

  // DecoupledLog3D
  else if ( (strcmp(argv[loc+2],"DecoupledLog3D") == 0) || (strcmp(argv[loc+2],"DecoupledLogarithmic3D") == 0) ) {
    
    double K_in = 0.0;
    double G_in = 0.0;
    double rho_in = 0.0;
    WEnergy  *wenergy =0;

    if (argc < 7) {
      opserr << "Warning: DecoupledLogarithmic3D -insufficient number of arguments\n";
      exit (-1);
    }

    if (Tcl_GetDouble(interp, argv[loc+3], &K_in) != TCL_OK) {
      opserr << "Warning: nDMaterial DecoupledLogarithmic3D - invalid K " << argv[loc+3] << "\n";
      exit (-1);
    }

    if (Tcl_GetDouble(interp, argv[loc+4], &G_in) != TCL_OK) {
      opserr << "Warning: nDMaterial DecoupledLogarithmic3D - invalid G " << argv[loc+4] << "\n";
      exit (-1);
    }

    if (Tcl_GetDouble(interp, argv[loc+5], &rho_in) != TCL_OK) {
      opserr << "Warning: nDMaterial DecoupledLogarithmic3D - invalid rho " << argv[loc+5] << "\n";
      exit (-1);
    }
    
    wenergy = new LogWEnergy(K_in, G_in);
    
    if ( wenergy != 0 ) {
      theMaterial = new FDdecoupledElastic3D(tag, wenergy, rho_in);
    }
    else {
      opserr << "Error: nDMaterial DecoupledLogarithmic3D -invalid material model\n";
      exit (-1);
    }
  
  }
  
  // DecoupledNeoHookean3D
  else if ( (strcmp(argv[loc+2],"DecoupledNeoHookean3D") == 0) || (strcmp(argv[loc+2],"DecoupledNH3D") == 0) ) {
    
    double K_in = 0.0;
    double G_in = 0.0;
    double rho_in = 0.0;
    WEnergy  *wenergy =0;

    if (argc < 7) {
      opserr << "Warning: DecoupledNeoHookean3D -insufficient number of arguments\n";
      exit (-1);
    }

    if (Tcl_GetDouble(interp, argv[loc+3], &K_in) != TCL_OK) {
      opserr << "Warning: nDMaterial DecoupledNeoHookean3D - invalid K " << argv[loc+3] << "\n";
      exit (-1);
    }

    if (Tcl_GetDouble(interp, argv[loc+4], &G_in) != TCL_OK) {
      opserr << "Warning: nDMaterial DecoupledNeoHookean3D - invalid G " << argv[loc+4] << "\n";
      exit (-1);
    }

    if (Tcl_GetDouble(interp, argv[loc+5], &rho_in) != TCL_OK) {
      opserr << "Warning: nDMaterial DecoupledNeoHookean3D - invalid rho " << argv[loc+5] << "\n";
      exit (-1);
    }
    
    wenergy = new NeoHookeanWEnergy(K_in, G_in);
    
    if ( wenergy != 0 ) {
      theMaterial = new FDdecoupledElastic3D(tag, wenergy, rho_in);
    }
    else {
      opserr << "Error: nDMaterial DecoupledNeoHookean3D -invalid material model\n";
      exit (-1);
    }
  
  }

  // DecoupledMooneyRivlinSimo3D
  else if ( (strcmp(argv[loc+2],"DecoupledMooneyRivlinSimo3D") == 0) || (strcmp(argv[loc+2],"DecoupledMRS3D") == 0) ) {
    
    double c1_in = 0.0;
    double c2_in = 0.0;
    double K_in = 0.0;
    double rho_in = 0.0;
    WEnergy  *wenergy =0;

    if (argc < 8) {
      opserr << "Warning: DecoupledNeoHookean3D -insufficient number of arguments\n";
      exit (-1);
    }

    if (Tcl_GetDouble(interp, argv[loc+3], &c1_in) != TCL_OK) {
      opserr << "Warning: nDMaterial DecoupledNeoHookean3D - invalid c1 " << argv[loc+3] << "\n";
      exit (-1);
    }

    if (Tcl_GetDouble(interp, argv[loc+4], &c2_in) != TCL_OK) {
      opserr << "Warning: nDMaterial DecoupledNeoHookean3D - invalid c2 " << argv[loc+4] << "\n";
      exit (-1);
    }

    if (Tcl_GetDouble(interp, argv[loc+5], &K_in) != TCL_OK) {
      opserr << "Warning: nDMaterial DecoupledNeoHookean3D - invalid K " << argv[loc+5] << "\n";
      exit (-1);
    }

    if (Tcl_GetDouble(interp, argv[loc+6], &rho_in) != TCL_OK) {
      opserr << "Warning: nDMaterial DecoupledNeoHookean3D - invalid rho " << argv[loc+6] << "\n";
      exit (-1);
    }
    
    wenergy = new MooneyRivlinSimoWEnergy(c1_in, c2_in, K_in);
    
    if ( wenergy != 0 ) {
      theMaterial = new FDdecoupledElastic3D(tag, wenergy, rho_in);
    }
    else {
      opserr << "Error: nDMaterial DecoupledMooneyRivlinSimo3D -invalid material model\n";
      exit (-1);
    }
  
  }

  // DecoupleOgdenSimo3D
  else if ( (strcmp(argv[loc+2],"DecoupledOgdenSimo3D") == 0) || (strcmp(argv[loc+2],"DecoupledOS3D") == 0) ) {
    
    int N_in = 0;
    double K_in = 0.0;
    double rho_in = 0.0;
    WEnergy  *wenergy =0;

    if (argc > 2*N_in+5) {
      if (Tcl_GetInt(interp, argv[loc+3], &N_in) != TCL_OK) {
        opserr << "Warning: invalid vector parameter number for Ogden Strain Energy Function\n";
        exit (-1);
      }
      
      double *cr_in = new double[N_in];
      double *mur_in = new double[N_in];

      for (int i=0; i<N_in; i++) {
        if (Tcl_GetDouble(interp, argv[loc+4+i], &cr_in[i]) != TCL_OK) {
          opserr << "Warning: invalid parameter for Ogden Strain Energy Function\n";
          exit (-1);
        }
        if (Tcl_GetDouble(interp, argv[loc+4+N_in+i], &mur_in[i]) != TCL_OK) {
          opserr << "Warning: invalid parameter for Ogden Strain Energy Function\n";
          exit (-1);
        }
      }

      if (Tcl_GetDouble(interp, argv[loc+2*N_in+4], &K_in) != TCL_OK) {
        opserr << "Warning: invalid Bulk Modulus number for Ogden Strain Energy Function\n";
        exit (-1);
      }

      if (Tcl_GetDouble(interp, argv[loc+2*N_in+5], &rho_in) != TCL_OK) {
        opserr << "Warning: nDMaterial DecoupledSimoPister3D - invalid rho\n";
        exit (-1);
      }
    
      wenergy = new OgdenSimoWEnergy(N_in, cr_in, mur_in, K_in);
    
    if ( wenergy != 0 ) {
      theMaterial = new FDdecoupledElastic3D(tag, wenergy, rho_in);
      //if ( cr_in  ) delete cr_in;
      //if ( mur_in ) delete mur_in;
    }
    else {
      opserr << "Error: nDMaterial DecoupledOgdenSimo3D -invalid material model\n";
      exit (-1);
    }
  
    }
  }

  // DecoupledMooneyRivlin3D
  else if ( (strcmp(argv[loc+2],"DecoupledMooneyRivlin3D") == 0) || (strcmp(argv[loc+2],"DecoupledMR3D") == 0) ) {
    
    double c1_in = 0.0;
    double c2_in = 0.0;
    double rho_in = 0.0;
    WEnergy  *wenergy =0;

    if (argc < 7) {
      opserr << "Warning: DecoupledNeoHookean3D -insufficient number of arguments\n";
      exit (-1);
    }

    if (Tcl_GetDouble(interp, argv[loc+3], &c1_in) != TCL_OK) {
      opserr << "Warning: nDMaterial DecoupledNeoHookean3D - invalid c1 " << argv[loc+3] << "\n";
      exit (-1);
    }

    if (Tcl_GetDouble(interp, argv[loc+4], &c2_in) != TCL_OK) {
      opserr << "Warning: nDMaterial DecoupledNeoHookean3D - invalid c2 " << argv[loc+4] << "\n";
      exit (-1);
    }

    if (Tcl_GetDouble(interp, argv[loc+5], &rho_in) != TCL_OK) {
      opserr << "Warning: nDMaterial DecoupledNeoHookean3D - invalid rho " << argv[loc+5] << "\n";
      exit (-1);
    }
    
    wenergy = new MooneyRivlinWEnergy(c1_in, c2_in);
    
    if ( wenergy != 0 ) {
      theMaterial = new FDdecoupledElastic3D(tag, wenergy, rho_in);
    }
    else {
      opserr << "Error: nDMaterial DecoupledMooneyRivlin3D -invalid material model\n";
      exit (-1);
    }
  
  }

  // DecoupleOgden3D
  else if ( (strcmp(argv[loc+2],"DecoupledOgden3D") == 0) ) {
    
    int N_in = 0;
    double rho_in = 0.0;    
    WEnergy  *wenergy =0;

    if (argc > 2*N_in+4) {
      if (Tcl_GetInt(interp, argv[loc+3], &N_in) != TCL_OK) {
        opserr << "Warning: invalid vector parameter number for Ogden Strain Energy Function\n";
        exit (-1);
      }
      
      double *cr_in = new double[N_in];
      double *mur_in = new double[N_in];

      for (int i=0; i<N_in; i++) {
        if (Tcl_GetDouble(interp, argv[loc+4+i], &cr_in[i]) != TCL_OK) {
          opserr << "Warning: invalid parameter for Ogden Strain Energy Function\n";
          exit (-1);
        }
        if (Tcl_GetDouble(interp, argv[loc+4+N_in+i], &mur_in[i]) != TCL_OK) {
          opserr << "Warning: invalid parameter for Ogden Strain Energy Function\n";
          exit (-1);
        }
      }

      if (Tcl_GetDouble(interp, argv[loc+2*N_in+4], &rho_in) != TCL_OK) {
        opserr << "Warning: nDMaterial Ogden - invalid rho\n";
        exit (-1);
      }
    
      wenergy = new OgdenWEnergy(N_in, cr_in, mur_in);
    
    if ( wenergy != 0 ) {
      theMaterial = new FDdecoupledElastic3D(tag, wenergy, rho_in);
      //if ( cr_in  ) delete cr_in;
      //if ( mur_in ) delete mur_in;
    }
    else {
      opserr << "Error: nDMaterial DecoupledOgden3D -invalid material model\n";
      exit (-1);
    }
  
    }
  }

  // DecoupledSimoPister3D
  else if ( (strcmp(argv[loc+2],"DecoupledSimoPister3D") == 0) || (strcmp(argv[loc+2],"DecoupledSP3D") == 0) ) {
    
    double K_in = 0.0;
    double rho_in = 0.0;
    WEnergy  *wenergy =0;

    if (argc < 6) {
      opserr << "Warning: DecoupledSimoPister3D -insufficient number of arguments\n";
      exit (-1);
    }

    if (Tcl_GetDouble(interp, argv[loc+3], &K_in) != TCL_OK) {
      opserr << "Warning: nDMaterial DecoupledSimoPister3D - invalid K " << argv[loc+3] << "\n";
      exit (-1);
    }

    if (Tcl_GetDouble(interp, argv[loc+4], &rho_in) != TCL_OK) {
      opserr << "Warning: nDMaterial DecoupledSimoPister3D - invalid rho " << argv[loc+4] << "\n";
      exit (-1);
    }
    
    wenergy = new SimoPisterWEnergy(K_in);
    
    if ( wenergy != 0 ) {
      theMaterial = new FDdecoupledElastic3D(tag, wenergy, rho_in);
    }
    else {
      opserr << "Error: nDMaterial DecoupledSimoPister3D -invalid material model\n";
      exit (-1);
    }
  
  }
  
  // Else
  else {
    opserr << "Error: nDMaterial - unknown FiniteDeformationElastic3D material model\n";
    exit (-1);
  }

  //cleanup(argv);
  
  return theMaterial;

}




