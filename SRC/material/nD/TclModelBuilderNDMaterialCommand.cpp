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
** With a lot of additions by                                         **
**   Boris Jeremic    (jeremic@ucdavis.edu)                           **
**   Zaohui Yang      (zhyang@ucdavis.edu)                            **
**   Zhao Cheng       (zcheng@ucdavis.edu)                            **
**                                                                    **
**                                                                    **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.44 $
// $Date: 2010-02-05 00:08:35 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/TclModelBuilderNDMaterialCommand.cpp,v $


// Description: This file contains the function invoked when the user invokes
// the nDMaterial command in the interpreter.
//
// What: "@(#) TclModelBuilderNDMaterialCommand.C, revA"

#include <TclModelBuilder.h>
#include <elementAPI.h>
#include <packages.h>

#include <PressureDependentElastic3D.h>
#include <J2Plasticity.h>
#include <MultiaxialCyclicPlasticity.h> //Gang Wang

#include <PlaneStressMaterial.h>
#include <PlaneStrainMaterial.h>  // Antonios Vytiniotis:
#include <PlateFiberMaterial.h>
#include <UniaxialFiber2d.h>
#include <UniaxialFiber3d.h>

//start Yuli Huang & Xinzheng Lu 
#include <PlateRebarMaterial.h>
#include <PlateFromPlaneStressMaterial.h>
#include <ConcreteS.h>
#include <PlaneStressUserMaterial.h>
//end Yuli Huang & Xinzheng Lu 

#include <CapPlasticity.h>          // Quan Gu & ZhiJian Qiu  2013
#include <SimplifiedJ2.h>           // Quan Gu & ZhiJian Qiu 2013-6-26
#include <PlaneStressSimplifiedJ2.h>// Quan Gu & ZhiJian Qiu 2013-6-26 

#include <BeamFiberMaterial.h>
#include <ConcreteMcftNonLinear5.h>
#include <ConcreteMcftNonLinear7.h>

#include <PressureIndependMultiYield.h>
#include <PressureDependMultiYield.h>
#include <PressureDependMultiYield02.h>
#include <PressureDependMultiYield03.h>
#include <FluidSolidPorousMaterial.h>

#include <J2PlasticityThermal.h>   //added by L.Jiang [SIF]
#include <PlateFiberMaterialThermal.h>//L.Jiang [SIF]
#include <PlateFromPlaneStressMaterialThermal.h>   //Liming Jiang [SIF]
#include <PlateRebarMaterialThermal.h> //Liming Jiang [SIF]

#include <MultiYieldSurfaceClay.h>
#include <string.h>

extern NDMaterial *
Tcl_addWrapperNDMaterial(matObj *, ClientData, Tcl_Interp *,  int, TCL_Char **, TclModelBuilder *);

extern  void *OPS_ReinforcedConcretePlaneStressMaterial(void);
extern  void *OPS_FAReinforcedConcretePlaneStressMaterial(void);
extern  void *OPS_FAFourSteelRCPlaneStressMaterial(void);
extern  void *OPS_RAFourSteelRCPlaneStressMaterial(void);
extern  void *OPS_PrestressedConcretePlaneStressMaterial(void);
extern  void *OPS_FAPrestressedConcretePlaneStressMaterial(void);
extern  void *OPS_FAFourSteelPCPlaneStressMaterial(void);
extern  void *OPS_RAFourSteelPCPlaneStressMaterial(void);
extern  void *OPS_MaterialCMM(void);
extern  void *OPS_NewMaterialCMM(void);
extern  void *  OPS_NewPlasticDamageConcrete3d(void);
extern  void *  OPS_NewPlasticDamageConcretePlaneStress(void);
extern  void *OPS_ElasticIsotropicMaterial(void);
extern  void *OPS_ElasticIsotropic3D(void);
extern  void *OPS_ElasticOrthotropicMaterial(void);
extern  void *OPS_DruckerPragerMaterial(void);
extern  void *OPS_BoundingCamClayMaterial(void);
extern  void *OPS_ContactMaterial2DMaterial(void);
extern  void *OPS_ContactMaterial3DMaterial(void);
extern  void *OPS_InitialStateAnalysisWrapperMaterial(void);
extern  void *OPS_ManzariDafaliasMaterial(void);
extern  void *OPS_ManzariDafaliasMaterialRO(void);
extern  void *OPS_PM4SandMaterial(void);
extern  void *OPS_PM4SiltMaterial(void);
extern  void *OPS_J2CyclicBoundingSurfaceMaterial(void);
extern  void *OPS_CycLiqCPMaterial(void);
extern  void *OPS_CycLiqCPSPMaterial(void);
extern  void *OPS_InitStressNDMaterial(void);
extern  void *OPS_StressDensityMaterial(void);
extern  void *OPS_J2Plasticity(void);
extern  void *OPS_J2BeamFiber2dMaterial(void);
extern  void *OPS_J2BeamFiber3dMaterial(void);
extern  void *OPS_J2PlateFibreMaterial(void);
extern  void *OPS_PlaneStressLayeredMaterial(void);
extern  void *OPS_PlaneStressRebarMaterial(void);
extern  void *OPS_PlateFiberMaterial(void);
extern  void *OPS_BeamFiberMaterial(void);
extern  void *OPS_BeamFiberMaterial2d(void);
extern  void *OPS_BeamFiberMaterial2dPS(void);
extern void *OPS_LinearCap(void);
extern void *OPS_AcousticMedium(void);
extern void* OPS_UVCmultiaxial(void);
extern void* OPS_UVCplanestress(void);

extern  void *OPS_ElasticIsotropicMaterialThermal(void);  //L.Jiang [SIF]
extern  void *OPS_DruckerPragerMaterialThermal(void);//L.Jiang [SIF]
//extern  void *OPS_PlasticDamageConcretePlaneStressThermal(void);//L.Jiang [SIF]

#ifdef _HAVE_Faria1998
extern void *OPS_NewFaria1998Material(void);
extern void *OPS_NewConcreteMaterial(void);
#endif

extern  void *OPS_FSAMMaterial(void); // K Kolozvari      

#ifdef _HAVE_Damage2p
extern void *OPS_Damage2p(void);
#endif


NDMaterial *
TclModelBuilder_addFeapMaterial(ClientData clientData, Tcl_Interp *interp,
				int argc, TCL_Char **argv,
				TclModelBuilder *theTclBuilder);


extern int OPS_ResetInput(ClientData clientData, 
			  Tcl_Interp *interp,  
			  int cArg, 
			  int mArg, 
			  TCL_Char **argv, 
			  Domain *domain,
			  TclModelBuilder *builder);

typedef struct ndMaterialPackageCommand {
  char *funcName;
  void * (*funcPtr)(); 
  struct ndMaterialPackageCommand *next;
} NDMaterialPackageCommand;

static NDMaterialPackageCommand *theNDMaterialPackageCommands = NULL;

static void printCommand(int argc, TCL_Char **argv)
{
    opserr << "Input command: ";
    for (int i=0; i<argc; i++)
	opserr << argv[i] << " ";
    opserr << endln;
}

int
TclModelBuilderNDMaterialCommand (ClientData clientData, Tcl_Interp *interp, int argc,
				  TCL_Char **argv, TclModelBuilder *theTclBuilder)
{
    // Make sure there is a minimum number of arguments
    if (argc < 3) {
	opserr << "WARNING insufficient number of ND material arguments\n";
	opserr << "Want: nDMaterial type? tag? <specific material args>" << endln;
	return TCL_ERROR;
    }

    OPS_ResetInput(clientData, interp, 2, argc, argv, 0, theTclBuilder);	  

    // Pointer to an ND material that will be added to the model builder
    NDMaterial *theMaterial = 0;

    // Check argv[1] for ND material type

    if ((strcmp(argv[1],"ReinforcedConcretePlaneStress") == 0) || (strcmp(argv[1],"ReinforceConcretePlaneStress") == 0)) {

      void *theMat = OPS_ReinforcedConcretePlaneStressMaterial();
      if (theMat != 0) 
	theMaterial = (NDMaterial *)theMat;
      else 
	return TCL_ERROR;
    }

    else if ((strcmp(argv[1],"PlasticDamageConcrete") == 0) || (strcmp(argv[1],"PlasticDamageConcrete3d") == 0)) {

      void *theMat = OPS_NewPlasticDamageConcrete3d();
      if (theMat != 0)  {
	theMaterial = (NDMaterial *)theMat;
      }
      else 
	return TCL_ERROR;
    }

    else if ((strcmp(argv[1],"PlasticDamageConcretePlaneStress") == 0)) {
      void *theMat = OPS_NewPlasticDamageConcretePlaneStress();
      if (theMat != 0) 
	theMaterial = (NDMaterial *)theMat;
      else 
	return TCL_ERROR;
    }

    else if ((strcmp(argv[1],"InitStressMaterial") == 0) || (strcmp(argv[1],"InitStress") == 0)) {
      void *theMat = OPS_InitStressNDMaterial();
      if (theMat != 0) 
        theMaterial = (NDMaterial *)theMat;
      else 
        return TCL_ERROR;
    }

    else if (strcmp(argv[1],"PlaneStressLayeredMaterial") == 0) {
      void *theMat = OPS_PlaneStressLayeredMaterial();
      if (theMat != 0) 
        theMaterial = (NDMaterial *)theMat;
      else 
        return TCL_ERROR;
    }

    else if (strcmp(argv[1],"PlaneStressRebarMaterial") == 0) {
      void *theMat = OPS_PlaneStressRebarMaterial();
      if (theMat != 0) 
        theMaterial = (NDMaterial *)theMat;
      else 
        return TCL_ERROR;
    }
    
    
    else if (strcmp(argv[1],"J2BeamFiber") == 0) {
      void *theMat = 0;
      if (theTclBuilder->getNDM() == 2)
	theMat = OPS_J2BeamFiber2dMaterial();
      else
	theMat = OPS_J2BeamFiber3dMaterial();

      if (theMat != 0) 
        theMaterial = (NDMaterial *)theMat;
      else 
        return TCL_ERROR;
    }

    else if (strcmp(argv[1],"J2PlateFibre") == 0) {
      void *theMat = OPS_J2PlateFibreMaterial();
      if (theMat != 0) 
        theMaterial = (NDMaterial *)theMat;
      else 
        return TCL_ERROR;
    }

#ifdef _HAVE_Faria1998
    else if (strcmp(argv[1],"Faria1998") == 0) {
      void *theMat = OPS_NewFaria1998Material();
      if (theMat != 0) 
	theMaterial = (NDMaterial *)theMat;
      else 
       	return TCL_ERROR;
    }
    else if (strcmp(argv[1],"Concrete") == 0) {
      void *theMat = OPS_NewConcreteMaterial();
      if (theMat != 0) 
	theMaterial = (NDMaterial *)theMat;
      else 
       	return TCL_ERROR;
    }
#endif

    else if ((strcmp(argv[1],"FAReinforceConcretePlaneStress") == 0) || (strcmp(argv[1],"FAReinforcedConcretePlaneStress") == 0)) {

      void *theMat = OPS_FAReinforcedConcretePlaneStressMaterial();
      if (theMat != 0) 
	theMaterial = (NDMaterial *)theMat;
      else 
	return TCL_ERROR;
    }

    else if ((strcmp(argv[1],"RAFourSteelRCPlaneStress") == 0)){

      void *theMat = OPS_RAFourSteelRCPlaneStressMaterial();
      if (theMat != 0) 
	theMaterial = (NDMaterial *)theMat;
      else 
	return TCL_ERROR;
    }

    else if ((strcmp(argv[1],"FAFourSteelRCPlaneStress") == 0)){

      void *theMat = OPS_FAFourSteelRCPlaneStressMaterial();
      if (theMat != 0) 
	theMaterial = (NDMaterial *)theMat;
      else 
	return TCL_ERROR;
    }

#ifdef _HAVE_Damage2p
    else if ((strcmp(argv[1],"Damage2p") == 0)){

      void *theMat = OPS_Damage2p();
      if (theMat != 0) 
	theMaterial = (NDMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
#else
    else if ((strcmp(argv[1],"Damage2p") == 0)){
      opserr << "SORRY - Damage2p source code not available in this version\n";
      return TCL_ERROR;
    }
#endif

    else if ((strcmp(argv[1],"PrestressedConcretePlaneStress") == 0)){

      void *theMat = OPS_PrestressedConcretePlaneStressMaterial();
      if (theMat != 0) 
	theMaterial = (NDMaterial *)theMat;
      else 
	return TCL_ERROR;
    }

    else if ((strcmp(argv[1],"FAPrestressedConcretePlaneStress") == 0)){

      void *theMat = OPS_FAPrestressedConcretePlaneStressMaterial();
      if (theMat != 0) 
	theMaterial = (NDMaterial *)theMat;
      else 
	return TCL_ERROR;
    }

    else if ((strcmp(argv[1],"RAFourSteetPCPlaneStress") == 0)){

      void *theMat = OPS_RAFourSteelPCPlaneStressMaterial();
      if (theMat != 0) 
	theMaterial = (NDMaterial *)theMat;
      else 
	return TCL_ERROR;
    }

    else if ((strcmp(argv[1],"FAFourSteelPCPlaneStress") == 0)){

      void *theMat = OPS_FAFourSteelPCPlaneStressMaterial();
      if (theMat != 0) 
	theMaterial = (NDMaterial *)theMat;
      else 
	return TCL_ERROR;
    }


    else if ((strcmp(argv[1],"DruckerPrager") == 0)){

      void *theMat = OPS_DruckerPragerMaterial();
      if (theMat != 0) 
	theMaterial = (NDMaterial *)theMat;
      else 
	return TCL_ERROR;
    }



    else if ((strcmp(argv[1],"TruncatedDP") == 0)){

      void *theMat = OPS_LinearCap();
      if (theMat != 0) 
	theMaterial = (NDMaterial *)theMat;
      else 
	return TCL_ERROR;
    }


    // K Kolozvari                                                            
    else if ((strcmp(argv[1], "FSAM") == 0) || 
	     (strcmp(argv[1], "FSAM") == 0)) {
      
      void *theMat = OPS_FSAMMaterial();
      if (theMat != 0)
	theMaterial = (NDMaterial *)theMat;
      else
	return TCL_ERROR;
    }
    

    else if ((strcmp(argv[1],"AcousticMedium") == 0)){

      void *theMat = OPS_AcousticMedium();
      if (theMat != 0) 
	theMaterial = (NDMaterial *)theMat;
      else 
	return TCL_ERROR;
    }

    else if ((strcmp(argv[1],"UVCplanestress") == 0)){

      void *theMat = OPS_UVCplanestress();
      if (theMat != 0) 
	theMaterial = (NDMaterial *)theMat;
      else 
	return TCL_ERROR;
    }

    else if ((strcmp(argv[1],"UVCmultiaxial") == 0)){

      void *theMat = OPS_UVCmultiaxial();
      if (theMat != 0) 
	theMaterial = (NDMaterial *)theMat;
      else 
	return TCL_ERROR;
    }

	  else if ((strcmp(argv[1],"MaterialCMM") == 0)){

      void *theMat = OPS_MaterialCMM();
      if (theMat != 0) 
	theMaterial = (NDMaterial *)theMat;
      else 
	return TCL_ERROR;
    }

    else if ((strcmp(argv[1],"CycLiqCP") == 0)){

      void *theMat = OPS_CycLiqCPMaterial();
      if (theMat != 0) 
        theMaterial = (NDMaterial *)theMat;
      else 
        return TCL_ERROR;
    }

    else if ((strcmp(argv[1],"CycLiqCPSP") == 0)){

      void *theMat = OPS_CycLiqCPSPMaterial();
      if (theMat != 0) 
        theMaterial = (NDMaterial *)theMat;
      else 
        return TCL_ERROR;
    }

    else if ((strcmp(argv[1],"BoundingCamClay") == 0)){

      void *theMat = OPS_BoundingCamClayMaterial();
      if (theMat != 0) 
	theMaterial = (NDMaterial *)theMat;
      else 
	return TCL_ERROR;
    }

    else if ((strcmp(argv[1],"ManzariDafalias") == 0)){

      void *theMat = OPS_ManzariDafaliasMaterial();
      if (theMat != 0) 
	theMaterial = (NDMaterial *)theMat;
      else 
	return TCL_ERROR;
    }

    else if ((strcmp(argv[1],"ManzariDafaliasRO") == 0)){

      void *theMat = OPS_ManzariDafaliasMaterialRO();
      if (theMat != 0) 
	theMaterial = (NDMaterial *)theMat;
      else 
	return TCL_ERROR;
    }	

    else if ((strcmp(argv[1],"PM4Sand") == 0)){

      void *theMat = OPS_PM4SandMaterial();
      if (theMat != 0) 
	theMaterial = (NDMaterial *)theMat;
      else 
	return TCL_ERROR;
    }

	else if ((strcmp(argv[1], "J2CyclicBoundingSurface") == 0)) {

		void *theMat = OPS_J2CyclicBoundingSurfaceMaterial();
		if (theMat != 0)
			theMaterial = (NDMaterial *)theMat;
		else
			return TCL_ERROR;
	}
	
	else if ((strcmp(argv[1], "PM4Silt") == 0)) {

		void *theMat = OPS_PM4SiltMaterial();
		if (theMat != 0)
			theMaterial = (NDMaterial *)theMat;
		else
			return TCL_ERROR;
	}

    else if ((strcmp(argv[1],"ContactMaterial2D") == 0)){

      void *theMat = OPS_ContactMaterial2DMaterial();
      if (theMat != 0)
    theMaterial = (NDMaterial *)theMat;
      else
    return TCL_ERROR;
    }

    else if ((strcmp(argv[1],"ContactMaterial3D") == 0)){

      void *theMat = OPS_ContactMaterial3DMaterial();
      if (theMat != 0)
	theMaterial = (NDMaterial *)theMat;
      else
	return TCL_ERROR;
    }

    else if ((strcmp(argv[1],"InitialStateAnalysisWrapper") == 0)){
      
      void *theMat = OPS_InitialStateAnalysisWrapperMaterial();
      if (theMat != 0)
	theMaterial = (NDMaterial *)theMat;
      else
	return TCL_ERROR;
    }

#if !_DLL
    else if ((strcmp(argv[1],"stressDensity") == 0) || (strcmp(argv[1],"StressDensity") == 0)) {
      
      void *theMat = OPS_StressDensityMaterial();
      if (theMat != 0)
	theMaterial = (NDMaterial *)theMat;
      else
	return TCL_ERROR;
    }
#endif
    else if ((strcmp(argv[1],"ElasticIsotropic3D") == 0)) {

      void *theMat = OPS_ElasticIsotropic3D();
      if (theMat != 0)
	theMaterial = (NDMaterial *)theMat;
      else
	return TCL_ERROR;
    }

    else if ((strcmp(argv[1],"ElasticIsotropic") == 0)) {

      void *theMat = OPS_ElasticIsotropicMaterial();
      if (theMat != 0)
	theMaterial = (NDMaterial *)theMat;
      else
	return TCL_ERROR;
    }

    else if ((strcmp(argv[1],"ElasticOrthotropic3D") == 0) || (strcmp(argv[1],"ElasticOrthotropic") == 0)) {

      void *theMat = OPS_ElasticOrthotropicMaterial();
      if (theMat != 0)
	theMaterial = (NDMaterial *)theMat;
      else
	return TCL_ERROR;
    }

    else if (strcmp(argv[1],"PressureDependentElastic3D") == 0) {
	if (argc < 6) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: nDMaterial PressureDependentElastic3D tag? E? v? rho?" << endln;
	    return TCL_ERROR;
	}

	int tag = 0;
	double E = 0.0;
	double v = 0.0;
	double rho = 0.0;
	double expp = 0.0;
	double prp = 0.0;
	double pop = 0.0;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid PressureDependentElastic3D tag" << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
	    opserr << "WARNING invalid E\n";
	    opserr << "nDMaterial PressureDependentElastic3D: E" << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[4], &v) != TCL_OK) {
	    opserr << "WARNING invalid v\n";
	    opserr << "nDMaterial PressureDependentElastic3D: v" << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[5], &rho) != TCL_OK) {
	    opserr << "WARNING invalid v\n";
	    opserr << "nDMaterial PressureDependentElastic3D: rho" << tag << endln;
	    return TCL_ERROR;
	}

//////////////////////////////////////////////////////////////////////////////////
	if( argc == 6 )
	{
     	   theMaterial = new PressureDependentElastic3D (tag, E, v, rho);
	   //opserr << "nDMaterial PressureDependentElastic3D: expp =" << expp << endln;
	}
//////////////////////////////////////////////////////////////////////////////////
	else if( argc == 7 )
	{
	   //get the exponent of the pressure sensitive elastic material)
	   if (Tcl_GetDouble(interp, argv[6], &expp) != TCL_OK) {
	       opserr << "WARNING invalid v\n";
	       opserr << "nDMaterial PressureDependentElastic3D: " << tag << endln;
	       return TCL_ERROR;
	   }
     	   theMaterial = new PressureDependentElastic3D (tag, E, v, rho, expp);
	   //opserr << "nDMaterial PressureDependentElastic3D: expp =" << expp << endln;
	}
//////////////////////////////////////////////////////////////////////////////////
	else if (argc == 8 )
	{
	   //get the exponent pressure of the pressure sensitive elastic material)
	   if (Tcl_GetDouble(interp, argv[6], &expp) != TCL_OK) {
	       opserr << "WARNING invalid v\n";
	       opserr << "nDMaterial PressureDependentElastic3D: expp" << tag << endln;
	       return TCL_ERROR;
	   }
	   //get the reference pressure of the pressure sensitive elastic material)
	   if (Tcl_GetDouble(interp, argv[7], &prp) != TCL_OK) {
	       opserr << "WARNING invalid v\n";
	       opserr << "nDMaterial PressureDependentElastic3D: prp " << tag << endln;
	       return TCL_ERROR;
	   }
	   //opserr << "nDMaterial ElasticIsotropic3D: prp =" << prp << endln;
     	   theMaterial = new PressureDependentElastic3D (tag, E, v, rho, expp, prp);
	}
//////////////////////////////////////////////////////////////////////////////////
	else if (argc >= 9 )
	{
	   //get the exponent of the pressure sensitive elastic material)
	   if (Tcl_GetDouble(interp, argv[6], &expp) != TCL_OK) {
	       opserr << "WARNING invalid v\n";
	       opserr << "nDMaterial PressureDependentElastic3D: expp" << tag << endln;
	       return TCL_ERROR;
	   }
	   //get the reference pressure of the pressure sensitive elastic material)
	   if (Tcl_GetDouble(interp, argv[7], &prp) != TCL_OK) {
	       opserr << "WARNING invalid v\n";
	       opserr << "nDMaterial PressureDependentElastic3D: prp" << tag << endln;
	       return TCL_ERROR;
	   }
	   //get the cutoff pressure po of the pressure sensitive elastic material)
	   if (Tcl_GetDouble(interp, argv[8], &pop) != TCL_OK) {
	       opserr << "WARNING invalid v\n";
	       opserr << "nDMaterial PressureDependentElastic3D: pop" << tag << endln;
	       return TCL_ERROR;
	   }
	   //opserr << "nDMaterial PressureDependentElastic3D: pop =" << pop << endln;
     	   theMaterial = new PressureDependentElastic3D (tag, E, v, rho, expp, prp, pop);
	}

    }

    // Check argv[1] for J2PlaneStrain material type
    else if ((strcmp(argv[1],"J2Plasticity") == 0)  ||
	     (strcmp(argv[1],"J2") == 0)) {

      void *theMat = OPS_J2Plasticity();
      if (theMat != 0) 
	theMaterial = (NDMaterial *)theMat;
      else 
	return TCL_ERROR;

    }

	/////////////////////////////////////////////////////////////////
/*
   nDmaterial PlaneStressJ2  $matTag  $G  $K  $sig0  $H_kin  $H_iso 


     PlaneStress (int tag, 
				 int nd,
				 NDMaterial &the3DMaterial);

*/


    else if ((strcmp(argv[1],"PlaneStressSimplifiedJ2") == 0)) {
	if (argc < 8) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: nDmaterial Simplified3DJ2  $matTag  $G  $K  $sig0  $H_kin  $H_iso" << endln;
	    return TCL_ERROR;
	}

	int tag;
	double K, G, sig0, H_kin, H_iso;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid SimplifiedJ2 tag" << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[3], &G) != TCL_OK) {
	    opserr << "WARNING invalid G\n";
	    opserr << "nDMaterial SimplifiedJ2: " << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[4], &K) != TCL_OK) {
	    opserr << "WARNING invalid K\n";
	    opserr << "nDMaterial SimplifiedJ2: " << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[5], &sig0) != TCL_OK) {
	    opserr << "WARNING invalid sig0\n";
	    opserr << "nDMaterial SimplifiedJ2: " << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[6], &H_kin) != TCL_OK) {
	    opserr << "WARNING invalid H_kin\n";
	    opserr << "nDMaterial SimplifiedJ2: " << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[7], &H_iso) != TCL_OK) {
	    opserr << "WARNING invalid H_iso\n";
	    opserr << "nDMaterial SimplifiedJ2: " << tag << endln;
	    return TCL_ERROR;
	}
	


	   NDMaterial *theMaterial2 =  new SimplifiedJ2 (tag, 
				 3,
				 G,   
				 K,
				 sig0,
				 H_kin,
				 H_iso);


	theMaterial = new PlaneStressSimplifiedJ2 (tag, 
				 2,
				 *theMaterial2);

//	delete theMaterial2;


    }
/////////////////////////////////////////////////////////////////

    //
    //  MultiAxialCyclicPlasticity Model   by Gang Wang
    //
    //  nDMaterial MultiaxialCyclicPlasticity $tag, $rho, $K, $G,
    //      $Su , $Ho , $h, $m, $beta, $KCoeff
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // Check argv[1] for MultiaxialCyclicPlasticity material type
    else if ((strcmp(argv[1],"MultiaxialCyclicPlasticity") == 0)  ||
	     (strcmp(argv[1],"MCP") == 0)) {
      if (argc < 12) {
	opserr << "WARNING insufficient arguments\n";
	printCommand(argc,argv);
	opserr << "Want: nDMaterial MultiaxialCyclicPlasticity tag? rho? K? G? Su? Ho? h? m? beta? KCoeff? <eta?>" << endln;
	return TCL_ERROR;
      }

      int tag;
      double K, G, rho, Su, Ho, h, m, beta, Kcoeff;
      double eta = 0.0;

      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING invalid MultiaxialCyclicPlasticity tag" << endln;
	return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[3], &rho) != TCL_OK) {
	opserr << "WARNING invalid rho\n";
	opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << endln;
	return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[4], &K) != TCL_OK) {
	opserr << "WARNING invalid K\n";
	opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << endln;
	return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[5], &G) != TCL_OK) {
	opserr << "WARNING invalid G\n";
	opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << endln;
	return TCL_ERROR;
      }


      if (Tcl_GetDouble(interp, argv[6], &Su) != TCL_OK) {
	opserr << "WARNING invalid alpha1\n";
	opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << endln;
	return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[7], &Ho) != TCL_OK) {
	opserr << "WARNING invalid Ho\n";
	opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << endln;
	return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[8], &h) != TCL_OK) {
	opserr << "WARNING invalid h\n";
	opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << endln;
	return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[9], &m) != TCL_OK) {
	opserr << "WARNING invalid m\n";
	opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << endln;
	return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[10], &beta) != TCL_OK) {
	opserr << "WARNING invalid beta\n";
	opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[11], &Kcoeff) != TCL_OK) {
	opserr << "WARNING invalid Kcoeff\n";
	opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << endln;
	return TCL_ERROR;
      }


      if (argc > 12 && Tcl_GetDouble(interp, argv[12], &eta) != TCL_OK) {
	opserr << "WARNING invalid eta\n";
	opserr << "nDMaterial MultiaxialCyclicPlasticity: " << tag << endln;
	return TCL_ERROR;
      }

      theMaterial = new MultiaxialCyclicPlasticity (tag, 0, rho, K, G, Su, Ho, h,m,
						    beta, Kcoeff, eta);
    }


    // Pressure Independent Multi-yield, by ZHY
    else if (strcmp(argv[1],"PressureIndependMultiYield") == 0) {
	const int numParam = 6;
	const int totParam = 10;
	int tag;  double param[totParam];
	param[6] = 0.0;
	param[7] = 100.;
	param[8] = 0.0;
	param[9] = 20;

	char * arg[] = {"nd", "rho", "refShearModul", "refBulkModul",
			"cohesi", "peakShearStra",
			"frictionAng (=0)", "refPress (=100)", "pressDependCoe (=0.0)",
	    "numberOfYieldSurf (=20)"};
	if (argc < (3+numParam)) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: nDMaterial PressureIndependMultiYield tag? " << arg[0];
	    opserr << "? "<< "\n";
	    opserr << arg[1] << "? "<< arg[2] << "? "<< arg[3] << "? "<< "\n";
	    opserr << arg[4] << "? "<< arg[5] << "? "<< arg[6] << "? "<< "\n";
	    opserr << arg[7] << "? "<< arg[8] << "? "<< arg[9] << "? "<<endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid PressureIndependMultiYield tag" << endln;
	    return TCL_ERROR;
	}

	for (int i=3; (i<argc && i<13); i++)
	    if (Tcl_GetDouble(interp, argv[i], &param[i-3]) != TCL_OK) {
		    opserr << "WARNING invalid " << arg[i-3] << "\n";
		    opserr << "nDMaterial PressureIndependMultiYield: " << tag << endln;
		    return TCL_ERROR;
	    }

	static double * gredu = 0;
	// user defined yield surfaces
	if (param[9] < 0 && param[9] > -40) {
     param[9] = -int(param[9]);
     gredu = new double[int(2*param[9])];
		 for (int i=0; i<2*param[9]; i++)
	      if (Tcl_GetDouble(interp, argv[i+13], &gredu[i]) != TCL_OK) {
		      opserr << "WARNING invalid " << arg[i-3] << "\n";
		      opserr << "nDMaterial PressureIndependMultiYield: " << tag << endln;
		      return TCL_ERROR;
				}
  }

	PressureIndependMultiYield * temp =
	    new PressureIndependMultiYield (tag, param[0], param[1], param[2],
					    param[3], param[4], param[5], param[6],
					    param[7], param[8], param[9], gredu);
	theMaterial = temp;

	   if (gredu != 0) {
		   delete [] gredu;
		   gredu = 0;
	   }
    }

    // Pressure Independent Multi-yield, by Quan Gu
    else if (strcmp(argv[1],"MultiYieldSurfaceClay") == 0) {
		const int numParam = 6;
		const int totParam = 10;
		int tag;  double param[totParam];
		param[6] = 0.0;
		param[7] = 100.;
		param[8] = 0.0;
		param[9] = 20;

		char * arg[] = {"nd", "rho", "refShearModul", "refBulkModul",
				"cohesi", "peakShearStra",
				"frictionAng (=0)", "refPress (=100)", "pressDependCoe (=0.0)",
			"numberOfYieldSurf (=20)"};
		if (argc < (3+numParam)) {
			opserr << "WARNING insufficient arguments\n";
			printCommand(argc,argv);
			opserr << "Want: nDMaterial MultiYieldSurfaceClay tag? " << arg[0];
			opserr << "? "<< "\n";
			opserr << arg[1] << "? "<< arg[2] << "? "<< arg[3] << "? "<< "\n";
			opserr << arg[4] << "? "<< arg[5] << "? "<< arg[6] << "? "<< "\n";
			opserr << arg[7] << "? "<< arg[8] << "? "<< arg[9] << "? "<<endln;
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
			opserr << "WARNING invalid MultiYieldSurfaceClay tag" << endln;
			return TCL_ERROR;
		}

		for (int i=3; (i<argc && i<13); i++)
			if (Tcl_GetDouble(interp, argv[i], &param[i-3]) != TCL_OK) {
				opserr << "WARNING invalid " << arg[i-3] << "\n";
				opserr << "nDMaterial MultiYieldSurfaceClay: " << tag << endln;
				return TCL_ERROR;
			}

		static double * gredu = 0;
		// user defined yield surfaces
		if (param[9] < 0 && param[9] > -40) {
		 param[9] = -int(param[9]);
		 gredu = new double[int(2*param[9])];
			 for (int i=0; i<2*param[9]; i++)
			  if (Tcl_GetDouble(interp, argv[i+13], &gredu[i]) != TCL_OK) {
				  opserr << "WARNING invalid " << arg[i-3] << "\n";
				  opserr << "nDMaterial MultiYieldSurfaceClay: " << tag << endln;
				  return TCL_ERROR;
					}
	  }

		MultiYieldSurfaceClay * temp =
			new MultiYieldSurfaceClay (tag, param[0], param[1], param[2],
							param[3], param[4], param[5], param[6],
							param[7], param[8], param[9], gredu);
		theMaterial = temp;

	   if (gredu != 0) {
		   delete [] gredu;
		   gredu = 0;
	   }
    }
	// ============

    // Pressure Dependent Multi-yield, by ZHY
    else if (strcmp(argv[1],"PressureDependMultiYield") == 0) {
	const int numParam = 15;
	const int totParam = 24;
	int tag;
	double param[totParam];
 	param[15] = 20;
 	param[16] = 0.6;
	param[17] = 0.9;
	param[18] = 0.02;
	param[19] = 0.7;
	param[20] = 101.;
	param[21] = .3;
	param[22] = 0.;
	param[23] = 1.;

	char * arg[] = {"nd", "rho", "refShearModul",
		  "refBulkModul", "frictionAng",
			"peakShearStra", "refPress", "pressDependCoe",
			"phaseTransformAngle", "contractionParam1",
			"dilationParam1", "dilationParam2",
			"liquefactionParam1", "liquefactionParam2",
			"liquefactionParam4", "numberOfYieldSurf (=20)",
			"e (=0.6)", "volLimit1 (=0.9)", "volLimit2 (=0.02)",
			"volLimit3 (=0.7)", "Atmospheric pressure (=101)", "cohesi (=.5)",
	        "Hv (=0)", "Pv (=1.)" };
	if (argc < (3+numParam)) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: nDMaterial PressureDependMultiYield tag? "<< arg[0];
	    opserr << "? "<< "\n";
	    opserr << arg[1] << "? "<< arg[2] << "? "<< arg[3] << "? "<< "\n";
	    opserr << arg[4] << "? "<< arg[5] << "? "<< arg[6] << "? "<< "\n";
	    opserr << arg[7] << "? "<< arg[8] << "? "<< arg[9] << "? "<< "\n";
	    opserr << arg[10] << "? "<< arg[11] << "? "<< arg[12] << "? "<< "\n";
	    opserr << arg[13] << "? "<< arg[14] << "? "<< arg[15] << "? "<< "\n";
	    opserr << arg[16] << "? "<< arg[17] << "? "<< arg[18] << "? "<< "\n";
	    opserr << arg[19] << "? "<< arg[20] << "? "<< arg[21] << "? "<< endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid PressureDependMultiYield tag" << endln;
	    return TCL_ERROR;
	}

	for (int i=3; (i<argc && i<19); i++)
	  if (Tcl_GetDouble(interp, argv[i], &param[i-3]) != TCL_OK) {
		    opserr << "WARNING invalid " << arg[i-3] << "\n";
		    opserr << "nDMaterial PressureDependMultiYield: " << tag << endln;
		    return TCL_ERROR;
	  }

	static double * gredu = 0;
	// user defined yield surfaces
	if (param[15] < 0 && param[15] > -40) {
     param[15] = -int(param[15]);
     gredu = new double[int(2*param[15])];

		 for (int i=0; i<2*param[15]; i++)
	      if (Tcl_GetDouble(interp, argv[i+19], &gredu[i]) != TCL_OK) {
		      opserr << "WARNING invalid " << arg[i-3] << "\n";
		      opserr << "nDMaterial PressureIndependMultiYield: " << tag << endln;
		      return TCL_ERROR;
		  }
	}

	if (gredu != 0) {
	  for (int i=19+int(2*param[15]); i<argc; i++)
	    if (Tcl_GetDouble(interp, argv[i], &param[i-3-int(2*param[15])]) != TCL_OK) {
		      opserr << "WARNING invalid " << arg[i-3-int(2*param[15])] << "\n";
		      opserr << "nDMaterial PressureDependMultiYield: " << tag << endln;
		      return TCL_ERROR;
			}
	} else {
	  for (int i=19; i<argc; i++)
	    if (Tcl_GetDouble(interp, argv[i], &param[i-3]) != TCL_OK) {
		      opserr << "WARNING invalid " << arg[i-3-int(2*param[15])] << "\n";
		      opserr << "nDMaterial PressureDependMultiYield: " << tag << endln;
		      return TCL_ERROR;
		}
	}

	PressureDependMultiYield * temp =
	    new PressureDependMultiYield (tag, param[0], param[1], param[2],
					  param[3], param[4], param[5],
					  param[6], param[7], param[8],
					  param[9], param[10], param[11],
					  param[12], param[13], param[14],
					  param[15], gredu, param[16], param[17],
					  param[18], param[19], param[20], param[21], param[22], param[23]);

	   theMaterial = temp;
	   if (gredu != 0) {
		   delete [] gredu;
		   gredu = 0;
	   }
	}

    // Pressure Dependent Multi-yield, by ZHY
    else if (strcmp(argv[1],"PressureDependMultiYield02") == 0) {
	const int numParam = 13;
	const int totParam = 26;
	int tag;
	double param[totParam];
	param[numParam] = 20;
 	param[numParam+1] = 5.0;
 	param[numParam+2] = 3.;
 	param[numParam+3] = 1.;
	param[numParam+4] = 0.;
 	param[numParam+5] = 0.6;
	param[numParam+6] = 0.9;
	param[numParam+7] = 0.02;
	param[numParam+8] = 0.7;
	param[numParam+9] = 101.;
	param[numParam+10] = 0.1;
	param[numParam+11] = 0.;
	param[numParam+12] = 1.;

	char * arg[] = {"nd", "rho", "refShearModul",
		  "refBulkModul", "frictionAng",
			"peakShearStra", "refPress", "pressDependCoe",
			"phaseTransformAngle", "contractionParam1",
			"contractionParam3","dilationParam1","dilationParam3",
			"numberOfYieldSurf (=20)",
			"contractionParam2=5.0", "dilationParam2=3.0",
			"liquefactionParam1=1.0", "liquefactionParam2=0.0",
			"e (=0.6)", "volLimit1 (=0.9)", "volLimit2 (=0.02)",
			"volLimit3 (=0.7)", "Atmospheric pressure (=101)", "cohesi (=.1)",
	        "Hv (=0)", "Pv (=1.)" };
	if (argc < (3+numParam)) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: nDMaterial PressureDependMultiYield02 tag? "<< arg[0];
	    opserr << "? "<< "\n";
	    opserr << arg[1] << "? "<< arg[2] << "? "<< arg[3] << "? "<< "\n";
	    opserr << arg[4] << "? "<< arg[5] << "? "<< arg[6] << "? "<< "\n";
	    opserr << arg[7] << "? "<< arg[8] << "? "<< arg[9] << "? "<< "\n";
	    opserr << arg[10] << "? "<< arg[11] << "? "<< arg[12] << "? "<< "\n";
	    opserr << arg[13] << "? "<< arg[14] << "? "<< arg[15] << "? "<< "\n";
	    opserr << arg[16] << "? "<< arg[17] << "? "<< arg[18] << "? "<< "\n";
	    opserr << arg[19] << "? "<< arg[20] << "? "<< arg[21] << "? "<< "\n";
		opserr << arg[22] << "? "<< arg[23] << "? " << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid PressureDependMultiYield02 tag" << endln;
	    return TCL_ERROR;
	}

	int in = 17;
	for (int i=3; (i<argc && i<in); i++)
	  if (Tcl_GetDouble(interp, argv[i], &param[i-3]) != TCL_OK) {
		    opserr << "WARNING invalid " << arg[i-3] << "\n";
		    opserr << "nDMaterial PressureDependMultiYield02: " << tag << endln;
		    return TCL_ERROR;
	  }

	static double * gredu = 0;

	// user defined yield surfaces
	if (param[numParam] < 0 && param[numParam] > -100) {
     param[numParam] = -int(param[numParam]);
     gredu = new double[int(2*param[numParam])];

		 for (int i=0; i<2*param[numParam]; i++)
	      if (Tcl_GetDouble(interp, argv[i+in], &gredu[i]) != TCL_OK) {
		      opserr << "WARNING invalid " << arg[i-3] << "\n";
		      opserr << "nDMaterial PressureIndependMultiYield: " << tag << endln;
		      return TCL_ERROR;
		  }
	}

	if (gredu != 0) {
	  for (int i=in+int(2*param[numParam]); i<argc; i++)
	    if (Tcl_GetDouble(interp, argv[i], &param[i-3-int(2*param[numParam])]) != TCL_OK) {
		      opserr << "WARNING invalid " << arg[i-3-int(2*param[numParam])] << "\n";
		      opserr << "nDMaterial PressureDependMultiYield02: " << tag << endln;
		      return TCL_ERROR;
			}
	} else {
	  for (int i=in; i<argc; i++)
	    if (Tcl_GetDouble(interp, argv[i], &param[i-3]) != TCL_OK) {
		      opserr << "WARNING invalid " << arg[i-3-int(2*param[numParam])] << "\n";
		      opserr << "nDMaterial PressureDependMultiYield02: " << tag << endln;
		      return TCL_ERROR;
		}
	}


	PressureDependMultiYield02 * temp =
	    new PressureDependMultiYield02 (tag, param[0], param[1], param[2],
					  param[3], param[4], param[5],
					  param[6], param[7], param[8],
					  param[9], param[10], param[11],
					  param[12], param[13], gredu, param[14],
					  param[15], param[16], param[17],
					  param[18], param[19], param[20], param[21],
					  param[22], param[23], param[24], param[25]);

	   theMaterial = temp;
	   if (gredu != 0) {
		   delete [] gredu;
		   gredu = 0;
	   }
  }

  // nDMaterial PressureDependMultiYield03  $tag  $nd  $rho  $refShearModul  $refBulkModul  
  // $frictionAng  $peakShearStra  $refPress  $pressDependCoe  $PTAng  
  // $mType $ca  $cb $cc $cd $ce $da $db $dc <$noYieldSurf=20 
  // <$r1 $Gs1 …>  $liquefac1=1. $liquefac2=0. $pa=101 <$c=1.73>>
  // PressureDependMultiYield03 (based on PressureDependMultiYield02). 
	else if (strcmp(argv[1], "PressureDependMultiYield03") == 0) {
		const int numParam = 18; 
		const int totParam = 23; 
		int tag;
		double param[totParam];
		param[numParam] = 20;
		param[numParam + 1] = 1.;
		param[numParam + 2] = 0.;
		param[numParam + 3] = 101.;
		param[numParam + 4] = 1.73;

		char * arg[] = { "nd", "rho", "refShearModul","refBulkModul", "frictionAng",
			"peakShearStra", "refPress", "pressDependCoe", "phaseTransformAngle", 
			"mType","ca", "cb", "cc", "cd", "ce", "da", "db", "dc",
			"numberOfYieldSurf (=20)", "liquefactionParam1=1.0", "liquefactionParam2=0.0",
			"Atmospheric pressure (=101)", "cohesi (=1.73)" };

		if (argc < (3 + numParam)) { // 3 refers to "nDMaterial PressureDependMultiYield03  $tag"
			opserr << "WARNING insufficient arguments\n";
			printCommand(argc, argv);
			opserr << "Want: nDMaterial PressureDependMultiYield03 tag? " << arg[0];
			opserr << "? " << "\n";
			opserr << arg[1] << "? " << arg[2] << "? " << arg[3] << "? " << "\n";
			opserr << arg[4] << "? " << arg[5] << "? " << arg[6] << "? " << "\n";
			opserr << arg[7] << "? " << arg[8] << "? " << arg[9] << "? " << "\n";
			opserr << arg[10] << "? " << arg[11] << "? " << arg[12] << "? " << "\n";
			opserr << arg[13] << "? " << arg[14] << "? " << arg[15] << "? " << "\n";
			opserr << arg[16] << "? " << arg[17] << "? " << arg[18] << "? " << "\n";
			opserr << arg[19] << "? " << arg[20] << "? " << arg[21] << "? " << arg[22] << "? " << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
			opserr << "WARNING invalid PressureDependMultiYield03 tag" << endln;
			return TCL_ERROR;
		}

		int in = 22;
		for (int i = 3; (i<argc && i<in); i++)
			if (Tcl_GetDouble(interp, argv[i], &param[i - 3]) != TCL_OK) {
				opserr << "WARNING invalid " << arg[i - 3] << "\n";
				opserr << "nDMaterial PressureDependMultiYield03: " << tag << endln;
				return TCL_ERROR;
			}

		static double * gredu = 0;

		// user defined yield surfaces
		if (param[numParam] < 0 && param[numParam] > -100) {
			param[numParam] = -int(param[numParam]);
			gredu = new double[int(2 * param[numParam])];

			for (int i = 0; i<2 * param[numParam]; i++)
				if (Tcl_GetDouble(interp, argv[i + in], &gredu[i]) != TCL_OK) {
					opserr << "WARNING invalid " << arg[i - 3] << "\n";
					opserr << "nDMaterial PressureDependMultiYield03: " << tag << endln;
					return TCL_ERROR;
				}
		}

		if (gredu != 0) {
			for (int i = in + int(2 * param[numParam]); i<argc; i++)
				if (Tcl_GetDouble(interp, argv[i], &param[i - 3 - int(2 * param[numParam])]) != TCL_OK) {
					opserr << "WARNING invalid " << arg[i - 3 - int(2 * param[numParam])] << "\n";
					opserr << "nDMaterial PressureDependMultiYield03: " << tag << endln;
					return TCL_ERROR;
				}
		}
		else {
			for (int i = in; i<argc; i++)
				if (Tcl_GetDouble(interp, argv[i], &param[i - 3]) != TCL_OK) {
					opserr << "WARNING invalid " << arg[i - 3 - int(2 * param[numParam])] << "\n";
					opserr << "nDMaterial PressureDependMultiYield03: " << tag << endln;
					return TCL_ERROR;
				}
		}


		PressureDependMultiYield03 * temp =
			new PressureDependMultiYield03(tag, param[0], param[1], param[2],
				param[3], param[4], param[5],
				param[6], param[7], param[8],
				param[9], param[10], param[11],
				param[12], param[13], param[14],
				param[15], param[16], param[17], param[18], gredu,
				param[19], param[20], param[21], param[22]);

		theMaterial = temp;
		if (gredu != 0) {
			delete[] gredu;
			gredu = 0;
		}
	}
	
    // Fluid Solid Porous, by ZHY
    else if (strcmp(argv[1],"FluidSolidPorous") == 0) {

	int tag;  double param[4];
	char * arg[] = {"nd", "soilMatTag", "combinedBulkModul", "Atmospheric pressure"};
	if (argc < 6) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: nDMaterial FluidSolidPorous tag? "<< arg[0];
	    opserr << "? "<< "\n";
	    opserr << arg[1] << "? "<< arg[2] << "? "<< endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid FluidSolidPorous tag" << endln;
	    return TCL_ERROR;
	}

	for (int i=3; i<6; i++)
	  if (Tcl_GetDouble(interp, argv[i], &param[i-3] ) != TCL_OK) {
	      opserr << "WARNING invalid " << arg[i-3] << "\n";
	      opserr << "nDMaterial FluidSolidPorous: " << tag << endln;
	      return TCL_ERROR;
	  }

	NDMaterial *soil = OPS_getNDMaterial(param[1]);
	if (soil == 0) {
	      opserr << "WARNING FluidSolidPorous: couldn't get soil material ";
	      opserr << "tagged: " << param[1] << "\n";
	      return TCL_ERROR;
	}

	param[3] = 101.;
	if (argc == 7) {
	  if (Tcl_GetDouble(interp, argv[6], &param[3] ) != TCL_OK) {
	      opserr << "WARNING invalid " << arg[3] << "\n";
	      opserr << "nDMaterial FluidSolidPorous: " << tag << endln;
	      return TCL_ERROR;
	  }
	}

	theMaterial = new FluidSolidPorousMaterial (tag, param[0], *soil,
						    param[2],param[3]);
    }

    else if (strcmp(argv[1],"PlaneStressMaterial") == 0 ||
 	     strcmp(argv[1],"PlaneStress") == 0) {
      if (argc < 4) {
	opserr << "WARNING insufficient arguments\n";
	printCommand(argc,argv);
	opserr << "Want: nDMaterial PlaneStress tag? matTag?" << endln;
	return TCL_ERROR;
      }
      
      int tag, matTag;
      
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING invalid nDMaterial PlaneStress tag" << endln;
	return TCL_ERROR;
 	}
      
      if (Tcl_GetInt (interp, argv[3], &matTag) != TCL_OK) {
	opserr << "WARNING invalid matTag" << endln;
	opserr << "PlaneStress: " << matTag << endln;
	return TCL_ERROR;
      }
      
      NDMaterial *threeDMaterial = OPS_getNDMaterial(matTag);
      if (threeDMaterial == 0) {
	opserr << "WARNING nD material does not exist\n";
	opserr << "nD material: " << matTag;
	opserr << "\nPlaneStress nDMaterial: " << tag << endln;
	return TCL_ERROR;
 	}
      
      theMaterial = new PlaneStressMaterial( tag, *threeDMaterial );
    }
    
    // PlaneStrainMaterial
    else if (strcmp(argv[1],"PlaneStrainMaterial") == 0 ||
	     strcmp(argv[1],"PlaneStrain") == 0) {
      if (argc < 4) {
	opserr << "WARNING insufficient arguments\n";
	printCommand(argc,argv);
	opserr << "Want: nDMaterial PlaneStrain tag? matTag?" << endln;
	return TCL_ERROR;
      }
      
      int tag, matTag;
       
       if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	 opserr << "WARNING invalid nDMaterial PlaneStrain tag" << endln;
	 return TCL_ERROR;
       }
       
       if (Tcl_GetInt (interp, argv[3], &matTag) != TCL_OK) {
	 opserr << "WARNING invalid matTag" << endln;
	 opserr << "PlaneStrain: " << matTag << endln;
	 return TCL_ERROR;
       }
       
       NDMaterial *threeDMaterial = OPS_getNDMaterial(matTag);
       if (threeDMaterial == 0) {
	 opserr << "WARNING nD material does not exist\n";
	 opserr << "nD material: " << matTag;
	 opserr << "\nPlaneStrain nDMaterial: " << tag << endln;
	 return TCL_ERROR;
       }
       
       theMaterial = new PlaneStrainMaterial( tag, *threeDMaterial );
     }
    
     else if (strcmp(argv[1],"PlateFiberMaterial") == 0 ||
	      strcmp(argv[1],"PlateFiber") == 0) {

       void *theMat = OPS_PlateFiberMaterial();
       if (theMat != 0) 
	 theMaterial = (NDMaterial *)theMat;
       else 
	 return TCL_ERROR;
     }

    // ----- Cap plasticity model ------    // Quan Gu & ZhiJian Qiu  2013

    // format nDmaterial CapPlasticity $tag $ndm $rho $G $K $X $D $W $R $lambda $theta $beta $alpha $T $tol
     else if (strcmp(argv[1],"CapPlasticity") == 0) {
       
       int tag;
       int ndm =3;
       double rho = 0.0;
       double G = 1.0e10;
       double K = 1.1e10;
       double X = 1.1032e8;
       double D = 4.6412e-10;
       double W = 0.42;
       double R = 4.43;
       double lambda = 7.9979e6;
       double theta = 0.11;
       double beta = 6.3816e-8;
       double alpha = 2.6614e7;
       double T = -2.0684e6; 
       double tol = 1.0e-10;
       
       if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	 opserr << "WARNING invalid CapPlasticity tag" << endln;
	 return TCL_ERROR;
       }
       
       if (Tcl_GetInt(interp, argv[3], &ndm) != TCL_OK) {
	 opserr << "WARNING invalid CapPlasticity nd" << endln;
	 return TCL_ERROR;
       }
       
       if (Tcl_GetDouble(interp, argv[4], &rho) != TCL_OK) {
	 opserr << "WARNING invalid CapPlasticity rho" << endln;
	 return TCL_ERROR;
       }
       
       if (Tcl_GetDouble(interp, argv[5], &G) != TCL_OK) {
	 opserr << "WARNING invalid CapPlasticity G" << endln;
	 return TCL_ERROR;
       }
       
       if (Tcl_GetDouble(interp, argv[6], &K) != TCL_OK) {
	 opserr << "WARNING invalid CapPlasticity K" << endln;
	 return TCL_ERROR;
       }
       
       if (argc > 7) {
	 
	 if (Tcl_GetDouble(interp, argv[7], &X) != TCL_OK) {
	   opserr << "WARNING invalid CapPlasticity X" << endln;
	   return TCL_ERROR;
	 }
	 
	 if (Tcl_GetDouble(interp, argv[8], &D) != TCL_OK) {
	   opserr << "WARNING invalid CapPlasticity D" << endln;
	   return TCL_ERROR;
	 }
	 
	 if (Tcl_GetDouble(interp, argv[9], &W) != TCL_OK) {
	   opserr << "WARNING invalid CapPlasticity W" << endln;
	   return TCL_ERROR;
	 }
	 
	 if (Tcl_GetDouble(interp, argv[10], &R) != TCL_OK) {
	   opserr << "WARNING invalid CapPlasticity R" << endln;
	   return TCL_ERROR;
	 }
	 
	 if (Tcl_GetDouble(interp, argv[11], &lambda) != TCL_OK) {
	   opserr << "WARNING invalid CapPlasticity lambda" << endln;
	   return TCL_ERROR;
	 }
	 
	 if (Tcl_GetDouble(interp, argv[12], &theta) != TCL_OK) {
	   opserr << "WARNING invalid CapPlasticity theta" << endln;
	   return TCL_ERROR;
	 }
	 
	 if (Tcl_GetDouble(interp, argv[13], &beta) != TCL_OK) {
	   opserr << "WARNING invalid CapPlasticity beta" << endln;
	   return TCL_ERROR;
	 }
	 
	 if (Tcl_GetDouble(interp, argv[14], &alpha) != TCL_OK) {
	   opserr << "WARNING invalid CapPlasticity alpha" << endln;
	   return TCL_ERROR;
	 }
	 
	 if (Tcl_GetDouble(interp, argv[15], &T) != TCL_OK) {
	   opserr << "WARNING invalid CapPlasticity T" << endln;
	   return TCL_ERROR;
	 }
	 
	 if (Tcl_GetDouble(interp, argv[16], &tol) != TCL_OK) {
	   opserr << "WARNING invalid CapPlasticity tol" << endln;
	   return TCL_ERROR;
	 }
	 
       } //end if
       
       theMaterial = 	new CapPlasticity(  tag,
					    G,
					    K,
					    rho,
					    X,
					    D,
					    W,
					    R,
					    lambda,
					    theta,
					    beta,
					    alpha,
					    T, 
					    ndm,
					    tol
					    ) ;
       
       
     }
    
    
/////////////////////////////////////////////////////////////////
/*
   nDmaterial Simplified3DJ2  $matTag  $G  $K  $sig0  $H_kin  $H_iso 


    SimplifiedJ2 (int tag, 
				 int nd,
				 double G,
				 double K,
				 double sigmaY0,
				 double H_kin,
				 double H_iso);

*/
 

    // Check argv[1] for J2PlaneStrain material type
    else if ((strcmp(argv[1],"Simplified3DJ2") == 0)  || (strcmp(argv[1],"3DJ2") == 0)) {
	if (argc < 8) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: nDmaterial Simplified3DJ2  $matTag  $G  $K  $sig0  $H_kin  $H_iso" << endln;
	    return TCL_ERROR;
	}

	int tag;
	double K, G, sig0, H_kin, H_iso;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid SimplifiedJ2 tag" << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[3], &G) != TCL_OK) {
	    opserr << "WARNING invalid G\n";
	    opserr << "nDMaterial SimplifiedJ2: " << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[4], &K) != TCL_OK) {
	    opserr << "WARNING invalid K\n";
	    opserr << "nDMaterial SimplifiedJ2: " << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[5], &sig0) != TCL_OK) {
	    opserr << "WARNING invalid sig0\n";
	    opserr << "nDMaterial SimplifiedJ2: " << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[6], &H_kin) != TCL_OK) {
	    opserr << "WARNING invalid H_kin\n";
	    opserr << "nDMaterial SimplifiedJ2: " << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[7], &H_iso) != TCL_OK) {
	    opserr << "WARNING invalid H_iso\n";
	    opserr << "nDMaterial SimplifiedJ2: " << tag << endln;
	    return TCL_ERROR;
	}
	
	
	theMaterial = new SimplifiedJ2 (tag, 
				 3,
				 G,
				 K,
				 sig0,
				 H_kin,
				 H_iso);
    }


     else if (strcmp(argv[1],"PlateRebarMaterial") == 0 ||
	      strcmp(argv[1],"PlateRebar") == 0) {
 	if (argc < 5) {
 	    opserr << "WARNING insufficient arguments\n";
 	    printCommand(argc,argv);
 	    opserr << "Want: nDMaterial PlateRebar tag? matTag? angle?" << endln;
 	    return TCL_ERROR;
 	}

 	int tag, matTag;
 	double angle;

 	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
 	    opserr << "WARNING invalid nDMaterial PlateRebar tag" << endln;
 	    return TCL_ERROR;
 	}

 	if (Tcl_GetInt (interp, argv[3], &matTag) != TCL_OK) {
 	    opserr << "WARNING invalid matTag" << endln;
 	    opserr << "PlateRebar: " << tag << endln;
 	    return TCL_ERROR;
 	}

 	UniaxialMaterial *theMat = OPS_getUniaxialMaterial(matTag);
 	if (theMat == 0) {
 	    opserr << "WARNING uniaxialmaterial does not exist\n";
 	    opserr << "UniaxialMaterial: " << matTag;
 	    opserr << "\nPlateRebar nDMaterial: " << tag << endln;
 	    return TCL_ERROR;
 	}

 	if (Tcl_GetDouble (interp, argv[4], &angle) != TCL_OK) {
 	    opserr << "WARNING invalid angle" << endln;
 	    opserr << "PlateRebar: " << tag << endln;
 	    return TCL_ERROR;
 	}

 	theMaterial = new PlateRebarMaterial( tag, *theMat, angle );
     }

    //start Yuli Huang & Xinzheng Lu PlateFromPlaneStressMaterial
     else if (strcmp(argv[1],"PlateFromPlaneStressMaterial") == 0 ||
	      strcmp(argv[1],"PlateFromPlaneStress") == 0) {
 	if (argc < 5) {
 	    opserr << "WARNING insufficient arguments\n";
 	    printCommand(argc,argv);
 	    opserr << "Want: nDMaterial PlateFromPlaneStress tag? matTag? gmod?" << endln;
 	    return TCL_ERROR;
 	}

 	int tag, matTag;
 	double gmod;

 	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
 	    opserr << "WARNING invalid nDMaterial PlateFromPlaneStress tag" << endln;
 	    return TCL_ERROR;
 	}

 	if (Tcl_GetInt (interp, argv[3], &matTag) != TCL_OK) {
 	    opserr << "WARNING invalid matTag" << endln;
 	    opserr << "PlateFromPlaneStress: " << tag << endln;
 	    return TCL_ERROR;
 	}

 	NDMaterial *theMat = OPS_getNDMaterial(matTag);
 	if (theMat == 0) {
 	    opserr << "WARNING ndMaterial does not exist\n";
 	    opserr << "ndMaterial: " << matTag;
 	    opserr << "\nPlateFromPlaneStress nDMaterial: " << tag << endln;
 	    return TCL_ERROR;
 	}

 	if (Tcl_GetDouble (interp, argv[4], &gmod) != TCL_OK) {
 	    opserr << "WARNING invalid gmod" << endln;
 	    opserr << "PlateFromPlaneStress: " << tag << endln;
 	    return TCL_ERROR;
 	}

 	theMaterial = new PlateFromPlaneStressMaterial( tag, *theMat, gmod );
     }


    //start Yuli Huang & Xinzheng Lu ConcreteS
     else if (strcmp(argv[1],"ConcreteS") == 0) {
 	if (argc < 8) {
 	    opserr << "WARNING insufficient arguments\n";
 	    printCommand(argc,argv);
 	    opserr << "Want: nDMaterial ConcreteS tag? E? nu? fc? ft? Es?" << endln;
 	    return TCL_ERROR;
 	}

 	int tag;
 	double E, nu, fc, ft, Es;

 	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
 	    opserr << "WARNING invalid nDMaterial ConcreteS tag" << endln;
 	    return TCL_ERROR;
 	}

 	if (Tcl_GetDouble (interp, argv[3], &E) != TCL_OK) {
 	    opserr << "WARNING invalid E" << endln;
 	    opserr << "ConcreteS: " << tag << endln;
 	    return TCL_ERROR;
 	}

 	if (Tcl_GetDouble (interp, argv[4], &nu) != TCL_OK) {
 	    opserr << "WARNING invalid nu" << endln;
 	    opserr << "ConcreteS: " << tag << endln;
 	    return TCL_ERROR;
 	}

 	if (Tcl_GetDouble (interp, argv[5], &fc) != TCL_OK) {
 	    opserr << "WARNING invalid fc" << endln;
 	    opserr << "ConcreteS: " << tag << endln;
 	    return TCL_ERROR;
 	}

 	if (Tcl_GetDouble (interp, argv[6], &ft) != TCL_OK) {
 	    opserr << "WARNING invalid ft" << endln;
 	    opserr << "ConcreteS: " << tag << endln;
 	    return TCL_ERROR;
 	}

 	if (Tcl_GetDouble (interp, argv[7], &Es) != TCL_OK) {
 	    opserr << "WARNING invalid Es" << endln;
 	    opserr << "ConcreteS: " << tag << endln;
 	    return TCL_ERROR;
 	}

 	theMaterial = new ConcreteS( tag, E, nu, fc, ft, Es);
     }
    //end Yuli Huang & Xinzheng Lu ConcreteS

    //start Yuli Huang & Xinzheng Lu PlaneStressUserMaterial
     else if (strcmp(argv[1],"PlaneStressUserMaterial") == 0) {
 	if (argc < 6) {
 	    opserr << "WARNING insufficient arguments\n";
 	    printCommand(argc,argv);
 	    opserr << "Want: nDMaterial PlaneStressUserMaterial tag? nstatevs? nprops? prop1? ... propn?" << endln;
 	    return TCL_ERROR;
 	}

 	int tag, nstatevs, nprops;
 	double *props, p;

 	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
 	    opserr << "WARNING invalid nDMaterial PlaneStressUserMaterial tag" << endln;
 	    return TCL_ERROR;
 	}

 	if (Tcl_GetInt(interp, argv[3], &nstatevs) != TCL_OK) {
 	    opserr << "WARNING invalid nDMaterial PlaneStressUserMaterial nstatevs" << endln;
 	    return TCL_ERROR;
 	}
 	
 	if (nstatevs < 1) nstatevs = 1;

 	if (Tcl_GetInt(interp, argv[4], &nprops) != TCL_OK) {
 	    opserr << "WARNING invalid nDMaterial PlaneStressUserMaterial nprops" << endln;
 	    return TCL_ERROR;
 	}
 	
 	if (nprops < 1) nprops = 1;

 	props = new double[nprops];
 	for (int i = 0; i < nprops; i++) {
 	  if (Tcl_GetDouble (interp, argv[5+i], &p) != TCL_OK) {
 	    opserr << "WARNING invalid prop" << endln;
 	    opserr << "PlaneStressUserMaterial: " << tag << endln;
 	    return TCL_ERROR;
 	  }
 	  props[i] = p;
 	}

 	theMaterial = new PlaneStressUserMaterial( tag, nstatevs, nprops, props);
	if (props != 0) delete props;
     }
    //end Yuli Huang & Xinzheng Lu PlaneStressUserMaterial

     else if (strcmp(argv[1],"BeamFiberMaterial") == 0 ||
 	     strcmp(argv[1],"BeamFiber") == 0) {

       void *theMat = OPS_BeamFiberMaterial();
       if (theMat != 0) 
	 theMaterial = (NDMaterial *)theMat;
       else 
	 return TCL_ERROR;
     }

     else if (strcmp(argv[1],"BeamFiberMaterial2d") == 0 ||
 	     strcmp(argv[1],"BeamFiber2d") == 0) {

       void *theMat = OPS_BeamFiberMaterial2d();
       if (theMat != 0) 
	 theMaterial = (NDMaterial *)theMat;
       else 
	 return TCL_ERROR;
     }

     else if (strcmp(argv[1],"BeamFiberMaterial2dPS") == 0 ||
 	     strcmp(argv[1],"BeamFiber2dPS") == 0) {

       void *theMat = OPS_BeamFiberMaterial2dPS();
       if (theMat != 0) 
	 theMaterial = (NDMaterial *)theMat;
       else 
	 return TCL_ERROR;
     }    

     else if (strcmp(argv[1],"ConcreteMcftNonLinear7") == 0 || strcmp(argv[1],"ConcreteMcftNonLinear5") == 0) {
       if (argc < 11) {
	 opserr << "WARNING insufficient arguments\n";
	 printCommand(argc,argv);
	 opserr << "Want: nDMaterial ConcreteMcftNonlinear7 tag? fcu? ecu? Ec? fcr? Esv? fyv? alphaV? RoV?" << endln;
	 return TCL_ERROR;
       }
       
       int tag = 0;
       double fcu = 0.0;
       double ecu = 0.0;
       double Ec = 0.0;
       double fcr = 0.0;
       double Esv = 0.0;
       double fyv = 0.0;
       double alphaV = 0.0;
       double RoV = 0.0;
       
       if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	 opserr << "WARNING invalid ConcreteMcftNonlinear7: tag" << endln;
	 return TCL_ERROR;
       }
       
       if (Tcl_GetDouble(interp, argv[3], &fcu) != TCL_OK) {
	 opserr << "WARNING invalid fcu\n";
	 opserr << "nDMaterial ConcreteMcftNonLinearNonLinear5: fcu" << tag << endln;
	 return TCL_ERROR;
       }
       
       if (Tcl_GetDouble(interp, argv[4], &ecu) != TCL_OK) {
	 opserr << "WARNING invalid ecu\n";
	 opserr << "nDMaterial ConcreteMcftNonLinearNonLinear5: ecu" << tag << endln;
	 return TCL_ERROR;
       }
       
       if (Tcl_GetDouble(interp, argv[5], &Ec) != TCL_OK) {
	 opserr << "WARNING invalid Ec\n";
	 opserr << "nDMaterial ConcreteMcftNonlinear7: Ec" << tag << endln;
	 return TCL_ERROR;
       }
       
       if (Tcl_GetDouble(interp, argv[6], &fcr) != TCL_OK) {
	 opserr << "WARNING invalid fcr\n";
	 opserr << "nDMaterial ConcreteMcftNonlinear7: fcr" << tag << endln;
	 return TCL_ERROR;
       }
       
       if (Tcl_GetDouble(interp, argv[7], &Esv) != TCL_OK) {
	 opserr << "WARNING invalid Esv\n";
	 opserr << "nDMaterial ConcreteMcftNonlinear7: Esv" << tag << endln;
	 return TCL_ERROR;
       }
       
       if (Tcl_GetDouble(interp, argv[8], &fyv) != TCL_OK) {
	 opserr << "WARNING invalid fyv\n";
	 opserr << "nDMaterial ConcreteMcftNonlinear7: fyv" << tag << endln;
	 return TCL_ERROR;
       }
       
       if (Tcl_GetDouble(interp, argv[9], &alphaV) != TCL_OK) {
	 opserr << "WARNING invalid alphaV\n";
	 opserr << "nDMaterial ConcreteMcftNonlinear7: alphaV" << tag << endln;
	 return TCL_ERROR;
       }
       
       if (Tcl_GetDouble(interp, argv[10], &RoV) != TCL_OK) {
	 opserr << "WARNING invalid RoV\n";
	 opserr << "nDMaterial ConcreteMcftNonlinear7: RoV" << tag << endln;
	 return TCL_ERROR;
       }
       
       if (strcmp(argv[1],"ConcreteMcftNonLinear7") == 0) 
	 theMaterial = new ConcreteMcftNonLinear7 (tag, fcu, ecu, Ec, fcr, Esv, fyv, alphaV, RoV);
       else 
	 theMaterial = new ConcreteMcftNonLinear5 (tag, fcu, ecu, Ec, fcr, Esv, fyv, alphaV, RoV);
     }

    else if (strcmp(argv[1],"Bidirectional") == 0) {
      opserr << "nDMaterial Bidirectional is now a section model, please "
	   << "change to \'section Bidirectional\'" << endln;
      return TCL_ERROR;
    }

	//-------nD materials for thermo-mechanical analysis---Added by L.Jiang[SIF]
	else if ((strcmp(argv[1], "DruckerPragerThermal") == 0)) {

		void *theMat = OPS_DruckerPragerMaterialThermal();
		if (theMat != 0)
			theMaterial = (NDMaterial *)theMat;
		else
			return TCL_ERROR;
	}
	//-------------------------------------------------------------
    /*
	else if ((strcmp(argv[1], "CDPPlaneStressThermal") == 0)) {
		void *theMat = OPS_PlasticDamageConcretePlaneStressThermal();
		if (theMat != 0)
			theMaterial = (NDMaterial *)theMat;
		else
			return TCL_ERROR;
	}
    */
	//-------------------------------------------------------------
	else if (strcmp(argv[1], "PlateFromPlaneStressThermal") == 0 ) {
		if (argc < 5) {
			opserr << "WARNING insufficient arguments\n";
			printCommand(argc, argv);
			opserr << "Want: nDMaterial PlateFromPlaneStress tag? matTag? gmod?" << endln;
			return TCL_ERROR;
		}

		int tag, matTag;
		double gmod;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
			opserr << "WARNING invalid nDMaterial PlateFromPlaneStress tag" << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
			opserr << "WARNING invalid matTag" << endln;
			opserr << "PlateFromPlaneStress: " << tag << endln;
			return TCL_ERROR;
		}

		NDMaterial *theMat = OPS_getNDMaterial(matTag);
		if (theMat == 0) {
			opserr << "WARNING ndMaterial does not exist\n";
			opserr << "ndMaterial: " << matTag;
			opserr << "\nPlateFromPlaneStress nDMaterial: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[4], &gmod) != TCL_OK) {
			opserr << "WARNING invalid gmod" << endln;
			opserr << "PlateFromPlaneStress: " << tag << endln;
			return TCL_ERROR;
		}

		theMaterial = new PlateFromPlaneStressMaterialThermal(tag, *theMat, gmod);
	}
	else if (strcmp(argv[1], "PlateRebarMaterialThermal") == 0 ||
		strcmp(argv[1], "PlateRebarThermal") == 0) {
		if (argc < 5) {
			opserr << "WARNING insufficient arguments\n";
			printCommand(argc, argv);
			opserr << "Want: nDMaterial PlateRebar tag? matTag? angle?" << endln;
			return TCL_ERROR;
		}

		int tag, matTag;
		double angle;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
			opserr << "WARNING invalid nDMaterial PlateRebar tag" << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
			opserr << "WARNING invalid matTag" << endln;
			opserr << "PlateRebar: " << tag << endln;
			return TCL_ERROR;
		}

		UniaxialMaterial *theMat = OPS_getUniaxialMaterial(matTag);
		if (theMat == 0) {
			opserr << "WARNING uniaxialmaterial does not exist\n";
			opserr << "UniaxialMaterial: " << matTag;
			opserr << "\nPlateRebar nDMaterial: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[4], &angle) != TCL_OK) {
			opserr << "WARNING invalid angle" << endln;
			opserr << "PlateRebar: " << tag << endln;
			return TCL_ERROR;
		}

		theMaterial = new PlateRebarMaterialThermal(tag, *theMat, angle);
	}
	else if ((strcmp(argv[1], "J2PlasticityThermal") == 0) ||
		(strcmp(argv[1], "J2Thermal") == 0)) {
		if (argc < 9) {
			opserr << "WARNING insufficient arguments\n";
			printCommand(argc, argv);
			opserr << "Want: nDMaterial J2PlasticityThermal tag? K? G? sig0? sigInf? delta? H? <eta?>" << endln;
			return TCL_ERROR;
		}

		int tag;
		double K, G, sig0, sigInf, delta, H;
		double eta = 0.0;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
			opserr << "WARNING invalid J2PlasticityThermal tag" << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[3], &K) != TCL_OK) {
			opserr << "WARNING invalid K\n";
			opserr << "nDMaterial J2PlasticityThermal: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[4], &G) != TCL_OK) {
			opserr << "WARNING invalid G\n";
			opserr << "nDMaterial J2PlasticityThermal: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[5], &sig0) != TCL_OK) {
			opserr << "WARNING invalid sig0\n";
			opserr << "nDMaterial J2PlasticityThermal: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[6], &sigInf) != TCL_OK) {
			opserr << "WARNING invalid sigInf\n";
			opserr << "nDMaterial J2PlasticityThermal: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[7], &delta) != TCL_OK) {
			opserr << "WARNING invalid delta\n";
			opserr << "nDMaterial J2PlasticityThermal: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[8], &H) != TCL_OK) {
			opserr << "WARNING invalid H\n";
			opserr << "nDMaterial J2PlasticityThermal: " << tag << endln;
			return TCL_ERROR;
		}
		if (argc > 9 && Tcl_GetDouble(interp, argv[9], &eta) != TCL_OK) {
			opserr << "WARNING invalid eta\n";
			opserr << "nDMaterial J2PlasticityThermal: " << tag << endln;
			return TCL_ERROR;
		}

		theMaterial = new J2PlasticityThermal(tag, 0, K, G, sig0, sigInf,
			delta, H, eta);
	}
	else if (strcmp(argv[1], "PlateFiberMaterialThermal") == 0 ||
		strcmp(argv[1], "PlateFiberThermal") == 0) {
		if (argc < 4) {
			opserr << "WARNING insufficient arguments\n";
			printCommand(argc, argv);
			opserr << "Want: nDMaterial PlateFiberThermal tag? matTag?" << endln;
			return TCL_ERROR;
		}

		int tag, matTag;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
			opserr << "WARNING invalid nDMaterial PlateFiberThermal tag" << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
			opserr << "WARNING invalid matTag" << endln;
			opserr << "PlateFiberThermal: " << matTag << endln;
			return TCL_ERROR;
		}

		NDMaterial *threeDMaterial = OPS_getNDMaterial(matTag);
		if (threeDMaterial == 0) {
			opserr << "WARNING nD material does not exist\n";
			opserr << "nD material: " << matTag;
			opserr << "\nPlateFiberThermal nDMaterial: " << tag << endln;
			return TCL_ERROR;
		}

		theMaterial = new PlateFiberMaterialThermal(tag, *threeDMaterial);
	}
	//--------End of adding PlateFiberMaterialThermal
	else if ( (strcmp(argv[1], "ElasticIsotropicThermal") == 0) || (strcmp(argv[1], "ElasticIsotropic3DThermal") == 0)) {

		void *theMat = OPS_ElasticIsotropicMaterialThermal();
		if (theMat != 0)
			theMaterial = (NDMaterial *)theMat;
		else
			return TCL_ERROR;
	}

	//end of adding thermo-mechanical nd materials-L.Jiang[SIF]
    else {
      theMaterial = TclModelBuilder_addFeapMaterial(clientData,
						    interp,
						    argc,
						    argv,
						    theTclBuilder);
    }

    if (theMaterial == 0) {
      //
      // maybe element in a class package already loaded
      //  loop through linked list of loaded functions comparing names & if find call it
      //
      
      NDMaterialPackageCommand *matCommands = theNDMaterialPackageCommands;
      bool found = false;
      while (matCommands != NULL && found == false) {
	if (strcmp(argv[1], matCommands->funcName) == 0) {
	  theMaterial = (NDMaterial *)(*(matCommands->funcPtr))();
	  found = true;;
	} else
	  matCommands = matCommands->next;
      }
    }

    //
    // check to see if element is a procedure
    //   the proc may already have been loaded from a package or may exist in a package yet to be loaded
    //
    if (theMaterial == 0) {

      // maybe material in a routine
      //

      char *matType = new char[strlen(argv[1])+1];
      strcpy(matType, argv[1]);
      matObj *matObject = OPS_GetMaterialType(matType, strlen(matType));
      
      delete [] matType;

      if (matObject != 0) {
	
	theMaterial = Tcl_addWrapperNDMaterial(matObject, clientData, interp,
						     argc, argv, theTclBuilder);
	
	if (theMaterial == 0)
	  delete matObject;
      }
    }


    //
    // maybe material class exists in a package yet to be loaded
    //

    if (theMaterial == 0) {
	
      void *libHandle;
      void * (*funcPtr)();
      
      int matNameLength = strlen(argv[1]);
      char *tclFuncName = new char[matNameLength+12];
      strcpy(tclFuncName, "OPS_");
      strcpy(&tclFuncName[4], argv[1]);    
      int res = getLibraryFunction(argv[1], tclFuncName, &libHandle, (void **)&funcPtr);
      
      delete [] tclFuncName;
      
      if (res == 0) {
	
	//
	// add loaded function to list of functions
	//
	
	char *matName = new char[matNameLength+1];
	strcpy(matName, argv[1]);
	NDMaterialPackageCommand *theMatCommand = new NDMaterialPackageCommand;
	theMatCommand->funcPtr = funcPtr;
	theMatCommand->funcName = matName;	
	theMatCommand->next = theNDMaterialPackageCommands;
	theNDMaterialPackageCommands = theMatCommand;
	
	theMaterial = (NDMaterial *)(*funcPtr)();
      }
    }


    if (theMaterial == 0) {
	opserr << "WARNING could not create nDMaterial: " << argv[1];
	return TCL_ERROR;
    }

    // Now add the material to the modelBuilder
    if (OPS_addNDMaterial(theMaterial) == false) {
	opserr << "WARNING could not add material to the domain\n";
	opserr << *theMaterial << endln;
	delete theMaterial; // invoke the material objects destructor, otherwise mem leak
	return TCL_ERROR;
    }

    return TCL_OK;
}


