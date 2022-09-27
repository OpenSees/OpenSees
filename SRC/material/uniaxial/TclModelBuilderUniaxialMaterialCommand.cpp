/* ****************************************************************** **
**    Opensee - Open System for Earthquake Engineering Simulation    **
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
                                                                        
// $Revision: 1.72 $
// $Date: 2010-09-16 00:04:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/TclModelBuilderUniaxialMaterialCommand.cpp,v $
                                                                        
                                                                        
// Written: fmk, MHS 
// Created: 07/99
//
// Description: This file contains the function invoked when the user invokes
// the uniaxialMaterial command in the interpreter. 
//
// What: "@(#) TclModelBuilderUniaxialMaterialCommand.C, revA"

//#include <TclModelBuilder.h>

#include <tcl.h>
#include <elementAPI.h>
extern "C" int OPS_ResetInputNoBuilder(ClientData clientData, Tcl_Interp * interp, int cArg, int mArg, TCL_Char * *argv, Domain * domain);


#include <Pinching4Material.h>   // NM
#include <ShearPanelMaterial.h>  // NM
#include <BarSlipMaterial.h>     // NM

#include <Vector.h>
#include <string.h>

extern void *OPS_SPSW02(void);		// SAJalali
extern void *OPS_TDConcreteEXP(void); // ntosic
extern void *OPS_TDConcrete(void); // ntosic
extern void *OPS_TDConcreteMC10(void); //ntosic
extern void *OPS_TDConcreteMC10NL(void); //ntosic
extern void *OPS_ECC01(void);
extern void *OPS_ElasticMaterial(void);
extern void *OPS_Elastic2Material(void);
extern void *OPS_ElasticPPMaterial(void);
extern void *OPS_ENTMaterial(void);
extern void *OPS_EPPGapMaterial(void);
extern void *OPS_ParallelMaterial(void);
extern void *OPS_SeriesMaterial(void);
extern void *OPS_PathIndependentMaterial(void);
extern void *OPS_BackboneMaterial(void);
extern void *OPS_FatigueMaterial(void);
extern void *OPS_HardeningMaterial(void);
extern void *OPS_UniaxialJ2Plasticity(void);
extern void *OPS_SmoothPSConcrete(void);
extern void *OPS_HystereticMaterial(void);
extern void *OPS_CableMaterial(void);
extern void *OPS_Bilin(void);
extern void *OPS_Bilin02(void);
extern void *OPS_Steel01(void);
extern void *OPS_SteelMP(void);
extern void *OPS_FRPConfinedConcrete02(void);
extern void *OPS_Steel02(void);
extern void *OPS_Steel03(void);
extern void *OPS_SteelFractureDI(void); // galvisf
extern void *OPS_Steel02Fatigue(void);
extern void *OPS_RambergOsgoodSteel(void);
extern void *OPS_ReinforcingSteel(void);
extern void *OPS_SteelDRC(void); // R. Carreno
extern void *OPS_Concrete01(void);
extern void *OPS_Concrete01WithSITC(void);
extern void *OPS_Concrete02(void);
extern void *OPS_Concrete02IS(void);
extern void *OPS_Concrete04(void);
extern void *OPS_Concrete06(void);
extern void *OPS_Concrete07(void);
extern void *OPS_PinchingLimitStateMaterial(void);
extern void *OPS_SAWSMaterial(void);
extern void *OPS_ConcreteZ01Material(void);
extern void *OPS_ConcreteL01Material(void);
extern void *OPS_SteelZ01Material(void);
extern void *OPS_TendonL01Material(void);
extern void *OPS_ConfinedConcrete01Material(void);
extern void *OPS_ElasticBilin(void);
extern void *OPS_SelfCenteringMaterial(void);
extern void *OPS_MinMaxMaterial(void);
extern void *OPS_SimpleFractureMaterial(void);
extern void *OPS_HoehlerStanton(void);
extern void *OPS_InitStrainMaterial(void);
extern void *OPS_InitStressMaterial(void);
extern void *OPS_pyUCLA(void);
extern void *OPS_Maxwell(void);
extern void *OPS_ViscousDamper(void);
extern void *OPS_DamperMaterial(void);
extern void *OPS_BilinearOilDamper(void);
extern void *OPS_Cast(void);
extern void *OPS_Dodd_Restrepo(void);
extern void *OPS_DoddRestr(void);
extern void *OPS_ElasticMultiLinear(void);
extern void *OPS_ImpactMaterial(void);
extern void *OPS_SteelBRB(void);
extern void *OPS_MultiLinear(void);
extern void *OPS_HookGap(void);
extern void *OPS_HyperbolicGapMaterial(void);
extern void *OPS_FRPConfinedConcrete(void);
extern void *OPS_FRPConfinedConcrete02(void);
extern void *OPS_UVCuniaxial(void);
extern void *OPS_Steel01Thermal(void);
extern void *OPS_Steel02Thermal(void);
extern void *OPS_Concrete02Thermal(void);
extern void *OPS_StainlessECThermal(void); // L.Jiang [SIF]
extern void *OPS_SteelECThermal(void); // L.Jiang [SIF]
extern void *OPS_ConcreteECThermal(void);// L.Jiang [SIF]
extern void *OPS_ElasticMaterialThermal(void); //L.Jiang[SIF]
//extern void *OPS_PlateBearingConnectionThermal(void);
extern void* OPS_ASD_SMA_3K(void); // Luca Aceto
extern void *OPS_BWBN(void);
extern void *OPS_IMKPeakOriented(void);
extern void *OPS_IMKBilin(void);
extern void *OPS_IMKPinching(void);
extern void *OPS_ModIMKPeakOriented(void);
extern void *OPS_ModIMKPeakOriented02(void);
extern void *OPS_ModIMKPinching(void);
extern void *OPS_ModIMKPinching02(void);
extern void *OPS_ConcretewBeta(void);
extern void *OPS_ConcreteD(void);
extern void *OPS_PinchingLimitState(void);
extern void *OPS_OriginCentered(void);
extern void *OPS_Steel2(void);
extern void *OPS_ConcreteSakaiKawashima(void);
extern void *OPS_ResilienceMaterialHR(void);
extern void *OPS_CFSSSWP(void);
extern void *OPS_CFSWSWP(void);
extern void *OPS_ResilienceLow(void);
extern void *OPS_ViscousMaterial(void);
extern void *OPS_SteelMPF(void); // K Kolozvari                                
extern void *OPS_ConcreteCM(void); // K Kolozvari
extern void *OPS_Bond_SP01(void); // K Kolozvari
extern void *OPS_FRCC(void); // Feras Khlef + Andre Barbosa
extern void *OPS_Steel4(void);
extern void *OPS_PySimple3(void);
extern void *OPS_BoucWenMaterial(void);
extern void *OPS_BoucWenOriginal(void);
extern void *OPS_GNGMaterial(void);
extern void *OPS_OOHystereticMaterial(void);
extern void *OPS_ElasticPowerFunc(void);
extern void *OPS_UVCuniaxial(void);
extern void *OPS_DegradingPinchedBW(void);
extern void* OPS_BoucWenInfill(void);  // S. Sirotti  18-January-2022  e-mail: stefano.sirotti@unimore.it
extern void *OPS_SLModel(void);
extern void *OPS_SMAMaterial(void);
extern void* OPS_HystereticPoly(void); // Salvatore Sessa 14-Jan-2021 Mail: salvatore.sessa2@unina.it
extern void* OPS_HystereticSmooth(void); // Salvatore Sessa 02-May-2022 Mail: salvatore.sessa2@unina.it
extern void* OPS_HystereticAsym(void); // Salvatore Sessa 02-May-2022 Mail: salvatore.sessa2@unina.it
extern void *OPS_Masonry(void);
extern void *OPS_Trilinwp(void);
extern void *OPS_Trilinwp2(void);
extern void *OPS_Masonryt(void);
extern void *OPS_DowelType(void);
extern void *OPS_DuctileFracture(void); // Kuanshi Zhong
extern void *OPS_MultiplierMaterial(void);
extern void *OPS_AxialSp(void);
extern void *OPS_AxialSpHD(void);
extern void *OPS_KikuchiAikenHDR(void);
extern void *OPS_KikuchiAikenLRB(void);

extern UniaxialMaterial *
Tcl_AddLimitStateMaterial(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **arg);

extern UniaxialMaterial *Tcl_addWrapperUniaxialMaterial(matObj *, ClientData clientData, Tcl_Interp *interp,
							int argc, TCL_Char **argv);

#include <packages.h>

typedef struct uniaxialPackageCommand {
  char *funcName;
  void * (*funcPtr)(); 
  struct uniaxialPackageCommand *next;
} UniaxialPackageCommand;

static UniaxialPackageCommand *theUniaxialPackageCommands = NULL;

static void printCommand(int argc, TCL_Char **argv)
{
    opserr << "Input command: ";
    for (int i=0; i<argc; i++)
	opserr << argv[i] << " ";
    opserr << endln;
} 


// external functions


UniaxialMaterial *
TclModelBuilder_addFedeasMaterial(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
				  

UniaxialMaterial *
TclModelBuilder_addDrainMaterial(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
				 

UniaxialMaterial *
TclModelBuilder_addSnapMaterial(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
				

UniaxialMaterial *
TclModelBuilder_addPyTzQzMaterial(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv, Domain *theDomain);

UniaxialMaterial *
TclModelBuilder_FRPCnfinedConcrete(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv, Domain *theDomain);

UniaxialMaterial *
TclModelBuilder_addDegradingMaterial(ClientData, Tcl_Interp *, int , TCL_Char **);				  


int
TclModelBuilderUniaxialMaterialCommand (ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv, Domain *theDomain)
{
  
  // Make sure there is a minimum number of arguments
    if (argc < 3) {
	opserr << "WARNING insufficient number of uniaxial material arguments\n";
	opserr << "Want: uniaxialMaterial type? tag? <specific material args>" << endln;
	return TCL_ERROR;
    }

    OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, theDomain);	  

    // Pointer to a uniaxial material that will be added to the model builder
    UniaxialMaterial *theMaterial = 0;

    // Check argv[2] for uniaxial material type
    if (strcmp(argv[1],"Elastic") == 0) {

      void *theMat = OPS_ElasticMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
	
	// SAJalali
	if (strcmp(argv[1], "SPSW02") == 0) {
		void *theMat = OPS_SPSW02();
		if (theMat != 0)
			theMaterial = (UniaxialMaterial *)theMat;
		else
			return TCL_ERROR;
	}

	// ntosic
	if (strcmp(argv[1], "TDConcreteEXP") == 0) {
		void *theMat = OPS_TDConcreteEXP();
		if (theMat != 0)
			theMaterial = (UniaxialMaterial *)theMat;
		else
			return TCL_ERROR;
	}

	// ntosic
	if (strcmp(argv[1], "TDConcrete") == 0) {
		void *theMat = OPS_TDConcrete();
		if (theMat != 0)
			theMaterial = (UniaxialMaterial *)theMat;
		else
			return TCL_ERROR;
	}

	// ntosic
	if (strcmp(argv[1], "TDConcreteMC10") == 0) {
		void *theMat = OPS_TDConcreteMC10();
		if (theMat != 0)
			theMaterial = (UniaxialMaterial *)theMat;
		else
			return TCL_ERROR;
	}

	// ntosic
	if (strcmp(argv[1], "TDConcreteMC10NL") == 0) {
		void *theMat = OPS_TDConcreteMC10NL();
		if (theMat != 0)
			theMaterial = (UniaxialMaterial *)theMat;
		else
			return TCL_ERROR;
	}

	if (strcmp(argv[1],"Steel01") == 0) {
	  
	  void *theMat = OPS_Steel01();
	  if (theMat != 0) 
	    theMaterial = (UniaxialMaterial *)theMat;
	  else 
	    return TCL_ERROR;
	  
	}
	if (strcmp(argv[1],"Steel02") == 0) {
	  void *theMat = OPS_Steel02();
	  if (theMat != 0) 
	    theMaterial = (UniaxialMaterial *)theMat;
	  else 
	    return TCL_ERROR;
	  
	}
    if (strcmp(argv[1], "SteelFractureDI") == 0) {
        void* theMat = OPS_SteelFractureDI();
        if (theMat != 0)
            theMaterial = (UniaxialMaterial*)theMat;
        else
            return TCL_ERROR;

    }
    if (strcmp(argv[1],"Steel02Fatigue") == 0) {
      void *theMat = OPS_Steel02Fatigue();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if (strcmp(argv[1],"Steel4") == 0) {
      void *theMat = OPS_Steel4();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if (strcmp(argv[1],"UVCuniaxial") == 0) {
      void *theMat = OPS_UVCuniaxial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if (strcmp(argv[1],"GNG") == 0) {
      void *theMat = OPS_GNGMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if (strcmp(argv[1],"PySimple3") == 0) {
      void *theMat = OPS_PySimple3();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if (strcmp(argv[1],"Concrete01") == 0) {
      void *theMat = OPS_Concrete01();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
      /*
	    } 
      if (strcmp(argv[1],"HoehlerStanton") == 0) {
      void *theMat = OPS_HoehlerStanton();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
      */
    }
    if (strcmp(argv[1],"Concrete02") == 0) {
      void *theMat = OPS_Concrete02();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if (strcmp(argv[1],"Concrete02IS") == 0) {
      void *theMat = OPS_Concrete02IS();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if ((strcmp(argv[1],"ElasticBilin") == 0) || (strcmp(argv[1],"ElasticBilinear") == 0)) {
      void *theMat = OPS_ElasticBilin();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if ((strcmp(argv[1],"ImpactMaterial") == 0) || (strcmp(argv[1],"Impact") == 0)) {
      void *theMat = OPS_ImpactMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if ((strcmp(argv[1],"SteelBRB") == 0)) {
      void *theMat = OPS_SteelBRB();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if ((strcmp(argv[1],"MinMaxMaterial") == 0) || (strcmp(argv[1],"MinMax") == 0)) {
      void *theMat = OPS_MinMaxMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if ((strcmp(argv[1],"SimpleFractureMaterial") == 0) || (strcmp(argv[1],"SimpleFracture") == 0)) {

      void *theMat = OPS_SimpleFractureMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if ((strcmp(argv[1],"Maxwell") == 0) || (strcmp(argv[1],"MaxwellMaterial") == 0)) {
      void *theMat = OPS_Maxwell();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if ((strcmp(argv[1],"ViscousDamper") == 0)) {
      void *theMat = OPS_ViscousDamper();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if (strcmp(argv[1],"DamperMaterial") == 0 || strcmp(argv[1],"Damper") == 0) {
      void *theMat = OPS_DamperMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if ((strcmp(argv[1],"BilinearOilDamper") == 0)) {
      void *theMat = OPS_BilinearOilDamper();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if ((strcmp(argv[1],"Bond_SP01") == 0) || (strcmp(argv[1],"Bond") == 0)) { 
      void *theMat = OPS_Bond_SP01();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;      

    }
    if (strcmp(argv[1], "FRCC") == 0) {
      void *theMat = OPS_FRCC();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;      
    }
    if ((strcmp(argv[1],"Cast") == 0) || (strcmp(argv[1],"CastFuse") == 0)) {
      void *theMat = OPS_Cast();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if ((strcmp(argv[1],"Dodd_Restrepo") == 0) || 
	       (strcmp(argv[1],"DoddRestrepo") == 0) || 
	       (strcmp(argv[1],"Restrepo") == 0)) {

      void *theMat = OPS_Dodd_Restrepo();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

#ifndef _NO_NEW_RESTREPO
	}
    if ((strcmp(argv[1],"DoddRestr") == 0)) {

      void *theMat = OPS_DoddRestr();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
#endif
    }
    if (strcmp(argv[1],"ElasticMultiLinear") == 0) {
      void *theMat = OPS_ElasticMultiLinear();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }

    if (strcmp(argv[1], "ElasticPowerFunc") == 0) {
      void* theMat = OPS_ElasticPowerFunc();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial*)theMat;
      else
        return TCL_ERROR;
    }
    if (strcmp(argv[1], "UVCuniaxial") == 0) {
      void *theMat = OPS_UVCuniaxial();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;
    }
    if (strcmp(argv[1], "SLModel") == 0) {
      void *theMat = OPS_SLModel();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;
    }
    if ((strcmp(argv[1],"RambergOsgood") == 0) || (strcmp(argv[1],"RambergOsgoodSteel") == 0)) {
      void *theMat = OPS_RambergOsgoodSteel();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if (strcmp(argv[1],"ReinforcingSteel") == 0) {
      void *theMat = OPS_ReinforcingSteel();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if (strcmp(argv[1], "SteelDRC") == 0) {
		void *theMat = OPS_SteelDRC();
		if (theMat != 0)
			theMaterial = (UniaxialMaterial *)theMat;
		else
			return TCL_ERROR;
    }
    if (strcmp(argv[1],"Steel2") == 0) {
      void *theMat = OPS_Steel2();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if (strcmp(argv[1],"OriginCentered") == 0) {
      void *theMat = OPS_OriginCentered();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if (strcmp(argv[1],"HookGap") == 0) {
      void *theMat = OPS_HookGap();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if (strcmp(argv[1],"HyperbolicGapMaterial") == 0) {
      void *theMat = OPS_HyperbolicGapMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if (strcmp(argv[1],"FRPConfinedConcrete02") == 0) {
      void *theMat = OPS_FRPConfinedConcrete02();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if ((strcmp(argv[1],"PinchingLimitState") == 0) || (strcmp(argv[1],"PinchingLimitStateMaterial") == 0)) {
      void *theMat = OPS_PinchingLimitState();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if ((strcmp(argv[1],"InitStrainMaterial") == 0) || (strcmp(argv[1],"InitStrain") == 0)) {
      void *theMat = OPS_InitStrainMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if ((strcmp(argv[1],"InitStressMaterial") == 0) || (strcmp(argv[1],"InitStress") == 0)) {
      void *theMat = OPS_InitStressMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if ((strcmp(argv[1],"pyUCLA") == 0) || (strcmp(argv[1],"PYUCLA") == 0)) {
      void *theMat = OPS_pyUCLA();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if ((strcmp(argv[1],"MultiLinear") == 0)) {
      void *theMat = OPS_MultiLinear();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if (strcmp(argv[1],"ModIMKPinching") == 0) {
      void *theMat = OPS_ModIMKPinching();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if (strcmp(argv[1],"ModIMKPinching02") == 0) {
      void *theMat = OPS_ModIMKPinching02();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if (strcmp(argv[1], "BoucWenOriginal") == 0) {
        void *theMat = OPS_BoucWenOriginal();
        if (theMat != 0)
            theMaterial = (UniaxialMaterial *)theMat;
        else
            return TCL_ERROR;

    }
    if (strcmp(argv[1], "BWBN") == 0) {
        void *theMat = OPS_BWBN();
        if (theMat != 0)
            theMaterial = (UniaxialMaterial *)theMat;
        else
            return TCL_ERROR;

    }
    if (strcmp(argv[1], "DegradingPinchedBW") == 0) {
        void *theMat = OPS_DegradingPinchedBW();
        if (theMat != 0)
            theMaterial = (UniaxialMaterial *)theMat;
        else
            return TCL_ERROR;

    }
    if (strcmp(argv[1], "BoucWenInfill") == 0) {
    void *theMat = OPS_BoucWenInfill();
    if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
    else
        return TCL_ERROR;
    }
    if (strcmp(argv[1], "IMKBilin") == 0) {
      void *theMat = OPS_IMKBilin();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;

    }
    if (strcmp(argv[1], "IMKPeakOriented") == 0) {
      void *theMat = OPS_IMKPeakOriented();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;

    }
    if (strcmp(argv[1], "IMKPinching") == 0) {
      void *theMat = OPS_IMKPinching();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;
    }
    if (strcmp(argv[1], "ModIMKPeakOriented") == 0) {
      void *theMat = OPS_ModIMKPeakOriented();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    } if (strcmp(argv[1],"ModIMKPeakOriented02") == 0) {
      void *theMat = OPS_ModIMKPeakOriented02();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if (strcmp(argv[1],"Bilin02") == 0) {
      void *theMat = OPS_Bilin02();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    if (strcmp(argv[1],"PlateBearingConnectionThermal") == 0) {
      //void *theMat = OPS_PlateBearingConnectionThermal();
      void *theMat = 0;
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if (strcmp(argv[1],"Steel01Thermal") == 0) {
      void *theMat = OPS_Steel01Thermal();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;      
    }
    if (strcmp(argv[1],"Steel02Thermal") == 0) {
      void *theMat = OPS_Steel02Thermal();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if (strcmp(argv[1], "SteelECThermal") == 0) {
		void *theMat = OPS_SteelECThermal();
		if (theMat != 0)
			theMaterial = (UniaxialMaterial *)theMat;
		else
			return TCL_ERROR;
		//------End of adding identity for SteelEcThermal	
	}
    if (strcmp(argv[1], "StainlessECThermal") == 0) {
		void *theMat = OPS_StainlessECThermal();
		if (theMat != 0)
			theMaterial = (UniaxialMaterial *)theMat;
		else
			return TCL_ERROR;
		//------End of adding identity for StainlessECThermal
	}
    if (strcmp(argv[1], "ElasticThermal") == 0) {
		void *theMat = OPS_ElasticMaterialThermal();
		if (theMat != 0)
			theMaterial = (UniaxialMaterial *)theMat;
		else
			return TCL_ERROR;

	}
    if (strcmp(argv[1], "ConcreteECThermal") == 0) {
		void *theMat = OPS_ConcreteECThermal();
		if (theMat != 0)
			theMaterial = (UniaxialMaterial *)theMat;
		else
			return TCL_ERROR;
		// end of adding More thermo-mechanical uniaxial materials, L.Jiang[SIF]

    }
    if (strcmp(argv[1],"ConcretewBeta") == 0) {
      void *theMat = OPS_ConcretewBeta();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if (strcmp(argv[1],"ConcreteD") == 0) {
      void *theMat = OPS_ConcreteD();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if (strcmp(argv[1],"ConcreteSakaiKawashima") == 0) {
      void *theMat = OPS_ConcreteSakaiKawashima();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if (strcmp(argv[1],"Concrete02Thermal") == 0) {
      void *theMat = OPS_Concrete02Thermal();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if (strcmp(argv[1],"SteelMPF") == 0) { // K Kolozvari            
      void *theMat = OPS_SteelMPF();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;
    }
    if (strcmp(argv[1],"ConcreteCM") == 0) { // K Kolozvari          
      void *theMat = OPS_ConcreteCM();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;
    }
    if (strcmp(argv[1],"ResilienceLow") == 0) {
      void *theMat = OPS_ResilienceLow();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if (strcmp(argv[1],"ResilienceMaterialHR") == 0) {
      void *theMat = OPS_ResilienceMaterialHR();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if (strcmp(argv[1],"CFSWSWP") == 0) {
      void *theMat = OPS_CFSWSWP();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if (strcmp(argv[1],"CFSSSWP") == 0) {
      void *theMat = OPS_CFSSSWP();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if (strcmp(argv[1],"FRPConfinedConcrete") == 0) {
      void *theMat = OPS_FRPConfinedConcrete();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if (strcmp(argv[1], "Masonry") == 0) {
	  void *theMat = OPS_Masonry();
	  if (theMat != 0)
	    theMaterial = (UniaxialMaterial *)theMat;
	  else
	    return TCL_ERROR;
	}
    if (strcmp(argv[1], "Trilinwp") == 0) {
	  void *theMat = OPS_Trilinwp();
	  if (theMat != 0)
	    theMaterial = (UniaxialMaterial *)theMat;
	  else
	    return TCL_ERROR;
	}
    if (strcmp(argv[1], "Trilinwp2") == 0) {
	  void *theMat = OPS_Trilinwp2();
	  if (theMat != 0)
	    theMaterial = (UniaxialMaterial *)theMat;
	  else
	    return TCL_ERROR;
	}
    if (strcmp(argv[1], "Masonryt") == 0) {
	  void *theMat = OPS_Masonryt();
	  if (theMat != 0)
	    theMaterial = (UniaxialMaterial *)theMat;
	  else
	    return TCL_ERROR;
    }
    if (strcmp(argv[1],"Elastic2") == 0) {
      void *theMat = OPS_Elastic2Material();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if (strcmp(argv[1],"ENT") == 0) {
      void *theMat = OPS_ENTMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if (strcmp(argv[1],"ElasticPP") == 0) {
      void *theMat = OPS_ElasticPPMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if (strcmp(argv[1],"ElasticPPGap") == 0) {
      void *theMat = OPS_EPPGapMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if (strcmp(argv[1],"Hardening") == 0 || strcmp(argv[1],"Hardening2") == 0) {
      
      void *theMat = OPS_HardeningMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if (strcmp(argv[1],"BoucWen") == 0) {

      void *theMat = OPS_BoucWenMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;      
    }
    if (strcmp(argv[1],"Parallel") == 0) {

      void *theMat = OPS_ParallelMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if (strcmp(argv[1],"Series") == 0) {

      void *theMat = OPS_SeriesMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if (strcmp(argv[1],"Steel03") == 0) {
      
      void *theMat = OPS_Steel03();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if (strcmp(argv[1],"Hysteretic") == 0) {

      void *theMat = OPS_HystereticMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if (strcmp(argv[1],"OOHysteretic") == 0) {

      void *theMat = OPS_OOHystereticMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if (strcmp(argv[1],"Concrete04") == 0) {

      void *theMat = OPS_Concrete04();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;      
    }
    if (strcmp(argv[1],"Concrete06") == 0) {

      void *theMat = OPS_Concrete06();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;            
    }
    if (strcmp(argv[1], "Concrete07") == 0) {

      void *theMat = OPS_Concrete07();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;            
    }
    if (strcmp(argv[1],"Viscous") == 0) {

      void *theMat = OPS_ViscousMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if (strcmp(argv[1],"PathIndependent") == 0) {

      void *theMat = OPS_PathIndependentMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;      
    }
    if (strcmp(argv[1],"Backbone") == 0) {

      void *theMat = OPS_BackboneMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;            
    }
    if (strcmp(argv[1],"Fatigue") == 0) {

      void *theMat = OPS_FatigueMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;      
    }

    if ((strcmp(argv[1],"SAWSMaterial") == 0) || (strcmp(argv[1],"SAWS") == 0)) {

      void *theMat = OPS_SAWSMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if ((strcmp(argv[1],"BilinMaterial") == 0) || (strcmp(argv[1],"Bilin") == 0)) {

      void *theMat = OPS_Bilin();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if ((strcmp(argv[1],"ConcreteZ01Material") == 0) || (strcmp(argv[1],"ConcreteZ01") == 0)) {

      void *theMat = OPS_ConcreteZ01Material();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if ((strcmp(argv[1],"ConcreteL01Material") == 0) || (strcmp(argv[1],"ConcreteL01") == 0)) {

      void *theMat = OPS_ConcreteL01Material();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if ((strcmp(argv[1],"SteelZ01Material") == 0) || (strcmp(argv[1],"SteelZ01") == 0)) {

      void *theMat = OPS_SteelZ01Material();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if ((strcmp(argv[1],"TendonL01Material") == 0) || (strcmp(argv[1],"TendonL01") == 0)) {
      void *theMat = OPS_TendonL01Material();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if ((strcmp(argv[1],"ConfinedConcrete01") == 0) || (strcmp(argv[1],"ConfinedConcrete") == 0)) {

      void *theMat = OPS_ConfinedConcrete01Material();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if (strcmp(argv[1],"Cable") == 0) {

      void *theMat = OPS_CableMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    if (strcmp(argv[1], "UVCuniaxial") == 0) {
      void* theMat = OPS_UVCuniaxial();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial*)theMat;
      else
        return TCL_ERROR;
    }
    if (strcmp(argv[1],"Pinching4") == 0) {
		if (argc != 42 && argc != 31 ) {
			opserr << "WARNING insufficient arguments\n";
			printCommand(argc,argv);
			opserr << "Want: uniaxialMaterial Pinching4 tag? stress1p? strain1p? stress2p? strain2p? stress3p? strain3p? stress4p? strain4p? "
				<< "\n<stress1n? strain1n? stress2n? strain2n? stress3n? strain3n? stress4n? strain4n?> rDispP? rForceP? uForceP? "
				<< "\n<rDispN? rForceN? uForceN?> gammaK1? gammaK2? gammaK3? gammaK4? gammaKLimit? gammaD1? gammaD2? gammaD3? gammaD4? "
				<< "\ngammaDLimit? gammaF1? gammaF2? gammaF3? gammaF4? gammaFLimit? gammaE? CycleOrEnergyDamage? ";
			return TCL_ERROR;
		}

		int tag, tDmg;
		double stress1p, stress2p, stress3p, stress4p;
		double strain1p, strain2p, strain3p, strain4p;
		double stress1n, stress2n, stress3n, stress4n;
		double strain1n, strain2n, strain3n, strain4n;
		double rDispP, rForceP, uForceP, rDispN, rForceN, uForceN;
		double gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit;
		double gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit;
		double gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit;
		double gammaE;

		int i = 2;

		if (Tcl_GetInt(interp, argv[i++], &tag) != TCL_OK) {
			opserr << "WARNING invalid uniaxialMaterial Pinching4 tag" << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &stress1p) != TCL_OK) {
			opserr << "WARNING invalid stress1p\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &strain1p) != TCL_OK) {
			opserr << "WARNING invalid strain1p\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &stress2p) != TCL_OK) {
			opserr << "WARNING invalid stress2p\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &strain2p) != TCL_OK) {
			opserr << "WARNING invalid strain2p\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &stress3p) != TCL_OK) {
			opserr << "WARNING invalid stress3p\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &strain3p) != TCL_OK) {
			opserr << "WARNING invalid strain3p\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &stress4p) != TCL_OK) {
			opserr << "WARNING invalid stress4p\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &strain4p) != TCL_OK) {
			opserr << "WARNING invalid strain4p\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

		if (argc == 42) {
			if (Tcl_GetDouble(interp, argv[i++], &stress1n) != TCL_OK) {
				opserr << "WARNING invalid stress1n\n";
				opserr << "Pinching4 material: " << tag << endln;
				return TCL_ERROR;
			}

			if (Tcl_GetDouble(interp, argv[i++], &strain1n) != TCL_OK) {
				opserr << "WARNING invalid strain1n\n";
				opserr << "Pinching4 material: " << tag << endln;
				return TCL_ERROR;
			}

			if (Tcl_GetDouble(interp, argv[i++], &stress2n) != TCL_OK) {
				opserr << "WARNING invalid stress2n\n";
				opserr << "Pinching4 material: " << tag << endln;
				return TCL_ERROR;
			}

			if (Tcl_GetDouble(interp, argv[i++], &strain2n) != TCL_OK) {
				opserr << "WARNING invalid strain2n\n";
				opserr << "Pinching4 material: " << tag << endln;
				return TCL_ERROR;
			}

			if (Tcl_GetDouble(interp, argv[i++], &stress3n) != TCL_OK) {
				opserr << "WARNING invalid stress3n\n";
				opserr << "Pinching4 material: " << tag << endln;
				return TCL_ERROR;
			}

			if (Tcl_GetDouble(interp, argv[i++], &strain3n) != TCL_OK) {
				opserr << "WARNING invalid strain3n\n";
				opserr << "Pinching4 material: " << tag << endln;
				return TCL_ERROR;
			}

			if (Tcl_GetDouble(interp, argv[i++], &stress4n) != TCL_OK) {
				opserr << "WARNING invalid stress4n\n";
				opserr << "Pinching4 material: " << tag << endln;
				return TCL_ERROR;
			}

			if (Tcl_GetDouble(interp, argv[i++], &strain4n) != TCL_OK) {
				opserr << "WARNING invalid strain4n\n";
				opserr << "Pinching4 material: " << tag << endln;
				return TCL_ERROR;
			}
		
		}
		

		if (Tcl_GetDouble(interp, argv[i++], &rDispP) != TCL_OK) {
			opserr << "WARNING invalid rDispP\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &rForceP) != TCL_OK) {
			opserr << "WARNING invalid rForceP\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &uForceP) != TCL_OK) {
			opserr << "WARNING invalid uForceP\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

		if (argc == 42) {
			if (Tcl_GetDouble(interp, argv[i++], &rDispN) != TCL_OK) {
				opserr << "WARNING invalid rDispN\n";
				opserr << "Pinching4 material: " << tag << endln;
				return TCL_ERROR;
			}

			if (Tcl_GetDouble(interp, argv[i++], &rForceN) != TCL_OK) {
				opserr << "WARNING invalid rForceN\n";
				opserr << "Pinching4 material: " << tag << endln;
				return TCL_ERROR;
			}

			if (Tcl_GetDouble(interp, argv[i++], &uForceN) != TCL_OK) {
				opserr << "WARNING invalid uForceN\n";
				opserr << "Pinching4 material: " << tag << endln;
				return TCL_ERROR;
			}
		}

		if (Tcl_GetDouble(interp, argv[i++], &gammaK1) != TCL_OK) {
			opserr << "WARNING invalid gammaK1\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[i++], &gammaK2) != TCL_OK) {
			opserr << "WARNING invalid gammaK2\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[i++], &gammaK3) != TCL_OK) {
			opserr << "WARNING invalid gammaK3\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[i++], &gammaK4) != TCL_OK) {
			opserr << "WARNING invalid gammaK4\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[i++], &gammaKLimit) != TCL_OK) {
			opserr << "WARNING invalid gammaKLimit\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[i++], &gammaD1) != TCL_OK) {
			opserr << "WARNING invalid gammaD1\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}										   
		if (Tcl_GetDouble(interp, argv[i++], &gammaD2) != TCL_OK) {
			opserr << "WARNING invalid gammaD2\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[i++], &gammaD3) != TCL_OK) {
			opserr << "WARNING invalid gammaD3\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[i++], &gammaD4) != TCL_OK) {
			opserr << "WARNING invalid gammaD4\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[i++], &gammaDLimit) != TCL_OK) {
			opserr << "WARNING invalid gammaDLimit\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[i++], &gammaF1) != TCL_OK) {
			opserr << "WARNING invalid gammaF1\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[i++], &gammaF2) != TCL_OK) {
			opserr << "WARNING invalid gammaF2\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[i++], &gammaF3) != TCL_OK) {
			opserr << "WARNING invalid gammaF3\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[i++], &gammaF4) != TCL_OK) {
			opserr << "WARNING invalid gammaF4\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}
		if (Tcl_GetDouble(interp, argv[i++], &gammaFLimit) != TCL_OK) {
			opserr << "WARNING invalid gammaFLimit\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[i++], &gammaE) != TCL_OK) {
			opserr << "WARNING invalid gammaE\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

		int y; 
		y = i;

		if ((strcmp(argv[y],"cycle") == 0) || (strcmp(argv[y],"Cycle") == 0) || (strcmp(argv[y],"DamageCycle") == 0) || (strcmp(argv[y],"damageCycle") == 0))
		{ tDmg = 1; }
		else if ((strcmp(argv[y],"energy") == 0) || (strcmp(argv[y],"Energy") == 0) || (strcmp(argv[y],"DamageEnergy") == 0) || (strcmp(argv[y],"damageEnergy") == 0))
		{ tDmg = 0; }
		else
		{
			opserr << "WARNING invalid type of damage calculation specified\n";
			opserr << "Pinching4 material: " << tag << endln;
			return TCL_ERROR;
		}

	// allocate the pinching material
		if (argc == 42) {
		theMaterial = new Pinching4Material (tag,
			stress1p, strain1p, stress2p, strain2p, stress3p, strain3p, stress4p, strain4p,
			stress1n, strain1n, stress2n, strain2n, stress3n, strain3n, stress4n, strain4n,
			rDispP, rForceP, uForceP, rDispN, rForceN, uForceN, 
			gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit,
			gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit,
			gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit, gammaE, tDmg);
		}
		if (argc == 31) {
		theMaterial = new Pinching4Material (tag,
			stress1p, strain1p, stress2p, strain2p, stress3p, strain3p, stress4p, strain4p,
			rDispP, rForceP, uForceP,  
			gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit,
			gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit,
			gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit, gammaE, tDmg);		
		}
   }
  if (strcmp(argv[1],"BarSlip") == 0) {
     if (argc != 17 && argc != 15)
       {
	 opserr << "WARNING insufficient arguments\n";
	 printCommand(argc,argv);
	 opserr << "Want: uniaxialMaterial BarSlip tag? fc? fy? Es? fu? Eh? db? ld? nb? width? depth? bsflag? type? <damage? unit?>"  << endln;
	 return TCL_ERROR;
       }
     
     int tag, nb, bsf, typ, dmg, unt;
     double fc, fy, Es, fu, Eh, ld, width, depth, db;
     
     int argStart = 2;
     
     if (Tcl_GetInt(interp, argv[argStart++], &tag) != TCL_OK)
       {
	 opserr << "WARNING invalid tag\n";
	 opserr << "BarSlip: " << tag << endln;
	 return TCL_ERROR;
       }
     if (Tcl_GetDouble(interp, argv[argStart++], &fc) != TCL_OK)
       {
	 opserr << "WARNING invalid fc\n";
	 opserr << "BarSlip: " << tag << endln;
	 return TCL_ERROR;
       }
     if (Tcl_GetDouble(interp, argv[argStart++], &fy) != TCL_OK)
       {
	 opserr << "WARNING invalid fy\n";
	 opserr << "BarSlip: " << tag << endln;
	 return TCL_ERROR;
       }
     if (Tcl_GetDouble(interp, argv[argStart++], &Es) != TCL_OK)
       {
	 opserr << "WARNING invalid Es\n";
	 opserr << "BarSlip: " << tag << endln;
	 return TCL_ERROR;
       }
     if (Tcl_GetDouble(interp, argv[argStart++], &fu) != TCL_OK)
       {
	 opserr << "WARNING invalid fu\n";
	 opserr << "BarSlip: " << tag << endln;
	 return TCL_ERROR;
       }
     if (Tcl_GetDouble(interp, argv[argStart++], &Eh) != TCL_OK)
       {
	 opserr << "WARNING invalid Eh\n";
	 opserr << "BarSlip: " << tag << endln;
	 return TCL_ERROR;
       }
     if (Tcl_GetDouble(interp, argv[argStart++], &db) != TCL_OK)
       {
	 opserr << "WARNING invalid db\n";
	 opserr << "BarSlip: " << tag << endln;
	 return TCL_ERROR;
       }
     if (Tcl_GetDouble(interp, argv[argStart++], &ld) != TCL_OK)
       {
	 opserr << "WARNING invalid ld\n";
	 opserr << "BarSlip: " << tag << endln;
	 return TCL_ERROR;
       }
     if (Tcl_GetInt(interp, argv[argStart++], &nb) != TCL_OK)
       {
	 opserr << "WARNING invalid nbars\n";
	 opserr << "BarSlip: " << tag << endln;
	 return TCL_ERROR;
       }
     if (Tcl_GetDouble(interp, argv[argStart++], &width) != TCL_OK)
       {
	 opserr << "WARNING invalid width\n";
	 opserr << "BarSlip: " << tag << endln;
	 return TCL_ERROR;
       }
     if (Tcl_GetDouble(interp, argv[argStart++], &depth) != TCL_OK)
       {
	 opserr << "WARNING invalid depth\n";
	 opserr << "BarSlip: " << tag << endln;
	 return TCL_ERROR;
       }
     
     int y;
     y = argStart;
     
     
     if ((strcmp(argv[y],"strong") == 0) || (strcmp(argv[y],"Strong") == 0) || (strcmp(argv[y],"weak") == 0) || (strcmp(argv[y],"Weak") == 0))
       {
	 if ((strcmp(argv[y],"strong") == 0) || (strcmp(argv[y],"Strong") == 0))
	   {
	     bsf = 0;
	   }
	 
	 if ((strcmp(argv[y],"weak") == 0) || (strcmp(argv[y],"Weak") == 0))
	   {
	     bsf = 1;
	   }
       }
     else
       {
	 opserr << "WARNING invalid bond strength specified\n";
	 opserr << "BarSlip: " << tag << endln;
	 return TCL_ERROR;
       }
     y ++;
     
     if ((strcmp(argv[y],"beamtop") == 0) || (strcmp(argv[y],"beamTop") == 0) || 
	 (strcmp(argv[y],"beambot") == 0) || (strcmp(argv[y],"beamBot") == 0) || (strcmp(argv[y],"beambottom") == 0) || (strcmp(argv[y],"beamBottom") == 0) ||
	 (strcmp(argv[y],"beam") == 0) || (strcmp(argv[y],"Beam") == 0) || (strcmp(argv[y],"Column") == 0) || (strcmp(argv[y],"column") == 0))
       {
	 if ((strcmp(argv[y],"beamtop") == 0) || (strcmp(argv[y],"beamTop") == 0) || (strcmp(argv[y],"beam") == 0) || (strcmp(argv[y],"Beam") == 0))
	   {
	     typ = 0;
	   }
	 
	 if ((strcmp(argv[y],"beambot") == 0) || (strcmp(argv[y],"beamBot") == 0) || (strcmp(argv[y],"beambottom") == 0) || (strcmp(argv[y],"beamBottom") == 0))
	   {
	     typ = 1;
	   }
	 
	 if ((strcmp(argv[y],"column") == 0) || (strcmp(argv[y],"Column") == 0))
	   {
	     typ = 2;
	   }
       }
     else
       {
	 opserr << "WARNING invalid location of bar specified\n";
	 opserr << "BarSlip: " << tag << endln;
	 return TCL_ERROR;
       }
     if (argc == 17) {
       y ++;
       
       if ((strcmp(argv[y],"damage1") == 0) || (strcmp(argv[y],"Damage1") == 0) || (strcmp(argv[y],"damage2") == 0) || (strcmp(argv[y],"Damage2") == 0) || 
	   (strcmp(argv[y],"nodamage") == 0) || (strcmp(argv[y],"Nodamage") == 0) || (strcmp(argv[y],"NoDamage") == 0) || (strcmp(argv[y],"noDamage") == 0))
	 {
	   if ((strcmp(argv[y],"damage1") == 0) || (strcmp(argv[y],"Damage1") == 0))
	     {
	       dmg = 1;
	     }
	   else if ((strcmp(argv[y],"damage2") == 0) || (strcmp(argv[y],"Damage2") == 0))
	     {
	       dmg = 2;
	     }
	   else if ((strcmp(argv[y],"nodamage") == 0) || (strcmp(argv[y],"Nodamage") == 0) || (strcmp(argv[y],"NoDamage") == 0) || (strcmp(argv[y],"noDamage") == 0))
	     {
	       dmg = 0;
	     }
	   
	 }
       else
	 {
	   opserr << "WARNING invalid damage specified\n";
	   opserr << "BarSlip: " << tag << endln;
	   return TCL_ERROR;
	 }
       
       y ++;
       
       if ((strcmp(argv[y],"mpa") == 0) || (strcmp(argv[y],"MPa") == 0) || (strcmp(argv[y],"mPa") == 0) || (strcmp(argv[y],"Mpa") == 0) ||
	   (strcmp(argv[y],"psi") == 0) || (strcmp(argv[y],"Psi") == 0) || (strcmp(argv[y],"PSI") == 0) || (strcmp(argv[y],"Pa") == 0) ||
	   (strcmp(argv[y],"pa") == 0) ||  (strcmp(argv[y],"psf") == 0) || (strcmp(argv[y],"Psf") == 0) || (strcmp(argv[y],"PSF") == 0) ||
	   (strcmp(argv[y],"ksi") == 0) || (strcmp(argv[y],"Ksi") == 0) || (strcmp(argv[y],"KSI") == 0) || (strcmp(argv[y],"ksf") == 0) ||
	   (strcmp(argv[y],"Ksf") == 0) || (strcmp(argv[y],"KSF") == 0))
	 {
	   if ((strcmp(argv[y],"mpa") == 0) || (strcmp(argv[y],"MPa") == 0) || (strcmp(argv[y],"mPa") == 0) || (strcmp(argv[y],"Mpa") == 0))
	     {
	       unt = 1;
	     }
	   else if ((strcmp(argv[y],"psi") == 0) || (strcmp(argv[y],"Psi") == 0) || (strcmp(argv[y],"PSI") == 0))
	     {
	       unt = 2;
	     }
	   else if ((strcmp(argv[y],"Pa") == 0) || (strcmp(argv[y],"pa") == 0))
	     {
	       unt = 3;
	     }
	   else if ((strcmp(argv[y],"psf") == 0) || (strcmp(argv[y],"Psf") == 0) || (strcmp(argv[y],"PSF") == 0))
	     {
	       unt = 4;
	     }
	   else if ((strcmp(argv[y],"ksi") == 0) || (strcmp(argv[y],"Ksi") == 0) || (strcmp(argv[y],"KSI") == 0))
	     {
	       unt = 5;
	     }
	   else if ((strcmp(argv[y],"ksf") == 0) || (strcmp(argv[y],"Ksf") == 0) || (strcmp(argv[y],"KSF") == 0))
	     {
	       unt = 6;
	     }
	 }
       else
	 {
	   opserr << "WARNING invalid unit specified\n";
	   opserr << "BarSlip: " << tag << endln;
	   return TCL_ERROR;
	 }
     }
     
     // allocate the material
     if (argc == 15 ) {
       theMaterial = new BarSlipMaterial (tag, fc, fy, Es, fu, Eh, db, ld, nb, width, depth, bsf, typ);
     }
     
     if (argc == 17) {
       theMaterial = new BarSlipMaterial (tag, fc, fy, Es, fu, Eh, db, ld, nb, width, depth, bsf, typ, dmg, unt);
     }
     
   } 
    
    
    
  if (strcmp(argv[1],"ShearPanel") == 0) {
    if (argc != 42 && argc != 31 ) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc,argv);
      opserr << "Want: uniaxialMaterial ShearPanel tag? stress1p? strain1p? stress2p? strain2p? stress3p? strain3p? stress4p? strain4p? "
	     << "\n<stress1n? strain1n? stress2n? strain2n? stress3n? strain3n? stress4n? strain4n?> rDispP? rForceP? uForceP? "
	     << "\n<rDispN? rForceN? uForceN?> gammaK1? gammaK2? gammaK3? gammaK4? gammaKLimit? gammaD1? gammaD2? gammaD3? gammaD4? "
	     << "\ngammaDLimit? gammaF1? gammaF2? gammaF3? gammaF4? gammaFLimit? gammaE? YieldStress? ";
      return TCL_ERROR;
    }
    
    int tag;
    double stress1p, stress2p, stress3p, stress4p;
    double strain1p, strain2p, strain3p, strain4p;
    double stress1n, stress2n, stress3n, stress4n;
    double strain1n, strain2n, strain3n, strain4n;
    double rDispP, rForceP, uForceP, rDispN, rForceN, uForceN;
    double gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit;
    double gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit;
    double gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit;
    double gammaE, yStr;
    
    int i = 2;
    
    if (Tcl_GetInt(interp, argv[i++], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial ShearPanel tag" << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &stress1p) != TCL_OK) {
      opserr << "WARNING invalid stress1p\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &strain1p) != TCL_OK) {
      opserr << "WARNING invalid strain1p\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &stress2p) != TCL_OK) {
      opserr << "WARNING invalid stress2p\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &strain2p) != TCL_OK) {
      opserr << "WARNING invalid strain2p\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &stress3p) != TCL_OK) {
      opserr << "WARNING invalid stress3p\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &strain3p) != TCL_OK) {
      opserr << "WARNING invalid strain3p\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &stress4p) != TCL_OK) {
      opserr << "WARNING invalid stress4p\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &strain4p) != TCL_OK) {
      opserr << "WARNING invalid strain4p\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (argc == 42) {
      if (Tcl_GetDouble(interp, argv[i++], &stress1n) != TCL_OK) {
	opserr << "WARNING invalid stress1n\n";
	opserr << "ShearPanel material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[i++], &strain1n) != TCL_OK) {
	opserr << "WARNING invalid strain1n\n";
	opserr << "ShearPanel material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[i++], &stress2n) != TCL_OK) {
	opserr << "WARNING invalid stress2n\n";
	opserr << "ShearPanel material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[i++], &strain2n) != TCL_OK) {
	opserr << "WARNING invalid strain2n\n";
	opserr << "ShearPanel material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[i++], &stress3n) != TCL_OK) {
	opserr << "WARNING invalid stress3n\n";
	opserr << "ShearPanel material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[i++], &strain3n) != TCL_OK) {
	opserr << "WARNING invalid strain3n\n";
	opserr << "ShearPanel material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[i++], &stress4n) != TCL_OK) {
	opserr << "WARNING invalid stress4n\n";
	opserr << "ShearPanel material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[i++], &strain4n) != TCL_OK) {
	opserr << "WARNING invalid strain4n\n";
	opserr << "ShearPanel material: " << tag << endln;
	return TCL_ERROR;
      }
      
    }
    
    
    if (Tcl_GetDouble(interp, argv[i++], &rDispP) != TCL_OK) {
      opserr << "WARNING invalid rDispP\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &rForceP) != TCL_OK) {
      opserr << "WARNING invalid rForceP\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &uForceP) != TCL_OK) {
      opserr << "WARNING invalid uForceP\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (argc == 42) {
      if (Tcl_GetDouble(interp, argv[i++], &rDispN) != TCL_OK) {
	opserr << "WARNING invalid rDispN\n";
	opserr << "ShearPanel material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[i++], &rForceN) != TCL_OK) {
	opserr << "WARNING invalid rForceN\n";
	opserr << "ShearPanel material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[i++], &uForceN) != TCL_OK) {
	opserr << "WARNING invalid uForceN\n";
	opserr << "ShearPanel material: " << tag << endln;
	return TCL_ERROR;
      }
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &gammaK1) != TCL_OK) {
      opserr << "WARNING invalid gammaK1\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaK2) != TCL_OK) {
      opserr << "WARNING invalid gammaK2\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaK3) != TCL_OK) {
      opserr << "WARNING invalid gammaK3\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaK4) != TCL_OK) {
      opserr << "WARNING invalid gammaK4\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaKLimit) != TCL_OK) {
      opserr << "WARNING invalid gammaKLimit\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaD1) != TCL_OK) {
      opserr << "WARNING invalid gammaD1\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }										   
    if (Tcl_GetDouble(interp, argv[i++], &gammaD2) != TCL_OK) {
      opserr << "WARNING invalid gammaD2\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaD3) != TCL_OK) {
      opserr << "WARNING invalid gammaD3\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaD4) != TCL_OK) {
      opserr << "WARNING invalid gammaD4\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaDLimit) != TCL_OK) {
      opserr << "WARNING invalid gammaDLimit\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaF1) != TCL_OK) {
      opserr << "WARNING invalid gammaF1\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaF2) != TCL_OK) {
      opserr << "WARNING invalid gammaF2\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaF3) != TCL_OK) {
      opserr << "WARNING invalid gammaF3\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaF4) != TCL_OK) {
      opserr << "WARNING invalid gammaF4\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[i++], &gammaFLimit) != TCL_OK) {
      opserr << "WARNING invalid gammaFLimit\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &gammaE) != TCL_OK) {
      opserr << "WARNING invalid gammaE\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    if (Tcl_GetDouble(interp, argv[i++], &yStr) != TCL_OK) {
      opserr << "WARNING invalid yield stress\n";
      opserr << "ShearPanel material: " << tag << endln;
      return TCL_ERROR;
    }
    
    // allocate the pinching material
    if (argc == 42) {
      theMaterial = new ShearPanelMaterial (tag,
					    stress1p, strain1p, stress2p, strain2p, stress3p, strain3p, stress4p, strain4p,
					    stress1n, strain1n, stress2n, strain2n, stress3n, strain3n, stress4n, strain4n,
					    rDispP, rForceP, uForceP, rDispN, rForceN, uForceN, 
					    gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit,
					    gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit,
					    gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit, gammaE, yStr);
    }
    if (argc == 31) {
      theMaterial = new ShearPanelMaterial (tag,
					    stress1p, strain1p, stress2p, strain2p, stress3p, strain3p, stress4p, strain4p,
					    rDispP, rForceP, uForceP,  
					    gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit,
					    gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit,
					    gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit, gammaE, yStr);		
    }
  }
    if (strcmp(argv[1],"Concrete01WithSITC") == 0) {
      void *theMat = OPS_Concrete01WithSITC();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;      
    }

    if (strcmp(argv[1],"SMA") == 0) {
      void *theMat = OPS_SMAMaterial();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;      
    }

    if (strcmp(argv[1],"ECC01") == 0) {
      void *theMat = OPS_ECC01();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;      
    }


    if (strcmp(argv[1], "ASD_SMA_3K") == 0) {
      void *theMat = OPS_ASD_SMA_3K();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;
    }
    if (strcmp(argv[1],"SelfCentering") == 0) {
      void *theMat = OPS_SelfCenteringMaterial();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;
    }
    
    if (strcmp(argv[1],"SteelMP") == 0) {
      void *theMat = OPS_SteelMP();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;            
    }

    if (strcmp(argv[1],"SmoothPSConcrete") == 0) {
      void *theMat = OPS_SmoothPSConcrete();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;            
    }
    
    if (strcmp(argv[1],"UniaxialJ2Plasticity") == 0) {
      void *theMat = OPS_UniaxialJ2Plasticity();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;            
    }
    if (strcmp(argv[1],"KikuchiAikenHDR") == 0) { 
      void *theMat = OPS_KikuchiAikenHDR();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;            
    }
    if (strcmp(argv[1],"KikuchiAikenLRB") == 0) { 
      void *theMat = OPS_KikuchiAikenLRB();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;            
    }
    if (strcmp(argv[1],"AxialSp") == 0) { 
      void *theMat = OPS_AxialSp();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;            
    }
    if (strcmp(argv[1],"AxialSpHD") == 0) { 
      void *theMat = OPS_AxialSpHD();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;            
    }
    if (strcmp(argv[1], "HystereticPoly") == 0) {		// BEGIN Salvatore Sessa 14-Jan-2021 Mail: salvatore.sessa2@unina.it
	void* theMat = OPS_HystereticPoly();
	if (theMat != 0)
		theMaterial = (UniaxialMaterial*)theMat;
	else
		return TCL_ERROR;
    }								// END Salvatore Sessa 14-Jan-2021 Mail: salvatore.sessa2@unina.it
    if (strcmp(argv[1], "HystereticSmooth") == 0) {		// BEGIN Salvatore Sessa 02-May-2022 Mail: salvatore.sessa2@unina.it
	void* theMat = OPS_HystereticSmooth();
	if (theMat != 0)
		theMaterial = (UniaxialMaterial*)theMat;
	else
		return TCL_ERROR;
	}
    if (strcmp(argv[1], "HystereticAsym") == 0) {
	void* theMat = OPS_HystereticAsym();
	if (theMat != 0)
		theMaterial = (UniaxialMaterial*)theMat;
	else
		return TCL_ERROR;
	}														// END Salvatore Sessa 02-May-2022 Mail: salvatore.sessa2@unina.it
	if (strcmp(argv[1], "DowelType") == 0) {
        void* theMat = OPS_DowelType();
        if (theMat != 0)
            theMaterial = (UniaxialMaterial*)theMat;
        else
            return TCL_ERROR;
    }
    if (strcmp(argv[1], "DuctileFracture") == 0) {
        void* theMat = OPS_DuctileFracture();
        if (theMat != 0)
            theMaterial = (UniaxialMaterial*)theMat;
        else
            return TCL_ERROR;
    }
    if (strcmp(argv[1], "Multiplier") == 0) {
        void* theMat = OPS_MultiplierMaterial();
        if (theMat != 0)
            theMaterial = (UniaxialMaterial*)theMat;
        else
            return TCL_ERROR;
    }
      // Fedeas
 #if defined(_STEEL2) || defined(OPSDEF_UNIAXIAL_FEDEAS)
    if (theMaterial == 0)
      theMaterial = TclModelBuilder_addFedeasMaterial(clientData, interp, argc, argv);
 #endif
      // Drain
      if (theMaterial == 0)
	theMaterial = TclModelBuilder_addDrainMaterial(clientData, interp, argc, argv);
      
      // SNAP
      if (theMaterial == 0)
	theMaterial = TclModelBuilder_addSnapMaterial(clientData, interp, argc, argv);
      
      // Py, Tz, Qz models
      if (theMaterial == 0)
	theMaterial = TclModelBuilder_addPyTzQzMaterial(clientData, interp, argc, argv, theDomain);
      
      // LimitState
      if (theMaterial == 0)
	theMaterial = Tcl_AddLimitStateMaterial(clientData, interp, argc, argv);
    

#if defined(OPSDEF_DAMAGE_FEDEAS)
      if (theMaterial == 0)
        theMaterial = TclModelBuilder_addDegradingMaterial(clientData, interp, argc, argv);
#endif

    if (theMaterial == 0) {
      
      //
      // maybe element in a class package already loaded
      //  loop through linked list of loaded functions comparing names & if find call it
      //
      
      UniaxialPackageCommand *matCommands = theUniaxialPackageCommands;
      bool found = false;
      while (matCommands != NULL && found == false) {
	if (strcmp(argv[1], matCommands->funcName) == 0) {
	  theMaterial = (UniaxialMaterial *)(*(matCommands->funcPtr))();
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
	
	theMaterial = Tcl_addWrapperUniaxialMaterial(matObject, clientData, interp, argc, argv);
	
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
	UniaxialPackageCommand *theMatCommand = new UniaxialPackageCommand;
	theMatCommand->funcPtr = funcPtr;
	theMatCommand->funcName = matName;	
	theMatCommand->next = theUniaxialPackageCommands;
	theUniaxialPackageCommands = theMatCommand;
	
	theMaterial = (UniaxialMaterial *)(*funcPtr)();
      }
    }
    
    //
    // if still here the element command does not exist
    //
    
    if (theMaterial == 0) {
      opserr << "WARNING could not create uniaxialMaterial " << argv[1] << endln;
      return TCL_ERROR;
    }
    
    // Now add the material to the modelBuilder
    if (OPS_addUniaxialMaterial(theMaterial) == false) {
      opserr << "WARNING could not add uniaxialMaterial to the modelbuilder\n";
      opserr << *theMaterial << endln;
      delete theMaterial; // invoke the material objects destructor, otherwise mem leak
      return TCL_ERROR;
    }

    return TCL_OK;
}
