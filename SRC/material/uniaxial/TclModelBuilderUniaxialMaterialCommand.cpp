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


#include <Elastic2Material.h>	// ZHY
#include <HardeningMaterial.h>	// MHS
#include <HardeningMaterial2.h>	// MHS
#include <Steel03.h>			// KM
#include <Concrete01WithSITC.h>		// Won Lee
#include <ECC01.h>                      // Won Lee
#include <Concrete04.h>
#include <Concrete05.h>
#include <Concrete06.h>			// LMS
#include <Concrete07.h>			// JDW
#include <HystereticBackbone.h>	// MHS
#include <EPPGapMaterial.h>		// Mackie
#include <PathIndependentMaterial.h>	// MHS
#include <BackboneMaterial.h>	// MHS
#include <FatigueMaterial.h>	// Patxi
#include <ENTMaterial.h>		// MHS
#include <BoucWenMaterial.h>	// Terje
#include <Pinching4Material.h>   // NM
#include <ShearPanelMaterial.h>  // NM
#include <BarSlipMaterial.h>     // NM
#include <Bond_SP01.h>	// JZ

#include <SteelMP.h>             //Quan & Michele
#include <SteelBRB.h>             //Quan & Michele
#include <SmoothPSConcrete.h>      //Quan & Michele
#include <SelfCenteringMaterial.h> //JAE
#include <ASD_SMA_3K.h> //LA

#include <KikuchiAikenHDR.h>
#include <KikuchiAikenLRB.h>
#include <AxialSp.h>
#include <AxialSpHD.h>

// #include <SMAMaterial.h>     // Davide Fugazza

#include <Vector.h>
#include <string.h>

#include <UniaxialJ2Plasticity.h>   // Quan 

extern void *OPS_SPSW02(void);		// SAJalali
extern void *OPS_TDConcreteEXP(void); // ntosic
extern void *OPS_TDConcrete(void); // ntosic
extern void *OPS_TDConcreteMC10(void); //ntosic
extern void *OPS_TDConcreteMC10NL(void); //ntosic
extern void *OPS_ElasticMaterial(void);
extern void *OPS_ElasticPPMaterial(void);
extern void *OPS_EPPGapMaterial(void);
extern void *OPS_ParallelMaterial(void);
extern void *OPS_SeriesMaterial(void);
extern void *OPS_HardeningMaterial(void);
extern void *OPS_HystereticMaterial(void);
extern void *OPS_CableMaterial(void);
extern void *OPS_Bilin(void);
extern void *OPS_Bilin02(void);
extern void *OPS_Steel01(void);
extern void *OPS_FRPConfinedConcrete02(void);
extern void *OPS_Steel02(void);
extern void *OPS_Steel02Fatigue(void);
extern void *OPS_RambergOsgoodSteel(void);
extern void *OPS_ReinforcingSteel(void);
extern void *OPS_Concrete01(void);
extern void *OPS_Concrete02(void);
extern void *OPS_Concrete02IS(void);
extern void *OPS_PinchingLimitStateMaterial(void);
extern void *OPS_SAWSMaterial(void);
extern void *OPS_ConcreteZ01Material(void);
extern void *OPS_ConcreteL01Material(void);
extern void *OPS_SteelZ01Material(void);
extern void *OPS_TendonL01Material(void);
extern void *OPS_ConfinedConcrete01Material(void);
extern void *OPS_ElasticBilin(void);
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
extern void *OPS_Steel4(void);
extern void *OPS_PySimple3(void);
extern void *OPS_BoucWenOriginal(void);
extern void *OPS_GNGMaterial(void);
extern void *OPS_OOHystereticMaterial(void);
extern void *OPS_ElasticPowerFunc(void);
extern void *OPS_UVCuniaxial(void);
extern void *OPS_DegradingPinchedBW(void);
extern void *OPS_SLModel(void);
extern void* OPS_HystereticPoly(void); // Salvatore Sessa 14-Jan-2021 Mail: salvatore.sessa2@unina.it

//extern int TclCommand_ConfinedConcrete02(ClientData clientData, Tcl_Interp *interp, int argc, 
//					 TCL_Char **argv, TclModelBuilder *theTclBuilder);

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

int
TclCommand_KikuchiAikenHDR(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int
TclCommand_KikuchiAikenLRB(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int
TclCommand_AxialSp(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int
TclCommand_AxialSpHD(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);


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
	else if (strcmp(argv[1], "SPSW02") == 0) {
		void *theMat = OPS_SPSW02();
		if (theMat != 0)
			theMaterial = (UniaxialMaterial *)theMat;
		else
			return TCL_ERROR;
	}

	// ntosic
	else if (strcmp(argv[1], "TDConcreteEXP") == 0) {
		void *theMat = OPS_TDConcreteEXP();
		if (theMat != 0)
			theMaterial = (UniaxialMaterial *)theMat;
		else
			return TCL_ERROR;
	}

	// ntosic
	else if (strcmp(argv[1], "TDConcrete") == 0) {
		void *theMat = OPS_TDConcrete();
		if (theMat != 0)
			theMaterial = (UniaxialMaterial *)theMat;
		else
			return TCL_ERROR;
	}

	// ntosic
	else if (strcmp(argv[1], "TDConcreteMC10") == 0) {
		void *theMat = OPS_TDConcreteMC10();
		if (theMat != 0)
			theMaterial = (UniaxialMaterial *)theMat;
		else
			return TCL_ERROR;
	}

	// ntosic
	else if (strcmp(argv[1], "TDConcreteMC10NL") == 0) {
		void *theMat = OPS_TDConcreteMC10NL();
		if (theMat != 0)
			theMaterial = (UniaxialMaterial *)theMat;
		else
			return TCL_ERROR;
	}

	else if (strcmp(argv[1],"Steel01") == 0) {

      void *theMat = OPS_Steel01();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"Steel02") == 0) {
      void *theMat = OPS_Steel02();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    } else if (strcmp(argv[1],"Steel02Fatigue") == 0) {
      void *theMat = OPS_Steel02Fatigue();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"Steel4") == 0) {
      void *theMat = OPS_Steel4();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"UVCuniaxial") == 0) {
      void *theMat = OPS_UVCuniaxial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"GNG") == 0) {
      void *theMat = OPS_GNGMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"PySimple3") == 0) {
      void *theMat = OPS_PySimple3();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"Concrete01") == 0) {
      void *theMat = OPS_Concrete01();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
      /*
	    } else if (strcmp(argv[1],"HoehlerStanton") == 0) {
      void *theMat = OPS_HoehlerStanton();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
      */
    } else if (strcmp(argv[1],"Concrete02") == 0) {
      void *theMat = OPS_Concrete02();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    } else if (strcmp(argv[1],"Concrete02IS") == 0) {
      void *theMat = OPS_Concrete02IS();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if ((strcmp(argv[1],"ElasticBilin") == 0) || (strcmp(argv[1],"ElasticBilinear") == 0)) {
      void *theMat = OPS_ElasticBilin();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if ((strcmp(argv[1],"ImpactMaterial") == 0) || (strcmp(argv[1],"Impact") == 0)) {
      void *theMat = OPS_ImpactMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if ((strcmp(argv[1],"SteelBRB") == 0)) {
      void *theMat = OPS_SteelBRB();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if ((strcmp(argv[1],"MinMaxMaterial") == 0) || (strcmp(argv[1],"MinMax") == 0)) {
      void *theMat = OPS_MinMaxMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if ((strcmp(argv[1],"SimpleFractureMaterial") == 0) || (strcmp(argv[1],"SimpleFracture") == 0)) {

      void *theMat = OPS_SimpleFractureMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if ((strcmp(argv[1],"Maxwell") == 0) || (strcmp(argv[1],"MaxwellMaterial") == 0)) {
      void *theMat = OPS_Maxwell();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if ((strcmp(argv[1],"ViscousDamper") == 0)) {
      void *theMat = OPS_ViscousDamper();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if ((strcmp(argv[1],"DamperMaterial") == 0)) {
      void *theMat = OPS_DamperMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if ((strcmp(argv[1],"BilinearOilDamper") == 0)) {
      void *theMat = OPS_BilinearOilDamper();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if ((strcmp(argv[1],"Bond_SP01") == 0) || (strcmp(argv[1],"Bond") == 0)) { 
      void *theMat = OPS_Bond_SP01();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;      

    } else if ((strcmp(argv[1],"Cast") == 0) || (strcmp(argv[1],"CastFuse") == 0)) {
      void *theMat = OPS_Cast();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if ((strcmp(argv[1],"Dodd_Restrepo") == 0) || 
	       (strcmp(argv[1],"DoddRestrepo") == 0) || 
	       (strcmp(argv[1],"Restrepo") == 0)) {

      void *theMat = OPS_Dodd_Restrepo();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

#ifndef _NO_NEW_RESTREPO
	} else if ((strcmp(argv[1],"DoddRestr") == 0)) {

      void *theMat = OPS_DoddRestr();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
#endif
    } else if (strcmp(argv[1],"ElasticMultiLinear") == 0) {
      void *theMat = OPS_ElasticMultiLinear();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    else if (strcmp(argv[1], "ElasticPowerFunc") == 0) {
      void* theMat = OPS_ElasticPowerFunc();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial*)theMat;
      else
        return TCL_ERROR;
    }
    /*
	} else if (strcmp(argv[1],"HoehlerStanton") == 0) {
      void *theMat = OPS_HoehlerStanton();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
*/
    else if (strcmp(argv[1], "UVCuniaxial") == 0) {
      void *theMat = OPS_UVCuniaxial();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;
    }
    else if (strcmp(argv[1], "SLModel") == 0) {
      void *theMat = OPS_SLModel();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;
    }
    else if ((strcmp(argv[1],"RambergOsgood") == 0) || (strcmp(argv[1],"RambergOsgoodSteel") == 0)) {
      void *theMat = OPS_RambergOsgoodSteel();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"ReinforcingSteel") == 0) {
      void *theMat = OPS_ReinforcingSteel();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"Steel2") == 0) {
      void *theMat = OPS_Steel2();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"OriginCentered") == 0) {
      void *theMat = OPS_OriginCentered();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"HookGap") == 0) {
      void *theMat = OPS_HookGap();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"HyperbolicGapMaterial") == 0) {
      void *theMat = OPS_HyperbolicGapMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"FRPConfinedConcrete02") == 0) {
      void *theMat = OPS_FRPConfinedConcrete02();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if ((strcmp(argv[1],"PinchingLimitState") == 0) || (strcmp(argv[1],"PinchingLimitStateMaterial") == 0)) {
      void *theMat = OPS_PinchingLimitState();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if ((strcmp(argv[1],"InitStrainMaterial") == 0) || (strcmp(argv[1],"InitStrain") == 0)) {
      void *theMat = OPS_InitStrainMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if ((strcmp(argv[1],"InitStressMaterial") == 0) || (strcmp(argv[1],"InitStress") == 0)) {
      void *theMat = OPS_InitStressMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if ((strcmp(argv[1],"pyUCLA") == 0) || (strcmp(argv[1],"PYUCLA") == 0)) {
      void *theMat = OPS_pyUCLA();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if ((strcmp(argv[1],"MultiLinear") == 0)) {
      void *theMat = OPS_MultiLinear();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"ModIMKPinching") == 0) {
      void *theMat = OPS_ModIMKPinching();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"ModIMKPinching02") == 0) {
      void *theMat = OPS_ModIMKPinching02();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1], "BoucWenOriginal") == 0) {
        void *theMat = OPS_BoucWenOriginal();
        if (theMat != 0)
            theMaterial = (UniaxialMaterial *)theMat;
        else
            return TCL_ERROR;

    } else if (strcmp(argv[1], "BWBN") == 0) {
        void *theMat = OPS_BWBN();
        if (theMat != 0)
            theMaterial = (UniaxialMaterial *)theMat;
        else
            return TCL_ERROR;

    } else if (strcmp(argv[1], "DegradingPinchedBW") == 0) {
        void *theMat = OPS_DegradingPinchedBW();
        if (theMat != 0)
            theMaterial = (UniaxialMaterial *)theMat;
        else
            return TCL_ERROR;

    }
    else if (strcmp(argv[1], "IMKBilin") == 0) {
      void *theMat = OPS_IMKBilin();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;

    }
    else if (strcmp(argv[1], "IMKPeakOriented") == 0) {
      void *theMat = OPS_IMKPeakOriented();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;

    }
    else if (strcmp(argv[1], "IMKPinching") == 0) {
      void *theMat = OPS_IMKPinching();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;

    }
    else if (strcmp(argv[1], "ModIMKPeakOriented") == 0) {
      void *theMat = OPS_ModIMKPeakOriented();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"ModIMKPeakOriented02") == 0) {
      void *theMat = OPS_ModIMKPeakOriented02();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"Bilin02") == 0) {
      void *theMat = OPS_Bilin02();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"PlateBearingConnectionThermal") == 0) {
      //void *theMat = OPS_PlateBearingConnectionThermal();
      void *theMat = 0;
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"Steel01Thermal") == 0) {
      void *theMat = OPS_Steel01Thermal();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;      

    } else if (strcmp(argv[1],"Steel02Thermal") == 0) {
      void *theMat = OPS_Steel02Thermal();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

	  // More thermo-mechanical uniaxial materials, L.Jiang[SIF]
	}
	else if (strcmp(argv[1], "SteelECThermal") == 0) {
		void *theMat = OPS_SteelECThermal();
		if (theMat != 0)
			theMaterial = (UniaxialMaterial *)theMat;
		else
			return TCL_ERROR;
		//------End of adding identity for SteelEcThermal	
	}
	else if (strcmp(argv[1], "StainlessECThermal") == 0) {
		void *theMat = OPS_StainlessECThermal();
		if (theMat != 0)
			theMaterial = (UniaxialMaterial *)theMat;
		else
			return TCL_ERROR;
		//------End of adding identity for StainlessECThermal
	}
	else if (strcmp(argv[1], "ElasticThermal") == 0) {
		void *theMat = OPS_ElasticMaterialThermal();
		if (theMat != 0)
			theMaterial = (UniaxialMaterial *)theMat;
		else
			return TCL_ERROR;

	} else if (strcmp(argv[1], "ConcreteECThermal") == 0) {
		void *theMat = OPS_ConcreteECThermal();
		if (theMat != 0)
			theMaterial = (UniaxialMaterial *)theMat;
		else
			return TCL_ERROR;
		// end of adding More thermo-mechanical uniaxial materials, L.Jiang[SIF]

    } else if (strcmp(argv[1],"ConcretewBeta") == 0) {
      void *theMat = OPS_ConcretewBeta();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"ConcreteD") == 0) {
      void *theMat = OPS_ConcreteD();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"ConcreteSakaiKawashima") == 0) {
      void *theMat = OPS_ConcreteSakaiKawashima();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"Concrete02Thermal") == 0) {
      void *theMat = OPS_Concrete02Thermal();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"SteelMPF") == 0) { // K Kolozvari            
      void *theMat = OPS_SteelMPF();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;
      
    } else if (strcmp(argv[1],"ConcreteCM") == 0) { // K Kolozvari          
      void *theMat = OPS_ConcreteCM();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;
      
    } else if (strcmp(argv[1],"ResilienceLow") == 0) {
      void *theMat = OPS_ResilienceLow();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"ResilienceMaterialHR") == 0) {
      void *theMat = OPS_ResilienceMaterialHR();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"CFSWSWP") == 0) {
      void *theMat = OPS_CFSWSWP();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"CFSSSWP") == 0) {
      void *theMat = OPS_CFSSSWP();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"FRPConfinedConcrete") == 0) {
      void *theMat = OPS_FRPConfinedConcrete();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    } else if (strcmp(argv[1],"Elastic2") == 0) {
	if (argc < 4 || argc > 5) {
	    opserr << "WARNING invalid number of arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: uniaxialMaterial Elastic tag? E? <eta?>" << endln;
	    return TCL_ERROR;
	}    

	int tag;
	double E;
        double eta = 0.0;
	
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid uniaxialMaterial Elastic tag" << endln;
	    return TCL_ERROR;		
	}

	if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
	    opserr << "WARNING invalid E\n";
	    opserr << "uniaxiaMaterial Elastic: " << tag << endln;
	    return TCL_ERROR;	
	}

	if (argc == 5) {
	    if (Tcl_GetDouble(interp,argv[4], &eta) != TCL_OK) {
               opserr << "WARNING invalid eta\n";
               opserr << "uniaxialMaterial Elastic: " << tag << endln;
               return TCL_ERROR;
            }
	}

	// Parsing was successful, allocate the material
	theMaterial = new Elastic2Material(tag, E, eta);       
      
    } else if (strcmp(argv[1],"ENT") == 0) {
      if (argc < 4) {
	opserr << "WARNING invalid number of arguments\n";
	printCommand(argc,argv);
	opserr << "Want: uniaxialMaterial ENT tag? E?" << endln;
	return TCL_ERROR;
      }    
      
      int tag;
      double E;
      
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING invalid uniaxialMaterial ENT tag" << endln;
	return TCL_ERROR;		
      }
      
      if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
	opserr << "WARNING invalid E\n";
	opserr << "uniaxiaMaterial ENT: " << tag << endln;
	return TCL_ERROR;	
      }
      
      // Parsing was successful, allocate the material
      theMaterial = new ENTMaterial(tag, E);       
      
    }
    
    else if (strcmp(argv[1],"ElasticPP") == 0) {
      void *theMat = OPS_ElasticPPMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    
    else if (strcmp(argv[1],"ElasticPPGap") == 0) {
      void *theMat = OPS_EPPGapMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }

    else if (strcmp(argv[1],"Hardening") == 0 || strcmp(argv[1],"Hardening2") == 0) {
      
      void *theMat = OPS_HardeningMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }

    else if (strcmp(argv[1],"BoucWen") == 0) {
      if (argc < 12) {
	opserr << "WARNING insufficient arguments\n";
	printCommand(argc,argv);
	opserr << "Want: uniaxialMaterial BoucWen tag? alpha? ko? n? gamma?" << endln 
	       << " beta? Ao? deltaA? deltaNu? deltaEta?" << endln;
	return TCL_ERROR;
      }
      
      int tag;
      double alpha, ko, n, gamma, beta, Ao, deltaA, deltaNu, deltaEta;
      
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING invalid uniaxialMaterial BoucWen tag" << endln;
	return TCL_ERROR;		
      }
      if (Tcl_GetDouble(interp, argv[3], &alpha) != TCL_OK) {
	opserr << "WARNING invalid alpha\n";
	opserr << "uniaxialMaterial BoucWen: " << tag << endln;
	return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[4], &ko) != TCL_OK) {
	opserr << "WARNING invalid ko\n";
	opserr << "uniaxialMaterial BoucWen: " << tag << endln;
	return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[5], &n) != TCL_OK) {
	opserr << "WARNING invalid n\n";
	opserr << "uniaxialMaterial BoucWen: " << tag << endln;
	return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[6], &gamma) != TCL_OK) {
	opserr << "WARNING invalid gamma\n";
	opserr << "uniaxialMaterial BoucWen: " << tag << endln;
	return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[7], &beta) != TCL_OK) {
	opserr << "WARNING invalid beta\n";
	opserr << "uniaxialMaterial BoucWen: " << tag << endln;
	return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[8], &Ao) != TCL_OK) {
	opserr << "WARNING invalid Ao\n";
	opserr << "uniaxialMaterial BoucWen: " << tag << endln;
	return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[9], &deltaA) != TCL_OK) {
	opserr << "WARNING invalid deltaA\n";
	opserr << "uniaxialMaterial BoucWen: " << tag << endln;
	return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[10], &deltaNu) != TCL_OK) {
	opserr << "WARNING invalid deltaNu\n";
	opserr << "uniaxialMaterial BoucWen: " << tag << endln;
	return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[11], &deltaEta) != TCL_OK) {
	opserr << "WARNING invalid deltaEta\n";
	opserr << "uniaxialMaterial BoucWen: " << tag << endln;
	return TCL_ERROR;	
      }
      
      // Check if the user has given a tolerance for the Newton scheme
      double tolerance = 1.0e-8;
      if (argc > 12) {
	if (Tcl_GetDouble(interp, argv[12], &tolerance) != TCL_OK) {
	  opserr << "WARNING invalid tolerance\n";
	  opserr << "uniaxialMaterial BoucWen: " << tolerance << endln;
	  return TCL_ERROR;	
	}
      }
      
      // Check if the user has given a maxNumIter for the Newton scheme
      int maxNumIter = 20;
      if (argc > 13) {
	if (Tcl_GetInt(interp, argv[13], &maxNumIter) != TCL_OK) {
	  opserr << "WARNING invalid maxNumIter\n";
	  opserr << "uniaxialMaterial BoucWen: " << maxNumIter << endln;
	  return TCL_ERROR;	
	}
      }
      
      // Parsing was successful, allocate the material
      theMaterial = new BoucWenMaterial(tag, alpha, ko, n, gamma, beta, 
					Ao, deltaA, deltaNu, deltaEta,tolerance,maxNumIter);       
    }

    else if (strcmp(argv[1],"Parallel") == 0) {

      void *theMat = OPS_ParallelMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }

    else if (strcmp(argv[1],"Series") == 0) {

      void *theMat = OPS_SeriesMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }


    else if (strcmp(argv[1],"Steel03") == 0) {
      // Check that there is the minimum number of arguments
      if (argc < 9) {
        opserr << "WARNING insufficient arguments\n";
        printCommand(argc,argv);
        opserr << "Want: uniaxialMaterial Steel03 tag? fy? E0? b? r? cR1 cR2?";
        opserr << " <a1? a2? a3? a4?>" << endln;    
        return TCL_ERROR;
      }
      
      int tag;
      
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
        opserr << "WARNING invalid uniaxialMaterial Steel03 tag" << endln;
        return TCL_ERROR;
      }
      
      
      // Read required Steel01 material parameters
      double fy, E, b, r, cR1, cR2;
      
      if (Tcl_GetDouble(interp, argv[3], &fy) != TCL_OK) {
        opserr << "WARNING invalid fy\n";
        opserr << "uniaxialMaterial Steel03: " << tag << endln;
        return TCL_ERROR;
      }
      
      
      if (Tcl_GetDouble(interp, argv[4], &E) != TCL_OK) {
        opserr << "WARNING invalid E0\n";
        opserr << "uniaxialMaterial Steel03: " << tag << endln;
        return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[5], &b) != TCL_OK) {
        opserr << "WARNING invalid b\n";
        opserr << "uniaxialMaterial Steel03: " << tag << endln;
        return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[6], &r) != TCL_OK) {
	opserr << "WARNING invalid r\n";
	opserr << "uniaxialMaterial Steel03: " << tag << endln;
	return TCL_ERROR;
      }
       

      if (Tcl_GetDouble(interp, argv[7], &cR1) != TCL_OK) {
	opserr << "WARNING invalid cR1\n";
	opserr << "uniaxialMaterial Steel03: " << tag << endln;
	return TCL_ERROR;
      }
      
      
      if (Tcl_GetDouble(interp, argv[8], &cR2) != TCL_OK) {
	opserr << "WARNING invalid cR2\n";
	opserr << "uniaxialMaterial Steel03: " << tag << endln;
	return TCL_ERROR;
      }
      
      // Read optional Steel01 material parameters
      double a1, a2, a3, a4;
      if (argc > 9) {
	if (argc < 13) {
	  opserr << "WARNING insufficient number of hardening parameters\n";
	  opserr << "uniaxialMaterial Steel03: " << tag << endln;
	  return TCL_ERROR;
	}
	
	if (Tcl_GetDouble(interp, argv[9], &a1) != TCL_OK) {
	  opserr << "WARNING invalid a1\n";
	  opserr << "uniaxialMaterial Steel03: " << tag << endln;
	  return TCL_ERROR;
	}
	
	if (Tcl_GetDouble(interp, argv[10], &a2) != TCL_OK) {
	  opserr << "WARNING invalid a2\n";
	  opserr << "uniaxialMaterial Steel03: " << tag << endln;
	  return TCL_ERROR;
	}
	
	if (Tcl_GetDouble(interp, argv[11], &a3) != TCL_OK) {
	  opserr << "WARNING invalid a3\n";
	  opserr << "uniaxialMaterial Steel03: " << tag << endln;
	  return TCL_ERROR;
	}
	
	
	if (Tcl_GetDouble(interp, argv[12], &a4) != TCL_OK) {
	  opserr << "WARNING invalid a4\n";
	  opserr << "uniaxialMaterial Steel03: " << tag << endln;
	  return TCL_ERROR;
	}
	
	
        // Parsing was successful, allocate the material
	theMaterial = new Steel03 (tag, fy, E, b, r, cR1, cR2, a1, a2, a3, a4);
      }
      else
	// Parsing was successful, allocate the material
	theMaterial = new Steel03 (tag, fy, E, b, r, cR1, cR2);
      
    }
    
    
    else if (strcmp(argv[1],"Hysteretic") == 0) {

      void *theMat = OPS_HystereticMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }

    else if (strcmp(argv[1],"OOHysteretic") == 0) {

      void *theMat = OPS_OOHystereticMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    
    else if (strcmp(argv[1],"Concrete04") == 0) {
      //        opserr << argc << endln;
      if (argc != 10 && argc != 9 && argc != 7) {
	opserr << "WARNING insufficient arguments\n";
	printCommand(argc,argv);
	opserr << "Want: uniaxialMaterial Concrete04 tag? fpc? epsc0? epscu? Ec0? <ft? etu? <beta?> >" << endln;
	return TCL_ERROR;
      }
      
      int tag;
      
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING invalid uniaxialMaterial Concrete04 tag"  
	       << endln;
	return TCL_ERROR;
      }
      
      // Read required Concrete04 material parameters
      double fpc, epsc0, ft, epscu, Ec0, etu, beta;
      
      if (Tcl_GetDouble(interp, argv[3], &fpc) != TCL_OK) {
	opserr << "WARNING invalid fpc\n";
	opserr << "Concrete04 material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[4], &epsc0) != TCL_OK) {
	opserr << "WARNING invalid epsc0\n";
	opserr << "Concrete04 material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[5], &epscu) != TCL_OK) {
	opserr << "WARNING invalid epscu\n";
	opserr << "Concrete04 material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[6], &Ec0) != TCL_OK) {
	opserr << "WARNING invalid Ec0\n";
	opserr << "Concrete04 material: " << tag << endln;
	return TCL_ERROR;
      }
      if (argc == 9 || argc == 10) {
	if (Tcl_GetDouble(interp, argv[7], &ft) != TCL_OK) {
	  opserr << "WARNING invalid ft\n";
	  opserr << "Concrete04 material: " << tag << endln;
	  return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[8], &etu) != TCL_OK) {
	  opserr << "WARNING invalid etu\n";
	  opserr << "Concrete04 material: " << tag << endln;
	  return TCL_ERROR;
	}
      }
      if (argc == 10) {
	if (Tcl_GetDouble(interp, argv[9], &beta) != TCL_OK) {
	  opserr << "WARNING invalid beta\n";
	  opserr << "Concrete04 material: " << tag << endln;
	  return TCL_ERROR;
	}
      }
      
      
      // Parsing was successful, allocate the material
      if (argc == 10) {
	theMaterial = new Concrete04(tag, fpc, epsc0, epscu, Ec0,  
				     ft, etu, beta);
      }
      else if (argc == 9) {
	theMaterial = new Concrete04(tag, fpc, epsc0, epscu, Ec0,  
				     ft, etu);
      }
      else if (argc == 7) {
	theMaterial = new Concrete04(tag, fpc, epsc0, epscu, Ec0);
      }
    }



    else if (strcmp(argv[1],"Concrete06") == 0) {
	if (argc < 12) {
	    opserr << "WARNING insufficient arguments\n";
	    printCommand(argc,argv);
	    opserr << "Want: uniaxialMaterial Concrete06 tag? fc? eo? r? k? alphaC? fcr? ecr? b? alphaT?" << endln;
	    return TCL_ERROR;
	}

	int tag;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid uniaxialMaterial Concrete06 tag" << endln;
	    return TCL_ERROR;
	}

	// Read required Concrete01 material parameters
	double fc, eo, r, k, fcr, ecr, b, alphaC, alphaT;

	if (Tcl_GetDouble(interp, argv[3], &fc) != TCL_OK) {
	    opserr << "WARNING invalid fc\n";
	    opserr << "Concrete06 material: " << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[4], &eo) != TCL_OK) {
	    opserr << "WARNING invalid eo\n";
	    opserr << "Concrete06 material: " << tag << endln;
	    return TCL_ERROR;
	}
	
	if (Tcl_GetDouble(interp, argv[5], &r) != TCL_OK) {
	    opserr << "WARNING invalid r\n";
	    opserr << "Concrete06 material: " << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[6], &k) != TCL_OK) {
	    opserr << "WARNING invalid k\n";
	    opserr << "Concrete06 material: " << tag << endln;
	    return TCL_ERROR;
	}
	
	if (Tcl_GetDouble(interp, argv[7], &alphaC) != TCL_OK) {
	  opserr << "WARNING invalid alphaC\n";
	  opserr << "Concrete06 material: " << tag << endln;
	  return TCL_ERROR;
	}
	
	if (Tcl_GetDouble(interp, argv[8], &fcr) != TCL_OK) {
	  opserr << "WARNING invalid fcr\n";
	  opserr << "Concrete06 material: " << tag << endln;
	  return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[9], &ecr) != TCL_OK) {
	  opserr << "WARNING invalid ecr\n";
	  opserr << "Concrete06 material: " << tag << endln;
	  return TCL_ERROR;
	}
	
	if (Tcl_GetDouble(interp, argv[10], &b) != TCL_OK) {
	  opserr << "WARNING invalid b\n";
	  opserr << "Concrete06 material: " << tag << endln;
	  return TCL_ERROR;
	}
	
	if (Tcl_GetDouble(interp, argv[11], &alphaT) != TCL_OK) {
	  opserr << "WARNING invalid alphaT\n";
	  opserr << "Concrete06 material: " << tag << endln;
	  return TCL_ERROR;
	}
	
	// Parsing was successful, allocate the material
	theMaterial = new Concrete06(tag, fc, eo, r, k, alphaC, fcr, ecr, b, alphaT);
    }

    
    
    else if (strcmp(argv[1], "Concrete07") == 0) {
      // Check to see if there are enough arquements
      if (argc < 11) {
	opserr << "WARNING: Insufficient arguments\n";
	printCommand(argc, argv);
	opserr << "Want: uniaxialMaterial Concrete07 tag? fpc? epsc0? Ec? fpt? epst0? xcrp? xcrn? r?\n";
	return TCL_ERROR;
      }
      
      int tag;

      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING: Invalid uniaxial Concrete07 tag\n";
	return TCL_ERROR;
      }
      
      // Read in the faluves required for the model
      double fpc, epsc0, Ec, fpt, epst0, xcrp, xcrn, r;
      
      if (Tcl_GetDouble(interp, argv[3], &fpc) != TCL_OK) {
	opserr << "WARNING: Invalid peak compression stress\n";
	opserr << "uniaxialMaterial Concrete07: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[4], &epsc0) != TCL_OK) {
	opserr << "WARNING: Invalid peak compression strain\n";
	opserr << "uniaxialMaterial Concrete07: " << tag <<endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[5], &Ec) != TCL_OK) {
	opserr << "WARNING: Invalid Young's Modulus\n";
	opserr << "uniaxialMaterial Concrete07: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[6], &fpt) != TCL_OK) {
	opserr << "WARNING: Invalid peak tension stress\n";
	opserr << "uniaxialMaterial Concrete07: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[7], &epst0) != TCL_OK) {
	opserr << "WARNING: Invalid peak tension strain\n";
	opserr << "uniaxialMaterial Concrete07: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[8], &xcrp) != TCL_OK) {
	opserr << "WARNING: Invalid critical nondimensional strain in tension\n";
	opserr << "uniaxialMaterial Concrete07: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[9], &xcrn) != TCL_OK) {
	opserr << "WARNING: Invalid critical nondimensional strain in compression\n";
	opserr << "uniaxialMaterial Concrete07: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[10], &r) != TCL_OK) {
	opserr << "WARNING: Invalid value for r\n";
	opserr << "uniaxialMaterial Concrete07: " << tag << endln;
      }
      
      //		opserr << "fpc: " << fpc << endln << "epsc0: " << epsc0 << endln << "Ec: " << Ec << endln;
      //		opserr << "fpt: " << fpt << endln << "epst0: " << epst0 << endln << "xcrp: " << xcrp << endln;
      //		opserr << "xcrn: " << xcrn << endln << "r: " << r << endln;
      
      // Parsing was successful, allocate the material
      
      theMaterial = new Concrete07(tag, fpc, epsc0, Ec, fpt, epst0, xcrp, xcrn, r);
    }

    
    else if (strcmp(argv[1],"Viscous") == 0) {

      void *theMat = OPS_ViscousMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }
    
    else if (strcmp(argv[1],"PathIndependent") == 0) {
		if (argc < 4)
		{
			opserr << "WARNING insufficient arguments\n";
			printCommand(argc,argv);
			opserr << "Want: uniaxialMaterial PathIndependent tag? matTag?" << endln;
			return TCL_ERROR;
		}

		int tag, matTag;

		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
		{
			opserr << "WARNING invalid tag\n";
			opserr << "PathIndependent material: " << tag << endln;
			return TCL_ERROR;
		}
		
		if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK)
		{
			opserr << "WARNING invalid matTag\n";
			opserr << "PathIndependent material: " << tag << endln;
			return TCL_ERROR;
		}

		UniaxialMaterial *material = OPS_getUniaxialMaterial(matTag);
		
		if (material == 0) {
		    opserr << "WARNING material does not exist\n";
		    opserr << "material: " << matTag; 
		    opserr << "\nuniaxialMaterial PathIndependent: " << tag << endln;
		    return TCL_ERROR;
		}

		theMaterial = new PathIndependentMaterial (tag, *material);
	}

    else if (strcmp(argv[1],"Backbone") == 0) {
      if (argc < 4) {
	opserr << "WARNING insufficient arguments\n";
	printCommand(argc,argv);
	opserr << "Want: uniaxialMaterial Backbone tag? bbTag?" << endln;
	return TCL_ERROR;
      }
      
      int tag, bbTag;

      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING invalid tag\n";
	opserr << "Backbone material: " << tag << endln;
	return TCL_ERROR;
      }
		
      if (Tcl_GetInt(interp, argv[3], &bbTag) != TCL_OK) {
	opserr << "WARNING invalid bTag\n";
	opserr << "Backbone material: " << tag << endln;
	return TCL_ERROR;
      }

      HystereticBackbone *backbone = OPS_getHystereticBackbone(bbTag);
		
      if (backbone == 0) {
	opserr << "WARNING backbone does not exist\n";
	opserr << "backbone: " << bbTag; 
	opserr << "\nuniaxialMaterial Backbone: " << tag << endln;
	return TCL_ERROR;
      }
      
      theMaterial = new BackboneMaterial(tag, *backbone);
    }

    else if (strcmp(argv[1],"Fatigue") == 0) {
      if (argc < 4) {
	opserr << "WARNING insufficient arguments\n";
	printCommand(argc,argv);
	opserr << "Want: uniaxialMaterial Fatigue tag? matTag?";
	opserr << " <-D_max dmax?> <-e0 e0?> <-m m?>" << endln;
	opserr << " <-min min?> <-max max?>" << endln;
	return TCL_ERROR;
      }
      
      int tag, matTag;
      
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING invalid uniaxialMaterial Fatigue tag" << endln;
	return TCL_ERROR;		
      }
      
      if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
	opserr << "WARNING invalid component tag\n";
	opserr << "uniaxialMaterial Fatigue: " << tag << endln;
	return TCL_ERROR;
      }
      
      double Dmax      =  1.0;
      double E0        =  0.191;
      double m         = -0.458;
      double epsmin    = NEG_INF_STRAIN;
      double epsmax    = POS_INF_STRAIN;
      
      for (int j = 4; j < argc; j++) {
	if (strcmp(argv[j],"-Dmax") == 0) {
	  if ((j+1 >= argc) || 
	      (Tcl_GetDouble (interp, argv[j+1], &Dmax) != TCL_OK)) {
	    opserr << "WARNING invalid -Dmax";
	    opserr << "uniaxialMaterial Fatigue: " << tag << endln;
	    return TCL_ERROR;
	  }
	} else if (strcmp(argv[j],"-E0") == 0) {
	  if ((j+1 >= argc) || 
	      (Tcl_GetDouble (interp, argv[j+1], &E0) != TCL_OK)) {
	    opserr << "WARNING invalid -E0";	 
	    opserr << "uniaxialMaterial Fatigue: " << tag << endln;
	    return TCL_ERROR;
	  }
	} else if (strcmp(argv[j],"-m") == 0) {
	  if ((j+1 >= argc) || 
	      (Tcl_GetDouble (interp, argv[j+1], &m) != TCL_OK)) {
	    opserr << "WARNING invalid -m"; 
	    opserr << "uniaxialMaterial Fatigue: " << tag << endln;
	    return TCL_ERROR;
	  }
	} else if (strcmp(argv[j],"-min") == 0) {
	  if ((j+1 >= argc) || 
	      (Tcl_GetDouble (interp, argv[j+1], &epsmin) != TCL_OK)) {
	    opserr << "WARNING invalid -min ";
	    opserr << "uniaxialMaterial Fatigue: " << tag << endln;
	    return TCL_ERROR;
	  }
	} else if (strcmp(argv[j],"-max") == 0) {
	  if ((j+1 >= argc) || 
	      (Tcl_GetDouble (interp, argv[j+1], &epsmax) != TCL_OK)) {
	    opserr << "WARNING invalid -max";
	    opserr << "uniaxialMaterial Fatigue: " << tag << endln;
	    return TCL_ERROR;
	  }
	}
	j++;
      }
      
      
      UniaxialMaterial *theMat = OPS_getUniaxialMaterial(matTag);
      
      if (theMat == 0) {
	opserr << "WARNING component material does not exist\n";
	opserr << "Component material: " << matTag; 
	opserr << "\nuniaxialMaterial Fatigue: " << tag << endln;
	return TCL_ERROR;
      }
      
      // Parsing was successful, allocate the material
      theMaterial = new FatigueMaterial(tag, *theMat, Dmax, E0, 
					m, epsmin, epsmax);
      
    }

    else if ((strcmp(argv[1],"SAWSMaterial") == 0) || (strcmp(argv[1],"SAWS") == 0)) {

      void *theMat = OPS_SAWSMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }

    else if ((strcmp(argv[1],"BilinMaterial") == 0) || (strcmp(argv[1],"Bilin") == 0)) {

      void *theMat = OPS_Bilin();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }

    else if ((strcmp(argv[1],"ConcreteZ01Material") == 0) || (strcmp(argv[1],"ConcreteZ01") == 0)) {

      void *theMat = OPS_ConcreteZ01Material();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }

    else if ((strcmp(argv[1],"ConcreteL01Material") == 0) || (strcmp(argv[1],"ConcreteL01") == 0)) {

      void *theMat = OPS_ConcreteL01Material();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }

    else if ((strcmp(argv[1],"SteelZ01Material") == 0) || (strcmp(argv[1],"SteelZ01") == 0)) {

      void *theMat = OPS_SteelZ01Material();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }

    else if ((strcmp(argv[1],"TendonL01Material") == 0) || (strcmp(argv[1],"TendonL01") == 0)) {
      void *theMat = OPS_TendonL01Material();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }


    else if ((strcmp(argv[1],"ConfinedConcrete01") == 0) || (strcmp(argv[1],"ConfinedConcrete") == 0)) {

      void *theMat = OPS_ConfinedConcrete01Material();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;
    }

    else if (strcmp(argv[1],"Cable") == 0) {

      void *theMat = OPS_CableMaterial();
      if (theMat != 0) 
	theMaterial = (UniaxialMaterial *)theMat;
      else 
	return TCL_ERROR;

    }
    
    else if (strcmp(argv[1], "UVCuniaxial") == 0) {
        void* theMat = OPS_UVCuniaxial();
    if (theMat != 0)
        theMaterial = (UniaxialMaterial*)theMat;
    else
        return TCL_ERROR;
    }
    
    else if (strcmp(argv[1],"Pinching4") == 0) {
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
   
  else if (strcmp(argv[1],"BarSlip") == 0)
   {
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
    
    
    
  else if (strcmp(argv[1],"ShearPanel") == 0) {
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
    
    else if (strcmp(argv[1],"Concrete01WithSITC") == 0) {
      if (argc < 7) {
	  opserr << "WARNING insufficient arguments\n";
	  printCommand(argc,argv);
	  opserr << "Want: uniaxialMaterial Concrete01 tag? fpc? epsc0? fpcu? epscu? <endStrainSITC?>" << endln;
	  return TCL_ERROR;
	}

	int tag;

	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	    opserr << "WARNING invalid uniaxialMaterial Concrete01 tag" << endln;
	    return TCL_ERROR;
	}

	// Read required Concrete01 material parameters
	double fpc, epsc0, fpcu, epscu;

	if (Tcl_GetDouble(interp, argv[3], &fpc) != TCL_OK) {
	    opserr << "WARNING invalid fpc\n";
	    opserr << "Concrete01 material: " << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[4], &epsc0) != TCL_OK) {
	    opserr << "WARNING invalid epsc0\n";
	    opserr << "Concrete01 material: " << tag << endln;
	    return TCL_ERROR;
	}
	
	if (Tcl_GetDouble(interp, argv[5], &fpcu) != TCL_OK) {
	    opserr << "WARNING invalid fpcu\n";
	    opserr << "Concrete01 material: " << tag << endln;
	    return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[6], &epscu) != TCL_OK) {
	    opserr << "WARNING invalid epscu\n";
	    opserr << "Concrete01 material: " << tag << endln;
	    return TCL_ERROR;
	}

	if (argc == 7)
	  // Parsing was successful, allocate the material
	  theMaterial = new Concrete01WithSITC(tag, fpc, epsc0, fpcu, epscu);
	else {
	  double endStrainSITC;
	  if (Tcl_GetDouble(interp, argv[7], &endStrainSITC) != TCL_OK) {
	    opserr << "WARNING invalid epscu\n";
	    opserr << "Concrete01 material: " << tag << endln;
	    return TCL_ERROR;
	  }
	  theMaterial = new Concrete01WithSITC(tag, fpc, epsc0, fpcu, epscu, endStrainSITC);
	}
    }

    /*
    else if (strcmp(argv[1],"SMA") == 0) {
      
      if (argc < 9) {
	opserr << "Warning insufficient arguments\n";
	printCommand(argc,argv);
	opserr << "Want: uniaxialMaterialSMA  tag  E  eps_L  sig_AS_s  sig_AS_f  sig_SA_s  sig_SA_f" << endln;
	
	return TCL_ERROR;
      }
      
      int tag;
      double E, eps_L, sig_AS_s, sig_AS_f, sig_SA_s, sig_SA_f;
      
      if (Tcl_GetInt(interp,argv[2], &tag) != TCL_OK){
	opserr << "warning invalid uniaxialMaterial SMA tag" << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp,argv[3], &E) != TCL_OK){
	opserr << "warning invalid E\n";
	opserr << "uniaxialMaterial SMA: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp,argv[4], &eps_L) != TCL_OK){
	opserr << "warning invalid eps_L\n";
	opserr << "uniaxialMaterial SMA: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp,argv[5], &sig_AS_s) != TCL_OK){
	opserr << "warning invalid sig_AS_s\n";
	opserr << "uniaxialMaterial SMA: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp,argv[6], &sig_AS_f) != TCL_OK){
	opserr << "warning invalid sig_AS_f\n";
	opserr << "uniaxialMaterial SMA: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp,argv[7], &sig_SA_s) != TCL_OK){
	opserr << "warning invalid sig_SA_s\n";
	opserr << "uniaxialMaterial SMA: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp,argv[8], &sig_SA_f) != TCL_OK){
	opserr << "warning invalid sig_SA_f\n";
	opserr << "uniaxialMaterial SMA: " << tag << endln;
	return TCL_ERROR;
      }
      
      // Parsing was successful, allocate the material
      
      theMaterial = new SMAMaterial(tag, E, eps_L, sig_AS_s, sig_AS_f, sig_SA_s, sig_SA_f);
      
    }
    */

    else if (strcmp(argv[1],"ECC01") == 0) {
      if (argc < 16) {
	opserr << "WARNING insufficient arguments\n";
	printCommand(argc,argv);
	opserr << "Want: uniaxialMaterial ECC01 TAG? SIGT0? EPST0? SIGT1? EPST1? EPST2? SIGC0? EPSC0? EPSC1? ";
	opserr << "ALPHAT1? ALPHAT2? ALPHAC? ALPHACU? BETAT? BETAC\n";
	return TCL_ERROR;
      }
      
      int tag;
      double SIGT0, EPST0, SIGT1, EPST1, EPST2, SIGC0, EPSC0, EPSC1, ALPHAT1, ALPHAT2, ALPHAC, ALPHACU, BETAT, BETAC;
      
      
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING invalid uniaxialMaterial ECC01 tag" << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[3], &SIGT0) != TCL_OK) {
	opserr << "WARNING invalid SIGTO\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[4], &EPST0) != TCL_OK) {
	opserr << "WARNING invalid EPSTO\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[5], &SIGT1) != TCL_OK) {
	opserr << "WARNING invalid SIGT1\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[6], &EPST1) != TCL_OK) {
	opserr << "WARNING invalid EPST1\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[7], &EPST2) != TCL_OK) {
	opserr << "WARNING invalid epscu\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[8], &SIGC0) != TCL_OK) {
	opserr << "WARNING invalid epscu\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[9], &EPSC0) != TCL_OK) {
	opserr << "WARNING invalid epscu\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[10], &EPSC1) != TCL_OK) {
	opserr << "WARNING invalid epscu\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[11], &ALPHAT1) != TCL_OK) {
	opserr << "WARNING invalid epscu\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[12], &ALPHAT2) != TCL_OK) {
	opserr << "WARNING invalid epscu\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[13], &ALPHAC) != TCL_OK) {
	opserr << "WARNING invalid epscu\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[14], &ALPHACU) != TCL_OK) {
	opserr << "WARNING invalid epscu\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[15], &BETAT) != TCL_OK) {
	opserr << "WARNING invalid epscu\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[16], &BETAC) != TCL_OK) {
	opserr << "WARNING invalid epscu\n";
	opserr << "ECC01 material: " << tag << endln;
	return TCL_ERROR;
      }
      
      theMaterial = new ECC01(tag, SIGT0, EPST0, SIGT1, EPST1, EPST2, SIGC0, EPSC0, EPSC1, 
			      ALPHAT1, ALPHAT2, ALPHAC, ALPHACU, BETAT, BETAC);
    }


    else if (strcmp(argv[1], "ASD_SMA_3K") == 0) {
      void *theMat = OPS_ASD_SMA_3K();
      if (theMat != 0)
        theMaterial = (UniaxialMaterial *)theMat;
      else
        return TCL_ERROR;
    }


    else if (strcmp(argv[1],"SelfCentering") == 0) {
      if (argc < 7) {
	opserr << "WARNING insufficient arguments\n";
	printCommand(argc,argv);
	opserr << "Want: uniaxialMaterial SelfCentering tag? k1? k2? ActF? beta? <SlipDef? BearDef? rBear?>" << endln;
	return TCL_ERROR;
      }
      
      int tag;
      double k1, k2, ActF, beta, rBear, SlipDef, BearDef;
      
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING invalid uniaxialMaterial SelfCentering tag" << endln;
	return TCL_ERROR;		
      }
      
      if (Tcl_GetDouble(interp, argv[3], &k1) != TCL_OK) {
	opserr << "WARNING invalid k1\n";
	opserr << "uniaxialMaterial SelfCentering: " << tag << endln;
	return TCL_ERROR;	
      }
      
      if (Tcl_GetDouble(interp, argv[4], &k2) != TCL_OK) {
	opserr << "WARNING invalid k2\n";
	opserr << "uniaxialMaterial SelfCentering: " << tag << endln;
	return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[5], &ActF) != TCL_OK) {
	opserr << "WARNING invalid ActF\n";
	opserr << "uniaxialMaterial SelfCentering: " << tag << endln;
	return TCL_ERROR;	
      }
      
      if (Tcl_GetDouble(interp, argv[6], &beta) != TCL_OK) {
	opserr << "WARNING invalid beta\n";
	opserr << "uniaxialMaterial SelfCentering: " << tag << endln;
	return TCL_ERROR;	
      }
      
      if (argc == 8) {
	if (Tcl_GetDouble(interp, argv[7], &SlipDef) != TCL_OK) {
	  opserr << "WARNING invalid SlipDef\n";
	  opserr << "uniaxialMaterial SelfCentering: " << tag << endln;
	  return TCL_ERROR;	
	}
	// Parsing was successful, allocate the material
	theMaterial = new SelfCenteringMaterial (tag, k1, k2, ActF, beta, SlipDef, 0, 0);
      }
      
      else if (argc > 8) {
	if (Tcl_GetDouble(interp, argv[7], &SlipDef) != TCL_OK) {
	  opserr << "WARNING invalid SlipDef\n";
	  opserr << "uniaxialMaterial SelfCentering: " << tag << endln;
	  return TCL_ERROR;	
	}
	if (Tcl_GetDouble(interp, argv[8], &BearDef) != TCL_OK) {
	  opserr << "WARNING invalid BearDef\n";
	  opserr << "uniaxialMaterial SelfCentering: " << tag << endln;
	  return TCL_ERROR;	
	}
	if (Tcl_GetDouble(interp, argv[9], &rBear) != TCL_OK) {
	  opserr << "WARNING invalid rBear\n";
	  opserr << "uniaxialMaterial SelfCentering: " << tag << endln;
	  return TCL_ERROR;	
	}
	// Parsing was successful, allocate the material
	theMaterial = new SelfCenteringMaterial (tag, k1, k2, ActF, beta, SlipDef, BearDef, rBear);
      }
      
      else {
	// Parsing was successful, allocate the material
	theMaterial = new SelfCenteringMaterial (tag, k1, k2, ActF, beta, 0, 0, 0);
      }
    }
    
    else if (strcmp(argv[1],"SteelMP") == 0) {
      // Check that there is the minimum number of arguments
      if (argc < 4) {
        opserr << "WARNING insufficient arguments\n";
        printCommand(argc,argv);
        opserr << "Want: uniaxialMaterial SteelMP tag? fy? E0? b? ";
        opserr << " <coeffR1?  coeffR2? a1? a2?>" << endln;    
        return TCL_ERROR;
      }

      int tag;
      
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
        opserr << "WARNING invalid uniaxialMaterial SteelMP tag" << endln;
        return TCL_ERROR;
      }
      
      // Read required Steel01 material parameters
      double fy, E, b;
      
      if (Tcl_GetDouble(interp, argv[3], &fy) != TCL_OK) {
        opserr << "WARNING invalid fy\n";
        opserr << "uniaxialMaterial SteelMP: " << tag << endln;
        return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[4], &E) != TCL_OK) {
        opserr << "WARNING invalid E0\n";
        opserr << "uniaxialMaterial SteelMP: " << tag << endln;
        return TCL_ERROR;
      }
      
      if (Tcl_GetDouble(interp, argv[5], &b) != TCL_OK) {
        opserr << "WARNING invalid b\n";
        opserr << "uniaxialMaterial SteelMP: " << tag << endln;
        return TCL_ERROR;
      }
        
      if (argc < 5) {
	opserr << "WARNING insufficient number of hardening parameters\n";
	opserr << "uniaxialMaterial Steel03: " << tag << endln;
	return TCL_ERROR;
      } 

      // Read optional Steel01 material parameters
      double r, coeffR1, coeffR2, a1, a2;
      r=20.0;
      coeffR1 =18.5;
      coeffR2 =.15;	 
      a1=0;
      a2=0;
      
      if (argc >6) { 
	if (Tcl_GetDouble(interp, argv[6], &r) != TCL_OK) {
	  opserr << "WARNING invalid r\n";
	  opserr << "uniaxialMaterial SteelMP: " << tag << endln;
	  return TCL_ERROR;
	}
		
	if (Tcl_GetDouble(interp, argv[7], &coeffR1) != TCL_OK) {
	  opserr << "WARNING invalid CR1\n";
	  opserr << "uniaxialMaterial SteelMP: " << tag << endln;
	  return TCL_ERROR;
	}

	if (Tcl_GetDouble(interp, argv[8], &coeffR2) != TCL_OK) {
	  opserr << "WARNING invalid CR2\n";
	  opserr << "uniaxialMaterial SteelMP: " << tag << endln;
	  return TCL_ERROR;
	  
	}

	if (Tcl_GetDouble(interp, argv[9], &a1) != TCL_OK) {
	  opserr << "WARNING invalid a1\n";
	  opserr << "uniaxialMaterial SteelMP: " << tag << endln;
	  return TCL_ERROR;
	  
	}
	
	if (Tcl_GetDouble(interp, argv[10], &a2) != TCL_OK) {
	  opserr << "WARNING invalid a2\n";
	  opserr << "uniaxialMaterial SteelMP: " << tag << endln;
	  return TCL_ERROR;
	  
	}
      } //if

      theMaterial = new SteelMP (tag, fy, E, b, r, coeffR1,coeffR2, a1, a2);
    }

    else if (strcmp(argv[1],"SmoothPSConcrete") == 0) {
      if (argc < 6 || argc > 9) {
	opserr << "WARNING invalid number of arguments\n";
	printCommand(argc,argv);
	opserr << "Want: uniaxialMaterial SmoothPSConcrete tag? fc? fu? Ec? <eps0?> <epsu?> <eta?>" << endln;
	return TCL_ERROR;
      }    
      
      int tag;
      double fu, Ec, fc;
      double eps0=0.002;
      double epsu=0.005;
      double eta=0.2;
      
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING invalid uniaxialMaterial SmoothPSConcrete tag" << endln;
	return TCL_ERROR;		
      }
      
      if (Tcl_GetDouble(interp, argv[3], &fc) != TCL_OK) {
	opserr << "WARNING invalid fc\n";
	opserr << "uniaxiaMaterial SmoothPSConcrete: " << tag << endln;
	return TCL_ERROR;	
      }
      
      if (Tcl_GetDouble(interp, argv[4], &fu) != TCL_OK) {
	opserr << "WARNING invalid fu\n";
	opserr << "uniaxiaMaterial SmoothPSConcrete: " << tag << endln;
	return TCL_ERROR;	
      }
      
      if (Tcl_GetDouble(interp, argv[5], &Ec) != TCL_OK) {
	opserr << "WARNING invalid Ec\n";
	opserr << "uniaxiaMaterial SmoothPSConcrete: " << tag << endln;
	return TCL_ERROR;	
      }
      
      if (argc >= 7) 
	if (Tcl_GetDouble(interp,argv[6], &eps0) != TCL_OK) {
	  opserr << "WARNING invalid eps0\n";
	  opserr << "uniaxialMaterial SmoothPSConcrete: " << tag << endln;
	  return TCL_ERROR;
	}
      
      if (argc >= 8) 
	if (Tcl_GetDouble(interp,argv[7], &epsu) != TCL_OK) {
	  opserr << "WARNING invalid epsu\n";
	  opserr << "uniaxialMaterial SmoothPSConcrete: " << tag << endln;
	  return TCL_ERROR;
	}
      
      if (argc >= 9) 
	if (Tcl_GetDouble(interp,argv[8], &eta) != TCL_OK) {
	  opserr << "WARNING invalid eta\n";
	  opserr << "uniaxialMaterial SmoothPSConcrete: " << tag << endln;
	  return TCL_ERROR;
	}
      
      // Parsing was successful, allocate the material
      theMaterial = new SmoothPSConcrete( tag, fc, fu, Ec, eps0, epsu, eta);       
    }
    
	    // ----- 1D J2 Plasticity ----
	    else if (strcmp(argv[1],"UniaxialJ2Plasticity") == 0) {
      if (argc < 7) {
		opserr << "WARNING invalid number of arguments\n";
		printCommand(argc,argv);
		opserr << "Want: uniaxialMaterial UniaxialJ2Plasticity tag? E? sigmaY? Hkin? <Hiso?>" << endln;
		return TCL_ERROR;
      }    
      
      int tag;
      double E, sigmaY, Hkin, Hiso;
	  Hiso =0.0;
      
	  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
		opserr << "WARNING invalid uniaxialMaterial UniaxialJ2Plasticity tag" << endln;
		return TCL_ERROR;		
	  }
	      
	  if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
			opserr << "WARNING invalid E\n";
			opserr << "uniaxiaMaterial UniaxialJ2Plasticity: " << tag << endln;
			return TCL_ERROR;	
	  }
	      
	   if (Tcl_GetDouble(interp, argv[4], &sigmaY) != TCL_OK) {
			opserr << "WARNING invalid sigmaY\n";
			opserr << "uniaxiaMaterial UniaxialJ2Plasticity: " << tag << endln;
			return TCL_ERROR;	
		  }
	      
		  if (Tcl_GetDouble(interp, argv[5], &Hkin) != TCL_OK) {
			opserr << "WARNING invalid Hkin\n";
			opserr << "uniaxiaMaterial SmoothPSConcrete: " << tag << endln;
			return TCL_ERROR;	
			}
	      
		  if (argc >= 7) 
			if (Tcl_GetDouble(interp,argv[6], &Hiso) != TCL_OK) {
			  opserr << "WARNING invalid Hiso\n";
			  opserr << "uniaxialMaterial UniaxialJ2Plasticity: " << tag << endln;
			  return TCL_ERROR;
			}
	      
		       
      // Parsing was successful, allocate the material
      theMaterial = new UniaxialJ2Plasticity(tag, E, sigmaY, Hkin, Hiso);       
    }
    else if (strcmp(argv[1],"KikuchiAikenHDR") == 0) { 
      return TclCommand_KikuchiAikenHDR(clientData, interp, argc, argv);
    }
    else if (strcmp(argv[1],"KikuchiAikenLRB") == 0) { 
      return TclCommand_KikuchiAikenLRB(clientData, interp, argc, argv);
    }
    else if (strcmp(argv[1],"AxialSp") == 0) { 
      return TclCommand_AxialSp(clientData, interp, argc, argv);
    }
    else if (strcmp(argv[1],"AxialSpHD") == 0) { 
      return TclCommand_AxialSpHD(clientData, interp, argc, argv);
    }
    else if (strcmp(argv[1], "HystereticPoly") == 0) {		// BEGIN Salvatore Sessa 14-Jan-2021 Mail: salvatore.sessa2@unina.it
	void* theMat = OPS_HystereticPoly();
	if (theMat != 0)
		theMaterial = (UniaxialMaterial*)theMat;
	else
		return TCL_ERROR;
    }								// END Salvatore Sessa 14-Jan-2021 Mail: salvatore.sessa2@unina.it
    else {
      // Fedeas
      theMaterial = TclModelBuilder_addFedeasMaterial(clientData, interp, argc, argv);
      
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
    }

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

int
TclCommand_KikuchiAikenHDR(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  //arguments (necessary)
  int tag;
  int tp;
  double ar;
  double hr;

  //arguments (optional)
  double cg = 1.0;
  double ch = 1.0;
  double cu = 1.0;
  double rs = 1.0;
  double rf = 1.0;
      
  //
  UniaxialMaterial *theMaterial = 0;


  //error flag
  bool ifNoError = true;
  
  if (argc < 6) { // uniaxialMaterial KikuchiAikenHDR matTag? tp? ar? hr?

    opserr << "WARNING invalid number of arguments\n";
    ifNoError = false;

  } else {

    //argv[2~5]
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid KikuchiAikenHDR tag" << endln;
      ifNoError = false;
    }
  
    if ((strcmp(argv[3],"X0.6") == 0) || (strcmp(argv[3],"1") == 0)) {
      tp = 1;
    } else if ((strcmp(argv[3],"X0.6-0MPa") == 0) || (strcmp(argv[3],"2") == 0)) {
      tp = 2;
    } else if ((strcmp(argv[3],"X0.4") == 0) || (strcmp(argv[3],"3") == 0)) {
      tp = 3;
    } else if ((strcmp(argv[3],"X0.4-0MPa") == 0) || (strcmp(argv[3],"4") == 0)) {
      tp = 4;
    } else if ((strcmp(argv[3],"X0.3") == 0) || (strcmp(argv[3],"5") == 0)) {
      tp = 5;
    } else if ((strcmp(argv[3],"X0.3-0MPa") == 0) || (strcmp(argv[3],"6") == 0)) {
      tp = 6;
    } else {
      opserr << "WARNING invalid KikuchiAikenHDR tp" << endln;
      ifNoError = false;
    }
  
    if (Tcl_GetDouble(interp, argv[4], &ar) != TCL_OK || ar <= 0.0) {
      opserr << "WARNING invalid ar\n";
      ifNoError = false;
    }
    
    if (Tcl_GetDouble(interp, argv[5], &hr) != TCL_OK || hr <= 0.0) {
      opserr << "WARNING invalid hr\n";
      ifNoError = false;
    }
    
    
    //argv[6~]
    for (int i=6; i<=(argc-1); i++) {
      
      if (strcmp(argv[i],"-coGHU")==0 && (i+3)<=(argc-1)) { // <-coGHU cg? ch? cu?>

	if (Tcl_GetDouble(interp,argv[i+1], &cg) != TCL_OK || cg < 0.0) {
	  opserr << "WARNING invalid cg\n";
	  ifNoError = false;
	}
	
	if (Tcl_GetDouble(interp,argv[i+2], &ch) != TCL_OK || ch < 0.0) {
	  opserr << "WARNING invalid ch\n";
	  ifNoError = false;
	}
	
	if (Tcl_GetDouble(interp,argv[i+3], &cu) != TCL_OK || cu < 0.0) {
	  opserr << "WARNING invalid cu\n";
	  ifNoError = false;
	}
	
	i += 3;
	
      } else if (strcmp(argv[i],"-coMSS")==0 && (i+2)<=(argc-1)) { // <-coMSS rs? rf?>
	
	if (Tcl_GetDouble(interp,argv[i+1], &rs) != TCL_OK || rs < 0.0) {
	  opserr << "WARNING invalid rs\n";
	  ifNoError = false;
	}
	
	if (Tcl_GetDouble(interp,argv[i+2], &rf) != TCL_OK || rf < 0.0) {
	  opserr << "WARNING invalid rf\n";
	  ifNoError = false;
	} 

	i += 2;
	
      } else { // invalid option
	opserr << "WARNING invalid optional arguments \n";
	ifNoError = false;
	break;
      }

    }

  }

  //if error detected
  if (!ifNoError) {
    //input:
    opserr << "Input command: ";
    for (int i=0; i<argc; i++){
      opserr << argv[i] << " ";
    }
    opserr << endln;
    
    //want:
    opserr << "Want: uniaxialMaterial KikuchiAikenHDR matTag? tp? ar? hr? <-coGHU cg? ch? cu?> <-coMSS rs? rf?>" << endln;
    return TCL_ERROR;
  }

  //regard 0.0 input as mistake (substitute 1.0 for 0.0)
  if (cg == 0.0) cg = 1.0;
  if (ch == 0.0) ch = 1.0;
  if (cu == 0.0) cu = 1.0;
  if (rs == 0.0) rs = 1.0;
  if (rf == 0.0) rf = 1.0;


  // Parsing was successful, allocate the material
  theMaterial = new KikuchiAikenHDR(tag, tp, ar, hr, cg, ch, cu, rs, rf);

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

  // succeeded
  return TCL_OK;

}

int
TclCommand_KikuchiAikenLRB(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  //arguments (necessary)
  int tag;
  int type = 1;
  double ar = 0.0;
  double hr = 0.0;
  double gr = 0.392e6;
  double ap = 0.0;
  double tp = 8.33e6;
  double alph = 0.588e6;
  double beta = 13.0;

  //arguments (optional)
  double temp = 15.0;
  double rk = 1.0;
  double rq = 1.0;
  double rs = 1.0;
  double rf = 1.0;

  //
  UniaxialMaterial *theMaterial = 0;


  //error flag
  bool ifNoError = true;

  if (argc < 11) { // uniaxialMaterial KikuchiAikenLRB matTag? type? ar? hr? gr? ap? tp? alph? beta?

    opserr << "WARNING KikuchiAikenLRB invalid number of arguments\n";
    ifNoError = false;

  } else {

    //argv[2~10]
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING KikuchiAikenLRB invalid tag" << endln;
      ifNoError = false;
    }

    if (Tcl_GetInt(interp, argv[3], &type) != TCL_OK) {
      opserr << "WARNING KikuchiAikenLRB invalid type" << endln;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[4], &ar) != TCL_OK || ar <= 0.0) {
      opserr << "WARNING KikuchiAikenLRB invalid ar" << endln;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[5], &hr) != TCL_OK || ar <= 0.0) {
      opserr << "WARNING KikuchiAikenLRB invalid hr" << endln;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[6], &gr) != TCL_OK || gr <= 0.0) {
      opserr << "WARNING KikuchiAikenLRB invalid gr" << endln;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[7], &ap) != TCL_OK || ap <= 0.0) {
      opserr << "WARNING KikuchiAikenLRB invalid ap" << endln;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[8], &tp) != TCL_OK || tp <= 0.0) {
      opserr << "WARNING KikuchiAikenLRB invalid tp" << endln;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[9], &alph) != TCL_OK || alph <= 0.0) {
      opserr << "WARNING KikuchiAikenLRB invalid alph" << endln;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[10], &beta) != TCL_OK || beta <= 0.0) {
      opserr << "WARNING KikuchiAikenLRB invalid beta" << endln;
      ifNoError = false;
    }


    //argv[11~]
    for (int i=11; i<=(argc-1); i++) {

      if (strcmp(argv[i],"-T")==0 && (i+1)<=(argc-1)) { // <-T temp?>

	if (Tcl_GetDouble(interp,argv[i+1], &temp) != TCL_OK) {
	  opserr << "WARNING KikuchiAikenLRB invalid temp" << endln;
	  ifNoError = false;
	}
	
	i += 1;
	
      } else if (strcmp(argv[i],"-coKQ")==0 && (i+2)<=(argc-1)) { // <-coKQ rk? rq?>
	
	if (Tcl_GetDouble(interp,argv[i+1], &rk) != TCL_OK || rk < 0.0) {
	  opserr << "WARNING KikuchiAikenLRB invalid rk" << endln;
	  ifNoError = false;
	}
	
	if (Tcl_GetDouble(interp,argv[i+2], &rq) != TCL_OK || rq < 0.0) {
	  opserr << "WARNING KikuchiAikenLRB invalid rq" << endln;
	  ifNoError = false;
	} 

	i += 2;
      } else if (strcmp(argv[i],"-coMSS")==0 && (i+2)<=(argc-1)) { // <-coMSS rs? rf?>
	
	if (Tcl_GetDouble(interp,argv[i+1], &rs) != TCL_OK || rs < 0.0) {
	  opserr << "WARNING KikuchiAikenLRB invalid rs" << endln;
	  ifNoError = false;
	}
	
	if (Tcl_GetDouble(interp,argv[i+2], &rf) != TCL_OK || rf < 0.0) {
	  opserr << "WARNING KikuchiAikenLRB invalid rf" << endln;
	  ifNoError = false;
	} 

	i += 2;
	
      } else { // invalid option
	opserr << "WAINING KikuchiAikenLRB invalid optional arguments" << endln;
	ifNoError = false;
	break;
      }

    }

  }

  //if error detected
  if (!ifNoError) {
    //input:
    opserr << "Input command: ";
    for (int i=0; i<argc; i++){
      opserr << argv[i] << " ";
    }
    opserr << endln;
    
    //want:
    opserr << "Want: uniaxialMaterial KikuchiAikenLRB matTag? type? ar? hr? gr? ap? tp? alph? beta? <-T temp? > <-coKQ rk? rq?> <-coMSS rs? rf?>" << endln;
//
    return TCL_ERROR;
  }

  //regard 0.0 input as misteke (substitute 1.0 for 0.0)
  if (rk == 0.0) rk = 1.0;
  if (rq == 0.0) rq = 1.0;
  if (rs == 0.0) rs = 1.0;
  if (rf == 0.0) rf = 1.0;

  // Parsing was successful, allocate the material
  theMaterial = new KikuchiAikenLRB(tag, type, ar, hr, gr, ap, tp, alph, beta, temp, rk, rq, rs, rf);


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
  } else {
    return TCL_OK;
  }

}

int
TclCommand_AxialSp(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  //arguments (necessary)
  int tag;
  double sce;
  double fty;
  double fcy;

  //arguments (optional)
  double bte = 0.0;
  double bty = 0.0;
  double bcy = 0.0;
  double fcr = 0.0;
      
  //
  UniaxialMaterial *theMaterial = 0;


  //error flag
  bool ifNoError = true;
  
  if (argc < 6 || argc > 10) { // uniaxialMaterial AxialSp matTag? sce? fty? fcy?

    opserr << "WARNING invalid number of arguments\n";
    ifNoError = false;

  }

    //argv[2~5]
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid AxialSp tag" << endln;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[3], &sce) != TCL_OK) {
      opserr << "WARNING invalid sce\n";
      opserr << "AxialSp: " << tag << endln;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[4], &fty) != TCL_OK) {
      opserr << "WARNING invalid fty\n";
      opserr << "AxialSp: " << tag << endln;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[5], &fcy) != TCL_OK) {
      opserr << "WARNING invalid fcy\n";
      opserr << "AxialSp: " << tag << endln;
      ifNoError = false;
    }
    
    //argv[6~]
    if (argc >= 7) {
      if (Tcl_GetDouble(interp, argv[6], &bte) != TCL_OK) {
	opserr << "WARNING invalid bte\n";
	opserr << "AxialSp: " << tag << endln;
	ifNoError = false;
      }
    }

    if (argc >= 8) {
      if (Tcl_GetDouble(interp, argv[7], &bty) != TCL_OK) {
	opserr << "WARNING invalid bty\n";
	opserr << "AxialSp: " << tag << endln;
	ifNoError = false;
      }
    }

    if (argc >= 9) {
      if (Tcl_GetDouble(interp, argv[8], &bcy) != TCL_OK) {
	opserr << "WARNING invalid bcy\n";
	opserr << "AxialSp: " << tag << endln;
	ifNoError = false;
      }
    }

    if (argc == 10) {
      if (Tcl_GetDouble(interp, argv[9], &fcr) != TCL_OK) {
	opserr << "WARNING invalid fcr\n";
	opserr << "AxialSp: " << tag << endln;
	ifNoError = false;
      }
    }


  //if error detected
  if (!ifNoError) {
    //input:
    opserr << "Input command: ";
    for (int i=0; i<argc; i++){
      opserr << argv[i] << " ";
    }
    opserr << endln;
    
    //want:
    opserr << "WANT: AxialSp tag? sce? fty? fcy? <bte?> <bty?> <bcy?> <fcr?>" << endln;

    return TCL_ERROR;
  }

 
  
  // Parsing was successful, allocate the material
  theMaterial = new AxialSp(tag, sce, fty, fcy, bte, bty, bcy, fcr);

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
  } else {
    return TCL_OK;
  }

}

int
TclCommand_AxialSpHD(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  //arguments (necessary)
  int tag;
  double sce;
  double fty;
  double fcy;

  //arguments (optional)
  double bte = 1.0;
  double bty = 1.0;
  double bth = 1.0;
  double bcy = 1.0;
  double fcr = 0.0;
  double ath = 1.0;
      
  //
  UniaxialMaterial *theMaterial = 0;


  //error flag
  bool ifNoError = true;
  
  if (argc < 6 || argc > 12) { // uniaxialMaterial AxialSpHD matTag? sce? fty? fcy?

    opserr << "WARNING invalid number of arguments\n";
    ifNoError = false;

  }

    //argv[2~5]
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid AxialSpHD tag" << endln;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[3], &sce) != TCL_OK) {
      opserr << "WARNING invalid sce\n";
      opserr << "AxialSpHD: " << tag << endln;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[4], &fty) != TCL_OK) {
      opserr << "WARNING invalid fty\n";
      opserr << "AxialSpHD: " << tag << endln;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[5], &fcy) != TCL_OK) {
      opserr << "WARNING invalid fcy\n";
      opserr << "AxialSpHD: " << tag << endln;
      ifNoError = false;
    }
    
    //argv[6~]
    if (argc >= 7) {
      if (Tcl_GetDouble(interp, argv[6], &bte) != TCL_OK) {
	opserr << "WARNING invalid bte\n";
	opserr << "AxialSpHD: " << tag << endln;
	ifNoError = false;
      }
    }

    if (argc >= 8) {
      if (Tcl_GetDouble(interp, argv[7], &bty) != TCL_OK) {
	opserr << "WARNING invalid bty\n";
	opserr << "AxialSpHD: " << tag << endln;
	ifNoError = false;
      }
    }

    if (argc >= 9) {
      if (Tcl_GetDouble(interp, argv[8], &bth) != TCL_OK) {
	opserr << "WARNING invalid bth\n";
	opserr << "AxialSpHD: " << tag << endln;
	ifNoError = false;
      }
    }

    if (argc >= 10) {
      if (Tcl_GetDouble(interp, argv[9], &bcy) != TCL_OK) {
	opserr << "WARNING invalid bcy\n";
	opserr << "AxialSpHD: " << tag << endln;
	ifNoError = false;
      }
    }

    if (argc >= 11) {
      if (Tcl_GetDouble(interp, argv[10], &fcr) != TCL_OK) {
	opserr << "WARNING invalid fcr\n";
	opserr << "AxialSpHD: " << tag << endln;
	ifNoError = false;
      }
    }

    if (argc == 12) {
      if (Tcl_GetDouble(interp, argv[11], &ath) != TCL_OK) {
	opserr << "WARNING invalid ath\n";
	opserr << "AxialSpHD: " << tag << endln;
	ifNoError = false;
      }
    }

  //if error detected
  if (!ifNoError) {
    //input:
    opserr << "Input command: ";
    for (int i=0; i<argc; i++){
      opserr << argv[i] << " ";
    }
    opserr << endln;
    
    //wand:
    opserr << "WANT: AxialSpHD tag? sce? fty? fcy? <bte?> <bty?> <bth?> <bcy?> <fcr?> <ath?>" << endln;

    return TCL_ERROR;
  }

  
  // Parsing was successful, allocate the material
  theMaterial = new AxialSpHD(tag, sce, fty, fcy, bte, bty, bth, bcy, fcr, ath);

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
  } else {
    return TCL_OK;
  }

}
