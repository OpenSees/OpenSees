// Code written/implemented by:	Shaohui Zhang
//								Tsinghua University, Beijing, China
//								Xiaodong Ji (jixd@mail.tsinghua.edu.cn)
//								Tsinghua University, Beijing, China
//								Yue Yu
//								Tsinghua University, Beijing, China           
//
// Created: 12/2023
//
// Description: This file contains the definition for
// FAM_CS
// For Detailed explanation of the model, please refer to the references.
//
// References:
// 1) Zhang S., Ji X., Sun L., et al. New OpenSees material model for simulating reinforced concrete shear walls subjected to 
//    coupled axial tension and cyclic lateral loads. Eng Struct 2024;318:118774. https://doi.org/10.1016/j.engstruct.2024.118774.
// 2) Zhang S. Seismic shear force and shear design of reinforced concrete core walls [Ph.D. Thesis].
//    Beijing: Tsinghua University; 2026. (in Chinese)

// modified from FSAM.cpp    author: Kristijan Kolozvari, Kutay Orakcal, John Wallace

#include "FAM_CS.h"
#include <UniaxialMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include <math.h>
#include <float.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <stdlib.h>
#include <MaterialResponse.h>
#include <DummyStream.h>
#include <elementAPI.h>
#define OPS_Export 

#include <string.h>

#ifndef fmin
#define fmin(a,b) ( ((a)<(b))?(a):(b) )
#endif

#ifndef fmax
#define fmax(a,b) ( ((a)>(b))?(a):(b) )
#endif

#include "ConcreteCM.h" // for creating ConcreteCM inside the panel element


static int numFAM_CSMaterials = 0;

// Read input parameters and build the material
OPS_Export void *OPS_FAM_CSMaterial()
{	// print out a message about who wrote this element & any copyright info wanted
	if (numFAM_CSMaterials == 0) {
		opserr << " \n";
		opserr << "                          FAM_CS element v1.0\n";
		opserr << "                    by Shaohui Zhang, Xiaodong Ji, Yue Yu\n";
		opserr << "                           Tsinghua University\n";
		opserr << " \n";
		numFAM_CSMaterials++;
	}

	// Pointer to a uniaxial material that will be returned
	NDMaterial *theMaterial = 0;

	int numRemainingArgs = OPS_GetNumRemainingInputArgs();

	// Parse the script for material parameters
	if (numRemainingArgs != 11) { // total # of input parameters
		opserr << "Invalid #Args want: NDMaterial FAM_CS $mattag $Rho $Tag_UniaxialSteelX $Tag_UniaxialSteelY $Tag_UniaxialConcrete $rouX $rouY $dY $Gamax $lm0 $sh\n";
		return 0;	
	}

	int tag;				// nDMaterial tag
	double rho;				// nDMaterial density
	int    iData[3];		// # of uniaxial materials
	double dData[6];		// # of material parameters
	int numData = 0;

	// nDMaterial tag
	numData = 1;
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "WARNING invalid uDMaterial FAM_CS tag" << endln;
		return 0;
	}

	// nDMaterial density
	numData = 1;
	if (OPS_GetDouble(&numData, &rho) != 0) {
		opserr << "Invalid Arg rho: nDMaterial FAM_CS $mattag $rho $sX $sY $conc $rouX $rouY $dY $Gamax $lm0 $sh" << endln;
		return 0;	
	}

	// Material tags of 2 steel (in X and Y directions) and 1 concrete materials
	numData = 3;
	if (OPS_GetInt(&numData, iData) != 0) {
		opserr << "WARNING invalid uniaxialMaterial FAM_CS tag" << endln;
		return 0;
	}

	// Other FAM_CS material parameters
	numData = 6;
	if (OPS_GetDouble(&numData, dData) != 0) {
		opserr << "WARNING invalid uniaxialMaterial FAM_CS tag" << endln;
		return 0;
	}

	// Get pointers to Uniaxial materials
	// Steel X
	UniaxialMaterial *theUniaxialMaterial1 = OPS_GetUniaxialMaterial(iData[0]);
	if (theUniaxialMaterial1 == 0) {
		opserr << "WARNING material not found\n";
		opserr << "Material: " << iData[0];
		opserr << "\nFAM_CS: " << tag << endln;
		return 0;
	}

	// Steel Y
	UniaxialMaterial *theUniaxialMaterial2 = OPS_GetUniaxialMaterial(iData[1]);
	if (theUniaxialMaterial2 == 0) {
		opserr << "WARNING material not found\n";
		opserr << "Material: " << iData[1];
		opserr << "\nFAM_CS: " << tag << endln;
		return 0;
	}

	// Concrete c1.1 - uncracked (dummy)
	UniaxialMaterial *theUniaxialMaterial3 = OPS_GetUniaxialMaterial(iData[2]);
	if (theUniaxialMaterial3 == 0) {
		opserr << "WARNING material not found\n";
		opserr << "Material: " << iData[2];
		opserr << "\nFAM_CS: " << tag << endln;
		return 0;
	}

	// Concrete c1.2 - uncracked (dummy)
	UniaxialMaterial *theUniaxialMaterial4 = OPS_GetUniaxialMaterial(iData[2]);  
	if (theUniaxialMaterial4 == 0) {
		opserr << "WARNING material not found\n";
		opserr << "Material: " << iData[2];
		opserr << "\nFAM_CS: " << tag << endln;
		return 0;
	}

	// Concrete cA1 - 1st, 2nd crack
	UniaxialMaterial *theUniaxialMaterial5 = OPS_GetUniaxialMaterial(iData[2]);
	if (theUniaxialMaterial5 == 0) {
		opserr << "WARNING material not found\n";
		opserr << "Material: " << iData[2];
		opserr << "\nFAM_CS: " << tag << endln;
		return 0;
	}

	// Concrete cA2 - 1st, 2nd crack
	UniaxialMaterial *theUniaxialMaterial6 = OPS_GetUniaxialMaterial(iData[2]);  
	if (theUniaxialMaterial6 == 0) {
		opserr << "WARNING material not found\n";
		opserr << "Material: " << iData[2];
		opserr << "\nFAM_CS: " << tag << endln;
		return 0;
	}

	// Concrete cB1 - 2nd crack
	UniaxialMaterial *theUniaxialMaterial7 = OPS_GetUniaxialMaterial(iData[2]);
	if (theUniaxialMaterial7 == 0) {
		opserr << "WARNING material not found\n";
		opserr << "Material: " << iData[2];
		opserr << "\nFAM_CS: " << tag << endln;
		return 0;
	}

	// Concrete cB2 - 2nd crack
	UniaxialMaterial *theUniaxialMaterial8 = OPS_GetUniaxialMaterial(iData[2]);  
	if (theUniaxialMaterial8 == 0) {
		opserr << "WARNING material not found\n";
		opserr << "Material: " << iData[2];
		opserr << "\nFAM_CS: " << tag << endln;
		return 0;
	}

	//Create the FAM_CS
	theMaterial = new FAM_CS (tag, rho, 
		theUniaxialMaterial1, theUniaxialMaterial2, 
		theUniaxialMaterial3, theUniaxialMaterial4, 
		theUniaxialMaterial5, theUniaxialMaterial6,
		theUniaxialMaterial7, theUniaxialMaterial8,
		dData[0], dData[1], dData[2], dData[3], dData[4], dData[5]);

	if (theMaterial == 0) {
		opserr << "WARNING ran out of memory creating material\n";
		opserr << "FAM_CS: " << tag << endln;
		return 0;
	}

	return theMaterial;
}

// Typical Constructor
FAM_CS::FAM_CS (int tag, 
	double RHO,				// density
	UniaxialMaterial *s1,	// steel X
	UniaxialMaterial *s2,	// steel Y
	UniaxialMaterial *c1,	// concrete 1.1 - uncracked
	UniaxialMaterial *c2,	// concrete 1.2 - uncracked
	UniaxialMaterial *cA1,	// concrete A1 - 1st crack, 2nd crack
	UniaxialMaterial *cA2,	// concrete A2 - 1st crack, 2nd crack
	UniaxialMaterial *cB1,	// concrete B1 - 2nd crack
	UniaxialMaterial *cB2,	// concrete B2 - 2nd crack	    
	double ROUX,            // Reinforcing ratio of Steel X
	double ROUY,			// Reinforcing ratio of Steel X
	double DY,                    // diameter of Steel 2
	double GAMAX,            // Maximum size of aggregate
	double LM0,              //Average spacing of X-direction cracks
	double SH)		    //Average spacing of Y-direction cracks

	: NDMaterial(tag, ND_TAG_FAM_CS), 
	rho(RHO), roux(ROUX), rouy(ROUY), 
	dY(DY), Gamax(GAMAX), lm0(LM0), sh(SH),
	strain_vec(3), stress_vec(3), tangent_matrix(3,3), 
	CStress(3), CStrain(3), 
	CPanelConcStress(3), CPanelSteelStress(3), TPanelConcStress(3), TPanelSteelStress(3), 
	CStrainStressSteel1(2), CStrainStressSteel2(2), TStrainStressSteel1(2), TStrainStressSteel2(2),
	CStrainStressConc1(2), CStrainStressConc2(2), TStrainStressConc1(2), TStrainStressConc2(2), 
	CStrainStressInterlock1(2), CStrainStressInterlock2(2), TStrainStressInterlock1(2), TStrainStressInterlock2(2),
	CCrackingAngles(2),  pi(3.1415926535)
{

	TeTaSt = 0.0; // Direction of horizontal reinforcement (fixed for now)

	// Material parameters
	E0x = 0.0;
	E0y = 0.0;

	Ec = 0.0;
	fpc = 0.0;
	epcc = 0.0;
	et = 0.0;

	// Principal strains
	Tprstrain1 = 0.0;		// Temp
	Tprstrain2 = 0.0;		// Temp
	Cprstrain1 = 0.0;		// Committed
	Cprstrain2 = 0.0;		// Committed

	alpha_strain = 10.0;	// Principal strain direction
	alfa_crackA = 10.0;		// Direction of 1st strut
	alfa_crackB = 10.0;		// Direction of 2nd strut

	// State Variables
	crackA = 0;				// Crack/Strut 1
	crackB = 0;				// Crack/Strut 2

	// Softening parameters
	// Temp
	beta = 0.0;
	delbeta = 0.0;
	epsiloncmax = 0.0;

	// Uncracked
	Tepscmax1 = 0.0;
	Tepscmax2 = 0.0;
	Cepscmax1 = 0.0;
	Cepscmax2 = 0.0;

	// 1st concrete strut
	TepscmaxA1 = 0.0;
	TepscmaxA2 = 0.0;
	CepscmaxA1 = 0.0;
	CepscmaxA2 = 0.0;
	
	// 2nd concrete strut
	TepscmaxB1 = 0.0;
	TepscmaxB2 = 0.0;
	CepscmaxB1 = 0.0;
	CepscmaxB2 = 0.0;
	
	// Shear Aggregate Interlock History Variables
	// "transform" variables
	gamma_cr = 0.0;
	Tau_Interlock = 0.0;
	Tau_Interlock0 = 0.0;
	dTau_de12_cr = 0.0;
	dTau_denormal = 0.0;

	//Maekawa model add by Z
	Beta = 0.0;

	// Crack Slip
	TeA12 = 0.0;
	TeB12 = 0.0;

	//Maekawa model add by Z
	Tbeta_A = 0.0;
	Tbeta_pmaxA = 0.0;	
	Cbeta_pmaxA = 0.0;
	Tbeta_nmaxA = 0.0;
	Cbeta_nmaxA = 0.0;
	Tbeta_B = 0.0;
	Tbeta_pmaxB = 0.0;
	Cbeta_pmaxB = 0.0;
	Tbeta_nmaxB = 0.0;
	Cbeta_nmaxB = 0.0;

	// Shear Stress
	Tgamma_cr_A = 0.0;
	Tgamma_cr_B = 0.0;
	Ttau_Interlock_A = 0.0;
	Ttau_Interlock_B = 0.0;	

	//Maekawa model add by Z
	Ctau0_pmaxA = 0.0;
	Ttau0_pmaxA = 0.0;
	Ctau0_nmaxA = 0.0;
	Ttau0_nmaxA = 0.0;
	Ctau0_pmaxB = 0.0;
	Ttau0_pmaxB = 0.0;
	Ctau0_nmaxB = 0.0;
	Ttau0_nmaxB = 0.0;

	// History variables for 2nd cracking criterium (cyclic cracking strain)
	TepsA2 = 0.0;
	CepsA2 = 0.0;

	//dowel action
	CDI_A_p = 0.0;
	TDI_A_p = 0.0;
	CDI_A_n = 0.0;
	TDI_A_n = 0.0;
	CDI_B_p = 0.0;
	TDI_B_p = 0.0;
	CDI_B_n = 0.0;
	TDI_B_n = 0.0;

	// Committed Stress
	CStress(0) = 0.0;
	CStress(1) = 0.0;
	CStress(2) = 0.0;

	// Committed Strain
	CStrain(0) = 0.0;
	CStrain(1) = 0.0;
	CStrain(2) = 0.0;

	// Fill arrays with zeros
	for (int i = 0; i<3; i++) {
	TPanelConcStress(i) = 0.0;
	TPanelSteelStress(i) = 0.0;
	
	CPanelConcStress(i) = 0.0;
	CPanelSteelStress(i) = 0.0;
	}

	for (int i =0; i<2; i++) {
	TStrainStressSteel1(i) = 0.0; 
	TStrainStressSteel2(i) = 0.0; 
	TStrainStressConc1(i) = 0.0;
	TStrainStressConc2(i) = 0.0; 
	TStrainStressInterlock1(i) = 0.0;
	TStrainStressInterlock2(i) = 0.0; 

	CStrainStressSteel1(i) = 0.0; 
	CStrainStressSteel2(i) = 0.0; 
	CStrainStressConc1(i) = 0.0;
	CStrainStressConc2(i) = 0.0; 
	CStrainStressInterlock1(i) = 0.0;
	CStrainStressInterlock2(i) = 0.0;
	CCrackingAngles(i) = 0.0;
	}

	// Allocate pointers for uniaxial materials ...................................................
	theMaterial = new UniaxialMaterial *[8];
	if ( theMaterial == 0 ) {
		opserr << " FAM_CS::FAM_CS - failed allocate material array\n";
		exit(-1);
	}

	// Get the copy for SteelX
	theMaterial[0] = s1->getCopy();
	// Check allocation    
	if ( theMaterial[0] == 0 ) {
		opserr << " FAM_CS::FAM_CS - failed to get a copy for Steel1\n";
		exit(-1);
	}

	// Get the copy for SteelY
	theMaterial[1] = s2->getCopy();	
	// Check allocation    
	if ( theMaterial[1] == 0 ) {
		opserr << " FAM_CS::FAM_CS - failed to get a copy for Steel2\n";
		exit(-1);
	}

	// Get the copy for Concrete A1
	theMaterial[4] = cA1->getCopy();	
	// Check allocation    
	if ( theMaterial[4] == 0 ) {
		opserr << " FAM_CS::FAM_CS - failed to get a copy for Concrete A1\n";
		exit(-1);
	}

	// Get the copy for Concrete A2
	theMaterial[5] = cA2->getCopy();	
	// Check allocation    
	if ( theMaterial[5] == 0 ) {
		opserr << " FAM_CS::FAM_CS - failed to get a copy for Concrete A2\n";
		exit(-1);
	}

	// Get the copy for Concrete B1
	theMaterial[6] = cB1->getCopy();	
	// Check allocation    
	if ( theMaterial[6] == 0 ) {
		opserr << " FAM_CS::FAM_CS - failed to get a copy for Concrete B1\n";
		exit(-1);
	}

	// Get the copy for Concrete B2
	theMaterial[7] = cB2->getCopy();	
	// Check allocation    
	if ( theMaterial[7] == 0 ) {
		opserr << " FAM_CS::FAM_CS - failed to get a copy for Concrete B2\n";
		exit(-1);
	}

	// get/set responses
	theResponses = new Response *[4];
	if ( theResponses == 0) {
		opserr << " FAM_CS::FAM_CS - failed allocate responses array\n";
		exit(-1);
	}

	//OPS_Stream *theDummyStream = new DummyStream();
	DummyStream theDummyStream;
	//const char **argv = new const char *[1];
	//argv[0] = "getCommittedCyclicCrackingConcreteStrain"; // to get committed concrete cyclic cracking strain from strut A2
	char aa[80] = "getCommittedCyclicCrackingConcreteStrain";
	const char *argv[1];
	argv[0] = aa;
	theResponses[0] = theMaterial[5]->setResponse(argv, 1, theDummyStream);

	if (theResponses[0] == 0) {
			opserr << " FAM_CS::FAM_CS - failed to get cracking strain for material with tag: " << tag << "\n";
			exit(-1);
	}

	//argv[0] = "getInputParameters"; // to get input parameters from ConcreteCM
	char bb[80] = "getInputParameters";
	argv[0] = bb;
	
	theResponses[1] = theMaterial[4]->setResponse(argv, 1, theDummyStream);

	if (theResponses[1] == 0) {
			opserr << " FAM_CS::FAM_CS - failed to get input parameters for material with tag: " << tag << "\n";
			exit(-1);
	}

	char cc[80] = "getCommittedNewOrigin";
	argv[0] = cc;

	theResponses[2] = theMaterial[4]->setResponse(argv, 1, theDummyStream);

	if (theResponses[2] == 0) {
			opserr << " FAM_CS::FAM_CS - failed to get committed new origin for material with tag: " << tag << "\n";
			exit(-1);
	}

	theResponses[3] = theMaterial[6]->setResponse(argv, 1, theDummyStream);

	if (theResponses[3] == 0) {
			opserr << " FAM_CS::FAM_CS - failed to get committed new origin for material with tag: " << tag << "\n";
			exit(-1);
	}

	//delete theDummyStream;

	// Get ConcreteCM material input variables
	theResponses[1]->getResponse();
	Information &theInfoInput = theResponses[1]->getInformation();
	const Vector &ConcreteInput = theInfoInput.getData();

	// Now create monotonic concrete materials for uncracked stage of behavior
	// Concrete 1.1
	// Instead of: theMaterial[2] = c1->getCopy(); we are creating monotonic ConcreteCM	
	theMaterial[2] = new ConcreteCM(-1111, ConcreteInput[1], ConcreteInput[2], ConcreteInput[3], 
		ConcreteInput[4], ConcreteInput[5], ConcreteInput[6], ConcreteInput[7], ConcreteInput[8], ConcreteInput[9], 1); // create monotonic concrete

	// Check allocation    
	if ( theMaterial[2] == 0 ) {
		opserr << " FAM_CS::FAM_CS - failed to get a copy for Concrete 1\n";
		exit(-1);
	}

	// Concrete 1.2
	//Instead of: theMaterial[3] = c2->getCopy();  we are creating monotonic ConcreteCM	
	theMaterial[3] = new ConcreteCM(-2222, ConcreteInput[1], ConcreteInput[2], ConcreteInput[3], 
		ConcreteInput[4], ConcreteInput[5], ConcreteInput[6], ConcreteInput[7], ConcreteInput[8], ConcreteInput[9], 1); // create monotonic concrete

	// Check allocation    
	if ( theMaterial[3] == 0 ) {
		opserr << " FAM_CS::FAM_CS - failed to get a copy for Concrete 2\n";
		exit(-1);
	}

	// Obtain some material properties used later in the FAM_CS model
	// Young's modulus for concrete
	Ec = theMaterial[4] -> getInitialTangent();

	// Strain at peak compressive stress for concrete
	epcc = ConcreteInput[2];
	
	// Peak compressive stress for concrete
	fpc = ConcreteInput[1];

	// Cracking strain for concrete
	et = ConcreteInput[7];

	// Young's modulus for steel
	E0x = theMaterial[0] -> getInitialTangent(); // Horizontal reinforcement
	E0y = theMaterial[1] -> getInitialTangent(); // Vertical reinforcement

	this->revertToStart();
}

// Blank constructor
FAM_CS::FAM_CS():NDMaterial(0, ND_TAG_FAM_CS), 
	strain_vec(3), stress_vec(3), tangent_matrix(3,3),
	pi(3.1415926535)

{
	theMaterial = 0;
	theResponses = 0;

	this->revertToStart();
}

// Destructor
FAM_CS::~FAM_CS()
{
	// Delete the pointers
	if (theMaterial != 0) {
		for (int i=0; i<8; i++)
		{
			if (theMaterial[i])
				delete theMaterial[i];
		}
		delete [] theMaterial;
	}

	if (theResponses != 0) {
		for (int j=0; j<4; j++)
		{
			if (theResponses[j] != 0)
				delete theResponses[j];
		}
		delete [] theResponses;
	}

}

// get copy
NDMaterial* FAM_CS::getCopy(void) 
{

	FAM_CS* theCopy =
		new FAM_CS( this->getTag(), 
		rho,
		theMaterial[0], 
		theMaterial[1], 
		theMaterial[2], 
		theMaterial[3], 
		theMaterial[4], 
		theMaterial[5],
		theMaterial[6],
		theMaterial[7], 
		roux, 
		rouy, 
		dY,
		Gamax,
		lm0,
		sh);
	
	return theCopy;
}

// get copy
NDMaterial* FAM_CS::getCopy(const char *type)
{

	FAM_CS* theModel =
		new FAM_CS( this->getTag(), 
		rho, 
		theMaterial[0], 
		theMaterial[1], 
		theMaterial[2], 
		theMaterial[3], 
		theMaterial[4], 
		theMaterial[5], 
		theMaterial[6], 
		theMaterial[7], 
		roux, 
		rouy,
		dY,
		Gamax,
		lm0,
		sh);

	return theModel;
}

// Print 
void FAM_CS::Print(OPS_Stream &s, int flag)
{
	s << "\nFAM_CS, nDMaterial tag: " << this->getTag() << endln;

	// Input values
	s << "density: " << rho << endln;
	s << "roux: " << roux << ", rouy: " << rouy << endln;
	s << ", dY: " << dY << ", Gamax: " << Gamax << endln;
	s << "lm0: " << lm0 << ", sh: " << sh << endln;

	// Output values
	// Strain and stress of the uniaxial materials
	s << "Strain and stress of the uniaxial materials:"<<endln;
	s<< " Steel X: Strain = "<<theMaterial[0]->getStrain()<<", Stress = "<<theMaterial[0]->getStress()<< endln;
	s<< " Steel Y: Strain = "<<theMaterial[1]->getStrain()<<", Stress = "<<theMaterial[1]->getStress()<< endln;
	s<< " Concrete 1-A1: Strain = "<<theMaterial[2]->getStrain()<<", Stress = "<<theMaterial[2]->getStress()<< endln;
	s<< " Concrete 1-A2: Strain = "<<theMaterial[3]->getStrain()<<", Stress = "<<theMaterial[3]->getStress()<< endln;
	s<< " Concrete 2-A1: Strain = "<<theMaterial[4]->getStrain()<<", Stress = "<<theMaterial[4]->getStress()<< endln;
	s<< " Concrete 2-A2: Strain = "<<theMaterial[5]->getStrain()<<", Stress = "<<theMaterial[5]->getStress()<< endln;
	s<< " Concrete 2-B1: Strain = "<<theMaterial[6]->getStrain()<<", Stress = "<<theMaterial[6]->getStress()<< endln;
	s<< " Concrete 2-B2: Strain = "<<theMaterial[7]->getStrain()<<", Stress = "<<theMaterial[7]->getStress()<< endln;
	s << " Crack Angle 1 = " << CCrackingAngles(0) << endln;
	s << " Crack Angle 2 = " << CCrackingAngles(1) << endln;

	// Strain and Stress of the RC panel
	s << "Panel strains:"<<endln;
	s<< " EpsX = "<< CStrain(0) << ", EpsY = " << CStrain(1) << ", GammaXY = " << CStrain(2) << endln;
	s << "Panel stresses:"<<endln;
	s<< " SigX = "<< CStress(0) << ", SigY = " << CStress(1) << ", TauXY = " << CStress(2) << endln;

}

int FAM_CS::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	int dataTag = this->getDbTag();

	// Packs its data into a Vector and sends this to theChannel
	static Vector data(8);

	data(0) = this->getTag();
	data(1) = rho;
	data(2) = roux;
	data(3) = rouy;
	data(4) = dY;
	data(5) = Gamax;
	data(6) = lm0;
	data(7) = sh;

	res += theChannel.sendVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING FAM_CS::sendSelf() - " << this->getTag() << " failed to send Vector\n";
		return res;
	}	      

	// Sends the IDs of its materials
	int matDbTag;

	static ID idData(16); // 2 x # of materials

	// NOTE: to ensure that the material has a database tag if sending to a database channel.

	int i;
	for (i=0; i<8; i++) // i < # of materials
	{
		idData(i) = theMaterial[i]->getClassTag();
		matDbTag = theMaterial[i]->getDbTag();
		if (matDbTag == 0) {
			matDbTag = theChannel.getDbTag();
			if (matDbTag != 0)
				theMaterial[i]->setDbTag(matDbTag);
		}
		idData(i+8) = matDbTag; // i + # of materials
	}

	res += theChannel.sendID(dataTag, commitTag, idData);
	if (res < 0) {
		opserr << "WARNING FAM_CS::sendSelf() - " << this->getTag() << " failed to send ID\n";
		return res;
	}

	// Quad asks its material objects to send themselves
	for (i = 0; i < 8; i++) { // i < # of materials
		res += theMaterial[i]->sendSelf(commitTag, theChannel);
		if (res < 0) {
			opserr << "FAM_CS::sendSelf() - " << this->getTag() << " failed to send its Material\n";
			return res;
		}
	}	

	return res;
}

int FAM_CS::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res = 0;

	int dataTag = this->getDbTag();

	// Quad creates a Vector, receives the Vector and then sets the internal data with the data in the Vector
	static Vector data(8);
	res += theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING FAM_CS::recvSelf() - failed to receive Vector\n";
		return res;
	}

	this->setTag((int)data(0));
	rho	  = data(1);
	roux  = data(2);
	rouy  = data(3);
	dY    = data(4);
	Gamax = data(5);
	lm0   = data(6);
	sh    = data(7);

	static ID idData(16); // idData(2 x # of uniaxial materials)

	// Receives the tags of its materials
	res += theChannel.recvID(dataTag, commitTag, idData);
	if (res < 0) {
		opserr << "WARNING FAM_CS::recvSelf() - " << this->getTag() << " failed to receive ID\n";
		return res;
	}

	if (theMaterial == 0) {
		// Allocate new materials
		theMaterial = new UniaxialMaterial *[8]; // *[# of uniaxial materials]
		if (theMaterial == 0) {
			opserr << "FAM_CS::recvSelf() - Could not allocate UniaxialMaterial* array\n";
			return -1;
		}
		for (int i = 0; i < 8; i++) { // i < # of uniaxial materials
			int matClassTag = idData(i);
			int matDbTag = idData(i+8); // (i + # of uniaxial materials)
			// Allocate new material with the sent class tag
			theMaterial[i] = theBroker.getNewUniaxialMaterial(matClassTag);
			if (theMaterial[i] == 0) {
				opserr << "FAM_CS::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
				return -1;
			}
			// Now receive materials into the newly allocated space
			theMaterial[i]->setDbTag(matDbTag);
			res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
			if (res < 0) {
				opserr << "FAM_CS::recvSelf() - material " << i << "failed to recv itself\n";
				return res;
			}
		}
	}

	// materials exist - ensure materials of correct type and recvSelf on them
	else {
		for (int i = 0; i < 8; i++) { // i < # of uniaxial materials
			int matClassTag = idData(i);
			int matDbTag = idData(i+8); // (i+ # of uniaxial materials)
			
			// Check that material is of the right type; if not,
			// delete it and create a new one of the right type
			if (theMaterial[i]->getClassTag() != matClassTag) {
				delete theMaterial[i];
				theMaterial[i] = theBroker.getNewUniaxialMaterial(matClassTag);
				if (theMaterial[i] == 0) {
					opserr << "FAM_CS::recvSelf() - material " << i << "failed to create\n";
					return -1;
				}
			}
			// Receive the material
			theMaterial[i]->setDbTag(matDbTag);
			res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
			if (res < 0) {
				opserr << "FAM_CS::recvSelf() - material " << i << "failed to recv itself\n";
				return res;
			}
		}
	}

	return res;
}

// Get density
double FAM_CS::getRho(void)
{
	return rho;
}

// Load strain from the element
int FAM_CS::setTrialStrain(const Vector &v)
{
	// Set values for strain_vec
	strain_vec(0) = v(0);
	strain_vec(1) = v(1);
	strain_vec(2) = v(2);

	// Calculate trial stress and tangent partial stiffness
	this -> determineTrialStressAndTangent();

	return 0;
}

// Load strain from the element
int FAM_CS::setTrialStrain(const Vector &v, const Vector &r)
{
	// Set values for strain_vec
	strain_vec(0) = v(0);
	strain_vec(1) = v(1);
	strain_vec(2) = v(2);

	// Calculate trial stress and tangent partial stiffness
	determineTrialStressAndTangent();

	return 0;
}

// Load strain from the element
int FAM_CS::setTrialStrainIncr(const Vector &v)
{
	// Set values for strain_vec
	strain_vec(0) = v(0);
	strain_vec(1) = v(1);
	strain_vec(2) = v(2);

	// Calculate trial stress and tangent partial stiffness
	determineTrialStressAndTangent();

	return 0;
}

// Load strain from the element
int FAM_CS::setTrialStrainIncr(const Vector &v, const Vector &r)
{
	// Set values for strain_vec
	strain_vec(0) = v(0);
	strain_vec(1) = v(1);
	strain_vec(2) = v(2);

	// Calculate trial stress and tangent partial stiffness
	determineTrialStressAndTangent();

	return 0;
}

// Get partial tangent stiffness
const Matrix& FAM_CS::getTangent (void)
{
	return tangent_matrix;
}

// Get trial stress vector
const Vector& FAM_CS::getStress(void)
{
	return stress_vec;
}

// Get strain vector
const Vector& FAM_CS::getStrain()
{
	return strain_vec;
}

// Get committed stress
const Vector& FAM_CS::getCommittedStress(void)
{
	return CStress;
}

// Get committed strain
const Vector& FAM_CS::getCommittedStrain(void)
{
	return CStrain;
}

// Functions below return values for recorders
// Concrete stresses
Vector FAM_CS::getPanelStressConcrete(void) 
{
	return CPanelConcStress;
}

// Steel stresses (multiplied with reinforcing ratio)
Vector FAM_CS::getPanelStressSteel(void) 
{
	return CPanelSteelStress;
}

// Steel stresses in horizontal direction (from uniaxial material)
Vector FAM_CS::getStrainStressSteel1(void) 
{
	return CStrainStressSteel1;
}

// Steel stresses in vertical direction (from uniaxial material)
Vector FAM_CS::getStrainStressSteel2(void) 
{
	return CStrainStressSteel2;
}

// Concrete stresses along 1st strut
Vector FAM_CS::getStrainStressConcrete1(void) 
{
	return CStrainStressConc1;
}

// Concrete stresses along 2nd strut
Vector FAM_CS::getStrainStressConcrete2(void) 
{
	return CStrainStressConc2;
}

// Aggregate interlock stresses along 1st strut
Vector FAM_CS::getStrainStressInterlock1(void) 
{
	return CStrainStressInterlock1;
}

// Aggregate interlock stresses along 2nd strut
Vector FAM_CS::getStrainStressInterlock2(void) 
{
	return CStrainStressInterlock2; 
}

// Cracking angles
Vector FAM_CS::getCrackingAngles(void) 
{
	return CCrackingAngles;
}

// Revert to start
int FAM_CS::revertToStart(void)
{

	// revert all uniaxial materials to start
	for (int i=0; i < 8; i++) {
		theMaterial[i]->revertToStart();
	}

	// set stress, strain and tangent to zero
	strain_vec.Zero();
	stress_vec.Zero();
	tangent_matrix.Zero();

	return 0;
}

// Calculate initial partial stiffness matrix
const Matrix& FAM_CS::getInitialTangent (void)
{
	double TeTaStper;
	double Ecsafe = 1.00*Ec;

	if (TeTaSt >= 0.0) {
		TeTaStper=TeTaSt - 0.5*pi;
	}
	else {
		TeTaStper=TeTaSt + 0.5*pi;
	}

	double dsxdex =Ecsafe+(4.0*(E0x*roux - E0y*rouy)*cos(2.0*TeTaSt) + (E0x*roux + E0y*rouy)*(3.0 + cos(4*TeTaSt)))/8.0;
	double dsxdey =(E0x*roux + E0y*rouy)*pow(cos(TeTaSt),2.0)*pow(sin(TeTaSt),2.0);
	double dsxdgamma =((E0x*roux - E0y*rouy + (E0x*roux + E0y*rouy)*cos(2*TeTaSt))*sin(2.0*TeTaSt))/4.0;

	double dsydex =(E0x*roux + E0y*rouy)*pow(cos(TeTaSt),2.0)*pow(sin(TeTaSt),2.0);
	double dsydey =Ecsafe+((-4.0*E0x*roux + 4.0*E0y*rouy)*cos(2.0*TeTaSt) + (E0x*roux + E0y*rouy)*(3.0 + cos(4*TeTaSt)))/8.0;
	double dsydgamma =-((-(E0x*roux) + E0y*rouy + (E0x*roux + E0y*rouy)*cos(2.0*TeTaSt))*sin(2.0*TeTaSt))/4.0;

	double dtxydex =((E0x*roux - E0y*rouy + (E0x*roux + E0y*rouy)*cos(2.0*TeTaSt))*sin(2.0*TeTaSt))/4.0;
	double dtxydey =-((-(E0x*roux) + E0y*rouy + (E0x*roux + E0y*rouy)*cos(2.0*TeTaSt))*sin(2.0*TeTaSt))/4.0;
	double dtxydgamma =(Ecsafe/2.0)+(E0x*roux + E0y*rouy)*pow(cos(TeTaSt),2.0)*pow(sin(TeTaSt),2.0);

	tangent_matrix(0,0) = dsxdex;
	tangent_matrix(0,1) = dsxdey;
	tangent_matrix(0,2) = dsxdgamma;

	tangent_matrix(1,0) = dsydex;
	tangent_matrix(1,1) = dsydey;
	tangent_matrix(1,2) = dsydgamma;

	tangent_matrix(2,0) = dtxydex;
	tangent_matrix(2,1) = dtxydey;
	tangent_matrix(2,2) = dtxydgamma;

	return tangent_matrix;
}

// Commit State - determines state of concrete panel based on strain field and commits uniaxial materials
int FAM_CS::commitState(void)
{
	double Tstrain[3];			//ex, ey, gamma

	// Get strain values from strain of element
	Tstrain[0] = strain_vec(0); //ex
	Tstrain[1] = strain_vec(1); //ey
	Tstrain[2] = strain_vec(2); //gxy

	// Determine State
	// Uncracked panel behavior ......................................................
	if (crackA == 0 && crackB == 0) 
	{
		// Stage 1 - uncracked
		Stage1(Tstrain[0], Tstrain[1], Tstrain[2]);

		// Store committed values of epscmax
		Cepscmax1 = Tepscmax1;
		Cepscmax2 = Tepscmax2;

		// Store committed values of principal strains
		Cprstrain1 = Tprstrain1;
		Cprstrain2 = Tprstrain2;

		// cracking criterium for 1st crack 
		if ( fmax(Cprstrain1,Cprstrain2) >= et) {

			// Initiate 1st crack
			crackA = 1;

			// Calculate Direction of 1st Concrete Strut
			double extest = 0.5*(Tstrain[0] + Tstrain[1]) + 0.5*(Tstrain[0] - Tstrain[1] )*cos( 2.0*alpha_strain ) + 0.5*Tstrain[2]*sin(2.0*alpha_strain);

			if (fabs(extest-Cprstrain1) < fabs(extest-Cprstrain2))
			{
				if (Cprstrain1 >= Cprstrain2) {
					alfa_crackA = alpha_strain;
				} else {
					if ( alpha_strain < 0.0) {
						alfa_crackA = alpha_strain + 0.5*pi;
					} else {
						alfa_crackA = alpha_strain - 0.5*pi;
					}
				}
			}
			else
			{
				if (Cprstrain2 >= Cprstrain1) {
					alfa_crackA = alpha_strain;
				} else {
					if (alpha_strain < 0.0) {
						alfa_crackA = alpha_strain + 0.5*pi;
					} else {
						alfa_crackA = alpha_strain - 0.5*pi;
					}
				}
			}

		CepscmaxA1 = Cepscmax1;
		CepscmaxA2 = Cepscmax2;

		// Stage 2 - 1st crack (strut A)
		Stage2(Tstrain[0], Tstrain[1], Tstrain[2]);

		// Store committed values of epscmax
		CepscmaxA1 = TepscmaxA1;
		CepscmaxA2 = TepscmaxA2;

		// Store committed values of shear aggregate interlock
		Cbeta_pmaxA = Tbeta_pmaxA;
		Cbeta_nmaxA = Tbeta_nmaxA;
		Ctau0_pmaxA = Ttau0_pmaxA;
		Ctau0_nmaxA = Ttau0_nmaxA;
		
		// Store committed values compressive strain in perpendicular direction (for 2nd cracking criterium)
		CepsA2 = TepsA2;

		//dowel action
		CDI_A_p = TDI_A_p;
		CDI_A_n = TDI_A_n;
		}

		// Commit State of Uniaxial Materials
		for (int i=0; i < 8; i++)
		{
			theMaterial[i]->commitState();
		}

	} 
	
	// 1st crack ....................................................................
	else if (crackA == 1 && crackB == 0) { 

		// Stage 2 - 1st crack (strut A)
		Stage2(Tstrain[0], Tstrain[1], Tstrain[2]);

		// Store committed values of epscmax
		CepscmaxA1 = TepscmaxA1;
		CepscmaxA2 = TepscmaxA2;

		// Store committed values of shear aggregate interlock
		Cbeta_pmaxA = Tbeta_pmaxA;
		Cbeta_nmaxA = Tbeta_nmaxA;
		Ctau0_pmaxA = Ttau0_pmaxA;
		Ctau0_nmaxA = Ttau0_nmaxA;

		// Store committed values compressive strain in perpendicular direction (for 2nd cracking criterium)
		CepsA2 = TepsA2;

		//dowel action
		CDI_A_p = TDI_A_p;
		CDI_A_n = TDI_A_n;

		// Commit State of Uniaxial Materials
		for (int i=0; i < 8; i++)
		{
			theMaterial[i]->commitState();
		}

		// Get Committed value of cyclic cracking strain for Concrete
		theResponses[0]->getResponse();
		Information &theInfoA2 = theResponses[0]->getInformation();
		double eunpA2 = theInfoA2.theDouble;

		// cracking criterium for 2nd crack .....................................
		if (CepsA2 >= eunpA2) { //eA2 from the 1st crack

			// Initiate 2nd crack
			crackB = 1; 
			if (alfa_crackA < 0.0) {
				alfa_crackB = alfa_crackA + 0.5*pi;
			} else {
				alfa_crackB = alfa_crackA - 0.5*pi;
			}

			// Stage 3 - 2nd crack (strut B)
			Stage3(Tstrain[0], Tstrain[1], Tstrain[2]);

			// Store committed values of epscmax
			CepscmaxA1 = TepscmaxA1;
			CepscmaxA2 = TepscmaxA2;
			CepscmaxB1 = TepscmaxB1;
			CepscmaxB2 = TepscmaxB2;

			// Store committed values of shear aggregate interlock
			Cbeta_pmaxA = Tbeta_pmaxA;
			Cbeta_nmaxA = Tbeta_nmaxA;
			Ctau0_pmaxA = Ttau0_pmaxA;
			Ctau0_nmaxA = Ttau0_nmaxA;
			Cbeta_pmaxB = Tbeta_pmaxB;
			Cbeta_nmaxB = Tbeta_nmaxB;
			Ctau0_pmaxB = Ttau0_pmaxB;
			Ctau0_nmaxB = Ttau0_nmaxB;

			//dowel action
			CDI_A_p = TDI_A_p;
			CDI_A_n = TDI_A_n;
			CDI_B_p = TDI_B_p;
			CDI_B_n = TDI_B_n;

			// Commit State of Uniaxial Materials
			for (int i=0; i < 8; i++)
			{
				theMaterial[i]->commitState();
			}

		}
	}

	else { // 2nd crack .............................................................

		// Stage 3 - 2nd crack (strut B)
		Stage3(Tstrain[0], Tstrain[1], Tstrain[2]);

		// Store committed values of epscmax
		CepscmaxA1 = TepscmaxA1;
		CepscmaxA2 = TepscmaxA2;
		CepscmaxB1 = TepscmaxB1;
		CepscmaxB2 = TepscmaxB2;

		// Store committed values of shear aggregate interlock
		Cbeta_pmaxA = Tbeta_pmaxA;
		Cbeta_nmaxA = Tbeta_nmaxA;
		Ctau0_pmaxA = Ttau0_pmaxA;
		Ctau0_nmaxA = Ttau0_nmaxA;
		Cbeta_pmaxB = Tbeta_pmaxB;
		Cbeta_nmaxB = Tbeta_nmaxB;
		Ctau0_pmaxB = Ttau0_pmaxB;
		Ctau0_nmaxB = Ttau0_nmaxB;

		//dowel action
		CDI_A_p = TDI_A_p;
		CDI_A_n = TDI_A_n;
		CDI_B_p = TDI_B_p;
		CDI_B_n = TDI_B_n;

		// Commit State of Uniaxial Materials
		for (int i=0; i < 8; i++)
		{
			theMaterial[i]->commitState();
		}

	}

	// ...............................................................................

	// Store committed stresses
	CStress(0) = stress_vec(0);
	CStress(1) = stress_vec(1);
	CStress(2) = stress_vec(2);

	// Store committed strains
	CStrain(0) = strain_vec(0);
	CStrain(1) = strain_vec(1);
	CStrain(2) = strain_vec(2);

	// For recorders
	CPanelConcStress(0) = TPanelConcStress(0);
	CPanelConcStress(1) = TPanelConcStress(1);
	CPanelConcStress(2) = TPanelConcStress(2);

	CPanelSteelStress(0) = TPanelSteelStress(0);
	CPanelSteelStress(1) = TPanelSteelStress(1);
	CPanelSteelStress(2) = TPanelSteelStress(2);

	CStrainStressSteel1(0) = TStrainStressSteel1(0); 
	CStrainStressSteel1(1) = TStrainStressSteel1(1); 

	CStrainStressSteel2(0) = TStrainStressSteel2(0); 
	CStrainStressSteel2(1) = TStrainStressSteel2(1); 

	CStrainStressConc1(0) = TStrainStressConc1(0);
	CStrainStressConc1(1) = TStrainStressConc1(1);

	CStrainStressConc2(0) = TStrainStressConc2(0);
	CStrainStressConc2(1) = TStrainStressConc2(1);

	CStrainStressInterlock1(0) = TStrainStressInterlock1(0);
	CStrainStressInterlock1(1) = TStrainStressInterlock1(1);

	CStrainStressInterlock2(0) = TStrainStressInterlock2(0);
	CStrainStressInterlock2(1) = TStrainStressInterlock2(1);

	CCrackingAngles(0) = alfa_crackA;
	CCrackingAngles(1) = alfa_crackB;

	return 0;
}

// Revert to previously converged state
int FAM_CS::revertToLastCommit(void)
{
	for (int i=0; i < 8; i++) {
		theMaterial[i]->revertToLastCommit();
	}

	// Shear aggregate interlock
	Tbeta_pmaxA = Cbeta_pmaxA;
	Tbeta_nmaxA = Cbeta_nmaxA;
	Ttau0_pmaxA = Ctau0_pmaxA;
	Ttau0_nmaxA = Ctau0_nmaxA;
	Tbeta_pmaxB = Cbeta_pmaxB;
	Tbeta_nmaxB = Cbeta_nmaxB;
	Ttau0_pmaxB = Ctau0_pmaxB;
	Ttau0_nmaxB = Ctau0_nmaxB;
	
	// for 2nd cracking criteria
	TepsA2 = CepsA2;

	//dowel action
	TDI_A_p = CDI_A_p;
	TDI_A_n = CDI_A_n;
	TDI_B_p = CDI_B_p;
	TDI_B_n = CDI_B_n;

	return 0;
}

// Calculate trial stress and tangent stiffness based on state
int FAM_CS::determineTrialStressAndTangent(void)
{
	double Tstrain[3];			//ex, ey, gamma

	// Get strain values from strain of element
	Tstrain[0] = strain_vec(0); //ex
	Tstrain[1] = strain_vec(1); //ey
	Tstrain[2] = strain_vec(2); //gxy

	if (crackA == 0 && crackB == 0) {		// Stage 1 - uncracked concrete
		
		Stage1(Tstrain[0], Tstrain[1], Tstrain[2]);
		
	} else if (crackA == 1 && crackB == 0) { // Stage 2 - 1st crack (strut A)

		Stage2(Tstrain[0], Tstrain[1], Tstrain[2]);

	} else {								// Stage 3 - 2nd crack (strut B)

		Stage3(Tstrain[0], Tstrain[1], Tstrain[2]);

	}

	return 0;

}

// Stage 1 - Uncracked concrete
void FAM_CS::Stage1(double &ex, double &ey, double &gamma)
{
	// START - Declaration of temporary variables
	// Principal Strains
	double alfa; // Direction

	// Stress
	double fx; 
	double fy;
	double tauxy;

	// Steel Material ..............................................
	// Strain
	double estxp;
	double estyp;

	// Stress
	double fstx;
	double fsty;
	double taustxy;
	double fstxp;		// Stress of steel layers in X direction
	double fstyp;		// Stress of steel layers in Y direction

	// Stiffness
	double Estxp;
	double Estyp;
	double stifstypF;
	double stifstxpF;

	// Concrete Material .............................................
	// Strains
	double ec1;			// Principal strain in concrete
	double ec2;			// Principal strain in concrete
	double ec12;		// Principal strain in concrete
	// Stresses
	double fc1;			// Concrete stress in direction One
	double fc2;			// Concrete stress in direction Two
	double fcmod1;		// Softened concrete stress in One direction
	double fcmod2;		// Softened concrete stress in Two direction
	double fcx;			// Concrete axial stress in X direction
	double fcy;			// Concrete axial stress in Y direction
	double taucxy;		// Concrete shear stress in X direction
	// Stiffness
	double Ect1;		// Concrete stiffness in direction One
	double Ect2;		// Concrete stiffness in direction Two
	double stifc12;		// Stiffness Parameter
	double stifc11;		// Stiffness Parameter
	double stifc21;		// Stiffness Parameter
	double stifc22;		// Stiffness Parameter
	double stifcu11F;	// Stiffness Parameter
	double stifcu12F;	// Stiffness Parameter
	double stifcu22F;	// Stiffness Parameter
	double stifcu21F;	// Stiffness Parameter
	double Fcu1F;		// Stiffness Parameter
	double Fcu2F;		// Stiffness Parameter

	// Misc ...............................................................
	double TeTaStper;	// Reinforcement Direction

	// Elements of Partial Stiffness Matrix ..........................................
	double R; // Parameter

	double dsxdex; 
	double dsxdey; 
	double dsxdgamma; 
	double dsydex; 
	double dsydey; 
	double dsydgamma; 
	double dtxydex; 
	double dtxydey; 
	double dtxydgamma;

	// Softening/Damage Parameters ....................................
	double beta1;
	double delbeta1;
	double beta2;
	double delbeta2;
	// END - Declaration of temporary variables .......................

	int first_iteration = 0;

	if (ex == 0.0 && ey == 0.0 && gamma == 0.0)
	{ 
		first_iteration = 1; 
	}

	if (gamma == 0.0) {
		gamma = 1e-20;
	}

	// Principal Strain Direction
	alfa = 0.5*atan((gamma)/(ex - ey)); // careful

	// Principal Strains
	double epslon1 = 0.5*(ex+ey) + 0.5*gamma/sin(2.0*alfa); //  strain along concrete strut One - epslon1
	double epslon2 = 0.5*(ex+ey) - 0.5*gamma/sin(2.0*alfa); //  strain along concrete strut Two - epslon2

	// Storing principanl strains and direction
	alpha_strain = alfa;		// principal strain direction
	Tprstrain1 = epslon1;		// principal strain 1 - for cracking criteria
	Tprstrain2 = epslon2;		// principal strain 2 - for cracking criteria

	ec1 = epslon1;
	ec2 = epslon2;
	ec12 = 0.0;

	// Strain transformation for steel reinforcements
	estxp = 0.5*(ex+ey) + 0.5*(ex-ey)*cos(2.0*TeTaSt) + 0.5*gamma*sin(2.0*TeTaSt);

	if (TeTaSt >= 0.0) {
		TeTaStper = TeTaSt - 0.5*pi;
	}
	else {
		TeTaStper = TeTaSt + 0.5*pi;
	}

	estyp = 0.5*(ex+ey) + 0.5*(ex-ey)*cos(2.0*TeTaStper) + 0.5*gamma*sin(2.0*TeTaStper);

	// Get Principal Concrete Stresses
	// Direction One
	theMaterial[2]->setTrialStrain( ec1 );
	fc1 = theMaterial[2]->getStress();
	Ect1 = theMaterial[2]->getTangent();
	
	TStrainStressConc1(0) = ec1;
	TStrainStressConc1(1) = fc1;

	// Get softened concrete One 
	betaf4(ec2, epcc, fc1, Cepscmax2);
	Tepscmax2 = epsiloncmax;
	beta1 = beta;
	delbeta1 = delbeta;

	fcmod1 = beta1 * fc1;
	stifc12 = delbeta1 * fc1;
	stifc11 = beta1 * Ect1;

	// Interlock Stress Strain One (ZERO)
	TStrainStressInterlock1(0) = 0.0;
	TStrainStressInterlock1(1) = 0.0;

	// Direction Two
	theMaterial[3]->setTrialStrain( ec2 );
	fc2 = theMaterial[3]->getStress();
	Ect2 = theMaterial[3]->getTangent();

	TStrainStressConc2(0) = ec2;
	TStrainStressConc2(1) = fc2;

	// Get Softened Concrete Two
	betaf4(ec1, epcc, fc2, Cepscmax1);
	Tepscmax1 = epsiloncmax;
	beta2 = beta;
	delbeta2 = delbeta;

	fcmod2 = beta2 * fc2;
	stifc21 = delbeta2 * fc2;
	stifc22 = beta2 * Ect2;

	// Interlock Stress Strain Two (ZERO)
	TStrainStressInterlock2(0) = 0.0;
	TStrainStressInterlock2(1) = 0.0;

	// Back Transform concrete stresses in X-Y coordinate system
	fcx = 0.5*(fcmod1+fcmod2) + 0.5*(fcmod1-fcmod2)*cos(2.0*alfa);
	fcy = 0.5*(fcmod1+fcmod2) - 0.5*(fcmod1-fcmod2)*cos(2.0*alfa);
	taucxy = 0.5*(fcmod1-fcmod2) * sin(2.0*alfa);

	TPanelConcStress(0) = fcx;
	TPanelConcStress(1) = fcy;
	TPanelConcStress(2) = taucxy;

	stifcu11F = stifc11;
	stifcu12F = stifc12;
	stifcu22F = stifc22;
	stifcu21F = stifc21;
	Fcu1F = fcmod1;
	Fcu2F = fcmod2;

	// Get STEEL Stress and Tangent Stiffness - X Direction
	theMaterial[0]->setTrialStrain( estxp );
	Estxp = theMaterial[0]->getTangent(); // EsX
	fstxp = theMaterial[0]->getStress(); // SigSX
	stifstxpF = Estxp;

	TStrainStressSteel1(0) = estxp;
	TStrainStressSteel1(1) = fstxp;

	// Get STEEL Stress and Tangent Stiffness - Y Direction
	theMaterial[1]->setTrialStrain( estyp );
	Estyp = theMaterial[1]->getTangent(); // EsY
	fstyp = theMaterial[1]->getStress(); // SigSY
	stifstypF = Estyp;
	
	TStrainStressSteel2(0) = estyp;
	TStrainStressSteel2(1) = fstyp;

	// Back Transfor Steel Stresses in X - Y coordinate system
	taustxy = (roux*fstxp-rouy*fstyp)/2.0*sin(2.0*TeTaSt); // Dowell action "taust xp yp" =0

	fstx = (roux*fstxp+rouy*fstyp)/2.0 + ((roux*fstxp-rouy*fstyp)/2.0)*cos(2.0*TeTaSt);

	fsty = (roux*fstxp+rouy*fstyp)/2.0 - ((roux*fstxp-rouy*fstyp)/2.0)*cos(2.0*TeTaSt);

	TPanelSteelStress(0) = fstx;
	TPanelSteelStress(1) = fsty; 
	TPanelSteelStress(2) = taustxy; 

	// Stress Superposition | Concrete + Steel
	fx = fcx+fstx;
	fy = fcy+fsty;
	tauxy = taucxy+taustxy;

	stress_vec(0) = fx;
	stress_vec(1) = fy;
	stress_vec(2) = tauxy;

	// Calculate tangent_matrix
	if (ex == ey) {
		R = 1.0;
	}
	else {
		R = 1.0 + pow(gamma,2.0)/pow((ex - ey),2.0);
	}

	// Calculating partial stiffnesses
	dsxdex = ((pow((1.0 + sqrt(R)),2.0)*stifcu11F)/R + stifcu12F + stifcu21F + stifcu22F +
		((2.0*(Fcu1F - Fcu2F)*pow(gamma,2.0))/pow((ex - ey),3.0) - sqrt(R)*(stifcu12F + stifcu21F - stifcu22F) - 2.0*R*stifcu22F)/pow(R,1.5))/4.0 + 
		(4.0*(stifstxpF*roux - stifstypF*rouy)*cos(2.0*TeTaSt) + (stifstxpF*roux + stifstypF*rouy)*(3.0 + cos(4.0*TeTaSt)))/8.0;

	dsxdey = (((-1.0 + R)*stifcu11F)/R + stifcu12F + stifcu21F + ((2.0*(-Fcu1F + Fcu2F)*pow(gamma,2.0))/pow((ex - ey),3.0) + 2.0*R*(stifcu12F - stifcu21F) + 
		sqrt(R)*(stifcu12F + stifcu21F - stifcu22F))/pow(R,1.5) + stifcu22F)/4.0 + 
		(stifstxpF*roux + stifstypF*rouy)*pow(cos(TeTaSt),2.0)*pow(sin(TeTaSt),2.0);

	dsxdgamma = (gamma*((-2.0*Fcu1F + 2.0*Fcu2F)/((pow((ex - ey),2.0) + pow(gamma,2.0))*sqrt(R)) + stifcu11F/((ex - ey)*sqrt(R)) + stifcu12F/((-ex + ey)*sqrt(R)) +
		stifcu21F/((ex - ey)*sqrt(R)) + stifcu22F/((-ex + ey)*sqrt(R)) + ((ex - ey)*(stifcu11F - stifcu12F - stifcu21F + stifcu22F))/(pow((ex - ey),2.0) + pow(gamma,2.0))))/4.0 + 
		((stifstxpF*roux - stifstypF*rouy + (stifstxpF*roux + stifstypF*rouy)*cos(2.0*TeTaSt))*sin(2.0*TeTaSt))/4.0;

	dsydex = (((-1.0 + R)*stifcu11F)/R + stifcu12F + stifcu21F + ((2.0*(-Fcu1F + Fcu2F)*pow(gamma,2.0))/pow((ex - ey),3.0) - 2.0*R*(stifcu12F - stifcu21F) + 
		sqrt(R)*(stifcu12F + stifcu21F - stifcu22F))/pow(R,1.5) + stifcu22F)/4.0 + 
		(stifstxpF*roux + stifstypF*rouy)*pow(cos(TeTaSt),2.0)*pow(sin(TeTaSt),2.0);

	dsydey = ((pow((-1 + sqrt(R)),2.0)*stifcu11F)/R + stifcu12F + stifcu21F + stifcu22F + 
		((2.0*(Fcu1F - Fcu2F)*pow(gamma,2.0))/pow((ex - ey),3.0) - sqrt(R)*(stifcu12F + stifcu21F - stifcu22F) + 2.0*R*stifcu22F)/pow(R,1.5))/4.0 + 
		((-4.0*stifstxpF*roux + 4.0*stifstypF*rouy)*cos(2.0*TeTaSt) + (stifstxpF*roux + stifstypF*rouy)*(3.0 + cos(4.0*TeTaSt)))/8.0;

	dsydgamma = (gamma*((2.0*(Fcu1F - Fcu2F))/((pow((ex - ey),2.0) + pow(gamma,2.0))*sqrt(R)) + stifcu11F/((ex - ey)*sqrt(R)) + stifcu12F/((-ex + ey)*sqrt(R)) + 
		stifcu21F/((ex - ey)*sqrt(R)) + stifcu22F/((-ex + ey)*sqrt(R)) - ((ex - ey)*(stifcu11F - stifcu12F - stifcu21F + stifcu22F))/(pow((ex - ey),2.0) + pow(gamma,2.0))))/4.0 + 
		-((-(stifstxpF*roux) + stifstypF*rouy + (stifstxpF*roux + stifstypF*rouy)*cos(2.0*TeTaSt))*sin(2.0*TeTaSt))/4.0;

	dtxydex = (gamma*(2.0*(Fcu1F - Fcu2F)*pow(gamma,2.0) - 2.0*(Fcu1F - Fcu2F)*(pow((ex - ey),2.0) + pow(gamma,2.0)) +
		((ex - ey)*(pow((ex - ey),2.0) + pow(gamma,2.0))*(stifcu11F - stifcu12F - stifcu21F + sqrt(R)*(stifcu11F + stifcu12F - stifcu21F - stifcu22F) + 
		stifcu22F))/sqrt(R)))/(4.0*pow((ex - ey),4.0)*pow(R,1.5)) + ((stifstxpF*roux - stifstypF*rouy + (stifstxpF*roux + stifstypF*rouy)*cos(2.0*TeTaSt))*sin(2.0*TeTaSt))/4.0;

	dtxydey = (gamma*(-2.0*(Fcu1F - Fcu2F)*pow(gamma,2.0) + 2.0*(Fcu1F - Fcu2F)*(pow((ex - ey),2.0) + pow(gamma,2.0)) + 
		((ex - ey)*(pow((ex - ey),2.0) + pow(gamma,2.0))*(-stifcu11F + stifcu12F + stifcu21F + sqrt(R)*(stifcu11F + stifcu12F - stifcu21F - stifcu22F) - 
		stifcu22F))/sqrt(R)))/(4.0*pow((ex - ey),4.0)*pow(R,1.5)) + 
		-((-(stifstxpF*roux) + stifstypF*rouy + (stifstxpF*roux + stifstypF*rouy)*cos(2.0*TeTaSt))*sin(2.0*TeTaSt))/4.0;

	dtxydgamma=(2.0*(ex - ey)*(Fcu1F - Fcu2F) + pow(gamma,2.0)*sqrt(R)*(stifcu11F - stifcu12F - stifcu21F + stifcu22F))/
		(4.0*(pow((ex - ey),2.0) + pow(gamma,2.0))*sqrt(R))+ (stifstxpF*roux + stifstypF*rouy)*pow(cos(TeTaSt),2.0)*pow(sin(TeTaSt),2.0);

	// Assemble partial tangent stiffness matrix
	tangent_matrix(0,0) = dsxdex;
	tangent_matrix(0,1) = dsxdey;
	tangent_matrix(0,2) = dsxdgamma;

	tangent_matrix(1,0) = dsydex;
	tangent_matrix(1,1) = dsydey;
	tangent_matrix(1,2) = dsydgamma;

	tangent_matrix(2,0) = dtxydex;
	tangent_matrix(2,1) = dtxydey;
	tangent_matrix(2,2) = dtxydgamma;

	// if 1st iteration calculate initial stiffness

	if (first_iteration == 1) {
		tangent_matrix = this->getInitialTangent();
	}

}

// STAGE 2 - 1st CRACK PANEL MODEL 
void FAM_CS::Stage2(double &ex, double &ey, double &gamma)
{
	// START - Declaration of temporary variables .......................

	// Stress
	double fx;
	double fy;
	double tauxy;

	// Steel Material ..............................................
	// Strain
	double estxp;
	double estyp;

	// Stress
	double fstx;
	double fsty;
	double fstxp;		// Stress of steel layers in X direction
	double fstyp;		// Stress of steel layers in Y direction
	double taustxy;

	// Stiffness
	double Estxp;
	double Estyp;
	double stifstypF;
	double stifstxpF;

	// Concrete Material .............................................
	// Strains
	double ec1;			// Principal strain in concrete
	double ec2;			// Principal strain in concrete

	// Stresses
	double fc1;			// Concrete stress in direction One
	double fc2;			// Concrete stress in direction Two
	double fcmod1;		// Softened concrete stress in One direction
	double fcmod2;		// Softened concrete stress in Two direction
	double fcx;			// Concrete axial stress in X direction
	double fcy;			// Concrete axial stress in Y direction
	double taucxy;		// Concrete shear stress in X direction

	// Stiffness
	double Ect1;		// Concrete stiffness in direction One
	double Ect2;		// Concrete stiffness in direction Two
	double stifc12;		// Stiffness parameter
	double stifc11;		// Stiffness parameter
	double stifc21;		// Stiffness parameter
	double stifc22;		// Stiffness parameter
	double stifcu11F;	// Stiffness parameter
	double stifcu12F;	// Stiffness parameter
	double stifcu22F;	// Stiffness parameter
	double stifcu21F;	// Stiffness parameter
	double Fcu1F;
	double Fcu2F;

	// Softening/Damage Parameters
	double beta1;
	double delbeta1;
	double beta2;
	double delbeta2;

	// Misc ...............................................................
	double TeTaStper;

	// Elements of Partial Stiffness Matrix ..........................................
	double dsxdex; 
	double dsxdey; 
	double dsxdgamma; 
	double dsydex; 
	double dsydey; 
	double dsydgamma; 
	double dtxydex; 
	double dtxydey; 
	double dtxydgamma;

	//New origin for tension
	double Ce0_A;

	// END - Declaration of temporary variables .......................
	
	double alfacrA = alfa_crackA; // alpha_crack = alfacrackold 
	
	// Strain in direction of 1st strut
	double e1 = 0.5*(ex+ey) + 0.5*(ex-ey)*cos(2.0*alfacrA) + 0.5*gamma*sin(2.0*alfacrA);

	double alfacrper; // Direction perpendicular to the crack
	if (alfacrA >= 0.0) {
		alfacrper = alfacrA - 0.5*pi;
	} else {
		alfacrper = alfacrA + 0.5*pi;
	}

	// Strain perpendicular to direction of 1st strut
	double e2 = 0.5*(ex+ey) + 0.5*(ex-ey)*cos(2.0*alfacrper) + 0.5*gamma*sin(2.0*alfacrper);

	// Crack Slip
	double e12 = ( (-1.0*(ex-ey)) * sin(2.0*alfacrA) ) + ( gamma * (cos(2.0*alfacrA)) );

	ec1 = e1;
	ec2 = e2;

	TeA12 = e12;		// History variable for shear aggregate interlock
	TepsA2 = e2;		// Stored for checking if 2nd crack is initiated

	// Strain Transformation for Steel Reinforcements
	estxp = ex;

	if (TeTaSt >= 0.0) {
		TeTaStper = TeTaSt - 0.5*pi;
	}
	else {
		TeTaStper = TeTaSt + 0.5*pi;
	}

	estyp = ey;
	
	// Get Principal Concrete Stresses
	// Concrete Strut A1
	theMaterial[4]->setTrialStrain( ec1 );
	fc1 = theMaterial[4]->getStress();
	Ect1 = theMaterial[4]->getTangent();
	theResponses[2]->getResponse();
	Ce0_A = theResponses[2]->getInformation().theDouble;
	TStrainStressConc1(0) = ec1;
	TStrainStressConc1(1) = fc1;

	// Get Softened Concrete A1 
	betaf4(ec2, epcc, fc1, CepscmaxA2);
	TepscmaxA2 = epsiloncmax;
	beta1 = beta;
	delbeta1 = delbeta;
	
	fcmod1 = beta1 * fc1;
	stifc12 = delbeta1 * fc1;
	stifc11 = beta1 * Ect1;

	// Concrete Strut A2
	theMaterial[5]->setTrialStrain( ec2 );
	fc2 = theMaterial[5]->getStress();
	Ect2 = theMaterial[5]->getTangent();	

	TStrainStressConc2(0) = ec2;
	TStrainStressConc2(1) = fc2;

	// Get Softened Concrete A2
	betaf4(ec1, epcc, fc2, CepscmaxA1);
	TepscmaxA1 = epsiloncmax;
	beta2 = beta;
	delbeta2 = delbeta;

	fcmod2 = beta2 * fc2;
	stifc21 = delbeta2 * fc2;
	stifc22 = beta2 * Ect2;


	// Elastic shear stiffness
	double Gc = 0.4 * Ec;

	// Shear Aggregate Interlock
	double e1_ture = e1 - Ce0_A;

	//Calculate crack shear strain using bisection method
	double lm_A = lm0;   //crack spacing
	if (abs(alfacrA) / (pi / 2) > 8.0 / 9.0)
	{
		lm_A = fmin(sh,lm0);
	}

	Bisection(fpc, lm_A, e1_ture, TeA12, Gc, Cbeta_pmaxA, Ctau0_pmaxA, Cbeta_nmaxA, Ctau0_nmaxA);
	Tgamma_cr_A = gamma_cr;
	InterLocker_improved(fpc, lm_A, e1_ture, Tgamma_cr_A, Gc, Cbeta_pmaxA, Ctau0_pmaxA, Cbeta_nmaxA, Ctau0_nmaxA);
	
	Ttau_Interlock_A = Tau_Interlock;
	// crack -- uncracked concrete
	double dtAggde12 = 1 / (1 / dTau_de12_cr + 1 / Gc);
	double dtAggde1 = dTau_denormal;

	//Maekawa model add by Z
	Tbeta_A = Beta;
	if (Tbeta_A > Cbeta_pmaxA)
	{
		Tbeta_pmaxA = Tbeta_A;
		Ttau0_pmaxA = Tau_Interlock0;
	} else	{
		Tbeta_pmaxA = Cbeta_pmaxA;
		Ttau0_pmaxA = Ctau0_pmaxA;
	}
	if (Tbeta_A < Cbeta_nmaxA)
	{
		Tbeta_nmaxA = Tbeta_A;
		Ttau0_nmaxA = Tau_Interlock0;
	} else	{
		Tbeta_nmaxA = Cbeta_nmaxA;
		Ttau0_nmaxA = Ctau0_nmaxA;
	}
	// Interlock Stress Strain One
	TStrainStressInterlock1(0) = TeA12;
	TStrainStressInterlock1(1) = Ttau_Interlock_A;

	// Interlock Stress Strain Two
	TStrainStressInterlock2(0) = 0.0;
	TStrainStressInterlock2(1) = 0.0;

	// Backtransform concrete stresses
	fcx = (fcmod1+fcmod2)/2.0 + (fcmod1-fcmod2)/2.0*cos(2.0*alfacrA)-Ttau_Interlock_A*(sin(2.0*alfacrA));
	fcy = (fcmod1+fcmod2)/2.0 - (fcmod1-fcmod2)/2.0*cos(2.0*alfacrA)+Ttau_Interlock_A*(sin(2.0*alfacrA));
	taucxy = (fcmod1-fcmod2)/2.0*sin(2.0*alfacrA) + Ttau_Interlock_A*cos(2.0*alfacrA);

	TPanelConcStress(0) = fcx;
	TPanelConcStress(1) = fcy;
	TPanelConcStress(2) = taucxy;

	stifcu11F = stifc11;
	stifcu12F = stifc12;
	stifcu22F = stifc22;
	stifcu21F = stifc21;
	Fcu1F = fcmod1;
	Fcu2F = fcmod2;

	// Get STEEL Stress and Tangent Stiffness - X Direction
	theMaterial[0]->setTrialStrain( estxp );
	Estxp = theMaterial[0]->getTangent(); // EsX
	fstxp = theMaterial[0]->getStress(); // SigSX
	stifstxpF = Estxp;

	TStrainStressSteel1(0) = estxp;
	TStrainStressSteel1(1) = fstxp;

	// Get STEEL Stress and Tangent Stiffness - Y Direction
	theMaterial[1]->setTrialStrain( estyp );
	Estyp = theMaterial[1]->getTangent(); // EsY
	fstyp = theMaterial[1]->getStress(); // SigSY
	stifstypF = Estyp;

	TStrainStressSteel2(0) = estyp;
	TStrainStressSteel2(1) = fstyp;

	// Dowel Action
	// Strain Transformation
	double gamma_dow_A_y = gamma;

	dowel_action_0(gamma_dow_A_y, E0y);
	double Tau_Dowel_0_A_y = Tau_Dowel_0;
	double dTau_dgamma_0_A_y = dTau_dgamma_0;

	// Obtain dowel stresses and stiffness
	dowel_action(fpc, dY, E0y, lm_A, alfacrA, TeA12, e1_ture, CDI_A_p, CDI_A_n);

	double Tau_Dowel_A = Tau_Dowel;
	double dTau_dgamma_A = dTau_dgamma;
	double dTau_deps_A = dTau_deps;

	if (TDI > CDI_A_p)
	{
		TDI_A_p = TDI;
	}
	else {
		TDI_A_p = CDI_A_p;
	}
	if (TDI < -CDI_A_n)
	{
		TDI_A_n = -TDI;
	}
	else {
		TDI_A_n = CDI_A_n;
	}

	// Back Transfor Steel Stresses in X - Y coordinate system
	taustxy = rouy*Tau_Dowel_0_A_y + rouy * Tau_Dowel_A * cos(2.0 * alfacrA); // + dowel stress 
	fstx = roux*fstxp - rouy * Tau_Dowel_A * sin(2.0 * alfacrA); // + dowel stress 
	fsty = rouy*fstyp + rouy * Tau_Dowel_A * sin(2.0 * alfacrA); // + dowel stress

	TPanelSteelStress(0) = fstx;
	TPanelSteelStress(1) = fsty; 
	TPanelSteelStress(2) = taustxy;

	// Stress Super Position | CONCRETE + STEEL
	fx = fcx+fstx;
	fy = fcy+fsty;
	tauxy = taucxy+taustxy;

	stress_vec(0) = fx;
	stress_vec(1) = fy;
	stress_vec(2) = tauxy;

	// Calculate tangent_matrix

	// Calculating partial stiffnesses
	dsxdex = (4.0 * (dtAggde12+rouy*dTau_dgamma_A) + 3.0 * stifcu11F + stifcu12F + stifcu21F + 3.0 * stifcu22F + 4.0 * (stifcu11F - stifcu22F) * cos(2.0 * alfacrA) +
		(-4.0 * (dtAggde12+rouy*dTau_dgamma_A) + stifcu11F - stifcu12F - stifcu21F + stifcu22F) * cos(4.0 * alfacrA) -
		4.0 * (dtAggde1+rouy*dTau_deps_A) * (1.0  +  cos(2.0 * alfacrA)) * sin(2.0 * alfacrA)) / 8.0 +
		(4.0 * (stifstxpF * roux - stifstypF * rouy) * cos(2.0 * TeTaSt) + (stifstxpF * roux + stifstypF * rouy) * (3.0 + cos(4.0 * TeTaSt))) / 8.0;

	dsxdey = stifcu12F * pow(cos(alfacrA), 4.0) + (-4.0 * (dtAggde12+rouy*dTau_dgamma_A) + stifcu11F + stifcu22F) * pow(cos(alfacrA), 2.0) * pow(sin(alfacrA), 2.0) -
		2.0 * (dtAggde1+rouy*dTau_deps_A) * cos(alfacrA) * pow(sin(alfacrA), 3.0) + stifcu21F * pow(sin(alfacrA), 4.0);

	dsxdgamma = -(sin(2.0 * alfacrA) * (-stifcu11F + stifcu12F - stifcu21F + stifcu22F + (4.0 * (dtAggde12+rouy*dTau_dgamma_A) - stifcu11F + stifcu12F + stifcu21F - stifcu22F) * cos(2.0 * alfacrA) +
		2.0 * (dtAggde1+rouy*dTau_deps_A) * sin(2.0 * alfacrA))) / 4.0;

	dsydex = (4.0 * stifcu12F * pow(sin(alfacrA), 3.0) * sin(alfacrA) - 4.0 * (dtAggde12+rouy*dTau_dgamma_A) * pow(sin(2.0 * alfacrA), 2.0) +
		4.0 * pow(cos(alfacrA), 2.0) * (stifcu21F * pow(cos(alfacrA), 2.0) + (stifcu11F + stifcu22F) * pow(sin(alfacrA), 2.0) + (dtAggde1+rouy*dTau_deps_A) * sin(2.0 * alfacrA))) / 4.0;

	dsydey = (4.0 * (dtAggde12+rouy*dTau_dgamma_A) + 3.0 * stifcu11F + stifcu12F + stifcu21F + 3.0 * stifcu22F + 4.0 * (-stifcu11F + stifcu22F) * cos(2.0 * alfacrA) +
		(-4.0 * (dtAggde12+rouy*dTau_dgamma_A) + stifcu11F - stifcu12F - stifcu21F + stifcu22F) * cos(4.0 * alfacrA) + 4.0 * (dtAggde1+rouy*dTau_deps_A) * sin(2.0 * alfacrA) +
		2.0 * (dtAggde1+rouy*dTau_deps_A) * (-1.0) * sin(4.0 * alfacrA)) / 8.0 +
		((-4.0 * stifstxpF * roux + 4.0 * stifstypF * rouy) * cos(2.0 * TeTaSt) + (stifstxpF * roux + stifstypF * rouy) * (3.0 + cos(4.0 * TeTaSt))) / 8.0;

	dsydgamma = (sin(2.0 * alfacrA) * (stifcu11F - stifcu12F + stifcu21F - stifcu22F + (4.0 * (dtAggde12+rouy*dTau_dgamma_A) - stifcu11F + stifcu12F + stifcu21F - stifcu22F) * cos(2.0 * alfacrA) +
		2.0 * (dtAggde1+rouy*dTau_deps_A) * sin(2.0 * alfacrA))) / 4.0;

	dtxydex = (dtAggde1+rouy*dTau_deps_A) * pow(cos(alfacrA), 2.0) * cos(2.0 * alfacrA) + (stifcu11F - stifcu21F) * pow(cos(alfacrA), 3.0) * sin(alfacrA) +
		(stifcu12F - stifcu22F) * cos(alfacrA) * pow(sin(alfacrA), 3.0) - ((dtAggde12+rouy*dTau_dgamma_A) * sin(4.0 * alfacrA)) / 2.0;

	dtxydey = (stifcu12F - stifcu22F) * pow(cos(alfacrA), 3.0) * sin(alfacrA) + (dtAggde1+rouy*dTau_deps_A) * cos(2.0 * alfacrA) * pow(sin(alfacrA), 2.0) +
		(stifcu11F - stifcu21F) * cos(alfacrA) * pow(sin(alfacrA), 3.0) + ((dtAggde12+rouy*dTau_dgamma_A) * sin(4.0 * alfacrA)) / 2.0;

	dtxydgamma = (dtAggde12+rouy*dTau_dgamma_A) * pow(cos(2 * alfacrA), 2.0) + (stifcu11F - stifcu12F - stifcu21F + stifcu22F) * pow(cos(alfacrA), 2.0) * pow(sin(alfacrA), 2.0) + 
		((dtAggde1+rouy*dTau_deps_A) * sin(4.0 * alfacrA)) / 4.0 + rouy * dTau_dgamma_0_A_y; // + dowel stiffness (only for TeTaSt = 0) - if TeTaSt><0 all should be derived

	// Assemble tangent stiffness matrix
	tangent_matrix(0,0) = dsxdex;
	tangent_matrix(0,1) = dsxdey;
	tangent_matrix(0,2) = dsxdgamma;

	tangent_matrix(1,0) = dsydex;
	tangent_matrix(1,1) = dsydey;
	tangent_matrix(1,2) = dsydgamma;

	tangent_matrix(2,0) = dtxydex;
	tangent_matrix(2,1) = dtxydey;
	tangent_matrix(2,2) = dtxydgamma;

}

// STAGE 3 - 2nd CRACK PANEL MODEL
void FAM_CS::Stage3(double &ex, double &ey, double &gamma)
{
	// START - Declaration of temporary variables ....................
	// Stress
	double fx;
	double fy;
	double tauxy;

	// Steel Material ................................................
	// Strain
	double estxp;
	double estyp;

	// Stress
	double fstx;
	double fsty;
	double fstxp;		// Stress of steel layers in X direction
	double fstyp;		// Stress of steel layers in Y direction
	double taustxy;

	// Stiffness
	double Estxp;
	double Estyp;
	double stifstypF;
	double stifstxpF;

	// Concrete Material .............................................
	// Strains
	double ecA1;		// Principal strain in concrete
	double ecA2;		// Principal strain in concrete
	double ecB1;		// Principal strain in concrete
	double ecB2;		// Principal strain in concrete

	// Stresses
	double fcA1;		// Concrete stress in direction One
	double fcA2;		// Concrete stress in direction Two
	double fcB1;		// Concrete stress in direction One
	double fcB2;		// Concrete stress in direction Two
	double fcmodA1;		// Softened concrete stress in One direction
	double fcmodA2;		// Softened concrete stress in Two direction
	double fcmodB1;		// Softened concrete stress in One direction
	double fcmodB2;		// Softened concrete stress in Two direction
	double fcxA;		// Concrete axial stress in X direction
	double fcyA;		// Concrete axial stress in Y direction
	double taucxyA;		// Concrete shear stress in X direction
	double fcxB;		// Concrete axial stress in X direction
	double fcyB;		// Concrete axial stress in Y direction
	double taucxyB;		// Concrete shear stress in X direction
	double fcx;			// Concrete axial stress in X direction
	double fcy;			// Concrete axial stress in Y direction
	double taucxy;		// Concrete shear stress in X direction

	// Stiffness
	double EctA1;		// Concrete stiffness in direction One
	double EctA2;		// Concrete stiffness in direction Two
	double EctB1;		// Concrete stiffness in direction One
	double EctB2;		// Concrete stiffness in direction Two
	double stifcA12;	// Stiffness parameter
	double stifcA11;	// Stiffness parameter
	double stifcA21;	// Stiffness parameter
	double stifcA22;	// Stiffness parameter
	double stifcB12;	// Stiffness parameter
	double stifcB11;	// Stiffness parameter
	double stifcB21;	// Stiffness parameter
	double stifcB22;	// Stiffness parameter
	double stifcu11FA;	// Stiffness parameter
	double stifcu12FA;	// Stiffness parameter
	double stifcu22FA;	// Stiffness parameter
	double stifcu21FA;	// Stiffness parameter
	double Fcu1FA;
	double Fcu2FA;
	double stifcu11FB;  // Stiffness parameter
	double stifcu12FB;	// Stiffness parameter
	double stifcu22FB;	// Stiffness parameter
	double stifcu21FB;	// Stiffness parameter
	double Fcu1FB;
	double Fcu2FB;

	// Softening/Damage Parameters
	double betaA1;
	double delbetaA1;
	double betaA2;
	double delbetaA2;
	double betaB1;
	double delbetaB1;
	double betaB2;
	double delbetaB2;

	// Misc ...........................................................
	double TeTaStper;

	// Elements of Stiffness Matrix ...................................
	double R;

	// Elements of Partial Stiffness Matrix
	double dsxdexA; 
	double dsxdeyA; 
	double dsxdgammaA; 
	double dsydexA; 
	double dsydeyA; 
	double dsydgammaA; 
	double dtxydexA; 
	double dtxydeyA; 
	double dtxydgammaA;

	double dsxdexB; 
	double dsxdeyB; 
	double dsxdgammaB; 
	double dsydexB; 
	double dsydeyB; 
	double dsydgammaB; 
	double dtxydexB; 
	double dtxydeyB; 
	double dtxydgammaB;

	//New origin for tension
	double Ce0_A;
	double Ce0_B;

	// END - Declaration of temporary variables ......................

	double alfacrA = alfa_crackA; // alpha_crackA=alfacrackAold
	double alfacrB = alfa_crackB; // alpha_crackB=alfacrackBold

	// Strain in direction of 1st strut
	double eA1 = 0.5*(ex+ey) + 0.5*(ex-ey)*cos(2.0*alfacrA) + 0.5*gamma*sin(2.0*alfacrA);

	double alfacrperA;
	if (alfacrA >= 0.0) {
		alfacrperA = alfacrA - 0.5*pi;
	} else {
		alfacrperA = alfacrA + 0.5*pi;
	}

	double eA2 = 0.5*(ex+ey) + 0.5*(ex-ey)*cos(2.0*alfacrperA) + 0.5*gamma*sin(2.0*alfacrperA);

	// Strain in direction of 2nd strut
	double eB1 = 0.5*(ex+ey) + 0.5*(ex-ey)*cos(2.0*alfacrB) + 0.5*gamma*sin(2.0*alfacrB);

	double alfacrperB;
	if (alfacrB >= 0.0) {
		alfacrperB = alfacrB - 0.5*pi;
	} else {
		alfacrperB = alfacrB + 0.5*pi;
	}

	double eB2 = 0.5*(ex+ey) + 0.5*(ex-ey)*cos(2.0*alfacrperB) + 0.5*gamma*sin(2.0*alfacrperB);

	// Crack Slip
	double eA12 = ( (-1.0*(ex-ey)) * sin(2.0*alfacrA) ) + ( gamma * (cos(2.0*alfacrA)) );
	double eB12 = ( (-1.0*(ex-ey)) * sin(2.0*alfacrB) ) + ( gamma * (cos(2.0*alfacrB)) );

	ecA1 = eA1;
	ecA2 = eA2;
	TeA12 = eA12;
	
	ecB1 = eB1;
	ecB2 = eB2;
	TeB12 = eB12;

	// Strain Transformation for Steel Reinforcements
	estxp = ex;

	if (TeTaSt >= 0.0) {
		TeTaStper = TeTaSt - 0.5*pi;
	}
	else {
		TeTaStper = TeTaSt + 0.5*pi;
	}

	estyp = 0.5*(ex+ey) + 0.5*(ex-ey)*cos(2.0*TeTaStper) + 0.5*gamma*sin(2.0*TeTaStper);

	// Get Principal Concrete Stresses - Strut A
	// Concrete Strut A1
	theMaterial[4]->setTrialStrain( ecA1 );
	fcA1 = theMaterial[4]->getStress();
	EctA1 = theMaterial[4]->getTangent();
	theResponses[2]->getResponse();
	Ce0_A = theResponses[2]->getInformation().theDouble;

	TStrainStressConc1(0) = ecA1;
	TStrainStressConc1(1) = fcA1;

	// Get Softened Concrete A1 
	betaf4(ecA2, epcc, fcA1, CepscmaxA2);
	TepscmaxA2 = epsiloncmax;
	betaA1 = beta;
	delbetaA1 = delbeta;

	fcmodA1 = betaA1 * fcA1;
	stifcA12 = delbetaA1 * fcA1;
	stifcA11 = betaA1 * EctA1;

	// Concrete Strut A2
	theMaterial[5]->setTrialStrain( ecA2 );
	fcA2 = theMaterial[5]->getStress();
	EctA2 = theMaterial[5]->getTangent();

	TStrainStressConc2(0) = ecA2;
	TStrainStressConc2(1) = fcA2;

	// Get Softened Concrete A2
	betaf4(ecA1, epcc, fcA2, CepscmaxA1); 
	TepscmaxA1 = epsiloncmax;
	betaA2 = beta;
	delbetaA2 = delbeta;

	fcmodA2 = betaA2 * fcA2;
	stifcA21 = delbetaA2 * fcA2;
	stifcA22 = betaA2 * EctA2;

	// Get Principal Concrete Stresses - Strut B
	// Concrete Strut B1
	theMaterial[6]->setTrialStrain( ecB1 );
	fcB1 = theMaterial[6]->getStress();
	EctB1 = theMaterial[6]->getTangent();
	theResponses[3]->getResponse();
	Ce0_B = theResponses[3]->getInformation().theDouble;

	// Get Softened Concrete B1 
	betaf4(ecB2, epcc, fcB1, CepscmaxB2);
	TepscmaxB2 = epsiloncmax;
	betaB1 = beta;
	delbetaB1 = delbeta;

	fcmodB1 = betaB1 * fcB1;
	stifcB12 = delbetaB1 * fcB1;
	stifcB11 = betaB1 * EctB1;

	// Concrete Strut B2
	theMaterial[7]->setTrialStrain( ecB2 );
	fcB2 = theMaterial[7]->getStress();
	EctB2 = theMaterial[7]->getTangent();

	// Get Softened Concrete B2
	betaf4(ecB1, epcc, fcB2, CepscmaxB1); 
	TepscmaxB1 = epsiloncmax;
	betaB2 = beta;
	delbetaB2 = delbeta;

	fcmodB2 = betaB2 * fcB2;
	stifcB21 = delbetaB2 * fcB2;
	stifcB22 = betaB2 * EctB2;

	// Elastic shear stiffness
	double Gc = 0.4 * Ec;

	double eA1_ture = eA1 - Ce0_A;
	double eB1_ture = eB1 - Ce0_B;
	double lm_A = lm0;   //crack spacing
	double lm_B = lm0;   //crack spacing
	if (abs(alfacrA) / (pi / 2) > 8.0 / 9.0)
	{
		lm_A = fmin(sh, lm0);
	}
	if (abs(alfacrB) / (pi / 2) > 8.0 / 9.0)
	{
		lm_B = fmin(sh, lm0);
	}
	
	double dtAggde12A;
	double dtAggde1A;
	double dtAggde12B;
	double dtAggde1B;

	//the softer crack is used
	if (eB1_ture <= eA1_ture) {
		// Shear Aggregate Interlock - Strut A
		//Calculate crack shear strain using bisection method
		// crack -- uncracked concrete
		Bisection(fpc, lm_A, eA1_ture, TeA12, Gc, Cbeta_pmaxA, Ctau0_pmaxA, Cbeta_nmaxA, Ctau0_nmaxA);
		Tgamma_cr_A = gamma_cr;
		InterLocker_improved(fpc, lm_A, eA1_ture, Tgamma_cr_A, Gc, Cbeta_pmaxA, Ctau0_pmaxA, Cbeta_nmaxA, Ctau0_nmaxA);

		Ttau_Interlock_A = Tau_Interlock;
		// crack -- uncracked concrete
		dtAggde12A = 1 / (1 / dTau_de12_cr + 1 / Gc);
		dtAggde1A = dTau_denormal;

		Ttau_Interlock_B = 0;
		dtAggde12B = 0;
		dtAggde1B = 0;
		//Maekawa model add by Z
		Tbeta_A = Beta;
		if (Tbeta_A > Cbeta_pmaxA)
		{
			Tbeta_pmaxA = Tbeta_A;
			Ttau0_pmaxA = Tau_Interlock0;
		} else {
			Tbeta_pmaxA = Cbeta_pmaxA;
			Ttau0_pmaxA = Ctau0_pmaxA;
		}
		if (Tbeta_A < Cbeta_nmaxA)
		{
			Tbeta_nmaxA = Tbeta_A;
			Ttau0_nmaxA = Tau_Interlock0;
		} else {
			Tbeta_nmaxA = Cbeta_nmaxA;
			Ttau0_nmaxA = Ctau0_nmaxA;
		}
	}
	else if (eB1_ture > eA1_ture) {
		Bisection(fpc, lm_B, eB1_ture, TeB12, Gc, Cbeta_pmaxB, Ctau0_pmaxB, Cbeta_nmaxB, Ctau0_nmaxB);
		Tgamma_cr_B = gamma_cr;
		// Shear Aggregate Interlock - Strut B
		//Calculate crack shear strain using bisection method
		InterLocker_improved(fpc, lm_B, eB1_ture, Tgamma_cr_B, Gc, Cbeta_pmaxB, Ctau0_pmaxB, Cbeta_nmaxB, Ctau0_nmaxB);

		Ttau_Interlock_B = Tau_Interlock;
		// crack -- uncracked concrete
		dtAggde12B = 1 / (1 / dTau_de12_cr + 1 / Gc);
		dtAggde1B = dTau_denormal;

		Ttau_Interlock_A = 0;
		dtAggde12A = 0;
		dtAggde1A = 0;
		//Maekawa model add by Z
		Tbeta_B = Beta;
		if (Tbeta_B > Cbeta_pmaxB)
		{
			Tbeta_pmaxB = Tbeta_B;
			Ttau0_pmaxB = Tau_Interlock0;
		} else {
			Tbeta_pmaxB = Cbeta_pmaxB;
			Ttau0_pmaxB = Ctau0_pmaxB;
		}
		if (Tbeta_B < Cbeta_nmaxB)
		{
			Tbeta_nmaxB = Tbeta_B;
			Ttau0_nmaxB = Tau_Interlock0;
		} else {
			Tbeta_nmaxB = Cbeta_nmaxB;
			Ttau0_nmaxB = Ctau0_nmaxB;
		}
	}

	// Interlock Stress Strain - Strut A
	TStrainStressInterlock1(0) = TeA12;
	TStrainStressInterlock1(1) = Ttau_Interlock_A;

	// Backtransform concrete stresses - Strut A
	fcxA = (0.0 + fcmodA2) / 2.0 + (((0.0 - fcmodA2) / 2.0) * cos(2.0 * alfacrA)) - Ttau_Interlock_A * sin(2.0 * alfacrA);
	fcyA = (0.0 + fcmodA2) / 2.0 - (((0.0 - fcmodA2) / 2.0) * cos(2.0 * alfacrA)) + Ttau_Interlock_A * sin(2.0 * alfacrA);
	taucxyA = (0.0 - fcmodA2) / 2.0 * sin(2.0 * alfacrA) + Ttau_Interlock_A * cos(2.0 * alfacrA);

	stifcu11FA = 0.0;
	stifcu12FA = 0.0;
	stifcu21FA = stifcA21;
	stifcu22FA = stifcA22;
	Fcu1FA = 0.0;
	Fcu2FA = fcmodA2;
	dtAggde1A = 0.0;

	// Interlock Stress Strain - Strut B
	TStrainStressInterlock2(0) = TeB12;
	TStrainStressInterlock2(1) = Ttau_Interlock_B;

	// Backtransform concrete stresses  - Strut B
	fcxB = (0.0+fcmodB2)/2.0 + (0.0-fcmodB2)/2.0*cos(2.0*alfacrB)-Ttau_Interlock_B*sin(2.0*alfacrB);
	fcyB = (0.0+fcmodB2)/2.0 - (0.0-fcmodB2)/2.0*cos(2.0*alfacrB)+Ttau_Interlock_B*sin(2.0*alfacrB);
	taucxyB = (0.0-fcmodB2)/2.0*sin(2.0*alfacrB) + Ttau_Interlock_B*cos(2.0*alfacrB);

	stifcu11FB = 0.0;
	stifcu12FB = 0.0;
	stifcu21FB = stifcB21;
	stifcu22FB = stifcB22;
	Fcu1FB = 0.0;
	Fcu2FB = fcmodB2;
	dtAggde1B = 0.0;

	// CONCRETE Stress Superposition - Strut A + Strut B
	taucxy = taucxyA + taucxyB;
	fcx = fcxA + fcxB;
	fcy = fcyA + fcyB;

	TPanelConcStress(0) = fcx;
	TPanelConcStress(1) = fcy;
	TPanelConcStress(2) = taucxy;

	// Get STEEL Stress and Tangent Stiffness - X Direction
	theMaterial[0]->setTrialStrain( estxp );
	Estxp = theMaterial[0]->getTangent(); // EsX
	fstxp = theMaterial[0]->getStress(); // SigSX
	stifstxpF = Estxp;

	TStrainStressSteel1(0) = estxp;
	TStrainStressSteel1(1) = fstxp;

	// Get STEEL Stress and Tangent Stiffness - Y Direction
	theMaterial[1]->setTrialStrain( estyp );
	Estyp = theMaterial[1]->getTangent(); // EsY
	fstyp = theMaterial[1]->getStress(); // SigSY
	stifstypF = Estyp;

	TStrainStressSteel2(0) = estyp;
	TStrainStressSteel2(1) = fstyp;

	// Dowel Action - Crack A - only on the horizontal plane (for now)
	double gamma_dow_A_x = 0;
	double gamma_dow_A_y = gamma;

	// Obtain dowel stresses and stiffness
	dowel_action_0(gamma_dow_A_x, E0x); // =0 for TeTaSt = 0 - if TeTaSt><0 this should be figure out
	double Tau_Dowel_0_A_x = Tau_Dowel_0;
	double dTau_dgamma_0_A_x = dTau_dgamma_0;

	dowel_action_0(gamma_dow_A_y, E0y);
	double Tau_Dowel_0_A_y = Tau_Dowel_0;
	double dTau_dgamma_0_A_y = dTau_dgamma_0;
	
	// Crack B
	double gamma_dow_B_x = 0;
	double gamma_dow_B_y = gamma;

	// Obtain dowel stresses and stiffness
	dowel_action_0(gamma_dow_B_x, E0x); // =0 for TeTaSt = 0 - if TeTaSt><0 this should be figure out
	double Tau_Dowel_0_B_x = Tau_Dowel_0;
	double dTau_dgamma_0_B_x = dTau_dgamma_0;

	dowel_action_0(gamma_dow_B_y, E0y);
	double Tau_Dowel_0_B_y = Tau_Dowel_0;
	double dTau_dgamma_0_B_y = dTau_dgamma_0;

	double Tau_Dowel_A;
	double dTau_dgamma_A;
	double dTau_deps_A;
	double Tau_Dowel_B;
	double dTau_dgamma_B;
	double dTau_deps_B;

	// the softer dowel action is used
	if (eB1_ture <= eA1_ture) {
		// Obtain dowel stresses and stiffness
		dowel_action(fpc, dY, E0y, lm_A, alfacrA, TeA12, eA1_ture, CDI_A_p, CDI_A_n);
		Tau_Dowel_A = Tau_Dowel;
		dTau_dgamma_A = dTau_dgamma;
		dTau_deps_A = dTau_deps;

		//dTau_deps_A_y = 0;

		if (TDI > CDI_A_p)
		{
			TDI_A_p = TDI;
		}
		else {
			TDI_A_p = CDI_A_p;
		}
		if (TDI < -CDI_A_n)
		{
			TDI_A_n = -TDI;
		}
		else {
			TDI_A_n = CDI_A_n;
		}
		Tau_Dowel_B = 0.0;
		dTau_dgamma_B = 0.0;
		dTau_deps_B = 0.0;
	}
	else if (eB1_ture > eA1_ture) {
		// Obtain dowel stresses and stiffness
		dowel_action(fpc, dY, E0y, lm_B, alfacrB, TeB12, eB1_ture, CDI_B_p, CDI_B_n);
		Tau_Dowel_B = Tau_Dowel;
		dTau_dgamma_B = dTau_dgamma;
		dTau_deps_B = dTau_deps;

		if (TDI > CDI_B_p)
		{
			TDI_B_p = TDI;
		}
		else {
			TDI_B_p = CDI_B_p;
		}
		if (TDI < -CDI_B_n)
		{
			TDI_B_n = -TDI;
		}
		else {
			TDI_B_n = CDI_B_n;
		}
		Tau_Dowel_A = 0.0;
		dTau_dgamma_A = 0.0;
		dTau_deps_A = 0.0;
	}

	//  Back Transformed Steel Stresses 
	taustxy = rouy * (Tau_Dowel_0_A_y + Tau_Dowel_0_B_y)*cos(2.0 * TeTaSt) + 
		rouy * Tau_Dowel_A * cos(2.0 * alfacrA) + rouy * Tau_Dowel_B * cos(2.0 * alfacrB); // + dowel stress 
	fstx = roux * fstxp - rouy * Tau_Dowel_A * sin(2.0 * alfacrA) - rouy * Tau_Dowel_B * sin(2.0 * alfacrB); // + dowel stress 
	fsty = rouy * fstyp + rouy * Tau_Dowel_A * sin(2.0 * alfacrA) + rouy * Tau_Dowel_B * sin(2.0 * alfacrB); // + dowel stress

	TPanelSteelStress(0) = fstx;
	TPanelSteelStress(1) = fsty; 
	TPanelSteelStress(2) = taustxy;

	// Stress Super Position | CONCRETE + STEEL
	fx = fcx+fstx; // fx
	fy = fcy+fsty; // fy
	tauxy = taucxy+taustxy; // tauxy

	stress_vec(0) = fx;
	stress_vec(1) = fy;
	stress_vec(2) = tauxy;
	
	// Calculate tangent_matrix
	if (ex == ey) {
		R = 1.0;
	}
	else {
		R = 1.0 + pow(gamma,2.0)/pow((ex - ey),2.0);
	}

	// Calculating partial stiffnesses
	// Strut A
	dsxdexA = (4.0 * (dtAggde12A+rouy*dTau_dgamma_A) + 3.0 * stifcu11FA + stifcu12FA + stifcu21FA + 3.0 * stifcu22FA + 4.0 * (stifcu11FA - stifcu22FA) * cos(2.0 * alfacrA) +
		(-4.0 * (dtAggde12A+rouy*dTau_dgamma_A) + stifcu11FA - stifcu12FA - stifcu21FA + stifcu22FA) * cos(4.0 * alfacrA) -
		4.0 * (dtAggde1A+rouy*dTau_deps_A) * (1.0 + cos(2.0 * alfacrA)) * sin(2.0 * alfacrA)) / 8.0 +
		(4.0 * (stifstxpF * roux - stifstypF * rouy) * cos(2.0 * TeTaSt) + (stifstxpF * roux + stifstypF * rouy) * (3.0 + cos(4.0 * TeTaSt))) / 8.0;

	dsxdeyA = stifcu12FA * pow(cos(alfacrA), 4.0) + (-4.0 * (dtAggde12A+rouy*dTau_dgamma_A) + stifcu11FA + stifcu22FA) * pow(cos(alfacrA), 2.0) * pow(sin(alfacrA), 2.0) -
		2.0 * (dtAggde1A+rouy*dTau_deps_A) * cos(alfacrA) * pow(sin(alfacrA), 3.0) + stifcu21FA * pow(sin(alfacrA), 4.0);

	dsxdgammaA = -(sin(2.0 * alfacrA) * (-stifcu11FA + stifcu12FA - stifcu21FA + stifcu22FA + (4.0 * (dtAggde12A+rouy*dTau_dgamma_A) - stifcu11FA + stifcu12FA + stifcu21FA - stifcu22FA) * cos(2.0 * alfacrA) +
		2.0 * (dtAggde1A+rouy*dTau_deps_A) * sin(2.0 * alfacrA))) / 4.0;

	dsydexA = (4.0 * stifcu12FA * pow(sin(alfacrA), 3.0) * sin(alfacrA) - 4.0 * (dtAggde12A+rouy*dTau_dgamma_A) * pow(sin(2.0 * alfacrA), 2.0) +
		4.0 * pow(cos(alfacrA), 2.0) * (stifcu21FA * pow(cos(alfacrA), 2.0) + (stifcu11FA + stifcu22FA) * pow(sin(alfacrA), 2.0) + (dtAggde1A+rouy*dTau_deps_A) * sin(2.0 * alfacrA))) / 4.0;

	dsydeyA = (4.0 * (dtAggde12A+rouy*dTau_dgamma_A) + 3.0 * stifcu11FA + stifcu12FA + stifcu21FA + 3.0 * stifcu22FA + 4.0 * (-stifcu11FA + stifcu22FA) * cos(2.0 * alfacrA) +
		(-4.0 * (dtAggde12A+rouy*dTau_dgamma_A) + stifcu11FA - stifcu12FA - stifcu21FA + stifcu22FA) * cos(4.0 * alfacrA) + 4.0 * (dtAggde1A+rouy*dTau_deps_A) * sin(2.0 * alfacrA) -
		2.0 * (dtAggde1A+rouy*dTau_deps_A) * sin(4.0 * alfacrA)) / 8.0 +
		((-4.0 * stifstxpF * roux + 4.0 * stifstypF * rouy) * cos(2.0 * TeTaSt) + (stifstxpF * roux + stifstypF * rouy) * (3.0 + cos(4.0 * TeTaSt))) / 8.0;

	dsydgammaA = (sin(2.0 * alfacrA) * (stifcu11FA - stifcu12FA + stifcu21FA - stifcu22FA + (4.0 * (dtAggde12A+rouy*dTau_dgamma_A) - stifcu11FA + stifcu12FA + stifcu21FA - stifcu22FA) * cos(2.0 * alfacrA) +
		2.0 * (dtAggde1A+rouy*dTau_deps_A) * sin(2.0 * alfacrA))) / 4.0;

	dtxydexA = (dtAggde1A+rouy*dTau_deps_A) * pow(cos(alfacrA), 2.0) * cos(2.0 * alfacrA) + (stifcu11FA - stifcu21FA) * pow(cos(alfacrA), 3.0) * sin(alfacrA) +
		(stifcu12FA - stifcu22FA) * cos(alfacrA) * pow(sin(alfacrA), 3.0) - ((dtAggde12A+rouy*dTau_dgamma_A) * sin(4.0 * alfacrA)) / 2.0;

	dtxydeyA = (stifcu12FA - stifcu22FA) * pow(cos(alfacrA), 3.0) * sin(alfacrA) + (dtAggde1A+rouy*dTau_deps_A) * cos(2.0 * alfacrA) * pow(sin(alfacrA), 2.0) +
		(stifcu11FA - stifcu21FA) * cos(alfacrA) * pow(sin(alfacrA), 3.0) + ((dtAggde12A+rouy*dTau_dgamma_A) * sin(4.0 * alfacrA)) / 2.0; // + dowel stiffness (only for TeTaSt = 0) - if TeTaSt><0 all should be derived

	dtxydgammaA = (dtAggde12A+rouy*dTau_dgamma_A) * pow(cos(2.0 * alfacrB), 2.0) + (stifcu11FA - stifcu12FA - stifcu21FA + stifcu22FA) * pow(cos(alfacrA), 2.0) * pow(sin(alfacrA), 2.0) + 
		((dtAggde1A+rouy*dTau_deps_A) * sin(4.0 * alfacrA)) / 4.0 + rouy * dTau_dgamma_0_A_y; // + dowel stiffness (only for TeTaSt = 0) - if TeTaSt><0 all should be derived

	// Strut B
	dsxdexB = (4.0 * (dtAggde12B+rouy*dTau_dgamma_B) + 3.0 * stifcu11FB + stifcu12FB + stifcu21FB + 3.0 * stifcu22FB + 4.0 * (stifcu11FB - stifcu22FB) * cos(2.0 * alfacrB) +
		(-4.0 * (dtAggde12B+rouy*dTau_dgamma_B) + stifcu11FB - stifcu12FB - stifcu21FB + stifcu22FB) * cos(4.0 * alfacrB) -
		4.0 * (dtAggde1B+rouy*dTau_deps_B) * (1.0 + cos(2.0 * alfacrB)) * sin(2.0 * alfacrB)) / 8.0;

	dsxdeyB = stifcu12FB * pow(cos(alfacrB), 4.0) + (-4.0 * (dtAggde12B+rouy*dTau_dgamma_B) + stifcu11FB + stifcu22FB) * pow(cos(alfacrB), 2.0) * pow(sin(alfacrB), 2.0) -
		2.0 * (dtAggde1B+rouy*dTau_deps_B) * cos(alfacrB) * pow(sin(alfacrB), 3.0) + stifcu21FB * pow(sin(alfacrB), 4.0);

	dsxdgammaB = -(sin(2.0 * alfacrB) * (-stifcu11FB + stifcu12FB - stifcu21FB + stifcu22FB + (4.0 * (dtAggde12B+rouy*dTau_dgamma_B) - stifcu11FB + stifcu12FB + stifcu21FB - stifcu22FB) * cos(2.0 * alfacrB) +
		2.0 * (dtAggde1B+rouy*dTau_deps_B) * sin(2.0 * alfacrB))) / 4.0;

	dsydexB = (4.0 * stifcu12FB * pow(sin(alfacrB), 3.0) * sin(alfacrB) - 4.0 * (dtAggde12B+rouy*dTau_dgamma_B) * pow(sin(2.0 * alfacrB), 2.0) +
		4.0 * pow(cos(alfacrB), 2.0) * (stifcu21FB * pow(cos(alfacrB), 2.0) + (stifcu11FB + stifcu22FB) * pow(sin(alfacrB), 2.0) + (dtAggde1B+rouy*dTau_deps_B) * sin(2.0 * alfacrB))) / 4.0;

	dsydeyB = (4.0 * (dtAggde12B+rouy*dTau_dgamma_B) + 3.0 * stifcu11FB + stifcu12FB + stifcu21FB + 3.0 * stifcu22FB + 4.0 * (-stifcu11FB + stifcu22FB) * cos(2.0 * alfacrB) +
		(-4.0 * (dtAggde12B+rouy*dTau_dgamma_B) + stifcu11FB - stifcu12FB - stifcu21FB + stifcu22FB) * cos(4.0 * alfacrB) + 4.0 * (dtAggde1B+rouy*dTau_deps_B) * sin(2.0 * alfacrB) -
		2.0 * (dtAggde1B+rouy*dTau_deps_B) * sin(4.0 * alfacrB)) / 8.0;

	dsydgammaB = (sin(2.0 * alfacrB) * (stifcu11FB - stifcu12FB + stifcu21FB - stifcu22FB + (4.0 * (dtAggde12B+rouy*dTau_dgamma_B) - stifcu11FB + stifcu12FB + stifcu21FB - stifcu22FB) * cos(2.0 * alfacrB) +
		2.0 * (dtAggde1B+rouy*dTau_deps_B) * sin(2.0 * alfacrB))) / 4.0;

	dtxydexB = (dtAggde1B+rouy*dTau_deps_B) * pow(cos(alfacrB), 2.0) * cos(2.0 * alfacrB) + (stifcu11FB - stifcu21FB) * pow(cos(alfacrB), 3.0) * sin(alfacrB) +
		(stifcu12FB - stifcu22FB) * cos(alfacrB) * pow(sin(alfacrB), 3.0) - ((dtAggde12B+rouy*dTau_dgamma_B) * sin(4.0 * alfacrB)) / 2.0;

	dtxydeyB = (stifcu12FB - stifcu22FB) * pow(cos(alfacrB), 3.0) * sin(alfacrB) + (dtAggde1B+rouy*dTau_deps_B) * cos(2.0 * alfacrB) * pow(sin(alfacrB), 2.0) +
		(stifcu11FB - stifcu21FB) * cos(alfacrB) * pow(sin(alfacrB), 3.0) + ((dtAggde12B+rouy*dTau_dgamma_B) * sin(4.0 * alfacrB)) / 2.0;

	dtxydgammaB = (dtAggde12B+rouy*dTau_dgamma_B) * pow(cos(2.0 * alfacrB), 2.0) + (stifcu11FB - stifcu12FB - stifcu21FB + stifcu22FB) * pow(cos(alfacrB), 2.0) * pow(sin(alfacrB), 2.0) +
		((dtAggde1B+rouy*dTau_deps_B) * sin(4.0 * alfacrB)) / 4.0 + rouy * dTau_dgamma_0_B_y; // + dowel stiffness (only for TeTaSt = 0) - if TeTaSt><0 all should be derived

	// Final tangent Stiffness = Strut A + Strut B
	tangent_matrix(0,0) = dsxdexA + dsxdexB;
	tangent_matrix(0,1) = dsxdeyA + dsxdeyB;
	tangent_matrix(0,2) = dsxdgammaA + dsxdgammaB;

	tangent_matrix(1,0) = dsydexA + dsydexB;
	tangent_matrix(1,1) = dsydeyA + dsydeyB;
	tangent_matrix(1,2) = dsydgammaA + dsydgammaB;

	tangent_matrix(2,0) = dtxydexA + dtxydexB;
	tangent_matrix(2,1) = dtxydeyA + dtxydeyB;
	tangent_matrix(2,2) = dtxydgammaA + dtxydgammaB;
}

// Shear Aggregate Interlock Model
void FAM_CS::InterLocker_improved(double &fpc, double &crack_spacing, double &e_cr_normal_true, double &e_cr_parallel, double &Gc, double &beta_pmax, double &tau0_pmax, double &beta_nmax, double &tau0_nmax)
{
	//double slope = 0.4*Ec; // shear stiffness coefficient G (fixed value)	
	double a = 8.0E-5; // The value of e_cr_normal_transform when (e_cr_normal_true = 0)
	double e_cr_normal_transform = 0.5 * (e_cr_normal_true + pow(e_cr_normal_true * e_cr_normal_true + 4 * a * a,0.5));
	//double e_cr_normal_true = fmax(e_cr_normal, 1.0E-10);
	double beta_m = e_cr_parallel / e_cr_normal_transform;			//  Maekawa model
	double crack_width = e_cr_normal_transform * crack_spacing;  // crack width with transformed normal strain
	double K = fmax(0, 1 - exp(1-Gamax/2/crack_width));   //  factor to consider crack width

	Beta = beta_m;

	//if (tension_tan == 0.0)
	//{
	//	tension_tan = -1.0E-20;
	//}
	//if (Beta == 0.0)
	Tau_Interlock = 0.0;
	dTau_de12_cr = 0.0;
	dTau_denormal = 0.0;

	if (beta_m > 0.0)
	{			
		double f_st = 3.8 * pow(-fpc, 1.0 / 3);         //shear stress capacity
		
		if (beta_m >= beta_pmax)
		{
			Tau_Interlock0 = f_st * beta_m * beta_m / (1 + beta_m * beta_m);
			Tau_Interlock = K * f_st * beta_m * beta_m / (1 + beta_m * beta_m);
			dTau_de12_cr = K / e_cr_normal_transform * 2.0 * beta_m * f_st / ((1.0 + beta_m * beta_m) * (1.0 + beta_m * beta_m));
			dTau_denormal = (Gamax * (K - 1) / 4.0 / crack_width - K / (1.0 + beta_m * beta_m) ) *
				f_st * beta_m * beta_m / e_cr_normal_transform / (1.0 + beta_m * beta_m) *
				(1.0 + e_cr_normal_true / pow(e_cr_normal_true * e_cr_normal_true + 4.0 * a * a, 0.5));
			
			if (e_cr_parallel >= 0.004)
			{
				dTau_de12_cr = dTau_de12_cr * pow(0.004 / e_cr_parallel, 0.4) + Tau_Interlock * (-0.4 / 0.004) * pow(0.004 / e_cr_parallel, 1.4);
				Tau_Interlock0 = Tau_Interlock0 * pow(0.004 / e_cr_parallel, 0.4);
				Tau_Interlock = Tau_Interlock * pow(0.004 / e_cr_parallel, 0.4);
				dTau_denormal = dTau_denormal * pow(0.004 / e_cr_parallel, 0.4);
			}
		}
		else if (beta_m >= 0.9 * beta_pmax)
		{
			Tau_Interlock = K * tau0_pmax/(0.15* beta_pmax)*(beta_m - 0.85*beta_pmax);
			dTau_de12_cr = K / e_cr_normal_transform * tau0_pmax / (0.15 * beta_pmax);
			dTau_denormal = (Gamax * (K - 1) / crack_width * (beta_m - 0.85 * beta_pmax) - 2 * K * beta_m) *
				tau0_pmax / (0.6 * e_cr_normal_transform * beta_pmax) *
				(1.0 + e_cr_normal_true / pow(e_cr_normal_true * e_cr_normal_true + 4.0 * a * a, 0.5));
		}
		else
		{
			Tau_Interlock = K * tau0_pmax / 3 * pow(beta_m / (0.9 * beta_pmax), 9);
			dTau_de12_cr = K / e_cr_normal_transform * tau0_pmax / (0.3 * beta_pmax) * pow(beta_m / (0.9 * beta_pmax), 8);
			dTau_denormal = (Gamax * (K - 1) / (18 * crack_width) - K) * tau0_pmax / (0.6 * beta_pmax) *
				e_cr_parallel / (e_cr_normal_transform * e_cr_normal_transform) * pow(beta_m / (0.9 * beta_pmax), 8) *
				(1.0 + e_cr_normal_true / pow(e_cr_normal_true * e_cr_normal_true + 4.0 * a * a, 0.5));
		}
	}
	else if (beta_m < 0.0)
	{
		double f_st = -3.8 * pow(-fpc, 1.0 / 3);     //shear stress capacity
			
		if (beta_m <= beta_nmax)
		{
			Tau_Interlock0 = f_st * beta_m * beta_m / (1.0 + beta_m * beta_m);
			Tau_Interlock = K * f_st * beta_m * beta_m / (1.0 + beta_m * beta_m);
			dTau_de12_cr = K / e_cr_normal_transform * 2.0 * beta_m * f_st / ((1 + beta_m * beta_m) * (1.0 + beta_m * beta_m));
			dTau_denormal = (Gamax * (K - 1) / 4.0 / crack_width - K / (1.0 + beta_m * beta_m)) *
				f_st * beta_m * beta_m / e_cr_normal_transform / (1.0 + beta_m * beta_m) *
				(1.0 + e_cr_normal_true / pow(e_cr_normal_true * e_cr_normal_true + 4.0 * a * a, 0.5));

			if (e_cr_parallel <= -0.004)
			{
				dTau_de12_cr = dTau_de12_cr * pow(-0.004 / e_cr_parallel, 0.4) + Tau_Interlock * (-0.4 / 0.004) * pow(-0.004 / e_cr_parallel, 1.4);
				Tau_Interlock0 = Tau_Interlock0 * pow(-0.004 / e_cr_parallel, 0.4);
				Tau_Interlock = Tau_Interlock * pow(-0.004 / e_cr_parallel, 0.4);
				dTau_denormal = dTau_denormal * pow(-0.004 / e_cr_parallel, 0.4);
			}
		}
		else if (beta_m <= 0.9 * beta_nmax)
		{
			Tau_Interlock = K * tau0_nmax / (0.15 * beta_nmax) * (beta_m - 0.85 * beta_nmax);
			dTau_de12_cr = K / e_cr_normal_transform * tau0_nmax / (0.15 * beta_nmax);
			dTau_denormal = (Gamax * (K - 1.0) / crack_width * (beta_m - 0.85 * beta_nmax) - 2 * K * beta_m) *
				tau0_nmax / (0.6 * e_cr_normal_transform * beta_nmax) *
				(1.0 + e_cr_normal_true / pow(e_cr_normal_true * e_cr_normal_true + 4.0 * a * a, 0.5));
		}
		else
		{
			Tau_Interlock = K * tau0_nmax / 3 * pow(beta_m / (0.9 * beta_nmax), 9);
			dTau_de12_cr = K / e_cr_normal_transform * tau0_nmax / 0.3 / beta_nmax * pow(beta_m / (0.9 * beta_nmax), 8);
			dTau_denormal = (Gamax * (K - 1) / (18 * crack_width) - K) * tau0_nmax / (0.6 * beta_nmax) *
				e_cr_parallel / (e_cr_normal_transform * e_cr_normal_transform) * pow(beta_m / (0.9 * beta_nmax), 8) *
				(1.0 + e_cr_normal_true / pow(e_cr_normal_true * e_cr_normal_true + 4.0 * a * a, 0.5));
		}
	}
	if (crack_width > Gamax / 2.0)
	{
		dTau_denormal = 0.0;
	}
}

// Bisection method, gamma = gamma_cr + gamma_c, total shear strain --> crack shear strain and uncracked concrete shear strain
void FAM_CS::Bisection(double &fpc, double &crack_spacing, double &e_cr_normal_true, double &gamma, double &Gc, double &beta_pmax, double &tau0_pmax, double &beta_nmax, double &tau0_nmax)
{
	double gamma_cr_down = fmin(0.0, gamma);
	double gamma_cr_up = fmax(0.0, gamma);
    double tol = 1.0E-8;  // allowable error
	InterLocker_improved(fpc, crack_spacing, e_cr_normal_true, gamma_cr_down, Gc, beta_pmax, tau0_pmax, beta_nmax, tau0_nmax);
	double tau_down = Tau_Interlock;
	InterLocker_improved(fpc, crack_spacing, e_cr_normal_true, gamma_cr_up, Gc, beta_pmax, tau0_pmax, beta_nmax, tau0_nmax);
	double tau_up = Tau_Interlock;
	double gamma_down = tau_down / Gc + gamma_cr_down;
	double gamma_up = tau_up / Gc + gamma_cr_up;

	while ((((gamma_down - gamma) * (gamma_up - gamma)) < 0.0) && (abs(gamma_cr_up - gamma_cr_down) > tol)) 
	{
		double gamma_cr_middle = (gamma_cr_down + gamma_cr_up) / 2.0;
		InterLocker_improved(fpc, crack_spacing, e_cr_normal_true, gamma_cr_middle,Gc, beta_pmax, tau0_pmax, beta_nmax, tau0_nmax);
		double tau_middle = Tau_Interlock;
		double gamma_middle = tau_middle / Gc + gamma_cr_middle;
		if (((gamma_down - gamma) * (gamma_middle - gamma)) < 0.0)
		{
			gamma_cr_up = gamma_cr_middle;
		}
		else
		{
			gamma_cr_down = gamma_cr_middle;
		}
	}
	gamma_cr = (gamma_cr_down + gamma_cr_up) / 2.0;
}

// Dowel Action
void FAM_CS::dowel_action_0(double &gamma, double &Es)
{
	double slope = 0.001*Es; //a low slope, for robustness

	Tau_Dowel_0 = slope*(gamma); // Linear Elastic
	dTau_dgamma_0 = slope;

}

void FAM_CS::dowel_action(double &fpc, double &dY, double &Es, double &lm, double &crack_angle, double &gamma, double &e_cr_normal_true, double &DI_p, double &DI_n)
{
	// if crack_angle != 0: 
	double a = 8.0E-5; // The value of e_cr_normal_transform when (e_cr_normal_true = 0)
	double e_cr_normal_transform = 0.5 * (e_cr_normal_true + pow(e_cr_normal_true * e_cr_normal_true + 4.0 * a * a, 0.5));
	double As = 0.25 * pi * dY * dY;  //Areas of steel 2
	double Is = 1.0 / 64.0 * pi * dY * dY * dY * dY; //Inertia moment of steel 2
	double DI_m = 0.02;
	double DI;
	double s = 75.0;
	if (gamma != 0) {
		DI = (1.0 + s * e_cr_normal_transform * lm / dY) * gamma * lm / (2.0 * dY); //Non-dimensional damage index
	} else
	{
		DI = 0.0;
	}

	TDI = DI;  //+  -> p  ,- -> n
	DI = abs(DI);

	if (gamma == 0.0)
	{
		Tau_Dowel = 0.0;
		dTau_dgamma = pow(Es * Is, 0.25) * lm / As * pow(55.0 * pow(-fpc, 0.85), 0.75);
		dTau_deps = 0.0;
	}
	else if (gamma > 0.0)
	{
		if (DI >= DI_p)
		{
			if (DI <= DI_m)
			{
				Tau_Dowel = pow(Es * Is, 0.25) * lm / As * pow(55.0 * pow(-fpc, 0.85), 0.75) * gamma;
				dTau_dgamma = pow(Es * Is, 0.25) * lm / As * pow(55.0 * pow(-fpc, 0.85), 0.75);
				dTau_deps = 0.0;
			}
			else
			{
				Tau_Dowel = pow(Es * Is, 0.25) * lm / As * pow(55.0 * pow(-fpc, 0.85), 0.75) * gamma  /
					pow(1.0 + 3.0 * pow(DI - DI_m, 0.8), 3.0);
				dTau_dgamma = pow(Es * Is, 0.25) * lm / As * pow(55.0 * pow(-fpc, 0.85), 0.75) *
					(1 + 3 * pow(DI - DI_m, 0.8) - 7.2 * gamma * pow(DI - DI_m, -0.2) * (1.0 + s * e_cr_normal_transform * lm / dY) * lm / 2.0 / dY) /
					pow(1.0 + 3.0 * pow(DI - DI_m, 0.8), 4.0);
				dTau_deps = pow(Es * Is, 0.25) * lm / As * pow(55.0 * pow(-fpc, 0.85), 0.75) * gamma *
					(-7.2) * pow(DI - DI_m, -0.2) * lm / 2.0 / dY * s * gamma * lm / dY *
					(1.0 + e_cr_normal_true / pow(e_cr_normal_true * e_cr_normal_true + 4.0 * a * a, 0.5)) /
					pow(1.0 + 3.0 * pow(DI - DI_m, 0.8), 4.0);
			}
		}
		else
		{
			if (DI_p <= DI_m)
			{
				Tau_Dowel = pow(Es * Is, 0.25) * lm / As * pow(55.0 * pow(-fpc, 0.85), 0.75) * gamma;
				dTau_dgamma = pow(Es * Is, 0.25) * lm / As * pow(55.0 * pow(-fpc, 0.85), 0.75);
				dTau_deps = 0.0;
			}
			else
			{
				Tau_Dowel = pow(Es * Is, 0.25) * lm / As * pow(55.0 * pow(-fpc, 0.85), 0.75) * gamma /
					pow(1.0 + 3.0 * pow(DI_p - DI_m, 0.8), 3.0);
				dTau_dgamma = pow(Es * Is, 0.25) * lm / As * pow(55.0 * pow(-fpc, 0.85), 0.75) /
					pow(1.0 + 3.0 * pow(DI_p - DI_m, 0.8), 3.0);
				dTau_deps = 0.0;
			}
		}
	}
	else
	{
		if (DI >= DI_n)
		{
			if (DI <= DI_m)
			{
				Tau_Dowel = pow(Es * Is, 0.25) * lm / As * pow(55.0 * pow(-fpc, 0.85), 0.75) * gamma;
				dTau_dgamma = pow(Es * Is, 0.25) * lm / As * pow(55.0 * pow(-fpc, 0.85), 0.75);
				dTau_deps = 0.0;
			}
			else
			{
				Tau_Dowel = pow(Es * Is, 0.25) * lm / As * pow(55.0 * pow(-fpc, 0.85), 0.75) * gamma /
					pow(1.0 + 3.0 * pow(DI - DI_m, 0.8), 3.0);
				dTau_dgamma = pow(Es * Is, 0.25) * lm / As * pow(55.0 * pow(-fpc, 0.85), 0.75) *
					(1.0 + 3.0 * pow(DI - DI_m, 0.8) + 7.2 * gamma * pow(DI - DI_m, -0.2) * (1.0 + s * e_cr_normal_transform * lm / dY) * lm / 2.0 / dY) /
					pow(1.0 + 3.0 * pow(DI - DI_m, 0.8), 4.0);
				dTau_deps = pow(Es * Is, 0.25) * lm / As * pow(55.0 * pow(-fpc, 0.85), 0.75) * gamma *
					(-7.2) * pow(DI - DI_m, -0.2) * lm / 2.0 / dY * (- s * gamma * lm / dY) *
					(1.0 + e_cr_normal_true / pow(e_cr_normal_true * e_cr_normal_true + 4.0 * a * a, 0.5)) /
					pow(1.0 + 3.0 * pow(DI - DI_m, 0.8), 4.0);
			}
		}
		else
		{
			if (DI_n <= DI_m)
			{
				Tau_Dowel = pow(Es * Is, 0.25) * lm / As * pow(55.0 * pow(-fpc, 0.85), 0.75) * gamma;
				dTau_dgamma = pow(Es * Is, 0.25) * lm / As * pow(55.0 * pow(-fpc, 0.85), 0.75);
				dTau_deps = 0.0;
			}
			else
			{
				Tau_Dowel = pow(Es * Is, 0.25) * lm / As * pow(55.0 * pow(-fpc, 0.85), 0.75) * gamma /
					pow(1.0 + 3.0 * pow(DI_n - DI_m, 0.8), 3.0);
				dTau_dgamma = pow(Es * Is, 0.25) * lm / As * pow(55.0 * pow(-fpc, 0.85), 0.75) /
					pow(1.0 + 3.0 * pow(DI_n - DI_m, 0.8), 3.0);
				dTau_deps = 0.0;
			}
		}
	}
}

// Damage
void FAM_CS::betaf4(double &eo, double &epc, double &fc, double &epsmax)
{

	// for damage function
	double x_coeff;

	double Kc=(0.27*(-eo/epc-0.37)); // Equation of Vecchio Model B

	double beta_m=1/(1+Kc);

	double delbeta_m=pow(beta_m,2.0)*0.27/epc;

	if ((beta_m>1)||(eo<0)){   
		beta_m=1;
		delbeta_m=0;
	}

	if (fc>0){
		beta_m=1;
		delbeta_m=0;
	}

	//  Beta Damage
	if (eo < epsmax) {
		epsiloncmax = eo;
	}
	else {
		epsiloncmax = epsmax;
	}

	x_coeff=(epsiloncmax/epc);

	if (x_coeff>1) {
		x_coeff=1;
	}

	if (x_coeff < 0) {
		//x_coeff=1;
		x_coeff=0;
		opserr << " Damage Coefficient ErRoR !\n";
		exit(-1); // STOP
	}

	double beta_d=(1-(0.4*x_coeff));

	beta=(beta_m*beta_d);
	delbeta=delbeta_m;
}

// Set Response
Response* FAM_CS::setResponse(const char **argv, int argc, OPS_Stream &theOutput) 
{
	Response *theResponse = 0;

	if (strcmp(argv[0],"panel_strain") == 0 || strcmp(argv[0],"Panel_Strain") == 0) {
		Vector data1(3);
		data1.Zero();

		theResponse = new MaterialResponse(this, 101, data1);
	
	} else if (strcmp(argv[0],"panel_stress") == 0 || strcmp(argv[0],"Panel_Stress") == 0) {

		Vector data2(3);
		data2.Zero();

		theResponse = new MaterialResponse(this, 102, data2);
	
	} else if (strcmp(argv[0],"panel_stress_concrete") == 0 || strcmp(argv[0],"Panel_Stress_Concrete") == 0) {

		Vector data3(3);
		data3.Zero();

		theResponse = new MaterialResponse(this, 103, data3);
	
	} else if (strcmp(argv[0],"panel_stress_steel") == 0 || strcmp(argv[0],"Panel_Stress_Steel") == 0) {

		Vector data4(3);
		data4.Zero();

		theResponse = new MaterialResponse(this, 104, data4);
	
	} else if (strcmp(argv[0],"strain_stress_steelX") == 0 || strcmp(argv[0],"Strain_Stress_SteelX") == 0) {

		Vector data5(2);
		data5.Zero();

		theResponse = new MaterialResponse(this, 105, data5);
	
	} else if (strcmp(argv[0],"strain_stress_steelY") == 0 || strcmp(argv[0],"Strain_Stress_SteelY") == 0) {

		Vector data6(2);
		data6.Zero();

		theResponse = new MaterialResponse(this, 106, data6);
	
	} else if (strcmp(argv[0],"strain_stress_concrete1") == 0 || strcmp(argv[0],"Strain_Stress_Concrete1") == 0) {

		Vector data7(2);
		data7.Zero();

		theResponse = new MaterialResponse(this, 107, data7);
	
	} else if (strcmp(argv[0],"strain_stress_concrete2") == 0 || strcmp(argv[0],"Strain_Stress_Concrete2") == 0) {
		Vector data8(2);
		data8.Zero();

		theResponse = new MaterialResponse(this, 108, data8);

	} else if (strcmp(argv[0],"strain_stress_interlock1") == 0 || strcmp(argv[0],"Strain_Stress_Interlock1") == 0) {

		Vector data9(2);
		data9.Zero();

		theResponse = new MaterialResponse(this, 109, data9);

	} else if (strcmp(argv[0],"strain_stress_interlock2") == 0 || strcmp(argv[0],"Strain_Stress_Interlock2") == 0) {

		Vector data10(2);
		data10.Zero();

		theResponse = new MaterialResponse(this, 110, data10);

	} else if (strcmp(argv[0],"cracking_angles") == 0 || strcmp(argv[0],"Cracking_Angles") == 0) {

		Vector data11(2);
		data11.Zero();

		theResponse = new MaterialResponse(this, 111, data11);
	
	} else if (strcmp(argv[0], "getInputParameters") == 0) {

		Vector data12(12);
		data12.Zero();
		theResponse = new MaterialResponse(this, 112, data12);

	} else if (strcmp(argv[0], "NewOrigin") == 0) {

		double data13 = 0.0;		
		theResponse = new MaterialResponse(this, 113, data13);

	} 

	else

		return this->NDMaterial::setResponse(argv, argc, theOutput);

	return theResponse;
}

// Get Response
int FAM_CS::getResponse(int responseID, Information &matInfo)
{
	if (responseID == 101) {
		return matInfo.setVector(this->getCommittedStrain());

	} else if (responseID == 102){
		return matInfo.setVector(this->getCommittedStress());

	} else if (responseID == 103){
		return matInfo.setVector(this->getPanelStressConcrete()); 

	} else if (responseID == 104){
		return matInfo.setVector(this->getPanelStressSteel()); 

	} else if (responseID == 105){
		return matInfo.setVector(this->getStrainStressSteel1()); 

	} else if (responseID == 106){
		return matInfo.setVector(this->getStrainStressSteel2()); 

	} else if (responseID == 107){
		return matInfo.setVector(this->getStrainStressConcrete1());

	} else if (responseID == 108){
		return matInfo.setVector(this->getStrainStressConcrete2()); 

	} else if (responseID == 109){
		return matInfo.setVector(this->getStrainStressInterlock1()); 

	} else if (responseID == 110){
		return matInfo.setVector(this->getStrainStressInterlock2()); 

	} else if (responseID == 111){
		return matInfo.setVector(this->getCrackingAngles()); 

	} else if (responseID == 112) {
		return matInfo.setVector(this->getInputParameters());

	} else if (responseID == 113) {
		theResponses[2]->getResponse();
		return matInfo.setDouble(theResponses[2]->getInformation().theDouble);
	}
	else {

	return 0;

	}
}

// Function that returns input parameters
Vector FAM_CS::getInputParameters(void)
{
	Vector input_par(12); // size = max number of parameters (assigned + default)

	input_par.Zero();

	input_par(0) = this->getTag();
	input_par(1) = rho;
	input_par(2) = fpc;
	input_par(3) = roux;
	input_par(4) = rouy;
	input_par(5) = dY;
	input_par(6) = Gamax;
	input_par(7) = lm0;
	input_par(8) = sh;
	input_par(9) = Ec;

	return input_par;
}
