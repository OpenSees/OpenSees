// Code written/implemented by:	Kristijan Kolozvari (kkolozvari@fullerton.edu)
//								California State University, Fullerton 
//								Kutay Orakcal
//								Bogazici University, Istanbul, Turkey
//								John Wallace
//								University of California, Los Angeles
//
// Created: 07/2015
//
// Description: This file contains the class definition for Fixed-Strut Angle
// Model (FSAM, Ulugtekin, 2010; Orakcal et al., 2012) which is a plane-stress 
// constitutive model for simulating the behavior of RC panel elements under 
// generalized, in-plane, reversed-cyclic loading conditions. The model assumes
// perfect bond assumption between concrete and reinforcing steel bars. The 
// reinforcing steel bars develop uniaxial stresses under strains in their 
// longitudinal direction, the behavior of concrete is defined using stress–strain 
// relationships in biaxial directions, the orientation of which is governed by 
// the state of cracking in concrete, and also incorporates biaxial softening 
// effects including compression softening and biaxial damage. For transfer of 
// shear stresses across the cracks, a friction-based elasto-plastic shear aggregate
// interlock model is adopted, together with a linear elastic model for representing
// dowel action on the reinforcing steel bars (Kolozvari, 2013).
//
// References:
// 1) Orakcal, K., Massone L.M., Ulugtekin, D.,“Constitutive Modeling of Reinforced Concrete 
// Panel Behavior under Cyclic Loading”, Proceedings of the 15th World Conference on 
// Earthquake Engineering, Lisbon, Portugal, 2012.
// 2) Ulugtekin, D., “Analytical Modeling of Reinforced Concrete Panel Elements under 
// Reversed Cyclic Loadings”, M.S. Thesis, Bogazici University, Istanbul, Turkey, 2010.
// 3) Kolozvari K. (2013). “Analytical Modeling of Cyclic Shear-Flexure Interaction in 
// Reinforced Concrete Structural Walls”, PhD Dissertation, University of California, Los Angeles.
//
// Source: /usr/local/cvs/OpenSees/SRC/material/nD/reinforcedConcretePlaneStress/FSAM.h
//
// Rev: 1


#include "FSAM.h"
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


static int numFSAMMaterials = 0;

// Read input parameters and build the material
OPS_Export void *OPS_FSAMMaterial()
{
	if (numFSAMMaterials == 0) {
		numFSAMMaterials++;
	}

	// Pointer to a uniaxial material that will be returned
	NDMaterial *theMaterial = 0;

	int numRemainingArgs = OPS_GetNumRemainingInputArgs();

	// Parse the script for material parameters
	if (numRemainingArgs != 9) { // total # of input parameters
		opserr << "Invalid #Args want: NDMaterial FSAM $mattag $Rho $Tag_UniaxialSteelX $Tag_UniaxialSteelY $Tag_UniaxialConcrete $rouX $rouY $nu $alfadow\n";
		return 0;	
	}

	int tag;				// nDMaterial tag
	double rho;				// nDMaterial density
	int    iData[3];		// # of uniaxial materials
	double dData[4];		// # of material parameters
	int numData = 0;

	// nDMaterial tag
	numData = 1;
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "WARNING invalid uniaxialMaterial FSAM tag" << endln;
		return 0;
	}

	// nDMaterial density
	numData = 1;
	if (OPS_GetDouble(&numData, &rho) != 0) {
		opserr << "Invalid Arg rho: nDMaterial FSAM $mattag $rho $sX $sY $conc $rouX $rouY $nu $alfadow" << endln;
		return 0;	
	}

	// Material tags of 2 steel (in X and Y directions) and 1 concrete materials
	numData = 3;
	if (OPS_GetInt(&numData, iData) != 0) {
		opserr << "WARNING invalid uniaxialMaterial FSAM tag" << endln;
		return 0;
	}

	// Other FSAM material parameters
	numData = 4;
	if (OPS_GetDouble(&numData, dData) != 0) {
		opserr << "WARNING invalid uniaxialMaterial FSAM tag" << endln;
		return 0;
	}

	// Get pointers to Uniaxial materials
	// Steel X
	UniaxialMaterial *theUniaxialMaterial1 = OPS_GetUniaxialMaterial(iData[0]);
	if (theUniaxialMaterial1 == 0) {
		opserr << "WARNING material not found\n";
		opserr << "Material: " << iData[0];
		opserr << "\nFSAM: " << tag << endln;
		return 0;
	}

	// Steel Y
	UniaxialMaterial *theUniaxialMaterial2 = OPS_GetUniaxialMaterial(iData[1]);
	if (theUniaxialMaterial2 == 0) {
		opserr << "WARNING material not found\n";
		opserr << "Material: " << iData[1];
		opserr << "\nFSAM: " << tag << endln;
		return 0;
	}

	// Concrete c1.1 - uncracked (dummy)
	UniaxialMaterial *theUniaxialMaterial3 = OPS_GetUniaxialMaterial(iData[2]);
	if (theUniaxialMaterial3 == 0) {
		opserr << "WARNING material not found\n";
		opserr << "Material: " << iData[2];
		opserr << "\nFSAM: " << tag << endln;
		return 0;
	}

	// Concrete c1.2 - uncracked (dummy)
	UniaxialMaterial *theUniaxialMaterial4 = OPS_GetUniaxialMaterial(iData[2]);  
	if (theUniaxialMaterial4 == 0) {
		opserr << "WARNING material not found\n";
		opserr << "Material: " << iData[2];
		opserr << "\nFSAM: " << tag << endln;
		return 0;
	}

	// Concrete cA1 - 1st, 2nd crack
	UniaxialMaterial *theUniaxialMaterial5 = OPS_GetUniaxialMaterial(iData[2]);
	if (theUniaxialMaterial5 == 0) {
		opserr << "WARNING material not found\n";
		opserr << "Material: " << iData[2];
		opserr << "\nFSAM: " << tag << endln;
		return 0;
	}

	// Concrete cA2 - 1st, 2nd crack
	UniaxialMaterial *theUniaxialMaterial6 = OPS_GetUniaxialMaterial(iData[2]);  
	if (theUniaxialMaterial6 == 0) {
		opserr << "WARNING material not found\n";
		opserr << "Material: " << iData[2];
		opserr << "\nFSAM: " << tag << endln;
		return 0;
	}

	// Concrete cB1 - 2nd crack
	UniaxialMaterial *theUniaxialMaterial7 = OPS_GetUniaxialMaterial(iData[2]);
	if (theUniaxialMaterial7 == 0) {
		opserr << "WARNING material not found\n";
		opserr << "Material: " << iData[2];
		opserr << "\nFSAM: " << tag << endln;
		return 0;
	}

	// Concrete cB2 - 2nd crack
	UniaxialMaterial *theUniaxialMaterial8 = OPS_GetUniaxialMaterial(iData[2]);  
	if (theUniaxialMaterial8 == 0) {
		opserr << "WARNING material not found\n";
		opserr << "Material: " << iData[2];
		opserr << "\nFSAM: " << tag << endln;
		return 0;
	}

	//Create the FSAM
	theMaterial = new FSAM (tag, rho, 
		theUniaxialMaterial1, theUniaxialMaterial2, 
		theUniaxialMaterial3, theUniaxialMaterial4, 
		theUniaxialMaterial5, theUniaxialMaterial6,
		theUniaxialMaterial7, theUniaxialMaterial8,
		dData[0], dData[1], dData[2], dData[3]);

	if (theMaterial == 0) {
		opserr << "WARNING ran out of memory creating material\n";
		opserr << "FSAM: " << tag << endln;
		return 0;
	}

	return theMaterial;
}

// Typical Constructor
FSAM::FSAM (int tag, 
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
	double NU,				// Friction coefficient of shear aggregate interlock
	double ALFADOW)			// Stiffness parameter of dowel action

	: NDMaterial(tag, ND_TAG_FSAM), 
	rho(RHO), roux(ROUX), rouy(ROUY),	
	nu(NU), alfadow(ALFADOW),					
	strain_vec(3), stress_vec(3), tangent_matrix(3,3), 
	CStress(3), CStrain(3), 
	CPanelConcStress(3), CPanelSteelStress(3), TPanelConcStress(3), TPanelSteelStress(3), 
	CStrainStressSteel1(2), CStrainStressSteel2(2), TStrainStressSteel1(2), TStrainStressSteel2(2),
	CStrainStressConc1(2), CStrainStressConc2(2), TStrainStressConc1(2), TStrainStressConc2(2), 
	CStrainStressInterlock1(2), CStrainStressInterlock2(2), TStrainStressInterlock1(2), TStrainStressInterlock2(2),
	CCrackingAngles(2), pi(3.1415926535)
{

	TeTaSt = 0.0; // Direction of horizontal reinforcement (fixed for now)

	// Material parameters
	E0x = 0.0;
	E0y = 0.0;

	Ec = 0.0;
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

	// State Vairables
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
	// "Transfer" variables
	Tau_Interlock = 0.0;
	dTau_de12 = 0.0;
	dTau_dfcnormal = 0.0;

	// Crack Slip
	TeA12 = 0.0;
	CeA12 = 0.0;
	TeB12 = 0.0;
	CeB12 = 0.0;

	// Shear Stress
	Ctau_Interlock_A = 0.0;
	Ttau_Interlock_A = 0.0;
	Ctau_Interlock_B = 0.0;
	Ttau_Interlock_B = 0.0;

	// History variables for 2nd cracking criterium (cyclic cracking strain)
	TepsA2 = 0.0;
	CepsA2 = 0.0;

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
		opserr << " FSAM::FSAM - failed allocate material array\n";
		exit(-1);
	}

	// Get the copy for SteelX
	theMaterial[0] = s1->getCopy();
	// Check allocation    
	if ( theMaterial[0] == 0 ) {
		opserr << " FSAM::FSAM - failed to get a copy for Steel1\n";
		exit(-1);
	}

	// Get the copy for SteelY
	theMaterial[1] = s2->getCopy();	
	// Check allocation    
	if ( theMaterial[1] == 0 ) {
		opserr << " FSAM::FSAM - failed to get a copy for Steel2\n";
		exit(-1);
	}

	// Get the copy for Concrete A1
	theMaterial[4] = cA1->getCopy();	
	// Check allocation    
	if ( theMaterial[4] == 0 ) {
		opserr << " FSAM::FSAM - failed to get a copy for Concrete A1\n";
		exit(-1);
	}

	// Get the copy for Concrete A2
	theMaterial[5] = cA2->getCopy();	
	// Check allocation    
	if ( theMaterial[5] == 0 ) {
		opserr << " FSAM::FSAM - failed to get a copy for Concrete A2\n";
		exit(-1);
	}

	// Get the copy for Concrete B1
	theMaterial[6] = cB1->getCopy();	
	// Check allocation    
	if ( theMaterial[6] == 0 ) {
		opserr << " FSAM::FSAM - failed to get a copy for Concrete B1\n";
		exit(-1);
	}

	// Get the copy for Concrete B2
	theMaterial[7] = cB2->getCopy();	
	// Check allocation    
	if ( theMaterial[7] == 0 ) {
		opserr << " FSAM::FSAM - failed to get a copy for Concrete B2\n";
		exit(-1);
	}

	// get/set responses
	theResponses = new Response *[2];  
	if ( theResponses == 0) {
		opserr << " FSAM::FSAM - failed allocate responses array\n";
		exit(-1);
	}

	OPS_Stream *theDummyStream = new DummyStream();
	const char **argv = new const char *[1];
	argv[0] = "getCommittedCyclicCrackingConcreteStrain"; // to get committed concrete cyclic cracking strain from strut A2
	theResponses[0] = theMaterial[5]->setResponse(argv, 1, *theDummyStream);

	if (theResponses[0] == 0) {
			opserr << " FSAM::FSAM - failed to set appropriate materials tag: " << tag << "\n";
			exit(-1);
	}

	argv[0] = "getInputParameters"; // to get input parameters from ConcreteCM
	theResponses[1] = theMaterial[4]->setResponse(argv, 1, *theDummyStream);

	if (theResponses[1] == 0) {
			opserr << " FSAM::FSAM - failed to set appropriate materials tag: " << tag << "\n";
			exit(-1);
	}

	delete theDummyStream;

	// Get ConcreteCM material input variables
	theResponses[1]->getResponse();
	Information &theInfoInput = theResponses[1]->getInformation();
	const Vector InputConc = theInfoInput.getData();

	for (int i=0; i<InputConc.Size() ; i++)
	ConcreteInput[i] = InputConc[i];

	// Now create monotonic concrete materials for uncracked stage of behavior
	// Concrete 1.1
	// Instead of: theMaterial[2] = c1->getCopy(); we are creating monotonic ConcreteCM	
	theMaterial[2] = new ConcreteCM(-1111, ConcreteInput[1], ConcreteInput[2], ConcreteInput[3], 
		ConcreteInput[4], ConcreteInput[5], ConcreteInput[6], ConcreteInput[7], ConcreteInput[8], ConcreteInput[9], 1); // create monotonic concrete

	// Check allocation    
	if ( theMaterial[2] == 0 ) {
		opserr << " FSAM::FSAM - failed to get a copy for Concrete 1\n";
		exit(-1);
	}

	// Concrete 1.2
	//Instead of: theMaterial[3] = c2->getCopy();  we are creating monotonic ConcreteCM	
	theMaterial[3] = new ConcreteCM(-2222, ConcreteInput[1], ConcreteInput[2], ConcreteInput[3], 
		ConcreteInput[4], ConcreteInput[5], ConcreteInput[6], ConcreteInput[7], ConcreteInput[8], ConcreteInput[9], 1); // create monotonic concrete

	// Check allocation    
	if ( theMaterial[3] == 0 ) {
		opserr << " FSAM::FSAM - failed to get a copy for Concrete 2\n";
		exit(-1);
	}

	// Obtain some material properties used later in the FSAM model
	// Young's modulus for concrete
	Ec = theMaterial[4] -> getInitialTangent();

	// Strain at peak compressive stress for concrete
	epcc = InputConc[2];

	// Cracking strain for concrete
	et = InputConc[7];

	// Young's modulus for steel
	E0x = theMaterial[0] -> getInitialTangent(); // Horizontal reinforcement
	E0y = theMaterial[1] -> getInitialTangent(); // Vertical reinforcement

	this->revertToStart();
}

// Blank constructor
FSAM::FSAM():NDMaterial(0, ND_TAG_FSAM), 
	strain_vec(3), stress_vec(3), tangent_matrix(3,3),
	pi(3.1415926535)

{
	theMaterial = 0;
	theResponses = 0;

	this->revertToStart();
}

// Destructor
FSAM::~FSAM()
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
		for (int j=0; j<2; j++)
		{
			if (theResponses[j] != 0)
				delete theResponses[j];
		}
		delete [] theResponses;
	}

}

// get copy
NDMaterial* FSAM::getCopy(void) 
{

	FSAM* theCopy =
		new FSAM( this->getTag(), 
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
		nu,
		alfadow);
	
	return theCopy;
}

// get copy
NDMaterial* FSAM::getCopy(const char *type)
{

	FSAM* theModel =
		new FSAM( this->getTag(), 
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
		nu,
		alfadow);

	return theModel;
}

// Print 
void FSAM::Print(OPS_Stream &s, int flag)
{
	s << "\nFSAM, nDMaterial tag: " << this->getTag() << endln;

	// Input values
	s << "density: " << rho << endln;
	s << "roux: " << roux << ", rouy: " << rouy << endln;
	s << "nu: " << nu  << ", alphadow: " << alfadow << endln;	

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

int FSAM::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	int dataTag = this->getDbTag();

	// Packs its data into a Vector and sends this to theChannel
	static Vector data(6);

	data(0) = this->getTag();
	data(1) = rho;
	data(2) = roux;
	data(3) = rouy;
	data(4) = nu;
	data(5) = alfadow;

	res += theChannel.sendVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING FSAM::sendSelf() - " << this->getTag() << " failed to send Vector\n";
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
		opserr << "WARNING FSAM::sendSelf() - " << this->getTag() << " failed to send ID\n";
		return res;
	}

	// Quad asks its material objects to send themselves
	for (i = 0; i < 8; i++) { // i < # of materials
		res += theMaterial[i]->sendSelf(commitTag, theChannel);
		if (res < 0) {
			opserr << "FSAM::sendSelf() - " << this->getTag() << " failed to send its Material\n";
			return res;
		}
	}	

	return res;
}

int FSAM::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res = 0;

	int dataTag = this->getDbTag();

	// Quad creates a Vector, receives the Vector and then sets the internal data with the data in the Vector
	static Vector data(16);
	res += theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING FSAM::recvSelf() - failed to receive Vector\n";
		return res;
	}

	this->setTag((int)data(0));
	rho		= data(1);
	roux	= data(2);
	rouy	= data(3);
	nu		= data(4);
	alfadow = data(5);

	static ID idData(16); // idData(2 x # of uniaxial materials)

	// Receives the tags of its materials
	res += theChannel.recvID(dataTag, commitTag, idData);
	if (res < 0) {
		opserr << "WARNING FSAM::recvSelf() - " << this->getTag() << " failed to receive ID\n";
		return res;
	}

	if (theMaterial == 0) {
		// Allocate new materials
		theMaterial = new UniaxialMaterial *[8]; // *[# of uniaxial materials]
		if (theMaterial == 0) {
			opserr << "FSAM::recvSelf() - Could not allocate UniaxialMaterial* array\n";
			return -1;
		}
		for (int i = 0; i < 8; i++) { // i < # of uniaxial materials
			int matClassTag = idData(i);
			int matDbTag = idData(i+8); // (i + # of uniaxial materials)
			// Allocate new material with the sent class tag
			theMaterial[i] = theBroker.getNewUniaxialMaterial(matClassTag);
			if (theMaterial[i] == 0) {
				opserr << "FSAM::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
				return -1;
			}
			// Now receive materials into the newly allocated space
			theMaterial[i]->setDbTag(matDbTag);
			res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
			if (res < 0) {
				opserr << "FSAM::recvSelf() - material " << i << "failed to recv itself\n";
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
					opserr << "FSAM::recvSelf() - material " << i << "failed to create\n";
					return -1;
				}
			}
			// Receive the material
			theMaterial[i]->setDbTag(matDbTag);
			res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
			if (res < 0) {
				opserr << "FSAM::recvSelf() - material " << i << "failed to recv itself\n";
				return res;
			}
		}
	}

	return res;
}

// Get density
double FSAM::getRho(void)
{
	return rho;
}

// Load strain from the element
int FSAM::setTrialStrain(const Vector &v)
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
int FSAM::setTrialStrain(const Vector &v, const Vector &r)
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
int FSAM::setTrialStrainIncr(const Vector &v)
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
int FSAM::setTrialStrainIncr(const Vector &v, const Vector &r)
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
const Matrix& FSAM::getTangent (void)
{
	return tangent_matrix;
}

// Get trial stress vector
const Vector& FSAM::getStress(void)
{
	return stress_vec;
}

// Get strain vector
const Vector& FSAM::getStrain()
{
	return strain_vec;
}

// Get commited stress
const Vector& FSAM::getCommittedStress(void)
{
	return CStress;
}

// Get commited strain
const Vector& FSAM::getCommittedStrain(void)
{
	return CStrain;
}

// Functions below return values for recorders
// Concrete stresses
Vector FSAM::getPanelStressConcrete(void) 
{
	return CPanelConcStress;
}

// Steel stresses (multiplied with reinforcing ratio)
Vector FSAM::getPanelStressSteel(void) 
{
	return CPanelSteelStress;
}

// Steel stresses in horizontal direction (from uniaxial material)
Vector FSAM::getStrainStressSteel1(void) 
{
	return CStrainStressSteel1;
}

// Steel stresses in vertical direction (from uniaxial material)
Vector FSAM::getStrainStressSteel2(void) 
{
	return CStrainStressSteel2;
}

// Concrete stresses along 1st strut
Vector FSAM::getStrainStressConcrete1(void) 
{
	return CStrainStressConc1;
}

// Concrete stresses along 2nd strut
Vector FSAM::getStrainStressConcrete2(void) 
{
	return CStrainStressConc2;
}

// Aggregate interlock stresses along 1st strut
Vector FSAM::getStrainStressInterlock1(void) 
{
	return CStrainStressInterlock1;
}

// Aggregate interlock stresses along 2nd strut
Vector FSAM::getStrainStressInterlock2(void) 
{
	return CStrainStressInterlock2; 
}

// Cracking angles
Vector FSAM::getCrackingAngles(void) 
{
	return CCrackingAngles;
}

// Revert to start
int FSAM::revertToStart(void)
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
const Matrix& FSAM::getInitialTangent (void)
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
int FSAM::commitState(void)
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

		// Store committed values of principal strains
		Cprstrain1 = Tprstrain1;
		Cprstrain2 = Tprstrain2;

		// Store committed values of shear aggregate interlock
		CeA12 = TeA12;
		Ctau_Interlock_A = Ttau_Interlock_A;

		// Store committed values compressive strain in perpendicular direction (for 2nd cracking criterium)
		CepsA2 = TepsA2;

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

		// Store committed values of principal strains
		Cprstrain1 = Tprstrain1;
		Cprstrain2 = Tprstrain2;

		// Store committed values of shear aggregate interlock
		CeA12 = TeA12;
		Ctau_Interlock_A = Ttau_Interlock_A;

		// Store committed values compressive strain in perpendicular direction (for 2nd cracking criterium)
		CepsA2 = TepsA2;

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

			// Store committed values of principal strains
			Cprstrain1 = Tprstrain1;
			Cprstrain2 = Tprstrain2;

			// Store committed values of shear aggregate interlock
			CeA12 = TeA12;
			CeB12 = TeB12;
			Ctau_Interlock_A = Ttau_Interlock_A;
			Ctau_Interlock_B = Ttau_Interlock_B;

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

		// Store committed values of principal strains
		Cprstrain1 = Tprstrain1;
		Cprstrain2 = Tprstrain2;

		// Store committed values of shear aggregate interlock
		CeA12 = TeA12;
		CeB12 = TeB12;
		Ctau_Interlock_A = Ttau_Interlock_A;
		Ctau_Interlock_B = Ttau_Interlock_B;

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
int FSAM::revertToLastCommit(void)
{
	for (int i=0; i < 8; i++) {
		theMaterial[i]->revertToLastCommit();
	}

	// Shear aggregate interlock
	TeA12 = CeA12;
	TeB12 = CeB12;
	Ttau_Interlock_A = Ctau_Interlock_A;
	Ttau_Interlock_B = Ctau_Interlock_B;
	
	// for 2nd cracking criteria
	TepsA2 = CepsA2;

	return 0;
}

// Calculate trial stress and tangent stiffness based on state
int FSAM::determineTrialStressAndTangent(void)
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
void FSAM::Stage1(double &ex, double &ey, double &gamma)
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

	// if 1st interation calculate initial stiffness

	if (first_iteration == 1) {
		tangent_matrix = this->getInitialTangent();
	}

}

// STAGE 2 - 1st CRACK PANEL MODEL 
void FSAM::Stage2(double &ex, double &ey, double &gamma)
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
	// END - Declaration of temporary variables .......................
	
	double alfacrA = alfa_crackA; // alpha_crack = alfacrackold 

	if (gamma == 0.0) {
		gamma = 1e-20;
	}

	double alfa = 0.5*atan((gamma)/(ex - ey));

	// Principal strains
	double epslon1 = 0.5*(ex+ey) + 0.5*gamma/sin(2.0*alfa); //  strain along concrete strut One - epslon1
	double epslon2 = 0.5*(ex+ey) - 0.5*gamma/sin(2.0*alfa); //  strain along concrete strut Two - epslon2

	// Storing principanl strains and direction
	alpha_strain = alfa;	// principal strain direction
	Tprstrain1 = epslon1;	// principal strain 1 - for cracking criteria
	Tprstrain2 = epslon2;	// principal strain 2 - for cracking criteria

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
	estxp = 0.5*(ex+ey) + 0.5*(ex-ey)*cos(2.0*TeTaSt) + 0.5*gamma*sin(2.0*TeTaSt);

	if (TeTaSt >= 0.0) {
		TeTaStper = TeTaSt - 0.5*pi;
	}
	else {
		TeTaStper = TeTaSt + 0.5*pi;
	}

	estyp = 0.5*(ex+ey) + 0.5*(ex-ey)*cos(2.0*TeTaStper) + 0.5*gamma*sin(2.0*TeTaStper);

	// Get Principal Concrete Stresses
	// Concrete Strut A1
	theMaterial[4]->setTrialStrain( ec1 );
	fc1 = theMaterial[4]->getStress();
	Ect1 = theMaterial[4]->getTangent();

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

	// Shear Aggregate Interlock
	InterLocker_improved(e1, fcmod1, TeA12, CeA12, epcc, Ec, Ctau_Interlock_A);

	Ttau_Interlock_A = Tau_Interlock;
	double dtAggde12 = dTau_de12;
	double dtAggdfc1 = dTau_dfcnormal;

	// Interlock Stress Strain One
	TStrainStressInterlock1(0) = TeA12;
	TStrainStressInterlock1(1) = Ttau_Interlock_A;

	// Interlock Stress Strain Two
	TStrainStressInterlock2(0) = 0.0;
	TStrainStressInterlock2(1) = 0.0;

	// Backtransform concrete stresses
	fcx = ((fcmod1+fcmod2)/2.0) + (((fcmod1-fcmod2)/2.0)*cos(2.0*alfacrA))-(Ttau_Interlock_A*(sin(2.0*alfacrA)));
	fcy = ((fcmod1+fcmod2)/2.0) - (((fcmod1-fcmod2)/2.0)*cos(2.0*alfacrA))+(Ttau_Interlock_A*(sin(2.0*alfacrA)));
	taucxy = (((fcmod1-fcmod2)/2.0)*sin(2.0*alfacrA)) + (((Ttau_Interlock_A*cos(2.0*alfacrA))));

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

	// Dowel Action - only on the horizontal plane (for now)
	// Strain Transformation
	double gamma_dow_A_x = 0;
	double gamma_dow_A_y = -0.5*(ex-ey)*sin(2.0*TeTaSt) + gamma*cos(2.0*TeTaSt);

	// Obtain dowel stresses and stiffness
	dowel_action(gamma_dow_A_x, E0x);
	double Tau_Dowel_A_x = Tau_Dowel;
	double dTau_dgamma_A_x = dTau_dgamma;

	dowel_action(gamma_dow_A_y, E0y);
	double Tau_Dowel_A_y = Tau_Dowel;
	double dTau_dgamma_A_y = dTau_dgamma;

	// Back Transfor Steel Stresses in X - Y coordinate system
	taustxy=(((roux*fstxp-rouy*fstyp)/2.0)*sin(2.0*TeTaSt)) + (rouy*Tau_Dowel_A_y*cos(2.0*TeTaSt)); // + dowel stress 
	fstx=(((roux*fstxp+rouy*fstyp)/2.0)+(((roux*fstxp-rouy*fstyp)/2.0)*cos(2.0*TeTaSt)))-(rouy*Tau_Dowel_A_y*sin(2.0*TeTaSt)); // + dowel stress (present only if TeTaSt><0);
	fsty=(((roux*fstxp+rouy*fstyp)/2.0)-(((roux*fstxp-rouy*fstyp)/2.0)*cos(2.0*TeTaSt)))+(rouy*Tau_Dowel_A_y*sin(2.0*TeTaSt)); // + dowel stress (present only if TeTaSt><0);

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

	if (ex == ey) {
		R = 1.0;
	}

	else {
		R = 1.0 + pow(gamma,2.0)/pow((ex - ey),2.0);
	}

	// Calculating partial stiffnesses
	dsxdex = (4.0*dtAggde12 + 3.0*stifcu11F + stifcu12F + stifcu21F + 3.0*stifcu22F + 4.0*(stifcu11F - stifcu22F)*cos(2.0*alfacrA) +
		(-4.0*dtAggde12 + stifcu11F - stifcu12F - stifcu21F + stifcu22F)*cos(4.0*alfacrA) - 
		4.0*dtAggdfc1*(stifcu11F + stifcu12F + (stifcu11F - stifcu12F)*cos(2.0*alfacrA))*sin(2.0*alfacrA))/8.0 + 
		(4.0*(stifstxpF*roux - stifstypF*rouy)*cos(2.0*TeTaSt) + (stifstxpF*roux + stifstypF*rouy)*(3.0 + cos(4.0*TeTaSt)))/8.0;

	dsxdey =  stifcu12F * pow(cos(alfacrA),4.0) - 2.0*dtAggdfc1*stifcu12F*pow(cos(alfacrA),3.0)*sin(alfacrA) + (-4.0*dtAggde12 + stifcu11F + stifcu22F)*pow(cos(alfacrA),2.0)*pow(sin(alfacrA),2.0) - 
		2.0*dtAggdfc1*stifcu11F*cos(alfacrA)*pow(sin(alfacrA),3.0) + stifcu21F*pow(sin(alfacrA),4.0)+ 
		(stifstxpF*roux + stifstypF*rouy)*pow(cos(TeTaSt),2.0)*pow(sin(TeTaSt),2.0);

	dsxdgamma = -(sin(2.0*alfacrA)*(-stifcu11F + stifcu12F - stifcu21F + stifcu22F + (4.0*dtAggde12 - stifcu11F + stifcu12F + stifcu21F - stifcu22F)*cos(2.0*alfacrA) + 
		2.0*dtAggdfc1*(stifcu11F - stifcu12F)*sin(2.0*alfacrA)))/4.0 + 
		((stifstxpF*roux - stifstypF*rouy + (stifstxpF*roux + stifstypF*rouy)*cos(2.0*TeTaSt))*sin(2.0*TeTaSt))/4.0;

	dsydex = (4.0*stifcu12F*pow(sin(alfacrA),3.0)*(2.0*dtAggdfc1*cos(alfacrA) + sin(alfacrA)) - 4.0*dtAggde12*pow(sin(2.0*alfacrA),2.0) +
		4.0*pow(cos(alfacrA),2.0)*(stifcu21F*pow(cos(alfacrA),2.0) + (stifcu11F + stifcu22F)*pow(sin(alfacrA),2.0) + dtAggdfc1*stifcu11F*sin(2.0*alfacrA)))/4.0 + 
		(stifstxpF*roux + stifstypF*rouy)*pow(cos(TeTaSt),2.0)*pow(sin(TeTaSt),2.0);

	dsydey = (4.0*dtAggde12 + 3.0*stifcu11F + stifcu12F + stifcu21F + 3.0*stifcu22F + 4.0*(-stifcu11F + stifcu22F)*cos(2.0*alfacrA) + 
		(-4.0*dtAggde12 + stifcu11F - stifcu12F - stifcu21F + stifcu22F)*cos(4.0*alfacrA) + 4.0*dtAggdfc1*(stifcu11F + stifcu12F)*sin(2.0*alfacrA) + 
		2.0*dtAggdfc1*(-stifcu11F + stifcu12F)*sin(4.0*alfacrA))/8.0 + 
		((-4.0*stifstxpF*roux + 4.0*stifstypF*rouy)*cos(2.0*TeTaSt) + (stifstxpF*roux + stifstypF*rouy)*(3.0 + cos(4.0*TeTaSt)))/8.0;

	dsydgamma = (sin(2.0*alfacrA)*(stifcu11F - stifcu12F + stifcu21F - stifcu22F + (4.0*dtAggde12 - stifcu11F + stifcu12F + stifcu21F - stifcu22F)*cos(2.0*alfacrA) + 
		2.0*dtAggdfc1*(stifcu11F - stifcu12F)*sin(2.0*alfacrA)))/4.0 + 
		-((-(stifstxpF*roux) + stifstypF*rouy + (stifstxpF*roux + stifstypF*rouy)*cos(2.0*TeTaSt))*sin(2.0*TeTaSt))/4.0;

	dtxydex =  dtAggdfc1*stifcu11F*pow(cos(alfacrA),2.0)*cos(2.0*alfacrA) + (stifcu11F - stifcu21F)*pow(cos(alfacrA),3.0)*sin(alfacrA) + dtAggdfc1*stifcu12F*cos(2.0*alfacrA)*pow(sin(alfacrA),2.0) + 
		(stifcu12F - stifcu22F)*cos(alfacrA)*pow(sin(alfacrA),3.0) - (dtAggde12*sin(4.0*alfacrA))/2.0+ 
		((stifstxpF*roux - stifstypF*rouy + (stifstxpF*roux + stifstypF*rouy)*cos(2.0*TeTaSt))*sin(2.0*TeTaSt))/4.0;

	dtxydey =  dtAggdfc1*stifcu12F*pow(cos(alfacrA),2.0)*cos(2.0*alfacrA) + (stifcu12F - stifcu22F)*pow(cos(alfacrA),3.0)*sin(alfacrA) + dtAggdfc1*stifcu11F*cos(2.0*alfacrA)*pow(sin(alfacrA),2.0) + 
		(stifcu11F - stifcu21F)*cos(alfacrA)*pow(sin(alfacrA),3.0) + (dtAggde12*sin(4.0*alfacrA))/2.0+ 
		-((-(stifstxpF*roux) + stifstypF*rouy + (stifstxpF*roux + stifstypF*rouy)*cos(2.0*TeTaSt))*sin(2.0*TeTaSt))/4.0;

	dtxydgamma = dtAggde12*pow(cos(2*alfacrA),2.0) + (stifcu11F - stifcu12F - stifcu21F + stifcu22F)*pow(cos(alfacrA),2.0)*pow(sin(alfacrA),2.0) + (dtAggdfc1*(stifcu11F - stifcu12F)*sin(4.0*alfacrA))/4.0 + 
		(stifstxpF*roux + stifstypF*rouy)*pow(cos(TeTaSt),2.0)*pow(sin(TeTaSt),2.0)+
		rouy*dTau_dgamma_A_y; // + dowel stiffness (only for TeTaSt = 0) - if TeTaSt><0 all should be derived

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

// STAGE 3 - 2st CRACK PANEL MODEL
void FSAM::Stage3(double &ex, double &ey, double &gamma)
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
	// END - Declaration of temporary variables ......................

	double alfacrA = alfa_crackA; // alpha_crackA=alfacrackAold
	double alfacrB = alfa_crackB; // alpha_crackB=alfacrackBold

	if (gamma == 0.0) {
		gamma = 1e-20;
	}

	// Principal strain direction
	double alfa = 0.5*atan((gamma)/(ex - ey)); 

	// Principal Strains
	double epslon1 = 0.5*(ex+ey) + 0.5*gamma/sin(2.0*alfa);
	double epslon2 = 0.5*(ex+ey) - 0.5*gamma/sin(2.0*alfa);

	// Storing principanl strains and direction
	alpha_strain = alfa; // principal strain direction
	Tprstrain1 = epslon1; // principal strain 1 - for cracking criteria
	Tprstrain2 = epslon2; // principal strain 2 - for cracking criteria

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
	estxp = 0.5*(ex+ey) + 0.5*(ex-ey)*cos(2.0*TeTaSt) + 0.5*gamma*sin(2.0*TeTaSt);

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

	// Shear Aggregate Interlock - Strut A
	InterLocker_improved(eA1, fcmodA1, TeA12, CeA12, epcc, Ec, Ctau_Interlock_A);

	Ttau_Interlock_A = Tau_Interlock;
	double dtAggde12A = dTau_de12;
	double dtAggdfc1A = dTau_dfcnormal;

	// Interlock Stress Strain - Strut A
	TStrainStressInterlock1(0) = TeA12;
	TStrainStressInterlock1(1) = Ttau_Interlock_A;

	// Backtransform concrete stresses - Strut A
	fcxA = ((0.0+fcmodA2)/2.0) + (((0.0-fcmodA2)/2.0)*cos(2.0*alfacrA))-(Ttau_Interlock_A*(sin(2.0*alfacrA)));
	fcyA = ((0.0+fcmodA2)/2.0) - (((0.0-fcmodA2)/2.0)*cos(2.0*alfacrA))+(Ttau_Interlock_A*(sin(2.0*alfacrA)));
	taucxyA = (((0.0-fcmodA2)/2.0)*sin(2.0*alfacrA)) + (((Ttau_Interlock_A*cos(2.0*alfacrA))));

	stifcu11FA = 0.0;
	stifcu12FA = 0.0;
	stifcu21FA = stifcA21;
	stifcu22FA = stifcA22;
	Fcu1FA = 0.0;
	Fcu2FA = fcmodA2;

	// Shear Aggregate Interlock - Strut B
	InterLocker_improved(eB1, fcmodB1, TeB12, CeB12, epcc, Ec, Ctau_Interlock_B);

	Ttau_Interlock_B = Tau_Interlock;
	double dtAggde12B = dTau_de12;
	double dtAggdfc1B = dTau_dfcnormal;

	// Interlock Stress Strain - Strut B
	TStrainStressInterlock2(0) = TeB12;
	TStrainStressInterlock2(1) = Ttau_Interlock_B;

	// Backtransform concrete stresses  - Strut B
	fcxB = ((0.0+fcmodB2)/2.0) + (((0.0-fcmodB2)/2.0)*cos(2.0*alfacrB))-(Ttau_Interlock_B*(sin(2.0*alfacrB)));
	fcyB = ((0.0+fcmodB2)/2.0) - (((0.0-fcmodB2)/2.0)*cos(2.0*alfacrB))+(Ttau_Interlock_B*(sin(2.0*alfacrB)));
	taucxyB = (((0.0-fcmodB2)/2.0)*sin(2.0*alfacrB)) + (((Ttau_Interlock_B*cos(2.0*alfacrB))));

	stifcu11FB = 0.0;
	stifcu12FB = 0.0;
	stifcu21FB = stifcB21;
	stifcu22FB = stifcB22;
	Fcu1FB = 0.0;
	Fcu2FB = fcmodB2;

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
	double gamma_dow_A_y = -0.5*(ex-ey)*sin(2.0*TeTaSt) + gamma*cos(2.0*TeTaSt);

	// Obtain dowel stresses and stiffness
	dowel_action(gamma_dow_A_x, E0x); // =0 for TeTaSt = 0 - if TeTaSt><0 this should be figure out
	double Tau_Dowel_A_x = Tau_Dowel;
	double dTau_dgamma_A_x = dTau_dgamma;

	dowel_action(gamma_dow_A_y, E0y);
	double Tau_Dowel_A_y = Tau_Dowel;
	double dTau_dgamma_A_y = dTau_dgamma;

	// Crack B
	double gamma_dow_B_x = 0;
	double gamma_dow_B_y = -0.5*(ex-ey)*sin(2.0*TeTaSt) + gamma*cos(2.0*TeTaSt);

	// Obtain dowel stresses and stiffness
	dowel_action(gamma_dow_B_x, E0x); // =0 for TeTaSt = 0 - if TeTaSt><0 this should be figure out
	double Tau_Dowel_B_x = Tau_Dowel;
	double dTau_dgamma_B_x = dTau_dgamma;

	dowel_action(gamma_dow_B_y, E0y);
	double Tau_Dowel_B_y = Tau_Dowel;
	double dTau_dgamma_B_y = dTau_dgamma;

	//  Back Transformed Steel Stresses 
	taustxy=(((roux*fstxp-rouy*fstyp)/2.0)*sin(2.0*TeTaSt))+ //Dowell action "taust xp yp" =0
		(rouy*Tau_Dowel_A_y*cos(2.0*TeTaSt))+(rouy*Tau_Dowel_B_y*cos(2.0*TeTaSt)); // + dowel stress; 
	fstx=(((roux*fstxp+rouy*fstyp)/2.0)+(((roux*fstxp-rouy*fstyp)/2.0)*cos(2.0*TeTaSt))); // + dowel stress (present only if TeTaSt><0);
	fsty=(((roux*fstxp+rouy*fstyp)/2.0)-(((roux*fstxp-rouy*fstyp)/2.0)*cos(2.0*TeTaSt))); // + dowel stress (present only if TeTaSt><0);

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
	dsxdexA = (4.0*dtAggde12A + 3.0*stifcu11FA + stifcu12FA + stifcu21FA + 3.0*stifcu22FA + 4.0*(stifcu11FA - stifcu22FA)*cos(2.0*alfacrA) + 
		(-4.0*dtAggde12A + stifcu11FA - stifcu12FA - stifcu21FA + stifcu22FA)*cos(4.0*alfacrA) - 
		4.0*dtAggdfc1A*(stifcu11FA + stifcu12FA + (stifcu11FA - stifcu12FA)*cos(2.0*alfacrA))*sin(2.0*alfacrA))/8.0 + 
		(4.0*(stifstxpF*roux - stifstypF*rouy)*cos(2.0*TeTaSt) + (stifstxpF*roux + stifstypF*rouy)*(3.0 + cos(4.0*TeTaSt)))/8.0;

	dsxdeyA =  stifcu12FA*pow(cos(alfacrA),4.0) - 2.0*dtAggdfc1A*stifcu12FA*pow(cos(alfacrA),3.0)*sin(alfacrA) + (-4.0*dtAggde12A + stifcu11FA + stifcu22FA)*pow(cos(alfacrA),2.0)*pow(sin(alfacrA),2.0) - 
		2.0*dtAggdfc1A*stifcu11FA*cos(alfacrA)*pow(sin(alfacrA),3.0) + stifcu21FA*pow(sin(alfacrA),4.0) + 
		(stifstxpF*roux + stifstypF*rouy)*pow(cos(TeTaSt),2.0)*pow(sin(TeTaSt),2.0);

	dsxdgammaA = -(sin(2.0*alfacrA)*(-stifcu11FA + stifcu12FA - stifcu21FA + stifcu22FA + (4.0*dtAggde12A - stifcu11FA + stifcu12FA + stifcu21FA - stifcu22FA)*cos(2.0*alfacrA) + 
		2.0*dtAggdfc1A*(stifcu11FA - stifcu12FA)*sin(2.0*alfacrA)))/4.0 + 
		((stifstxpF*roux - stifstypF*rouy + (stifstxpF*roux + stifstypF*rouy)*cos(2.0*TeTaSt))*sin(2.0*TeTaSt))/4.0;

	dsydexA = (4.0*stifcu12FA*pow(sin(alfacrA),3.0)*(2.0*dtAggdfc1A*cos(alfacrA) + sin(alfacrA)) - 4.0*dtAggde12A*pow(sin(2.0*alfacrA),2.0) + 
		4.0*pow(cos(alfacrA),2.0)*(stifcu21FA*pow(cos(alfacrA),2.0) + (stifcu11FA + stifcu22FA)*pow(sin(alfacrA),2.0) + dtAggdfc1A*stifcu11FA*sin(2.0*alfacrA)))/4.0 + 
		(stifstxpF*roux + stifstypF*rouy)*pow(cos(TeTaSt),2.0)*pow(sin(TeTaSt),2.0);

	dsydeyA = (4.0*dtAggde12A + 3.0*stifcu11FA + stifcu12FA + stifcu21FA + 3.0*stifcu22FA + 4.0*(-stifcu11FA + stifcu22FA)*cos(2.0*alfacrA) + 
		(-4.0*dtAggde12A + stifcu11FA - stifcu12FA - stifcu21FA + stifcu22FA)*cos(4.0*alfacrA) + 4.0*dtAggdfc1A*(stifcu11FA + stifcu12FA)*sin(2.0*alfacrA) + 
		2.0*dtAggdfc1A*(-stifcu11FA + stifcu12FA)*sin(4.0*alfacrA))/8.0 + 
		((-4.0*stifstxpF*roux + 4.0*stifstypF*rouy)*cos(2.0*TeTaSt) + (stifstxpF*roux + stifstypF*rouy)*(3.0 + cos(4.0*TeTaSt)))/8.0;

	dsydgammaA = (sin(2.0*alfacrA)*(stifcu11FA - stifcu12FA + stifcu21FA - stifcu22FA + (4.0*dtAggde12A - stifcu11FA + stifcu12FA + stifcu21FA - stifcu22FA)*cos(2.0*alfacrA) + 
		2.0*dtAggdfc1A*(stifcu11FA - stifcu12FA)*sin(2.0*alfacrA)))/4.0 + 
		-((-(stifstxpF*roux) + stifstypF*rouy + (stifstxpF*roux + stifstypF*rouy)*cos(2.0*TeTaSt))*sin(2.0*TeTaSt))/4.0;

	dtxydexA =  dtAggdfc1A*stifcu11FA*pow(cos(alfacrA),2.0)*cos(2.0*alfacrA) + (stifcu11FA - stifcu21FA)*pow(cos(alfacrA),3.0)*sin(alfacrA) + dtAggdfc1A*stifcu12FA*cos(2.0*alfacrA)*pow(sin(alfacrA),2.0) + 
		(stifcu12FA - stifcu22FA)*cos(alfacrA)*pow(sin(alfacrA),3.0) - (dtAggde12A*sin(4.0*alfacrA))/2.0 + 
		((stifstxpF*roux - stifstypF*rouy + (stifstxpF*roux + stifstypF*rouy)*cos(2.0*TeTaSt))*sin(2.0*TeTaSt))/4.0;

	dtxydeyA =  dtAggdfc1A*stifcu12FA*pow(cos(alfacrA),2.0)*cos(2.0*alfacrA) + (stifcu12FA - stifcu22FA)*pow(cos(alfacrA),3.0)*sin(alfacrA) + dtAggdfc1A*stifcu11FA*cos(2.0*alfacrA)*pow(sin(alfacrA),2.0) + 
		(stifcu11FA - stifcu21FA)*cos(alfacrA)*pow(sin(alfacrA),3.0) + (dtAggde12A*sin(4.0*alfacrA))/2.0 + 
		-((-(stifstxpF*roux) + stifstypF*rouy + (stifstxpF*roux + stifstypF*rouy)*cos(2.0*TeTaSt))*sin(2.0*TeTaSt))/4.0;

	dtxydgammaA = dtAggde12A*pow(cos(2.0*alfacrB),2.0) + (stifcu11FA - stifcu12FA - stifcu21FA + stifcu22FA)*pow(cos(alfacrA),2.0)*pow(sin(alfacrA),2.0) + (dtAggdfc1A*(stifcu11FA - stifcu12FA)*sin(4.0*alfacrA))/4.0 + 
		(stifstxpF*roux + stifstypF*rouy)*pow(cos(TeTaSt),2.0)*pow(sin(TeTaSt),2.0) +
		rouy*dTau_dgamma_A_y; // + dowel stiffness (only for TeTaSt = 0) - if TeTaSt><0 all should be derived

	// Strut B
	dsxdexB = (4.0*dtAggde12B + 3.0*stifcu11FB + stifcu12FB + stifcu21FB + 3.0*stifcu22FB + 4.0*(stifcu11FB - stifcu22FB)*cos(2.0*alfacrB) + 
		(-4.0*dtAggde12B + stifcu11FB - stifcu12FB - stifcu21FB + stifcu22FB)*cos(4.0*alfacrB) - 
		4.0*dtAggdfc1B*(stifcu11FB + stifcu12FB + (stifcu11FB - stifcu12FB)*cos(2.0*alfacrB))*sin(2.0*alfacrB))/8.0;

	dsxdeyB =  stifcu12FB*pow(cos(alfacrB),4.0) - 2.0*dtAggdfc1B*stifcu12FB*pow(cos(alfacrB),3.0)*sin(alfacrB) + (-4.0*dtAggde12B + stifcu11FB + stifcu22FB)*pow(cos(alfacrB),2.0)*pow(sin(alfacrB),2.0) - 
		2.0*dtAggdfc1B*stifcu11FB*cos(alfacrB)*pow(sin(alfacrB),3.0) + stifcu21FB*pow(sin(alfacrB),4.0);

	dsxdgammaB = -(sin(2.0*alfacrB)*(-stifcu11FB + stifcu12FB - stifcu21FB + stifcu22FB + (4.0*dtAggde12B - stifcu11FB + stifcu12FB + stifcu21FB - stifcu22FB)*cos(2.0*alfacrB) +
		2.0*dtAggdfc1B*(stifcu11FB - stifcu12FB)*sin(2.0*alfacrB)))/4.0;

	dsydexB = (4.0*stifcu12FB*pow(sin(alfacrB),3.0)*(2.0*dtAggdfc1B*cos(alfacrB) + sin(alfacrB)) - 4.0*dtAggde12B*pow(sin(2.0*alfacrB),2.0) + 
		4.0*pow(cos(alfacrB),2.0)*(stifcu21FB*pow(cos(alfacrB),2.0) + (stifcu11FB + stifcu22FB)*pow(sin(alfacrB),2.0) + dtAggdfc1B*stifcu11FB*sin(2.0*alfacrB)))/4.0;

	dsydeyB = (4.0*dtAggde12B + 3.0*stifcu11FB + stifcu12FB + stifcu21FB + 3.0*stifcu22FB + 4.0*(-stifcu11FB + stifcu22FB)*cos(2.0*alfacrB) + 
		(-4.0*dtAggde12B + stifcu11FB - stifcu12FB - stifcu21FB + stifcu22FB)*cos(4.0*alfacrB) + 4.0*dtAggdfc1B*(stifcu11FB + stifcu12FB)*sin(2.0*alfacrB) + 
		2.0*dtAggdfc1B*(-stifcu11FB + stifcu12FB)*sin(4.0*alfacrB))/8.0;

	dsydgammaB = (sin(2.0*alfacrB)*(stifcu11FB - stifcu12FB + stifcu21FB - stifcu22FB + (4.0*dtAggde12B - stifcu11FB + stifcu12FB + stifcu21FB - stifcu22FB)*cos(2.0*alfacrB) + 
		2.0*dtAggdfc1B*(stifcu11FB - stifcu12FB)*sin(2.0*alfacrB)))/4.0;

	dtxydexB =  dtAggdfc1B*stifcu11FB*pow(cos(alfacrB),2.0)*cos(2.0*alfacrB) + (stifcu11FB - stifcu21FB)*pow(cos(alfacrB),3.0)*sin(alfacrB) + dtAggdfc1B*stifcu12FB*cos(2.0*alfacrB)*pow(sin(alfacrB),2.0) + 
		(stifcu12FB - stifcu22FB)*cos(alfacrB)*pow(sin(alfacrB),3.0) - (dtAggde12B*sin(4.0*alfacrB))/2.0;

	dtxydeyB =  dtAggdfc1B*stifcu12FB*pow(cos(alfacrB),2.0)*cos(2.0*alfacrB) + (stifcu12FB - stifcu22FB)*pow(cos(alfacrB),3.0)*sin(alfacrB) + dtAggdfc1B*stifcu11FB*cos(2.0*alfacrB)*pow(sin(alfacrB),2.0) + 
		(stifcu11FB - stifcu21FB)*cos(alfacrB)*pow(sin(alfacrB),3.0) + (dtAggde12B*sin(4.0*alfacrB))/2.0;

	dtxydgammaB = dtAggde12B*pow(cos(2.0*alfacrB),2.0) + (stifcu11FB - stifcu12FB - stifcu21FB + stifcu22FB)*pow(cos(alfacrB),2.0)*pow(sin(alfacrB),2.0) + (dtAggdfc1B*(stifcu11FB - stifcu12FB)*sin(4.0*alfacrB))/4.0 +
		rouy*dTau_dgamma_B_y; // + dowel stiffness (only for TeTaSt = 0) - if TeTaSt><0 all should be derived

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
void FSAM::InterLocker_improved(double &e_cr_normal, double &f_cr_normal, double &e_cr_parallel, double &e_cr_parallel_old,  double &epc, double &Ec, double &Tau_Interlock_old)
{

	double slope = 0.4*Ec; // shear stiffness coefficient G (fixed value)
 
	if (e_cr_parallel == e_cr_parallel_old)
	{
		if (f_cr_normal >= 0.0) // interlock stress is zero when crack under tension
		{
			Tau_Interlock = 0.0;
			dTau_de12 = 0.0;
			dTau_dfcnormal = 0.0;
		}

		else 
		{
			Tau_Interlock = Tau_Interlock_old+slope*(e_cr_parallel-e_cr_parallel_old);
			dTau_de12 = slope;
			dTau_dfcnormal = 0.0;

			if (Tau_Interlock > -(nu*f_cr_normal)) 
			{
				Tau_Interlock = -(nu*f_cr_normal);
				dTau_de12 = 0.0;
				dTau_dfcnormal = -nu;  // original
				//dTau_dfcnormal = nu;
			} 
			else if (Tau_Interlock < (nu*f_cr_normal)) 
			{
				Tau_Interlock = (nu*f_cr_normal);
				dTau_de12 = 0.0;
				dTau_dfcnormal = nu;   //original
				//dTau_dfcnormal = -nu;
			}
		}
	}

	else if (e_cr_parallel>e_cr_parallel_old)

	{
		if (f_cr_normal>=0)
		{
			Tau_Interlock=0;
			dTau_de12=0;
			dTau_dfcnormal=0;
		}
		else
		{
			Tau_Interlock = Tau_Interlock_old+slope*(e_cr_parallel-e_cr_parallel_old);
			dTau_de12 = slope;
			dTau_dfcnormal = 0.0;

			if (Tau_Interlock>-(nu*f_cr_normal))
			{
				Tau_Interlock = -(nu*f_cr_normal);
				dTau_de12 = 0.0;
				dTau_dfcnormal = -nu;  //original
				//dTau_dfcnormal=nu;
			}
		}
	}

	else if (e_cr_parallel < e_cr_parallel_old)
	{
		if (f_cr_normal > 0.0)
		{
			Tau_Interlock = 0.0;
			dTau_de12 = 0.0;
			dTau_dfcnormal = 0.0;
		}

		else

		{
			Tau_Interlock=Tau_Interlock_old+slope*(e_cr_parallel-e_cr_parallel_old);
			dTau_de12=slope;
			dTau_dfcnormal=0;

			if (Tau_Interlock<(nu*f_cr_normal))
			{
				Tau_Interlock=(nu*f_cr_normal);
				dTau_de12=0;
				dTau_dfcnormal=nu;   //original
				//dTau_dfcnormal=-nu;
			}
		}
	}

}

// Dowel Action
void FSAM::dowel_action(double &gama, double &Es)
{
	double slope = alfadow*Es;

	Tau_Dowel = slope*(gama); // Linear Elastic
	dTau_dgamma = slope;

}

// Damage
void FSAM::betaf4(double &eo, double &epc, double &fc, double &epsmax)
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
Response* FSAM::setResponse(const char **argv, int argc, OPS_Stream &theOutput) 
{
	Response *theResponse = 0;

	if (strcmp(argv[0],"panel_strain") == 0 || strcmp(argv[0],"Panel_strain") == 0) {

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
	
	} else

		return this->NDMaterial::setResponse(argv, argc, theOutput);

	return theResponse;
}

// Get Response
int FSAM::getResponse(int responseID, Information &matInfo)
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

	} else {

	return 0;

	}
}
