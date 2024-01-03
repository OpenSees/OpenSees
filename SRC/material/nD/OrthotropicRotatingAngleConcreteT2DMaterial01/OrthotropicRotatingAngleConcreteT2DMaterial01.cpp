// Code written/implemented by: Fabian Rojas B.
//								Maria Jose Nunez
//
// Created: 12/2022
// 
// Description: This file contains the OrthotropicRotatingAngleConcreteT2DMaterial01 class definition.
// An OrthotropicRotatingAngleConcreteT2DMaterial01 is a subclass of the class NDMaterial and corresponds to the abstract representation
// of an Orthotropic Concrete Layer (plane stress) 2D Material with the Rotating Angle and Tangent Formulation for Cycling or Reversed
// Loading with damage that is used in Finite Element Method or Structural Analysis.
//
// Reference:
// 1. Rojas, F., Anderson, J. C., Massones, L. M. (2016). A nonlinear quadrilateral layered membrane with drilling degrees of freedom for 
// the modeling of reinforced concrete walls. Engineering Structures, 124, 521-538.
//
// Source: \OpenSees\SRC\material\nD\OrthotropicRotatingAngleConcreteT2DMaterial01
//
// Rev: 1.0

#include "OrthotropicRotatingAngleConcreteT2DMaterial01.h"
#include <UniaxialMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include <cmath>				
#include <float.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <stdlib.h>
#include <MaterialResponse.h>
#include <DummyStream.h>
#include <elementAPI.h>
#include <algorithm>
#include <iostream>			
#define OPS_Export

#define dbl_Epsilon 2.22e-16	 // a very small double

#include <string.h>
using namespace std;

static int numOrthotropicRotatingAngleConcreteT2DMaterials = 0;

// Read input parameters and build the material
OPS_Export void* OPS_OrthotropicRotatingAngleConcreteT2DMaterial01()
{
	if (numOrthotropicRotatingAngleConcreteT2DMaterials == 0) {
		numOrthotropicRotatingAngleConcreteT2DMaterials++;
	}

	// Pointer to a ND material that will be returned
	NDMaterial* theMaterial = 0;

	int numRemainingArgs = OPS_GetNumRemainingInputArgs();

	// Parse the script for material parameters
	if (numRemainingArgs < 5) {	// total # of input parameters
		opserr << "Invalid #Args want: nDMaterial OrthotropicRotatingAngleConcreteT2DMaterial01 $matTag $Tag_UniaxialConcrete $epscr $epsc $rho <-damageCte1 $DamageCte1> <-damageCte2 $DamageCte2>\n";
		return 0;
	}

	int tag;					// nDMaterial tag
	int UniaxialMatTag;			// uniaxial material tag
	double dData[2];			// # of material parameters
	double rhoNDM;				// nDMaterial density
	double damageCte1 = 0.14;	// default value for damage constant 1
	double damageCte2 = 0.6;	// default value for damage constant 2
	
	// nDMaterial tag
	int numData = 1;
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "WARNING invalid integer tag for nDMaterial OrthotropicRotatingAngleConcreteT2DMaterial" << endln;
		return 0;
	}

	// Material tag of 1 concrete material
	numData = 1;
	if (OPS_GetInt(&numData, &UniaxialMatTag) != 0) {
		opserr << "WARNING invalid uniaxial OrthotropicRotatingAngleConcreteT2DMaterial01 tag" << endln;
		return 0;
	}

	// Other OrthotropicRotatingAngleConcreteT2DMaterial01 material parameters
	numData = 2;
	if (OPS_GetDouble(&numData, dData) != 0) {
		opserr << "WARNING invalid OrthotropicRotatingAngleConcreteT2DMaterial01 material parameters" << endln;
		return 0;
	}

	// nDMaterial density
	numData = 1;
	if (OPS_GetDouble(&numData, &rhoNDM) != 0) {
		opserr << "Invalid arg rho: nDMaterial OrthotropicRotatingAngleConcreteT2DMaterial01 $matTag $Tag_UniaxialConcrete $epscr $epsc $rho <-damageCte1 $DamageCte1> <-damageCte2 $DamageCte2>" << endln;
		return 0;
	}

	// Get pointer to Uniaxial material
	// Concrete 1
	UniaxialMaterial* theUniaxialMaterial1 = OPS_getUniaxialMaterial(UniaxialMatTag);
	if (theUniaxialMaterial1 == 0) {
		opserr << "WARNING material not found\n";
		opserr << "Material: " << UniaxialMatTag;
		opserr << "\nOrthotropicRotatingAngleConcreteT2DMaterial01: " << tag << endln;
		return 0;
	}

	// Concrete 2
	UniaxialMaterial* theUniaxialMaterial2 = OPS_getUniaxialMaterial(UniaxialMatTag);
	if (theUniaxialMaterial2 == 0) {
		opserr << "WARNING material not found\n";
		opserr << "Material: " << UniaxialMatTag;
		opserr << "\nOrthotropicRotatingAngleConcreteT2DMaterial01: " << tag << endln;
		return 0;
	}

	numRemainingArgs -= 5;
	while (numRemainingArgs > 1) {
		const char* str = OPS_GetString();
		if (strcmp(str, "-damageCte1") == 0) {
			numData = 1;
			if (OPS_GetDouble(&numData, &damageCte1) != 0) {
				opserr << "Invalid damageConstant1 for OrthotropicRotatingAngleConcreteT2DMaterial01 " << tag <<
					" $Tag_UniaxialConcrete $epscr $epsc $rho <-damageCte1 $DamageCte1> <-damageCte2 $DamageCte2>\n";
				return 0;
			}
		}
		else if (strcmp(str, "-damageCte2") == 0) {
			numData = 1;
			if (OPS_GetDouble(&numData, &damageCte2) != 0) {
				opserr << "Invalid damageConstant1 for OrthotropicRotatingAngleConcreteT2DMaterial01 " << tag <<
					" $Tag_UniaxialConcrete $epscr $epsc $rho <-damageCte1 $DamageCte1> <-damageCte2 $DamageCte2>\n";
				return 0;
			}
		}
		else {
			opserr << "WARNING: Invalid option " << str << " in OrthotropicRotatingAngleConcreteT2DMaterial01 " << tag << 
				" $Tag_UniaxialConcrete $epscr $epsc $rho <-damageCte1 $DamageCte1> <-damageCte2 $DamageCte2>\n";
			return 0;
		}
		numRemainingArgs -= 2;
	}

	// Create the OrthotropicRotatingAngleConcreteT2DMaterial01
	theMaterial = new OrthotropicRotatingAngleConcreteT2DMaterial01(tag, 
		theUniaxialMaterial1,theUniaxialMaterial2, 
		dData[0], dData[1], 
		rhoNDM, damageCte1, damageCte2);

	if (theMaterial == 0) {
		opserr << "WARNING ran out memory creating material\n";
		opserr << "OrthotropicRotatingAngleConcreteT2DMaterial01: " << tag << endln;
		return 0;
	}

	return theMaterial;
}

// Typical Constructor
OrthotropicRotatingAngleConcreteT2DMaterial01::OrthotropicRotatingAngleConcreteT2DMaterial01(int tag,			// nDMaterial tag
	UniaxialMaterial* concreteUniaxialMatObject1,																// concrete in a principal direction
	UniaxialMaterial* concreteUniaxialMatObject2,																// concrete in the other principal direction 
	double strainAtFcr,																							// strain at tension cracking of the concrete
	double strainAtFc,																							// strain at the compression strength of the concrete
	double rho,																									// density
	double damageCte1,																							// damage constant 1
	double damageCte2)																							// damage constant 2

	: NDMaterial(tag, ND_TAG_OrthotropicRotatingAngleConcreteT2DMaterial01),
	ecr(strainAtFcr), ec(strainAtFc), rhoNDM(rho), damageConstant1(damageCte1), damageConstant2(damageCte2),
	strain_vec(3), stress_vec(3), tangent_matrix(3,3), initialTangentNDM(3,3), strainPrincipalDirection(3), 
	poissonRatios(2), Cstrain(3), CMaxMinStrainRec(3), TMaxMinStrainRec(3), pi(3.1415926535)
{
	thetaPrincipalDirection = 0.0; 
	strainPrincipalDirection(0) = 0.0;
	strainPrincipalDirection(1) = 0.0;
	strainPrincipalDirection(2) = 0.0;
	isConcreteCracked = false;		

	poissonRatios(0) = 0.0;
	poissonRatios(1) = 0.0;

	Cstrain(0) = 0.0;
	Cstrain(1) = 0.0;
	Cstrain(2) = 0.0;

	CMaxMinStrainRec(0) = 1.0;
	TMaxMinStrainRec(0) = 1.0;
	for (int i = 1; i < 3; i++) {
		CMaxMinStrainRec(i) = 0.0;
		TMaxMinStrainRec(i) = 0.0;
	}

	// Allocate pointers for uniaxial materials........................................................
	theMaterial = new UniaxialMaterial * [2];
	if (theMaterial == 0) {
		opserr << " OrthotropicRotatingAngleConcreteT2DMaterial01::OrthotropicRotatingAngleConcreteT2DMaterial01 - failed allocate material array\n";
		exit(-1);
	}

	// We create the copies of the concrete uniaxial material objects that define the behavior at each principal direction 
	theMaterial[0] = concreteUniaxialMatObject1->getCopy();
	// Check allocation
	if (theMaterial[0] == 0) {
		opserr << "OrthotropicRotatingAngleConcreteT2DMaterial01::OrthotropicRotatingAngleConcreteT2DMaterial01 - failed to get a copy for Concrete 1\n";
		exit(-1);
	}

	theMaterial[1] = concreteUniaxialMatObject2->getCopy();
	// Check allocation
	if (theMaterial[1] == 0) {
		opserr << "OrthotropicRotatingAngleConcreteT2DMaterial01::OrthotropicRotatingAngleConcreteT2DMaterial01 - failed to get a copy for Concrete 2\n";
		exit(-1);
	}
		
	Eo = theMaterial[0]->getInitialTangent();
	Gmin = (1.0 / 100) * Eo / (2.0 * (1.0 + 0.2));

	this->revertToStart();
}

// Blank constructor
OrthotropicRotatingAngleConcreteT2DMaterial01::OrthotropicRotatingAngleConcreteT2DMaterial01()
	:NDMaterial(0, ND_TAG_OrthotropicRotatingAngleConcreteT2DMaterial01),
	strain_vec(3), stress_vec(3), tangent_matrix(3,3), initialTangentNDM(3,3), strainPrincipalDirection(3), 
	poissonRatios(2), Cstrain(3), CMaxMinStrainRec(3), TMaxMinStrainRec(3), pi(3.1415926535)	
{
	theMaterial = 0;

	this->revertToStart();
}

// Destructor
OrthotropicRotatingAngleConcreteT2DMaterial01::~OrthotropicRotatingAngleConcreteT2DMaterial01()
{
	// Delete the pointers
	if (theMaterial != 0) {
		for (int i = 0; i < 2; i++)
		{
			if (theMaterial[i])
				delete theMaterial[i];
		}
		delete[] theMaterial;
	}

}

// Get copy
NDMaterial* OrthotropicRotatingAngleConcreteT2DMaterial01::getCopy(void)
{
	OrthotropicRotatingAngleConcreteT2DMaterial01* theCopy =
		new OrthotropicRotatingAngleConcreteT2DMaterial01(this->getTag(),
			theMaterial[0],
			theMaterial[1],
			ecr,
			ec,
			rhoNDM,
			damageConstant1,
			damageConstant2);

	return theCopy;
}

// Get copy
NDMaterial* OrthotropicRotatingAngleConcreteT2DMaterial01::getCopy(const char* type)
{
	OrthotropicRotatingAngleConcreteT2DMaterial01* theModel =
		new OrthotropicRotatingAngleConcreteT2DMaterial01(this->getTag(),
			theMaterial[0],
			theMaterial[1],
			ecr,
			ec,
			rhoNDM,
			damageConstant1,
			damageConstant2);

	return theModel;
}

// Print
void OrthotropicRotatingAngleConcreteT2DMaterial01::Print(OPS_Stream& s, int flag)
{
	s << "\nOrthotropicRotatingAngleConcreteT2DMaterial01, nDMaterial tag: " << this->getTag() << endln;

	// Input values
	s << "Strain at Fcr: "<< ecr << endln;
	s << "Strain at Fc: " << ec << endln;
	s << "Density: " << rhoNDM << endln;

	s << "Hsu/Shu Poisson Ratios nu12: " << poissonRatios(0) << endln;
	s << "Hsu/Shu Poisson Ratios nu21: " << poissonRatios(1) << endln;

	if (isConcreteCracked) {
		s << "Concrete has already Cracked" << endln;
	}
	else {
		s << "Concrete has not Cracked" << endln;
	}
	//
	s << "Strains:" << endln;
	s << "EpsX = " << strain_vec(0) << ", EpsY = " << strain_vec(1) << ", GammaXY = " << strain_vec(2) << endln;
	s << "Damage Factor: " << CMaxMinStrainRec(0) << endln;
	//Strain and stress of the uniaxial materials
	s << "Strain and stress of Concrete Uniaxial Materials:" << endln;
	s << "Concrete 1: Strain = " << theMaterial[0]->getStrain() << ", Stress = " << theMaterial[0]->getStress() << endln;
	s << "Concrete 2: Strain = " << theMaterial[1]->getStrain() << ", Stress = " << theMaterial[1]->getStress() << endln;
}

int OrthotropicRotatingAngleConcreteT2DMaterial01::sendSelf(int commitTag, Channel& theChannel)
{	
	int res = 0;

	int dataTag = this->getDbTag();

	static Vector data(6);

	data(0) = this->getTag();
	data(1) = ecr;
	data(2) = ec;
	data(3) = rhoNDM;
	data(4) = damageConstant1;
	data(5) = damageConstant2;

	res += theChannel.sendVector(dataTag, commitTag,data);
	if (res < 0) {
		opserr << "WARNING OrthotropicRotatingAngleConcreteT2DMaterial01::sendself() - " << this->getTag() << " failed to send Vector\n";
		return res;
	}

	int matDbTag;

	static ID idData(4);
	int i;
	for (i = 0; i < 2; i++) {
		idData(i) = theMaterial[i]->getClassTag();
		matDbTag = theMaterial[i]->getDbTag();
		if (matDbTag == 0) {
			matDbTag = theChannel.getDbTag();
			if (matDbTag != 0)
				theMaterial[i]->setDbTag(matDbTag);
		}
		idData(i + 2) = matDbTag;
	}

	res += theChannel.sendID(dataTag, commitTag, idData);
	if (res < 0) {
		opserr << "WARNING OrthotropicRotatingAngleConcreteT2DMaterial01::sendself() - " << this->getTag() << " failed to send ID\n";
		return res;
	}

	for (i = 0; i < 2; i++) {
		res += theMaterial[i]->sendSelf(commitTag, theChannel);
		if (res < 0) {
			opserr << "OrthotropicRotatingAngleConcreteT2DMaterial01::sendself() - " << this->getTag() << " failed to send its Material\n";
			return res;
		}
	}

	return res;
}

int OrthotropicRotatingAngleConcreteT2DMaterial01::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	int res = 0;

	int dataTag = this->getDbTag();

	static Vector data(6);
	res += theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING OrthotropicRotatingAngleConcreteT2DMaterial01::recvSelf() - failed to receive Vector\n";
		return res;
	}

	this->setTag((int)data(0));
	ecr             = data(1);
	ec              = data(2);
	rhoNDM          = data(3);
	damageConstant1 = data(4);
	damageConstant2 = data(5);

	static ID idData(4);

	res += theChannel.recvID(dataTag, commitTag, idData);
	if (res < 0) {
		opserr << "WARNING OrthotropicRotatingAngleConcreteT2DMaterial01::recvSelf() - " << this->getTag() << " failed to receive ID\n";
		return res;
	}

	if (theMaterial == 0) {
		theMaterial = new UniaxialMaterial * [2];
		if (theMaterial == 0) {
			opserr << "OrthotropicRotatingAngleConcreteT2DMaterial01::recvSelf() - Could not allocate Uniaxial* array\n";
			return -1;
		}
		for (int i = 0; i < 2; i++) {
			int matClassTag = idData(i);
			int matDbTag = idData(i + 2);

			theMaterial[i] = theBroker.getNewUniaxialMaterial(matClassTag);
			if (theMaterial[i] == 0) {
				opserr << "OrthotropicRotatingAngleConcreteT2DMaterial01::recvself() - Broker could not create Uniaxial of class type " << matClassTag << endln;
				return -1;
			}
			theMaterial[i]->setDbTag(matDbTag);
			res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
			if (res < 0) {
				opserr << "OrthotropicRotatingAngleConcreteT2DMaterial01::recvSelf() - material " << i << " failed to recv itself\n";
				return res;
			}
		}
	}

	else {
		for (int i = 0; i < 2; i++) {
			int matClassTag = idData(i);
			int matDbTag = idData(i + 2);

			if (theMaterial[i]->getClassTag() != matClassTag) {
				delete theMaterial[i];
				theMaterial[i] = theBroker.getNewUniaxialMaterial(matClassTag);
				if (theMaterial[i] == 0) {
					opserr << "OrthotropicRotatingAngleConcreteT2DMaterial01::recvSelf() - material " << i << "failed to create\n";
					return -1;
				}
			}
			theMaterial[i]->setDbTag(matDbTag);
			res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
			if (res < 0) {
				opserr << "OrthotropicRotatingAngleConcreteT2DMaterial01::recvSelf() - material " << i << " failed to recv itself\n";
				return res;
			}
		}
	}

	return res;
}

// Get density
double OrthotropicRotatingAngleConcreteT2DMaterial01::getRho(void)
{
	return rhoNDM;
}

// Load strain from the element
int OrthotropicRotatingAngleConcreteT2DMaterial01::setTrialStrain(const Vector& v)
{
	// Calculate trial stress and tangent partial stiffness
	setTrialStrainPrincipalDirection(v);

	return 0;
}

int OrthotropicRotatingAngleConcreteT2DMaterial01::setTrialStrain(const Vector& v, const Vector& r)
{
	opserr << "YOU SHOULD NOT SEE THIS: OrthotropicRotatingAngleConcreteT2DMaterial01::setTrialStrain(const Vector& v, const Vector& r)\n";
	return this->setTrialStrain(v);
}

int OrthotropicRotatingAngleConcreteT2DMaterial01::setTrialStrianIncr(const Vector& v)
{
	// Add the increment trial strain to the trial strain
	Vector newTrialStrain(3);

	for (int i = 0; i < 3; i++) {
		newTrialStrain(i) = strain_vec(i) + v(i);
	}
	return this->setTrialStrain(newTrialStrain);
 }

int OrthotropicRotatingAngleConcreteT2DMaterial01::setTrialStrainIncr(const Vector& v, const Vector& r)
{
	opserr << "OrthotropicRotatingAngleConcreteT2DMaterial01::setTrialStrainIncr(const Vector& v, const Vector& r) -- should not be used!\n";
	return this->setTrialStrianIncr(v);
}

int OrthotropicRotatingAngleConcreteT2DMaterial01::setTrialStrainPrincipalDirection(const Vector& v)
{
	// Set values for strain_vec
	strain_vec(0) = v(0);
	strain_vec(1) = v(1);
	strain_vec(2) = v(2);

	// Calculate the Principal Direction of Strain
	this->calculateStrainPrincipalDirections01();

	// Calculate the Poisson Ratios
	this->calculatePoissonRatios(strainPrincipalDirection(0), strainPrincipalDirection(1));

	// Define the Poisson Ratio Matrix V = [1 / (1 - nu12 * nu21) nu12 / (1 - nu12 * nu21) 0; nu21 / (1 - nu12 * nu21) 1 / (1 - nu12 * nu21) 0; 0 0 1]
	double nu[2]; 
	nu[0] = poissonRatios(0); nu[1] = poissonRatios(1);
	double oneOverOneMinusNu12Nu21 = 1.0 / (1.0 - nu[0] * nu[1]);

	double V[3][3] = { {oneOverOneMinusNu12Nu21,nu[0] * oneOverOneMinusNu12Nu21,0.0},{nu[1] * oneOverOneMinusNu12Nu21,oneOverOneMinusNu12Nu21,0.0},{0.0,0.0,1.0} };

	// Calculate newUniaxialStrain in the orientation of the Principal Direction of Strain
	double newUniaxialStrainPD[3] = { 0.0,0.0,0.0 };

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			newUniaxialStrainPD[i] = newUniaxialStrainPD[i] + (strainPrincipalDirection(j) * V[i][j]);
		}
	}

	double e_11 = newUniaxialStrainPD[0];
	double e_22 = newUniaxialStrainPD[1];
	
	theMaterial[0]->setTrialStrain(e_11);
	theMaterial[1]->setTrialStrain(e_22);


	// Define the Max and Min Strains Recorded
	TMaxMinStrainRec(0) = CMaxMinStrainRec(0);
	TMaxMinStrainRec(1) = CMaxMinStrainRec(1);
	TMaxMinStrainRec(2) = CMaxMinStrainRec(2);

	double TMMSR[2];
	TMMSR[0] = TMaxMinStrainRec(1);
	TMMSR[1] = TMaxMinStrainRec(2);

	if (TMMSR[0] > min(newUniaxialStrainPD[0], newUniaxialStrainPD[1])) {
		TMMSR[0] = min(newUniaxialStrainPD[0], newUniaxialStrainPD[1]);
	}
	if (TMMSR[1] < max(newUniaxialStrainPD[0], newUniaxialStrainPD[1])) {
		TMMSR[1] = max(newUniaxialStrainPD[0], newUniaxialStrainPD[1]);
	}
	TMaxMinStrainRec(1) = TMMSR[0];
	TMaxMinStrainRec(2) = TMMSR[1];

	// Update the Damage Factor
	double TstrainRec = fabs(TMaxMinStrainRec(2) - TMaxMinStrainRec(1));
	double damageFactor = 1 / (1 + damageConstant1 * (pow(fabs(TstrainRec / ec), damageConstant2)));
	TMaxMinStrainRec(0) = damageFactor;

	// Calculate the Strain Transformation Matrix that goes from the Local Coord System to the orientation of the Principal Direction of Strain
	double TthetaPD[3][3], TthetaPDT[3][3], tangentNDM[3][3];;
	double* pTthetaPD = &TthetaPD[0][0];
	this->calculateStrainTransformationMatrix(pTthetaPD, thetaPrincipalDirection);

	// Extract the Stress in the Orientation of the Principal Direction
	double fc1 = theMaterial[0]->getStress();
	double fc2 = theMaterial[1]->getStress();

	// Aply the damage factor
	double s1 = damageFactor * fc1;
	double s2 = damageFactor * fc2;

	// Set the Stresses in the Principal Directions
	double stressNDMPD[3] = { s1,s2,0.0 }, stressNDM[3] = { 0.0,0.0,0.0 };

	// Calculate the traspose of TthetaPD
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			TthetaPDT[i][j] = TthetaPD[j][i];
		}
	}

	// Calculate the Stress in the Local Coord System 
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			stressNDM[i] += TthetaPDT[i][j] * stressNDMPD[j];
		}
	}
	for (int i = 0; i < 3; i++) {
		stress_vec(i) = stressNDM[i];
	}

	// Extract the Strain in the Orientation of the Principal Direction
	double e1 = theMaterial[0]->getStrain();
	double e2 = theMaterial[1]->getStrain();

	// Calculate the Tangent Matrix in the orientation of the Principal Direction of Strain in the Layer [E1 0 0; 0 E2 0; 0 0 G]

	// Extract the Tangent Stiffness in the Orientation of the Principal Direction 
	double Ect1 = theMaterial[0]->getTangent();
	double Ect2 = theMaterial[1]->getTangent();

	// Apply the damage factor
	double E1 = damageFactor * Ect1;
	double E2 = damageFactor * Ect2;

	// Calculate the Delta Strain
	double deltaStrain = e1 - e2;
	// Calculate the Delta Stress
	double deltaStress = s1 - s2;
	// Calculate the Shear Tangent for the Orthotropic Concrete
	double G;
	if (fabs(deltaStrain) <= dbl_Epsilon || fabs(deltaStress) <= dbl_Epsilon) {
		if (fabs(E1) < dbl_Epsilon && fabs(E2) < dbl_Epsilon) {
			G = 1.0e-5;
		}
		else {
			G = (E1 * E2) / (E1 * (1 + nu[0]) + E2 * (1 + nu[1]));
		}
	}
	else {
		G = 0.5 * deltaStress / deltaStrain;
	}
	if (G < Gmin) {
		G = Gmin;
	}
	// Calculate the Tangent Matrix in the Local Coord System [Tstress(-thetaPD)]*[E1 0  0 ; 0 E2 0 ; 0 0 G]*[V]*[Tstrain(thetaPD)], but Tstress(-thetaPD)] = [Tstrain(thetaPD)]^T
	double D[3][3] = { {E1,0,0},{0,E2,0},{0,0,G} };		// Material tangent matrix in the Principal Strain Direction

	// Inicializacion de la matriz
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			tangentNDM[i][j] = 0.0;
		}
	}
	// Calculate the Tangent Matrix in the Local Coord System
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				for (int m = 0; m < 3; m++) {
					for (int n = 0; n < 3; n++) {
						tangentNDM[m][n] += TthetaPDT[m][i] * (D[i][k] * V[k][j]) * TthetaPD[j][n];
					}
				}
			}
		}
	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			tangent_matrix(i, j) = tangentNDM[i][j];
		}
	}

	return 0;
}

const Matrix& OrthotropicRotatingAngleConcreteT2DMaterial01::getTangent(void)
{
	return tangent_matrix;
}

const Vector& OrthotropicRotatingAngleConcreteT2DMaterial01::getStress(void)
{
	return stress_vec;
}

const Vector& OrthotropicRotatingAngleConcreteT2DMaterial01::getStrain()
{
	return strain_vec;
}

int OrthotropicRotatingAngleConcreteT2DMaterial01::revertToStart(void)
{
	// Revert to start the values in the Concrete UniaxialMaterial Objects
	theMaterial[0]->revertToStart();
	theMaterial[1]->revertToStart();

	// Set stress, strain and tangent to zero
	strain_vec.Zero();
	stress_vec.Zero();
	tangent_matrix.Zero();
	initialTangentNDM.Zero();
	
	// Set to inicial some propierties of the NDMaterial
	poissonRatios.Zero();
	thetaPrincipalDirection = 0.0;
	isConcreteCracked = false;
	// Trial State Variables
	TMaxMinStrainRec.Zero();
	TMaxMinStrainRec(0) = 1.0;
	// Committed State Variables
	Cstrain.Zero();
	CMaxMinStrainRec.Zero();
	CMaxMinStrainRec(0) = 1.0;

	return 0;
}

// Calculate initial partial stiffness matrix
const Matrix& OrthotropicRotatingAngleConcreteT2DMaterial01::getInitialTangent(void)
{
	// Calculate the Strain Transformation Matrix
	double TthetaPD[3][3], TthetaPDT[3][3], tangentNDM[3][3];
	double* pTthetaPD = &TthetaPD[0][0];
	this->calculateStrainTransformationMatrix(pTthetaPD, thetaPrincipalDirection);

	// Calculate the Tangent Matrix in the orientation of the Principal Direction of Strain in the Layer [E1 0 0; 0 E2 0; 0 0 G]
	// Define the Poisson Ratio Matrix V = [ 1/(1-nu12*nu21) nu12/(1-nu12*nu21) 0 ; nu21/(1-nu12*nu21) 1/(1-nu12*nu21) 0 ; 0 0 1 ]
	double nu = 0.2;
	double oneOverOneMinusNu12Nu21 = 1.0 / (1.0 - nu * nu);
	double V[3][3] = { {oneOverOneMinusNu12Nu21,nu * oneOverOneMinusNu12Nu21,0.0},{nu * oneOverOneMinusNu12Nu21,oneOverOneMinusNu12Nu21,0.0},{0.0,0.0,1.0} };

	// Extract the Initial Tangent Stiffness in the Orientation of the Principal Direction 
	double E1 = theMaterial[0]->getInitialTangent();
	double E2 = theMaterial[1]->getInitialTangent();

	// Calculate the Initial Shear Tangent for the Orthotropic Concrete
	double G;
	G = Gmin * 100.0;
	// Calculate the Tangent Matrix in the Local Coord System [Tstress(-thetaPD)]*[E1 0  0 ; 0 E2 0 ; 0 0 G]*[V]*[Tstrain(thetaPD)], but Tstress(-thetaPD)] = [Tstrain(thetaPD)]^T
	double D[3][3] = { {E1,0.0,0.0},{0.0,E2,0.0},{0.0,0.0,G} };		// Material tangent matrix in the Principal Strain Direction

	// Traspose of TthetaPD
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			TthetaPDT[i][j] = TthetaPD[j][i];
		}
	}

	// Inicializacion de la matriz
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			tangentNDM[i][j] = 0.0;
		}
	}
	// Calculate the Tangent Matrix in the Local Coord System
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				for (int m = 0; m < 3; m++) {
					for (int n = 0; n < 3; n++) {
						tangentNDM[m][n] += TthetaPDT[m][i] * (D[i][k] * V[k][j]) * TthetaPD[j][n];
					}
				}
			}
		}
	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			initialTangentNDM(i, j) = tangentNDM[i][j];
		}
	}
	return initialTangentNDM;
}

Response* OrthotropicRotatingAngleConcreteT2DMaterial01::setResponse(const char** argv, int argc, OPS_Stream& theOutput)
{
	Response* theResponse = 0;
	if (strcmp(argv[0], "getEc") == 0) {
		theOutput.tag("NdMaterialOutput");
		theOutput.attr("matType", this->getClassType());
		theOutput.attr("matTag", this->getTag());
		theOutput.tag("ResponseType", "MaterialTag");
		theOutput.tag("ResponseType", "Ec");
		theOutput.endTag();

		Vector data1(2);
		data1.Zero();
		theResponse = new MaterialResponse(this, 101, data1);
	}
	else
		return this->NDMaterial::setResponse(argv, argc, theOutput);

	return theResponse;
}

int OrthotropicRotatingAngleConcreteT2DMaterial01::getResponse(int responseID, Information& matInformation)
{
	if (responseID == 101) {
		return matInformation.setVector(this->getEc());
	}
	else {
		return 0;
	}
}

// Function that returns input parameters (concrete Young's modulus) - added for MEFI3D by Maria Jose Nunez, UChile
Vector OrthotropicRotatingAngleConcreteT2DMaterial01::getEc(void)
{
	Vector input_par(2);

	input_par.Zero();

	input_par(0) = this->getTag();
	input_par(1) = Eo;

	return input_par;
}


int OrthotropicRotatingAngleConcreteT2DMaterial01::commitState(void)
{
	// commit the State in the Concrete UniaxialMaterial Objects
	theMaterial[0]->commitState();
	theMaterial[1]->commitState();

	// State Variable
	Cstrain(0) = strain_vec(0);
	Cstrain(1) = strain_vec(1);
	Cstrain(2) = strain_vec(2);

	// update the Damage Factor
	double TstrainRec = fabs(TMaxMinStrainRec(2) - TMaxMinStrainRec(1));
	double CstrainRec = fabs(CMaxMinStrainRec(2) - CMaxMinStrainRec(1));

	if (TstrainRec == CstrainRec) {
		CMaxMinStrainRec(0) = 1 / (1 + damageConstant1 * (pow(fabs(CstrainRec / ec), damageConstant2)));
		for (int i = 0; i < 3; i++) {
			TMaxMinStrainRec(i) = CMaxMinStrainRec(i);
		}
	}
	else {
		for (int i = 0; i < 3; i++) {
			CMaxMinStrainRec(i) = TMaxMinStrainRec(i);
		}
	}

	return 0;
}

// Revert to previously converged state
int OrthotropicRotatingAngleConcreteT2DMaterial01::revertToLastCommit(void)
{
	// revert to last committed state the values in the Concrete UniaxialMaterial Objects
	theMaterial[0]->revertToLastCommit();
	theMaterial[1]->revertToLastCommit();
	// State Variables
	for (int i = 0; i < 3; i++) {
		strain_vec(i) = Cstrain(i);
		TMaxMinStrainRec(i) = CMaxMinStrainRec(i);
	}
	return  0;
}

void OrthotropicRotatingAngleConcreteT2DMaterial01::calculateStrainPrincipalDirections01(void)
{
	double Tstrain[3];		//ex, ey, gamma
	double doubleThetaPD, cos2Theta, sin2Theta;
	
	// Get strain values from strain of element
	Tstrain[0] = strain_vec(0);		// ex
	Tstrain[1] = strain_vec(1);		// ey
	Tstrain[2] = strain_vec(2);		// gxy
	
	double averageStrains = 0.5 * (Tstrain[0] + Tstrain[1]);
	double deltaStrains = Tstrain[0] - Tstrain[1];
	double constantCS = pow(pow(deltaStrains, 2) + pow(Tstrain[2], 2), 0.5);
	if (deltaStrains == 0.0 && Tstrain[2] == 0.0) {
		doubleThetaPD = 0.0;
	}
	else {
		cos2Theta = deltaStrains / constantCS;
		sin2Theta = Tstrain[2] / constantCS;
		this->calculateAngle01(cos2Theta, sin2Theta, doubleThetaPD);

	}

	thetaPrincipalDirection = 0.5 * doubleThetaPD;
	strainPrincipalDirection(0) = averageStrains + 0.5 * constantCS;
	strainPrincipalDirection(1) = averageStrains - 0.5 * constantCS;
	strainPrincipalDirection(2) = 0.0;

	return;
}

void OrthotropicRotatingAngleConcreteT2DMaterial01::calculateAngle01(double cosTheta, double sinTheta, double &theta)
{
	double thetaFromCos, thetaFromSin;

	thetaFromCos = acos(cosTheta);
	thetaFromSin = asin(sinTheta);

	if (thetaFromCos <= pi / 2.0 + 1e-7) {
		theta = thetaFromSin;
	}
	else {
		if (thetaFromSin >= 0.0) {
			theta = thetaFromCos;
		}
		else {
			theta = -thetaFromCos;
		}
	}

	return;
}

void OrthotropicRotatingAngleConcreteT2DMaterial01::calculatePoissonRatios(double e1, double e2)
{
	double nu[2];

	if (e1 <= ec || e1 >= ecr || e2 <= ec || e2 >= ecr) {
		isConcreteCracked = true;
	}

	double Nuo = 0.2;

	if (isConcreteCracked) {
		nu[0] = 0.0;
		nu[1] = 0.0;
	}
	else {
		nu[0] = Nuo;
		nu[1] = Nuo;
	}
	if (!isConcreteCracked) {
		if (e2 <= ec / 2.0) {
			nu[0] = min(Nuo * (1.0 + 1.5 * (pow(2.0 * e2 / ec - 1.0, 2.0))), 0.5);
		}
		if (e1 <= ec / 2.0) {
			nu[1] = min(Nuo * (1.0 + 1.5 * (pow(2.0 * e1 / ec - 1.0, 2.0))), 0.5);
		}
	}

	poissonRatios(0) = nu[0];
	poissonRatios(1) = nu[1];

	return;
}

void OrthotropicRotatingAngleConcreteT2DMaterial01::calculateStrainTransformationMatrix(double* pTmatrixStrain, double theta)
{
	double cosTheta, sinTheta;
	int N = 3;
	cosTheta = cos(theta);
	sinTheta = sin(theta);

	// Uso de la formula de direccionamiento *(p+i*N+j): Accede al elemento C[i][j] de una matriz C de N columnas
	*(pTmatrixStrain + 0 * N + 0) = pow(cosTheta, 2);						//TmatrixStrain[0][0]
	*(pTmatrixStrain + 0 * N + 1) = pow(sinTheta, 2);						//TmatrixStrain[0][1]
	*(pTmatrixStrain + 0 * N + 2) = sinTheta * cosTheta;				    //TmatrixStrain[0][2]
	*(pTmatrixStrain + 1 * N + 0) = pow(sinTheta, 2);			    		//TmatrixStrain[1][0]
	*(pTmatrixStrain + 1 * N + 1) = pow(cosTheta, 2);						//TmatrixStrain[1][1]
	*(pTmatrixStrain + 1 * N + 2) = -sinTheta * cosTheta;				    //TmatrixStrain[1][2]
	*(pTmatrixStrain + 2 * N + 0) = -2.0 * sinTheta * cosTheta;				//TmatrixStrain[2][0]
	*(pTmatrixStrain + 2 * N + 1) = 2.0 * sinTheta * cosTheta;				//TmatrixStrain[2][1]
	*(pTmatrixStrain + 2 * N + 2) = pow(cosTheta, 2) - pow(sinTheta, 2);	//TmatrixStrain[2][2]
	return;
}

void OrthotropicRotatingAngleConcreteT2DMaterial01::calculateStressTransformationMatrix(double* pTmatrixStress, double theta)
{
	double cosTheta, sinTheta;
	int N = 3;
	cosTheta = cos(theta);
	sinTheta = sin(theta);		

	// Uso de la formula de direccionamiento *(p+i*N+j): Accede al elemento C[i][j] de una matriz C de N columnas
	*(pTmatrixStress + 0 * N + 0) = pow(cosTheta, 2);						//TmatrixStress[0][0]
	*(pTmatrixStress + 0 * N + 1) = pow(sinTheta, 2);						//TmatrixStress[0][1]
	*(pTmatrixStress + 0 * N + 2) = 2.0 * sinTheta * cosTheta;				//TmatrixStress[0][2]
	*(pTmatrixStress + 1 * N + 0) = pow(sinTheta, 2);			    		//TmatrixStress[1][0]
	*(pTmatrixStress + 1 * N + 1) = pow(cosTheta, 2);						//TmatrixStress[1][1]
	*(pTmatrixStress + 1 * N + 2) = -2.0 * sinTheta * cosTheta;				//TmatrixStress[1][2]
	*(pTmatrixStress + 2 * N + 0) = -sinTheta * cosTheta;				    //TmatrixStress[2][0]
	*(pTmatrixStress + 2 * N + 1) = sinTheta * cosTheta;				    //TmatrixStress[2][1]
	*(pTmatrixStress + 2 * N + 2) = pow(cosTheta, 2) - pow(sinTheta, 2);	//TmatrixStress[2][2]
	return;
}
