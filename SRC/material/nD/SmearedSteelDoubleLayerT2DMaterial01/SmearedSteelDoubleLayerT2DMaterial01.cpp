// Code written/implemented by: Fabian Rojas B.
//								Maria Jose Nunez
//
// Created: 03/2023
// 
// Description: This file contains the class definition for SmearedSteelDoubleLayerT2DMaterial01. 
// A SmearedSteelDoubleLayerT2DMaterial01 is a subclass of the class NDMaterial and corresponds to the abstract representation
// of a double perpendicular Smeared Steel layers (plane stress) 2D Material with a tangent formulation. Each layer works only 
// in the direction of the bars, so a uniaxial constitutive model is used to represent the behavior of reinforcing steel bars in each direction.
//
// Reference:
// 1. Rojas, F., Anderson, J. C., Massones, L. M. (2016). A nonlinear quadrilateral layered membrane with drilling degrees of freedom for 
// the modeling of reinforced concrete walls. Engineering Structures, 124, 521-538.
//
// Source: \OpenSees\SRC\material\nD\SmearedSteelDoubleLayerT2DMaterial01
//
// Rev: 1.0

#include <SmearedSteelDoubleLayerT2DMaterial01.h>
#include <UniaxialMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <FEM_ObjectBroker.h>
#include <DummyStream.h>
#include <elementAPI.h>
#include <cmath>

// Read input parameters and build the material
void* OPS_SmearedSteelDoubleLayerT2DMaterial01()
{
	// Pointer to a ND material that will be returned
	NDMaterial* theMaterial = 0;

	int numArgs = OPS_GetNumRemainingInputArgs();
	
	// Parse the script for material parameters
	if (numArgs != 6) {			// total # of input parameters
		opserr << "Invalid #Args want: nDMaterial SmearedSteelDoubleLayerT2DMaterial01 $tag $s1 $s2 $ratioSteelLayer1 $ratioSteelLayer2 $OrientationEmbeddedSteel" << endln;
		return 0;
	}

	int tag;					// nDMaterial tag
	int iData[2];				// # of uniaxial materials
	double dData[3];			// # of material parameters

	// nDMaterial tag
	int numData = 1;
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "WARNING invalid integer tag: nDMaterial SmearedSteelDoubleLayerT2DMaterial01" << endln; 
		return 0;
	}

	// Material tags of 2 steel layers
	numData = 2;
	if (OPS_GetInt(&numData, iData) != 0) {
		opserr << "WARNING: Invalid uniaxialMaterial tags: nDMaterial SmearedSteelDoubleLayerT2DMaterial01: " << tag << "\n";
		return 0;
	}

	// Material parameters
	numData = 3;
	if (OPS_GetDouble(&numData, dData) != 0) {
		opserr << "WARNING invalid data: nDMaterial SmearedSteelDoubleLayerT2DMaterial01: " << tag << "\n";
		return 0;
	}

	// Get pointers to Uniaxial material
	// Steel Layer 1
	UniaxialMaterial* theUniaxialMaterial1 = OPS_getUniaxialMaterial(iData[0]);
	if (theUniaxialMaterial1 == 0) {
		opserr << "WARNING material not found\n";
		opserr << "Material: " << iData[0];
		opserr << "\nSmearedSteelDoubleLayerT2DMaterial01: " << tag << endln;
		return 0;
	}

	// Steel Layer 2
	UniaxialMaterial* theUniaxialMaterial2 = OPS_getUniaxialMaterial(iData[1]);
	if (theUniaxialMaterial2 == 0) {
		opserr << "WARNING material not found\n";
		opserr << "Material: " << iData[1];
		opserr << "\nSmearedSteelDoubleLayerT2DMaterial01: " << tag << endln;
		return 0;
	}
	// Create the SmearedSteelDoubleLayerT2DMaterial01
	theMaterial = new SmearedSteelDoubleLayerT2DMaterial01(tag, theUniaxialMaterial1 , theUniaxialMaterial2 ,dData[0], dData[1], dData[2]);
	
	if (theMaterial == 0) {
		opserr << "WARNING ran out of memory creating material\n";
		opserr << "SmearedSteelDoubleLayerT2DMaterial01: " << tag << endln;
		return 0;
	}

	return theMaterial;
}

// Typical Constructor
SmearedSteelDoubleLayerT2DMaterial01::SmearedSteelDoubleLayerT2DMaterial01(int tag,			// nDMaterial tag
	UniaxialMaterial* steelULayer1Obj,														// steel Layer 1
	UniaxialMaterial* steelULayer2Obj,														// steel Layer 2
	double ratioSteelLayer1,																// ratio of smeared steel layer 1
	double ratioSteelLayer2,																// ratio of smeared steel layer 2
	double OrientationEmbeddedSteel)														// orientation of the smeared steel layer
	: NDMaterial(tag, ND_TAG_SmearedSteelDoubleLayerT2DMaterial01),
	ratioLayer1(ratioSteelLayer1), ratioLayer2(ratioSteelLayer2), thetaSmearedSteel(OrientationEmbeddedSteel),
	Tstrain(2), Tstress(2), Ttangent(2), TstrainLayer(3), strainPrincipalDirection(3),
	Cstrain(2), Cstress(2), Ctangent(2), CstrainLayer(3), stressNDM(3), tangentNDM(3,3), initialTangentNDM(3,3), pi(3.1415926535)
{

	// Allocate pointers for uniaxial materials 
	theMaterial = new UniaxialMaterial * [2];
	if (theMaterial == 0) {
		opserr << " SmearedSteelDoubleLayerT2DMaterial01::SmearedSteelDoubleLayerT2DMaterial01 - failed allocate material array\n";
		exit(-1);
	}

	// Get the copy for steelLayer1
	theMaterial[0] = steelULayer1Obj->getCopy();
	// Check allocation
	if (theMaterial[0] == 0) {
		opserr << " SmearedSteelDoubleLayerT2DMaterial01::SmearedSteelDoubleLayerT2DMaterial01 - failed to get a copy for SteelLayer1\n";
		exit(-1);
	}

	// Get the copy for steelLayer2
	theMaterial[1] = steelULayer2Obj->getCopy();
	// Check allocation
	if (theMaterial[1] == 0) {
		opserr << " SmearedSteelDoubleLayerT2DMaterial01::SmearedSteelDoubleLayerT2DMaterial01 - failed to get a copy for SteelLayer2\n";
		exit(-1);
	}

	// Set initial values
	thetaPrincipalDirection = 0.0;
	strainPrincipalDirection(0) = 0.0;
	strainPrincipalDirection(1) = 0.0;
	strainPrincipalDirection(2) = 0.0;

	for (int i = 0; i < 2; i++) {
		Ctangent(i) = theMaterial[i]->getInitialTangent();
	}

	for (int i = 0; i < 2; i++) {
		Cstrain(i) = 0.0;
		Cstress(i) = 0.0;
		Tstrain(i) = Cstrain(i);
		Tstress(i) = Cstress(i);
		Ttangent(i) = Ctangent(i);
	}

	for (int i = 0; i < 3; i++) {
		CstrainLayer(i) = 0.0;
		TstrainLayer(i) = CstrainLayer(i);
	}

	// Set stress, strain and tangent to zero
	stressNDM.Zero();
	tangentNDM.Zero();
	initialTangentNDM.Zero();
}

// Blanck constructor 
SmearedSteelDoubleLayerT2DMaterial01::SmearedSteelDoubleLayerT2DMaterial01()
	:NDMaterial(0, ND_TAG_SmearedSteelDoubleLayerT2DMaterial01),
	Tstrain(2), Tstress(2), Ttangent(2), TstrainLayer(3), strainPrincipalDirection(3),
	Cstrain(2), Cstress(2), Ctangent(2), CstrainLayer(3), stressNDM(3), tangentNDM(3, 3), initialTangentNDM(3, 3), pi(3.1415926535)
{
	theMaterial = 0;
	this->revertToStart();
}

// Destructor
SmearedSteelDoubleLayerT2DMaterial01::~SmearedSteelDoubleLayerT2DMaterial01()
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
NDMaterial* SmearedSteelDoubleLayerT2DMaterial01::getCopy(void)
{
	SmearedSteelDoubleLayerT2DMaterial01* theCopy =
		new SmearedSteelDoubleLayerT2DMaterial01(this->getTag(),
			theMaterial[0],
			theMaterial[1],
			ratioLayer1,
			ratioLayer2,
			thetaSmearedSteel);

	return theCopy;
}

// Get copy
NDMaterial* SmearedSteelDoubleLayerT2DMaterial01::getCopy(const char* type)
{
	SmearedSteelDoubleLayerT2DMaterial01* theCopy =
		new SmearedSteelDoubleLayerT2DMaterial01(this->getTag(),
			theMaterial[0],
			theMaterial[1],
			ratioLayer1,
			ratioLayer2,
			thetaSmearedSteel);

	return theCopy;
}

// Print
void SmearedSteelDoubleLayerT2DMaterial01::Print(OPS_Stream& s, int flag)
{
	s << "\n SmearedSteelDoubleLayerT2DMaterial01, nDMaterial tag: " << this->getTag() << endln;

	// Input values
	s << "Orientation of the smeared steel layers: " << thetaSmearedSteel << endln;
	s << "Reinforcing ratio of the smeared steel layer 1: " << ratioLayer1 << endln;
	s << "Reinforcing ratio of the smeared steel layer 2: " << ratioLayer2 << endln;

	s << "Strain:" << endln;
	s << "EpsX = " << Tstrain(0) << ", EpsY = " << Tstrain(1) << ", GammaXY = " << Tstrain(2) << endln;
}

int SmearedSteelDoubleLayerT2DMaterial01::sendSelf(int commitTag, Channel& theChannel)
{
	int res = 0;

	int dataTag = this->getDbTag();

	static Vector data(4);

	data(0) = this->getTag();
	data(1) = ratioLayer1;
	data(2) = ratioLayer2;
	data(3) = thetaSmearedSteel;

	res += theChannel.sendVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING SmearedSteelDoubleLayerT2DMaterial01::sendSelf() - " << this->getTag() << " failed to send Vector\n";
		return res;
	}

	int matDbTag;

	static ID idData(4);

	int i;
	for (i = 0; i < 2; i++)
	{
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
		opserr << "WARNING SmearedSteelDoubleLayerT2DMaterial01::sendSelf() - " << this - getTag() << " failed to send ID\n";
		return res;
	}

	for (i = 0; i < 2; i++) {
		res += theMaterial[i]->sendSelf(commitTag, theChannel);
		if (res < 0) {
			opserr << "SmearedSteelDoubleLayerT2DMaterial01::sendSelf() - " << this - getTag() << " failed to send its Material\n";
			return res;
		}
	}

	return res;
}

int SmearedSteelDoubleLayerT2DMaterial01::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	int res = 0;

	int dataTag = this->getDbTag();

	static Vector data(4);
	res += theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING SmearedSteelDoubleLayerT2DMaterial01::recvSelf() - failed to receive Vector\n";
		return res;
	}

	this->setTag((int)data(0));
	ratioLayer1       = data(1);
	ratioLayer2       = data(2);
	thetaSmearedSteel = data(3);

	static ID idData(4);
	res += theChannel.recvID(dataTag, commitTag, idData);
	if (res < 0) {
		opserr << "WARNING SmearedSteelDoubleLayerT2DMaterial01::recvSelf() - " << this->getTag() << " failed to receive ID\n";
		return res;
	}

	if (theMaterial == 0) {
		theMaterial = new UniaxialMaterial * [2];
		if (theMaterial == 0) {
			opserr << "SmearedSteelDoubleLayerT2DMaterial01::recvSelf() - Could not allocate Uniaxial* array\n";
			return -1;
		}
		for (int i = 0; i < 2; i++) {
			int matClassTag = idData(i);
			int matDbTag = idData(i + 2);
			theMaterial[i] = theBroker.getNewUniaxialMaterial(matClassTag);
			if (theMaterial[i] == 0) {
				opserr << "SmearedSteelDoubleLayerT2DMaterial01::recvSelf() - Broker could not create Uniaxial of class type " << matClassTag << endln;
				return -1;
			}
			theMaterial[i]->setDbTag(matDbTag);
			res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
			if (res < 0) {
				opserr << "SmearedSteelDoubleLayerT2DMaterial01::recvSelf() - material " << i << " failed to receive itself\n";
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
					opserr << "SmearedSteelDoubleLayerT2DMaterial01::recvself() - material " << i << " failed to create\n";
					return -1;
				}
			}
			theMaterial[i]->setDbTag(matDbTag);
			res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
			if (res < 0) {
				opserr << "SmearedSteelDoubleLayerT2DMaterial01::recvSelf() - material " << i << " failed to receive itself\n";
				return res;
			}
		}
	}

	return res;
}

int SmearedSteelDoubleLayerT2DMaterial01::setTrialStrain(const Vector& v)
{
	// Calculate trial stress and tangent partial stiffness
	setTrialStrainPrincipalDirection(v);

	return 0;
}

// Unused trial strain functions
int SmearedSteelDoubleLayerT2DMaterial01::setTrialStrain(const Vector& v, const Vector& r)
{
	opserr << "SmearedSteelDoubleLayerT2DMaterial01 :: setTrialStrain(const Vector& v, const Vector& r) -- should not be used! \n";
	return this->setTrialStrain(v);
}

int SmearedSteelDoubleLayerT2DMaterial01::setTrialStrainIncr(const Vector& v)
{
	// Add the increment trial strain to the trial strain
	Vector newTrialStrain(3);

	for (int i = 0; i < 3; i++) {
		newTrialStrain(i) = TstrainLayer(i) + v(i);
	}
	return this->setTrialStrain(newTrialStrain);
}

int SmearedSteelDoubleLayerT2DMaterial01::setTrialStrainIncr(const Vector& v, const Vector& r)
{
	opserr << "SmearedSteelDoubleLayerT2DMaterial01::setTrialStrainIncr(const Vector& v, const Vector& r) -- should not be used!\n";
	return this->setTrialStrainIncr(v);
}

int SmearedSteelDoubleLayerT2DMaterial01::setTrialStrainPrincipalDirection(const Vector& v)
{
	TstrainLayer(0) = v(0);
	TstrainLayer(1) = v(1);
	TstrainLayer(2) = v(2);

	// Define the Principal Direction of Strain
	this->calculateStrainPrincipalDirections01();

	// Calculate the Transformation Matrix that go from the Principal Direction of Strain to the orientation of the Steel Layer T(thetaES - thetaPD)
	double TthetaESMinusthetaPD[3][3];
	double* pTthetaESMinusthetaPD = &TthetaESMinusthetaPD[0][0];

	double theta = thetaSmearedSteel - thetaPrincipalDirection;

	this->calculateStrainTransformationMatrix(pTthetaESMinusthetaPD, theta);

	// Calculate the newUniaxialStrain in the in the Orientation of the Smeared Embedded Steel Layer 
	double newTrialStrain[3] = { 0,0,0 };

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			newTrialStrain[i] = newTrialStrain[i] + (TthetaESMinusthetaPD[i][j]* strainPrincipalDirection(j));
		}
	}

	Tstrain(0) = newTrialStrain[0];
	Tstrain(1) = newTrialStrain[1];

	theMaterial[0]->setTrialStrain(newTrialStrain[0]);
	theMaterial[1]->setTrialStrain(newTrialStrain[1]);

	// Calculate the Strain Transformation Matrix that go from the Local Coordinate System to the orientation of the Steel Layer T(thetaES)
	double TthetaES[3][3], TthetaEST[3][3];
	double* pTthetaES = &TthetaES[0][0];

	this->calculateStrainTransformationMatrix(pTthetaES, thetaSmearedSteel);

	// Calculate the Stress in the Local Coord System
	double fs1, fs2;
	fs1 = theMaterial[0]->getStress();
	fs2 = theMaterial[1]->getStress();

	double stressES[3] = { ratioLayer1 * fs1, ratioLayer2 * fs2, 0.0 };
	
	// traspose of TthetaES
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			TthetaEST[i][j] = TthetaES[j][i];
		}
	}

	stressNDM.Zero();

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			stressNDM(i) += TthetaEST[i][j] * stressES[j];
		}
	}

	// Calculate the Tangent Matrix in the Local Coord System
	double Es1, Es2;
	Es1 = theMaterial[0]->getTangent();
	Es2 = theMaterial[1]->getTangent();

	double Ds[3][3] = { {ratioLayer1 * Es1,0.0,0.0},{0.0,ratioLayer2 * Es2,0.0},{0.0,0.0,0.0} };

	tangentNDM.Zero();

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				for (int l = 0; l < 3; l++) {
					tangentNDM(i, l) += TthetaEST[i][j] * Ds[j][k] * TthetaES[k][l];
				}
			}
		}
	}

	return 0;
}

const Vector& SmearedSteelDoubleLayerT2DMaterial01::getStrain(void)
{
	return Tstrain;
}

const Vector& SmearedSteelDoubleLayerT2DMaterial01::getStress(void)
{
	return stressNDM;
}

const Matrix& SmearedSteelDoubleLayerT2DMaterial01::getTangent(void)
{
	return tangentNDM;
}

const Matrix& SmearedSteelDoubleLayerT2DMaterial01::getInitialTangent(void)
{
	// Calculate the Strain Transformation Matrix that go from the Local Coordinate System to the orientation of the Steel Layer T(thetaES)
	double TthetaES[3][3], TthetaEST[3][3];
	double* pTthetaES = &TthetaES[0][0];

	this->calculateStrainTransformationMatrix(pTthetaES, thetaSmearedSteel);

	// Calculate the Tangent Matrix in the Orientation of the Smeared Embedded Steel Layer 
	double E0s1, E0s2;
	E0s1 = theMaterial[0]->getInitialTangent();
	E0s2 = theMaterial[1]->getInitialTangent();

	double D0s[3][3] = { {ratioLayer1 * E0s1,0.0,0.0},{0.0,ratioLayer2 * E0s2,0.0},{0.0,0.0,0.0} };
	
	// Calculate the Tangent Matrix in the Local Coord System
	// Traspose of TthetaES
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			TthetaEST[i][j] = TthetaES[j][i];
		}
	}

	initialTangentNDM.Zero();
	
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				for (int l = 0; l < 3; l++) {
					initialTangentNDM(i, l) += TthetaEST[i][j] * D0s[j][k] * TthetaES[k][l];
				}
			}
		}
	}

	// Add a small value int he 3,3 component of the stiffness Matrix to avoid inestability in the analysis
	initialTangentNDM(2, 2) = initialTangentNDM(2, 2) + 1e-5;


	return initialTangentNDM;
}

int SmearedSteelDoubleLayerT2DMaterial01::commitState(void)
{
	theMaterial[0]->commitState();
	theMaterial[1]->commitState();

	// Commit the state values
	CstrainLayer(0) = TstrainLayer(0);
	CstrainLayer(1) = TstrainLayer(1);
	CstrainLayer(2) = TstrainLayer(2);

	Cstrain(0) = Tstrain(0);
	Cstrain(1) = Tstrain(1);

	Cstress(0) = theMaterial[0]->getStress();
	Cstress(1) = theMaterial[1]->getStress();

	Ctangent(0) = theMaterial[0]->getTangent();
	Ctangent(1) = theMaterial[1]->getTangent();

	return 0;
}

int SmearedSteelDoubleLayerT2DMaterial01::revertToLastCommit(void)
{
	// Revert to last committed values
	TstrainLayer(0) = CstrainLayer(0);
	TstrainLayer(1) = CstrainLayer(1);
	TstrainLayer(2) = CstrainLayer(2);

	for (int i = 0; i < 2; i++) {
		Tstrain(i) = Cstrain(i);
		Tstress(i) = Cstress(i);
		Ttangent(i) = Ctangent(i);
	}

	theMaterial[0]->revertToLastCommit();
	theMaterial[1]->revertToLastCommit();

	return 0;
}

int SmearedSteelDoubleLayerT2DMaterial01::revertToStart(void)
{
	// Set to the initial values some of the propierties
	thetaPrincipalDirection = 0.0;

	// State Variables
	CstrainLayer.Zero();
	Cstrain.Zero();
	Cstress.Zero();
	Ctangent(0) = theMaterial[0]->getInitialTangent();
	Ctangent(1) = theMaterial[1]->getInitialTangent();

	// Set stress, strain and tangent to zero
	stressNDM.Zero();
	tangentNDM.Zero();
	initialTangentNDM.Zero();

	// Revert all uniaxial materials to start
	theMaterial[0]->revertToStart();
	theMaterial[1]->revertToStart();

	// Set Initial Trial Values
	TstrainLayer(0) = CstrainLayer(0);
	TstrainLayer(1) = CstrainLayer(1);
	TstrainLayer(2) = CstrainLayer(2);

	for (int i = 0; i < 2; i++) {
		Tstrain(i) = Cstrain(i);
		Tstress(i) = Cstress(i);
		Ttangent(i) = Ctangent(i);
	}

	return 0;
}

Response* SmearedSteelDoubleLayerT2DMaterial01::setResponse(const char** argv, int argc, OPS_Stream& theOutput)
{
	return NDMaterial::setResponse(argv, argc, theOutput);
}

int SmearedSteelDoubleLayerT2DMaterial01::getResponse(int responseID, Information& matInformation)
{
	return NDMaterial::getResponse(responseID, matInformation);
}

void SmearedSteelDoubleLayerT2DMaterial01::calculateStrainPrincipalDirections01(void)
{
	double strain_vec[3];					// ex, ey, gamma
	double doubleThetaPD, cos2Theta, sin2Theta;

	// Get strain values from strain of element
	strain_vec[0] = TstrainLayer(0);		// ex
	strain_vec[1] = TstrainLayer(1);		// ey
	strain_vec[2] = TstrainLayer(2);		// gxy

	double averageStrains = 0.5 * (strain_vec[0] + strain_vec[1]);
	double deltaStrains = strain_vec[0] - strain_vec[1];
	double constantCS = pow(pow(deltaStrains, 2) + pow(strain_vec[2], 2), 0.5);
	if (deltaStrains == 0.0 && strain_vec[2] == 0.0) {
		doubleThetaPD = 0.0;
	}
	else {
		cos2Theta = deltaStrains / constantCS;
		sin2Theta = strain_vec[2] / constantCS;
		this->calculateAngle01(cos2Theta, sin2Theta, doubleThetaPD);

	}
	thetaPrincipalDirection = 0.5 * doubleThetaPD;

	strainPrincipalDirection(0) = averageStrains + 0.5 * constantCS;
	strainPrincipalDirection(1) = averageStrains - 0.5 * constantCS;
	strainPrincipalDirection(2) = 0.0;

	return;
}

void SmearedSteelDoubleLayerT2DMaterial01::calculateAngle01(double cosTheta, double sinTheta, double& theta)
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

void SmearedSteelDoubleLayerT2DMaterial01::calculateStrainTransformationMatrix(double* pTmatrixStrain, double theta)
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
