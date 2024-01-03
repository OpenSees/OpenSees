// Code written/implemented by: Fabian Rojas B.
//								Maria Jose Nunez
//
// Created: 07/2023
// 
// Description: This file contains the ReinforcedConcreteLayerMembraneSection02 class definition
// A ReinforcedConcreteLayerMembraneSection02 is a subclass of the sectionForceDeformation class and corresponds to the abstract representation
// for the stress-strain behavior for a Reinforced Concrete Layer Membrane Element in the Finite Element Method or Structural Analysis. 
//
// Reference:
// 1. Rojas, F., Anderson, J. C., Massones, L. M. (2016). A nonlinear quadrilateral layered membrane with drilling degrees of freedom for 
// the modeling of reinforced concrete walls. Engineering Structures, 124, 521-538.
//
// Source: \OpenSees\SRC\material\section\ReinforcedConcreteLayerMembraneSection
//
// Rev: 1.0

#include <ReinforcedConcreteLayerMembraneSection02.h>
#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <Parameter.h>
#include <MaterialResponse.h>
#include <algorithm>				/*min, max*/
#include <DummyStream.h>

using namespace std;

// Read input parameters and build the section
void* OPS_ReinforcedConcreteLayerMembraneSection02()
{
	int numArgs = OPS_GetNumRemainingInputArgs();

	// Parse the script for material parameters
	if (numArgs != 3) {
		opserr << "Want: ReinforcedConcreteLayerMembraneSection02 $secTag $matTag $Thickness" << endln;
		return 0;
	}

	int tag;						// section tag
	int matTag;						// nDMaterial tag
	double thickness;				// macro-fiber thickness

	// section tag
	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) {
		opserr << "WARNING invalid section ReinforcedConcreteLayerMembraneSection02 tag" << endln;
		return 0;
	}
	
	// nDMaterial tag
	numdata = 1;
	if (OPS_GetIntInput(&numdata, &matTag) < 0) {
		opserr << "WARNING invalid nDMaterial tag" << endln;
		return 0;
	}
	
	// thickness
	numdata = 1;
	if (OPS_GetDoubleInput(&numdata, &thickness) < 0) {
		opserr << "Invalid thick parameter for ReinforcedConcreteLayerMembraneSection02 " << tag << endln;
		return 0;
	}

	// Get pointer to nDMaterial
	NDMaterial* theNDMaterial = OPS_getNDMaterial(matTag);
	if (theNDMaterial == 0) {
		opserr << "WARNING material not found\n";
		opserr << "Material: " << matTag;
		opserr << "\nReinforcedConcreteLayerMembraneSection02: " << tag << endln;
		return 0;
	}

	SectionForceDeformation* theSection = new ReinforcedConcreteLayerMembraneSection02(tag, theNDMaterial, thickness);

	return theSection;
}

//static vector and matrices
ID ReinforcedConcreteLayerMembraneSection02::array(8);

// Full constructor
ReinforcedConcreteLayerMembraneSection02::ReinforcedConcreteLayerMembraneSection02(int tag,		// section tag
	NDMaterial* RCMaterialObject,																// reinforced concrete nDMaterial tag
	double Thickness)																			//macro-fiber thickness

	:SectionForceDeformation(tag, SEC_TAG_ReinforcedConcreteLayerMembraneSection02),
	t(Thickness), CSectionStrain(3), CSectionStress(3), CSectionTangent(3,3),
	TSectionStrain(3), TSectionStress(3), TSectionTangent(3,3), InitialTangent(3,3)
{
	// Set initial values
	for (int i = 0; i < 3; i++) {
		CSectionStrain(i) = 0.0;
		CSectionStress(i) = 0.0;
		TSectionStrain(i) = 0.0;
		TSectionStress(i) = 0.0;
	}

	// Get a copy of the RC material
	TheRCMaterial = RCMaterialObject->getCopy();
	// Check allocation
	if (TheRCMaterial == 0) {
		opserr << "ReinforcedConcreteLayerMembraneSection02::ReinforcedConcreteLayerMembraneSection02() - Failed to get a copy for RCMaterial\n";
		exit(-1);
	}

}

// Blank constructor (constructor for blank object that recvSelf needs to be invoked upon) (constructor which should be invoked by an FEM_ObjectBroker only)
ReinforcedConcreteLayerMembraneSection02::ReinforcedConcreteLayerMembraneSection02() :
	SectionForceDeformation(0, SEC_TAG_ReinforcedConcreteLayerMembraneSection02),
	CSectionStrain(3), CSectionStress(3), CSectionTangent(3, 3),
	TSectionStrain(3), TSectionStress(3), TSectionTangent(3, 3), InitialTangent(3, 3)
{
	this->revertToStart();
}

// Destructor (clean up memory the ReinforcedConcreteLayerMembraneSection02 objects allocates)
ReinforcedConcreteLayerMembraneSection02::~ReinforcedConcreteLayerMembraneSection02()
{
	if (TheRCMaterial != 0)
		delete TheRCMaterial;
}

int ReinforcedConcreteLayerMembraneSection02::setTrialSectionDeformation(const Vector& newTrialSectionStrain)
{
	// Set the membrane strain (e11,e22,e12) in the section
	TSectionStrain(0) = newTrialSectionStrain(0);
	TSectionStrain(1) = newTrialSectionStrain(1);
	TSectionStrain(2) = newTrialSectionStrain(2);

	double sectionResultantStress[3];		// section resultant stress: {N11 ; N22 ; N12 } = Sum(Stress_ic(zTop_ic - zBottom_ic)) + Sum(Stress_is(t_is))
	sectionResultantStress[0] = 0.0; sectionResultantStress[1] = 0.0; sectionResultantStress[2] = 0.0;

	double sectionTangent[3][3];			// membrane section stiffness: [Dm] = h*[MaterialTangent] is the membrane stiffness Tangent = Sum(Tangent_ic(zTop_ic - zBottom_ic)) + Sum(Tangent_is(t_is))
	sectionTangent[0][0] = 0.0; sectionTangent[0][1] = 0.0;	sectionTangent[0][2] = 0.0;
	sectionTangent[1][0] = 0.0; sectionTangent[1][1] = 0.0;	sectionTangent[1][2] = 0.0;
	sectionTangent[2][0] = 0.0; sectionTangent[2][1] = 0.0;	sectionTangent[2][2] = 0.0;

	// Set the strain
	TheRCMaterial->setTrialStrain(newTrialSectionStrain);
	const Vector& stressNDM = TheRCMaterial->getStress();
	const Matrix& tangentNDM = TheRCMaterial->getTangent();

	double fc11 = stressNDM(0);
	double fc22 = stressNDM(1);
	double fc12 = stressNDM(2);

	sectionResultantStress[0] = t * fc11;
	sectionResultantStress[1] = t * fc22;
	sectionResultantStress[2] = t * fc12;

	double Ec00 = tangentNDM(0, 0); double Ec01 = tangentNDM(0, 1); double Ec02 = tangentNDM(0, 2);
	double Ec10 = tangentNDM(1, 0); double Ec11 = tangentNDM(1, 1); double Ec12 = tangentNDM(1, 2);
	double Ec20 = tangentNDM(2, 0); double Ec21 = tangentNDM(2, 1); double Ec22 = tangentNDM(2, 2);

	sectionTangent[0][0] = t * Ec00; sectionTangent[0][1] = t * Ec01;	sectionTangent[0][2] = t * Ec02;
	sectionTangent[1][0] = t * Ec10; sectionTangent[1][1] = t * Ec11;	sectionTangent[1][2] = t * Ec12;
	sectionTangent[2][0] = t * Ec20; sectionTangent[2][1] = t * Ec21;	sectionTangent[2][2] = t * Ec22;
	
	// Set the stress and tangent
	TSectionStress(0) = sectionResultantStress[0]; TSectionStress(1) = sectionResultantStress[1]; TSectionStress(2) = sectionResultantStress[2];

	TSectionTangent(0, 0) = sectionTangent[0][0]; TSectionTangent(0, 1) = sectionTangent[0][1]; TSectionTangent(0, 2) = sectionTangent[0][2];
	TSectionTangent(1, 0) = sectionTangent[1][0]; TSectionTangent(1, 1) = sectionTangent[1][1]; TSectionTangent(1, 2) = sectionTangent[1][2];
	TSectionTangent(2, 0) = sectionTangent[2][0]; TSectionTangent(2, 1) = sectionTangent[2][1]; TSectionTangent(2, 2) = sectionTangent[2][2];

	return 0;
}

const Vector& ReinforcedConcreteLayerMembraneSection02::getSectionDeformation(void)
{
	return TSectionStrain;
}

const Vector& ReinforcedConcreteLayerMembraneSection02::getStressResultant(void)
{
	return TSectionStress;
}

const Matrix& ReinforcedConcreteLayerMembraneSection02::getSectionTangent(void)
{
	return TSectionTangent;
}

const Vector& ReinforcedConcreteLayerMembraneSection02::getCommittedStrain(void)
{
	return CSectionStrain;
}

const Vector& ReinforcedConcreteLayerMembraneSection02::getCommittedStress(void)
{
	return CSectionStress;
}

double ReinforcedConcreteLayerMembraneSection02::getRho(void)
{
	double rhoH = 0.0;

	rhoH = (TheRCMaterial->getRho()) * t;

	return rhoH;
}

double ReinforcedConcreteLayerMembraneSection02::getEcAvg(void)
{
	DummyStream theDummyStream;
	char aa[80] = "getInputParameters";
	const char* argv[1];
	argv[0] = aa;

	double EcAvg = 0.0;

	Response* theResponse = TheRCMaterial->setResponse(argv, 1, theDummyStream);

	if (theResponse == 0) {
		opserr << "ReinforcedConcreteLayerMembraneSection01::ReinforcedConcreteLayerMembraneSection01 - failed to get input parameters for OrthotropicRotatingAngleConcreteT2DMaterial01 with tag : " << this->getTag() << "\n";
		exit(-1);
	}
	theResponse->getResponse();
	Information& theInfoInput = theResponse->getInformation();
	const Vector& InputNDMat = theInfoInput.getData();

	EcAvg = InputNDMat[9];

	return EcAvg;
}

const Matrix& ReinforcedConcreteLayerMembraneSection02::getInitialTangent(void)
{
	// Get the initial membrane section stiffness
	double initialSectionTangent[3][3];			// membrane section stiffness: [Dm] = h*[MaterialTangent] is the membrane stiffness Tangent = Sum(Tangent_ic(zTop_ic - zBottom_ic)) + Sum(Tangent_is(t_is))
	initialSectionTangent[0][0] = 0.0; initialSectionTangent[0][1] = 0.0;	initialSectionTangent[0][2] = 0.0;
	initialSectionTangent[1][0] = 0.0; initialSectionTangent[1][1] = 0.0;	initialSectionTangent[1][2] = 0.0;
	initialSectionTangent[2][0] = 0.0; initialSectionTangent[2][1] = 0.0;	initialSectionTangent[2][2] = 0.0;

	InitialTangent.Zero();

	const Matrix& RCTangent = TheRCMaterial->getInitialTangent();

	double Eco00 = RCTangent(0, 0); double Eco01 = RCTangent(0, 1); double Eco02 = RCTangent(0, 2);
	double Eco10 = RCTangent(1, 0); double Eco11 = RCTangent(1, 1); double Eco12 = RCTangent(1, 2);
	double Eco20 = RCTangent(2, 0); double Eco21 = RCTangent(2, 1); double Eco22 = RCTangent(2, 2);

	initialSectionTangent[0][0] = t * Eco00; initialSectionTangent[0][1] = t * Eco01;	initialSectionTangent[0][2] = t * Eco02;
	initialSectionTangent[1][0] = t * Eco10; initialSectionTangent[1][1] = t * Eco11;	initialSectionTangent[1][2] = t * Eco12;
	initialSectionTangent[2][0] = t * Eco20; initialSectionTangent[2][1] = t * Eco21;	initialSectionTangent[2][2] = t * Eco22;
	
	InitialTangent(0, 0) = initialSectionTangent[0][0]; InitialTangent(0, 1) = initialSectionTangent[0][1]; InitialTangent(0, 2) = initialSectionTangent[0][2];
	InitialTangent(1, 0) = initialSectionTangent[1][0]; InitialTangent(1, 1) = initialSectionTangent[1][1]; InitialTangent(1, 2) = initialSectionTangent[1][2];
	InitialTangent(2, 0) = initialSectionTangent[2][0]; InitialTangent(2, 1) = initialSectionTangent[2][1]; InitialTangent(2, 2) = initialSectionTangent[2][2];

	return InitialTangent;
}

SectionForceDeformation* ReinforcedConcreteLayerMembraneSection02::getCopy(void)
{
	ReinforcedConcreteLayerMembraneSection02* theCopy = new ReinforcedConcreteLayerMembraneSection02(this->getTag(),
		TheRCMaterial,
		t);

	return theCopy;
}

const ID& ReinforcedConcreteLayerMembraneSection02::getType(void)
{
	static bool initialized = false;
	if (!initialized) {
		array(0) = SECTION_RESPONSE_FXX;
		array(1) = SECTION_RESPONSE_FYY;
		array(2) = SECTION_RESPONSE_FXY;
		array(3) = SECTION_RESPONSE_MXX;
		array(4) = SECTION_RESPONSE_MYY;
		array(5) = SECTION_RESPONSE_MXY;
		array(6) = SECTION_RESPONSE_VXZ;
		array(7) = SECTION_RESPONSE_VYZ;
		initialized = true;
	}
	return array;
}

int ReinforcedConcreteLayerMembraneSection02::getOrder(void) const
{
	return 3;
}

int ReinforcedConcreteLayerMembraneSection02::commitState(void)
{
	int success = 0;

	success = TheRCMaterial->commitState();

	// Commit the history variables
	CSectionStrain = TSectionStrain;
	CSectionStress = TSectionStress;
	CSectionTangent = TSectionTangent;

	return success;
}

int ReinforcedConcreteLayerMembraneSection02::revertToLastCommit(void)
{
	int success = 0;

	success = TheRCMaterial->revertToLastCommit();

	// Revert the history variables to last commit
	TSectionStrain = CSectionStrain;
	TSectionStress = CSectionStress;
	TSectionTangent = TSectionTangent;

	return success;
}

int ReinforcedConcreteLayerMembraneSection02::revertToStart(void)
{
	int success = 0;

	success = TheRCMaterial->revertToStart();

	// Revert the state variables to start
	//// Committed state variables
	CSectionStrain.Zero();
	CSectionStress.Zero();
	CSectionTangent = this->getInitialTangent();
	//// Trial state variables
	TSectionStrain.Zero();
	TSectionStress.Zero();
	TSectionTangent = CSectionTangent;

	return success;
}

int ReinforcedConcreteLayerMembraneSection02::sendSelf(int commitTag, Channel& theChannel)
{
	int res = 0;
	int dataTag = this->getDbTag();

	static Vector data(2);
	data(0) = this->getTag();
	data(1) = t;

	res += theChannel.sendVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING ReinforcedConcreteLayerMembraneSection02::sendself() - " << this->getTag() << " failed to send Vector\n";
		return res;
	}

	int matDbTag;
	static ID idData(2);
	idData(0) = TheRCMaterial->getClassTag();
	matDbTag = TheRCMaterial->getDbTag();
	if (matDbTag == 0) {
		matDbTag = theChannel.getDbTag();
		if (matDbTag != 0)
			TheRCMaterial->setDbTag(matDbTag);
	}
	idData(1) = matDbTag;

	res += theChannel.sendID(dataTag, commitTag, idData);
	if (res < 0) {
		opserr << "WARNING ReinforcedConcreteLayerMembraneSection02::sendself() - " << this->getTag() << " failed to send ID\n";
		return res;
	}

	res += TheRCMaterial->sendSelf(commitTag, theChannel);
	if (res < 0) {
		opserr << "WARNING ReinforcedConcreteLayerMembraneSection02::sendself() - " << this->getTag() << " failed to send its Material\n";
		return res;
	}

	return res;
}

int ReinforcedConcreteLayerMembraneSection02::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	int res = 0;

	int dataTag = this->getDbTag();

	static Vector data(2);
	res += theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING ReinforcedConcreteLayerMembraneSection02::recvSelf() - failed to receive Vector\n";
		return res;
	}
	this->setTag((int)data(0));
	t = data(1);

	static ID idData(2);

	res += theChannel.recvID(dataTag, commitTag, idData);
	if (res < 0) {
		opserr << "WARNING ReinforcedConcreteLayerMembraneSection02::recvSelf() - "<<this->getTag()<<" failed to receive ID\n";
		return res;
	}

	int matClassTag = idData(0);
	if (TheRCMaterial != 0) {
		if (matClassTag != TheRCMaterial->getClassTag()) {
			delete TheRCMaterial;
			TheRCMaterial = 0;
		}
	}

	if (TheRCMaterial == 0) {
		TheRCMaterial = theBroker.getNewNDMaterial(matClassTag);
		if (TheRCMaterial == 0) {
			opserr << "ReinforcedConcreteLayerMembraneSection02::recvSelf() - failed to get a NDMaterial of type " << matClassTag << endln;
			return -1;
		}
	}

	TheRCMaterial->setDbTag(idData(1));

	res += TheRCMaterial->recvSelf(commitTag, theChannel, theBroker);
	if (res < 0) {
		opserr << "ReinforcedConcreteLayerMembraneSection02::recvSelf() - the material failed in recvSelf()\n";
		return res;
	}

	return res;
}

void ReinforcedConcreteLayerMembraneSection02::Print(OPS_Stream& s, int flag)
{
	s << "ReinforcedConcreteLayerMembraneSection02 tag: " << this->getTag() << endln;
	s << "Thickness t = " << t << endln;
	
	s << "Reinforced Concrete section: " << endln;
	TheRCMaterial->Print(s, flag);
	s << endln;

}

Response* ReinforcedConcreteLayerMembraneSection02::setResponse(const char** argv, int argc, OPS_Stream& s)
{
	Response* theResponse = 0;
	if (strcmp(argv[0], "panel_strain") == 0 || strcmp(argv[0], "Panel_strain") == 0) {
		s.tag("SectionOutput");
		s.attr("secType", "ReinforcedConcreteLayerMembraneSection02");
		s.attr("secTag", this->getTag());
		s.tag("ResponseType", "eps11");
		s.tag("ResponseType", "eps22");
		s.tag("ResponseType", "eps12");
		s.endTag();

		Vector data1(3);
		data1.Zero();

		theResponse = new MaterialResponse(this, 1,data1);

	}
	else if (strcmp(argv[0], "panel_stress") == 0 || strcmp(argv[0], "Panel_Stress") == 0) {
		s.tag("SectionOutput");
		s.attr("secType", "ReinforcedConcreteLayerMembraneSection02");
		s.attr("secTag", this->getTag());
		s.tag("ResponseType", "sigma11");
		s.tag("ResponseType", "sigma22");
		s.tag("ResponseType", "sigma12");
		s.endTag();

		Vector data2(3);
		data2.Zero();

		theResponse = new MaterialResponse(this,2,data2);

	}
	else if (strcmp(argv[0], "getBendingParameters") == 0) {
		s.tag("SectionOutput");
		s.attr("secType", "ReinforcedConcreteLayerMembraneSection02");
		s.attr("secTag", this->getTag());
		s.tag("ResponseType", "Eave");
		s.tag("ResponseType", "Tave");
		s.endTag();

		Vector data3(2);
		data3.Zero();

		theResponse = new MaterialResponse(this, 3, data3);
	}
	else {

		return this->SectionForceDeformation::setResponse(argv, argc, s);
	}

	return theResponse;
}


int ReinforcedConcreteLayerMembraneSection02::getResponse(int responseID, Information& info)
{
	if (responseID == 1) {
		return info.setVector(this->getCommittedStrain());
	}
	else if (responseID == 2) {
		return info.setVector(this->getCommittedStress());
	}
	else if (responseID == 3) {
		return info.setVector(this->getBendingParameters());
	}
	else {
		return 0;
	}
}

// Function that returns bending parameters - added for MEFI3D by Maria Jose Nunez, UChile
Vector ReinforcedConcreteLayerMembraneSection02::getBendingParameters(void)
{
	Vector input_par(2);

	input_par.Zero();

	input_par(0) = this->getEcAvg();
	input_par(1) = t;

	return input_par;
}