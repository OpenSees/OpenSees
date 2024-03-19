// Code written/implemented by: Maria Jose Nunez
//
// Created: 03/2024
// 
// Description: This file contains the ElasticMembraneSection class definition
// A ElasticMembraneSection is a subclass of the sectionForceDeformation class and corresponds to the abstract representation
// for the stress-strain behavior for a elastic membrane section in the Finite Element Method or Structural Analysis. 
//
// Reference:
// 1. Rojas, F., Anderson, J. C., Massone, L. M. (2016). A nonlinear quadrilateral layered membrane element with drilling degrees of freedom for 
// the modeling of reinforced concrete walls. Engineering Structures, 124, 521-538.
//
// Source: \OpenSees\SRC\material\section\LayeredMembraneSection
//
// Rev: 1.0

#include <ElasticMembraneSection.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <Parameter.h>

using namespace std;

// Read input parameters and build the section
void* OPS_ElasticMembraneSection()
{
	int numArgs = OPS_GetNumRemainingInputArgs();

	// Parse the script for material parameters
	if (numArgs < 4) {
		opserr << "Want: ElasticMembraneSection $secTag $E $nu $t <$rho>" << endln;
		return 0;
	}

	int tag;						// section tag
	double dData[4];				// # of material parameters
	dData[3] = 0.0;

	// section tag
	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) {
		opserr << "WARNING invalid section ElasticMembraneSection tag" << endln;
		return 0;
	}
	
	// section parameters
	if (numArgs > 4)
		numdata = 4;
	else
		numdata = 3;

	if (OPS_GetDouble(&numdata, dData) < 0) {
		opserr << "Invalid section parameters for ElasticMembraneSection " << tag << endln;
		return 0;
	}

	SectionForceDeformation* theSection = new ElasticMembraneSection(tag, dData[0], dData[1], dData[2], dData[3]);

	return theSection;
}

//static vector and matrices
ID ElasticMembraneSection::array(3);

// Full constructor
ElasticMembraneSection::ElasticMembraneSection(int tag,		// section tag
	double elasticModulus,									// Young's Modulus
	double poisson,											// Poisson ratio
	double thickness,										// section thickness
	double r)							                    // mass density			

	:SectionForceDeformation(tag, SEC_TAG_ElasticMembraneSection),
	E(elasticModulus), nu(poisson), t(thickness), rho(r*thickness),
	CSectionStrain(3), CSectionStress(3), CSectionTangent(3,3),
	TSectionStrain(3), TSectionStress(3), TSectionTangent(3,3), InitialTangent(3,3)
{
	// Set initial values
	for (int i = 0; i < 3; i++) {
		CSectionStrain(i) = 0.0;
		CSectionStress(i) = 0.0;
		TSectionStrain(i) = 0.0;
		TSectionStress(i) = 0.0;
	}
}

// Blank constructor (constructor for blank object that recvSelf needs to be invoked upon) (constructor which should be invoked by an FEM_ObjectBroker only)
ElasticMembraneSection::ElasticMembraneSection() :
	SectionForceDeformation(0, SEC_TAG_ElasticMembraneSection),
	E(0), nu(0), t(0), rho(0),
	CSectionStrain(3), CSectionStress(3), CSectionTangent(3, 3),
	TSectionStrain(3), TSectionStress(3), TSectionTangent(3, 3), InitialTangent(3, 3)
{
	this->revertToStart();
}

// Destructor (clean up memory the ElasticMembraneSection objects allocates)
ElasticMembraneSection::~ElasticMembraneSection()
{

}

int ElasticMembraneSection::setTrialSectionDeformation(const Vector& newTrialSectionStrain)
{
	// Set the membrane strain (e11,e22,e12) in the section
	TSectionStrain(0) = newTrialSectionStrain(0);
	TSectionStrain(1) = newTrialSectionStrain(1);
	TSectionStrain(2) = newTrialSectionStrain(2);

	double sectionResultantStress[3];		// section resultant stress: {N11 ; N22 ; N12 }
	sectionResultantStress[0] = 0.0; sectionResultantStress[1] = 0.0; sectionResultantStress[2] = 0.0;

	double sectionTangent[3][3];			// membrane section stiffness: [Dm] = h*[MaterialTangent] is the membrane stiffness Tangent
	sectionTangent[0][0] = 0.0; sectionTangent[0][1] = 0.0;	sectionTangent[0][2] = 0.0;
	sectionTangent[1][0] = 0.0; sectionTangent[1][1] = 0.0;	sectionTangent[1][2] = 0.0;
	sectionTangent[2][0] = 0.0; sectionTangent[2][1] = 0.0;	sectionTangent[2][2] = 0.0;

	// elastic material
	double d00 = E / (1.0 - nu * nu);
	double d01 = nu * d00;
	double d22 = 0.5 * (d00 - d01);

	double fc11 = d00 * TSectionStrain(0) + d01 * TSectionStrain(1);
	double fc22 = d01 * TSectionStrain(0) + d00 * TSectionStrain(1);
	double fc12 = d22 * TSectionStrain(2);

	sectionResultantStress[0] = t * fc11;
	sectionResultantStress[1] = t * fc22;
	sectionResultantStress[2] = t * fc12;

	double Ec00 = d00; double Ec01 = d01; double Ec02 = 0.0;
	double Ec10 = d01; double Ec11 = d00; double Ec12 = 0.0;
	double Ec20 = 0.0; double Ec21 = 0.0; double Ec22 = d22;

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

const Vector& ElasticMembraneSection::getSectionDeformation(void)
{
	return TSectionStrain;
}

const Vector& ElasticMembraneSection::getStressResultant(void)
{
	return TSectionStress;
}

const Matrix& ElasticMembraneSection::getSectionTangent(void)
{
	return TSectionTangent;
}

const Vector& ElasticMembraneSection::getCommittedStrain(void)
{
	return CSectionStrain;
}

const Vector& ElasticMembraneSection::getCommittedStress(void)
{
	return CSectionStress;
}

// density per unit area
double ElasticMembraneSection::getRho(void)
{
	return rho;
}

const Matrix& ElasticMembraneSection::getInitialTangent(void)
{
	// Get the initial membrane section stiffness
	double initialSectionTangent[3][3];			// membrane section stiffness: [Dm] = h*[MaterialTangent] is the membrane stiffness Tangent 
	initialSectionTangent[0][0] = 0.0; initialSectionTangent[0][1] = 0.0;	initialSectionTangent[0][2] = 0.0;
	initialSectionTangent[1][0] = 0.0; initialSectionTangent[1][1] = 0.0;	initialSectionTangent[1][2] = 0.0;
	initialSectionTangent[2][0] = 0.0; initialSectionTangent[2][1] = 0.0;	initialSectionTangent[2][2] = 0.0;

	InitialTangent.Zero();

	double d00 = E / (1.0 - nu * nu);
	double d01 = nu * d00;
	double d22 = 0.5 * (d00 - d01);

	double Eco00 = d00; double Eco01 = d01; double Eco02 = 0.0;
	double Eco10 = d01; double Eco11 = d00; double Eco12 = 0.0;
	double Eco20 = 0.0; double Eco21 = 0.0; double Eco22 = d22;

	initialSectionTangent[0][0] = t * Eco00; initialSectionTangent[0][1] = t * Eco01;	initialSectionTangent[0][2] = t * Eco02;
	initialSectionTangent[1][0] = t * Eco10; initialSectionTangent[1][1] = t * Eco11;	initialSectionTangent[1][2] = t * Eco12;
	initialSectionTangent[2][0] = t * Eco20; initialSectionTangent[2][1] = t * Eco21;	initialSectionTangent[2][2] = t * Eco22;
	
	InitialTangent(0, 0) = initialSectionTangent[0][0]; InitialTangent(0, 1) = initialSectionTangent[0][1]; InitialTangent(0, 2) = initialSectionTangent[0][2];
	InitialTangent(1, 0) = initialSectionTangent[1][0]; InitialTangent(1, 1) = initialSectionTangent[1][1]; InitialTangent(1, 2) = initialSectionTangent[1][2];
	InitialTangent(2, 0) = initialSectionTangent[2][0]; InitialTangent(2, 1) = initialSectionTangent[2][1]; InitialTangent(2, 2) = initialSectionTangent[2][2];

	return InitialTangent;
}

SectionForceDeformation* ElasticMembraneSection::getCopy(void)
{
	ElasticMembraneSection* theCopy = new ElasticMembraneSection(this->getTag(),
		E,
		nu,
		t,
		rho);

	return theCopy;
}

const ID& ElasticMembraneSection::getType(void)
{
	static bool initialized = false;
	if (!initialized) {
		array(0) = SECTION_RESPONSE_FXX;
		array(1) = SECTION_RESPONSE_FYY;
		array(2) = SECTION_RESPONSE_FXY;
		initialized = true;
	}
	return array;
}

int ElasticMembraneSection::getOrder(void) const
{
	return 3;
}

int ElasticMembraneSection::commitState(void)
{
	int success = 0;

	// Commit the history variables
	CSectionStrain = TSectionStrain;
	CSectionStress = TSectionStress;
	CSectionTangent = TSectionTangent;

	return success;
}

int ElasticMembraneSection::revertToLastCommit(void)
{
	int success = 0;

	// Revert the history variables to last commit
	TSectionStrain = CSectionStrain;
	TSectionStress = CSectionStress;
	TSectionTangent = TSectionTangent;

	return success;
}

int ElasticMembraneSection::revertToStart(void)
{
	int success = 0;

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

int ElasticMembraneSection::sendSelf(int commitTag, Channel& theChannel)
{
	int res = 0;
	int dataTag = this->getDbTag();

	static Vector data(5);
	data(0) = this->getTag();
	data(1) = E;
	data(2) = nu;
	data(3) = t;
	data(4) = rho;

	res += theChannel.sendVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING ElasticMembraneSection::sendself() - " << this->getTag() << " failed to send Vector\n";
		return res;
	}

	return res;
}

int ElasticMembraneSection::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	int res = 0;
	int dataTag = this->getDbTag();

	static Vector data(5);

	res += theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING ElasticMembraneSection::recvSelf() - failed to receive Vector\n";
		return res;
	}

	this->setTag((int)data(0));
	E   = data(1);
	nu  = data(2);
	t   = data(3);
	rho = data(4);

	return res;
}

void ElasticMembraneSection::Print(OPS_Stream& s, int flag)
{
	s << "ElasticMembraneSection tag: " << this->getTag() << endln;
	s << "Young's Modulus E = " << E << endln;
	s << "Poisson's Ratio v = " << nu << endln;
	s << "Thickness t = " << t << endln;
	s << "Density rho = " << rho << endln;
}