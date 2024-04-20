// Code written/implemented by: Fabian Rojas B.
//								Maria Jose Nunez
//
// Created: 04/2023
// 
// Description: This file contains the ReinforcedConcreteLayeredMembraneSection class definition
// A ReinforcedConcreteLayeredMembraneSection is a subclass of the sectionForceDeformation class and corresponds to the abstract representation
// for the stress-strain behavior for a reinforced concrete layered membrane element in the Finite Element Method or Structural Analysis. 
//
// Reference:
// 1. Rojas, F., Anderson, J. C., Massone, L. M. (2016). A nonlinear quadrilateral layered membrane element with drilling degrees of freedom for 
// the modeling of reinforced concrete walls. Engineering Structures, 124, 521-538.
//
// Source: \OpenSees\SRC\material\section\LayeredMembraneSection
//
// Rev: 1.0

#include <ReinforcedConcreteLayeredMembraneSection.h>
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
void* OPS_ReinforcedConcreteLayeredMembraneSection()
{
	int numArgs = OPS_GetNumRemainingInputArgs();

	// Parse the script for material parameters
	if (numArgs < 9) {
		opserr << "Want: ReinforcedConcreteLayeredMembraneSection $secTag $nSteelLayers $nConcLayers -reinfSteel {$RSteelAtEachLayer} -conc {$concAtEachLayer} -concThick {$concThicknessAtEachLayer}" << endln;
		return 0;
	}

	int tag, nClayers, nSlayers;
	int* concMatTags, * steelMatTags;
	double* cThickness;
	const char* str;
	NDMaterial** theConcMats;
	NDMaterial** theSteelMats;

	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) {
		opserr << "WARNING invalid section ReinforcedConcreteLayeredMembraneSection tag" << endln;
		return 0;
	}
	
	if (OPS_GetIntInput(&numdata, &nSlayers) < 0) {
		opserr << "WARNING invalid nSlayers" << endln;
		opserr << "ReinforcedConcreteLayeredMembraneSection section: " << tag << endln;
		return 0;
	}
	
	if (OPS_GetIntInput(&numdata, &nClayers) < 0) {
		opserr << "WARNING invalid nClayers" << endln;
		opserr << "ReinforcedConcreteLayeredMembraneSection section: " << tag << endln;
		return 0;
	}

	theConcMats = new NDMaterial * [nClayers];
	theSteelMats = new NDMaterial * [nSlayers];
	cThickness = new double[nClayers];
	concMatTags = new int[nClayers];
	steelMatTags = new int[nSlayers];

	numArgs = OPS_GetNumRemainingInputArgs();

	while (numArgs > 0) {
		str = OPS_GetString();
		if (strcmp(str, "-reinfSteel") == 0) {
			numdata = nSlayers;
			if (OPS_GetIntInput(&numdata, steelMatTags) != 0) {
				opserr << "Invalid nDMaterial reinforced steel tag for ReinforcedConcreteLayeredMembraneSection " << tag << endln;
				return 0;
			}
			for (int i = 0; i < nSlayers; i++) {
				theSteelMats[i] = 0;
				theSteelMats[i] = OPS_getNDMaterial(steelMatTags[i]);
				if (theSteelMats[i] == 0) {
					opserr << "Invalid nDMaterial reinforced steel tag " << steelMatTags[i] << " for ReinforcedConcreteLayeredMembraneSection " << tag << endln;
					return 0;
				}
			}
		}
		else if (strcmp(str, "-conc") == 0) {
			numdata = nClayers;
			if (OPS_GetIntInput(&numdata, concMatTags) != 0) {
				opserr << "Invalid nDMaterial concrete tag for ReinforcedConcreteLayeredMembraneSection " << tag << endln;
				return 0;
			}
			for (int i = 0; i < nClayers; i++) {
				theConcMats[i] = 0;
				theConcMats[i] = OPS_getNDMaterial(concMatTags[i]);
				if (theConcMats[i] == 0) {
					opserr << "Invalid nDMaterial concrete tag " << concMatTags[i] << " for ReinforcedConcreteLayeredMembraneSection " << tag << endln;
					return 0;
				}
			}
		}
		else if (strcmp(str, "-concThick") == 0) {
			numdata = nClayers;
			if (OPS_GetDoubleInput(&numdata, cThickness) != 0) {
				opserr << "Invalid thick parameter for ReinforcedConcreteLayeredMembraneSection " << tag << endln;
				return 0;
			}
		}
		else {
			opserr << "WARNING: Invalid option " << str << " in ReinforcedConcreteLayeredMembraneSection " << tag << endln;
			return 0;
		}
		
		numArgs = OPS_GetNumRemainingInputArgs();
    
	}

	SectionForceDeformation* theSection = new ReinforcedConcreteLayeredMembraneSection(tag, nSlayers, nClayers, theSteelMats, theConcMats, cThickness);

	// Clean up dynamic memory
	if (theSteelMats != 0)
		delete[] theSteelMats;
	if (theConcMats != 0)
		delete[] theConcMats;
	if (steelMatTags != 0)
		delete[] steelMatTags;
	if (concMatTags != 0)
		delete[] concMatTags;
	if (cThickness != 0)
		delete[] cThickness;

	return theSection;
}

//static vector and matrices
ID ReinforcedConcreteLayeredMembraneSection::array(3);

// Full constructor
ReinforcedConcreteLayeredMembraneSection::ReinforcedConcreteLayeredMembraneSection(int tag,							// section tag
	int nSteelLayer,																								// number of reinforced steel layers
	int nConcLayer,																									// number of concrete layers
	NDMaterial** reinforcedSteelMaterialObjects,																	// array of nDMaterial reinforced steel tags for each layer
	NDMaterial** concrete2DMaterialObjects,																			// array of nDMaterial concrete tags for each layer
	double* concThickness)																							// array of concrete layers thickness

	:SectionForceDeformation(tag, SEC_TAG_ReinforcedConcreteLayeredMembraneSection),
	numberReinforcedSteelLayers(nSteelLayer), numberConcreteLayers(nConcLayer), 
	CSectionStrain(3), CSectionStress(3), CSectionTangent(3, 3),
	TSectionStrain(3), TSectionStress(3), TSectionTangent(3, 3), InitialTangent(3, 3),
	strainPrincipalDirection(3), poissonRatios(2), crackPattern(4), pi(3.1415926535)
{

	// Set initial values
	thetaPrincipalDirection = 0.0;
	strainPrincipalDirection(0) = 0.0;
	strainPrincipalDirection(1) = 0.0;
	strainPrincipalDirection(2) = 0.0;
	isConcreteCracked = false;

	poissonRatios(0) = 0.0;
	poissonRatios(1) = 0.0;

	ecr = 0.0;
	ec = 0.0;

	crackPattern.Zero();

	for (int i = 0; i < 3; i++) {
		CSectionStrain(i) = 0.0;
		CSectionStress(i) = 0.0;
		TSectionStrain(i) = 0.0;
		TSectionStress(i) = 0.0;
	}

	// Check concrete layers thickness input
	if (concThickness == 0) {
		opserr << "ReinforcedConcreteLayeredMembraneSection::ReinforcedConcreteLayeredMembraneSection - Null concrete layers thickness array passed.\n";
		exit(-1);
	}

	// Check materials input
	if (reinforcedSteelMaterialObjects == 0) {
		opserr << "ReinforcedConcreteLayeredMembraneSection::ReinforcedConcreteLayeredMembraneSection() - Null nDMaterial reinforced steel array passed.\n";
		exit(-1);
	}
	if (concrete2DMaterialObjects == 0) {
		opserr <<"ReinforcedConcreteLayeredMembraneSection::ReinforcedConcreteLayeredMembraneSection() - Null nDMaterial concrete array passed.\n";
		exit(-1);
	}

	// Allocate memory for concrete layers thickness
	t = new double[numberConcreteLayers];
	h = 0.0;

	for (int i = 0; i < numberConcreteLayers; i++) {
		t[i] = concThickness[i];
		h += t[i];		// Total thickness of the section
	}

	// Allocate memory for the NDMaterials
	TheReinforcedSteel2DMaterial = new NDMaterial * [numberReinforcedSteelLayers];
	TheConcrete2DMaterial = new NDMaterial * [numberConcreteLayers];

	if (TheReinforcedSteel2DMaterial == 0) {
		opserr << "ReinforcedConcreteLayeredMembraneSection::ReinforcedConcreteLayeredMembraneSection() - Failed to allocate pointers for ND reinforced steel materials.\n";
		exit(-1);
	}
	if (TheConcrete2DMaterial == 0) {
		opserr << "ReinforcedConcreteLayeredMembraneSection::ReinforcedConcreteLayeredMembraneSection() - Failed to allocate pointers for ND concrete materials.\n";
		exit(-1);
	}

	// Get copies of the ND reinforced steel materials 
	for (int i = 0; i < numberReinforcedSteelLayers; i++) {
		if (reinforcedSteelMaterialObjects[i] == 0) {
			opserr << "ReinforcedConcreteLayeredMembraneSection::ReinforcedConcreteLayeredMembraneSection() - Null ND reinforced steel material pointer passed.\n";
			exit(-1);
		}

		TheReinforcedSteel2DMaterial[i] = reinforcedSteelMaterialObjects[i]->getCopy("PlaneStress2D");

		if (TheReinforcedSteel2DMaterial[i] == 0) {
			opserr << "ReinforcedConcreteLayeredMembraneSection::ReinforcedConcreteLayeredMembraneSection() - Failed to copy ND reinforced steel material.\n";
			exit(-1);
		}
	}

	// Get copies of the ND concrete materials 
	for (int i = 0; i < numberConcreteLayers; i++) {
		if (concrete2DMaterialObjects[i] == 0) {
			opserr << "ReinforcedConcreteLayeredMembraneSection::ReinforcedConcreteLayeredMembraneSection() - Null ND concrete material pointer passed.\n";
			exit(-1);
		}

		TheConcrete2DMaterial[i] = concrete2DMaterialObjects[i]->getCopy("PlaneStress2D");

		if (TheConcrete2DMaterial[i] == 0) {
			opserr << "ReinforcedConcreteLayeredMembraneSection::ReinforcedConcreteLayeredMembraneSection() - Failed to copy ND concrete material.\n";
			exit(-1);
		}
	}

}

// Blank constructor (constructor for blank object that recvSelf needs to be invoked upon) (constructor which should be invoked by an FEM_ObjectBroker only)
	ReinforcedConcreteLayeredMembraneSection::ReinforcedConcreteLayeredMembraneSection() :
	SectionForceDeformation(0, SEC_TAG_ReinforcedConcreteLayeredMembraneSection),
	numberReinforcedSteelLayers(0), numberConcreteLayers(0),
	CSectionStrain(3), CSectionStress(3), CSectionTangent(3, 3),
	TSectionStrain(3), TSectionStress(3), TSectionTangent(3, 3), InitialTangent(3, 3),
	strainPrincipalDirection(3), poissonRatios(2), pi(3.1415926535)
{
	TheConcrete2DMaterial = 0;
	TheReinforcedSteel2DMaterial = 0;

	this->revertToStart();
}

// Destructor (clean up memory the ReinforcedConcreteLayeredMembraneSection objects allocates)
	ReinforcedConcreteLayeredMembraneSection::~ReinforcedConcreteLayeredMembraneSection()
{
	if (t != 0) delete []t;

	if (TheReinforcedSteel2DMaterial != 0) {
		for (int i = 0; i < numberReinforcedSteelLayers; i++)
			if (TheReinforcedSteel2DMaterial[i] != 0)
				delete TheReinforcedSteel2DMaterial[i];
		delete[] TheReinforcedSteel2DMaterial;
	}

	if (TheConcrete2DMaterial != 0) {
		for (int i = 0; i < numberConcreteLayers; i++)
			if (TheConcrete2DMaterial[i] != 0)
				delete TheConcrete2DMaterial[i];
		delete[] TheConcrete2DMaterial;
	}
}

int ReinforcedConcreteLayeredMembraneSection::setTrialSectionDeformation(const Vector& newTrialSectionStrain)
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

	// Set the strain at each concrete layer
	for (int ic = 0; ic < numberConcreteLayers; ic++) {
		// Get concrete stress and tangent stiffness
		TheConcrete2DMaterial[ic]->setTrialStrain(newTrialSectionStrain);
		const Vector& stressNDM = TheConcrete2DMaterial[ic]->getStress();
		const Matrix& tangentNDM = TheConcrete2DMaterial[ic]->getTangent();

		double fc11 = stressNDM(0);
		double fc22 = stressNDM(1);
		double fc12 = stressNDM(2);
		// Calculate the concrete contribution to section resultant stress
		sectionResultantStress[0] += t[ic] * fc11;
		sectionResultantStress[1] += t[ic] * fc22;
		sectionResultantStress[2] += t[ic] * fc12;

		double Ec00 = tangentNDM(0, 0); double Ec01 = tangentNDM(0, 1); double Ec02 = tangentNDM(0, 2);
		double Ec10 = tangentNDM(1, 0); double Ec11 = tangentNDM(1, 1); double Ec12 = tangentNDM(1, 2);
		double Ec20 = tangentNDM(2, 0); double Ec21 = tangentNDM(2, 1); double Ec22 = tangentNDM(2, 2);
		// Calculate the concrete contribution to section tangent
		sectionTangent[0][0] += t[ic] * Ec00; sectionTangent[0][1] += t[ic] * Ec01;	sectionTangent[0][2] += t[ic] * Ec02;
		sectionTangent[1][0] += t[ic] * Ec10; sectionTangent[1][1] += t[ic] * Ec11;	sectionTangent[1][2] += t[ic] * Ec12;
		sectionTangent[2][0] += t[ic] * Ec20; sectionTangent[2][1] += t[ic] * Ec21;	sectionTangent[2][2] += t[ic] * Ec22;
	}

	// Set the strain at each reinforced steel layer
	for (int iRs = 0; iRs < numberReinforcedSteelLayers; iRs++) {
		// Get reinforced steel stress and tangent stiffness
		TheReinforcedSteel2DMaterial[iRs]->setTrialStrain(newTrialSectionStrain);
		const Vector& stressNDM = TheReinforcedSteel2DMaterial[iRs]->getStress();
		const Matrix& tangentNDM = TheReinforcedSteel2DMaterial[iRs]->getTangent();

		double fs11 = stressNDM(0);
		double fs22 = stressNDM(1);
		double fs12 = stressNDM(2);
		// Calculate the reinforced steel contribution to section resultant stress
		sectionResultantStress[0] += h * fs11;
		sectionResultantStress[1] += h * fs22;
		sectionResultantStress[2] += h * fs12;

		double Es00 = tangentNDM(0, 0); double Es01 = tangentNDM(0, 1); double Es02 = tangentNDM(0, 2);
		double Es10 = tangentNDM(1, 0); double Es11 = tangentNDM(1, 1); double Es12 = tangentNDM(1, 2);
		double Es20 = tangentNDM(2, 0); double Es21 = tangentNDM(2, 1); double Es22 = tangentNDM(2, 2);
		// Calculate the reinforced steel contribution to section tangent
		sectionTangent[0][0] += h * Es00; sectionTangent[0][1] += h * Es01; sectionTangent[0][2] += h * Es02;
		sectionTangent[1][0] += h * Es10; sectionTangent[1][1] += h * Es11; sectionTangent[1][2] += h * Es12;
		sectionTangent[2][0] += h * Es20; sectionTangent[2][1] += h * Es21; sectionTangent[2][2] += h * Es22;
	}
	// Set the stress and tangent
	TSectionStress(0) = sectionResultantStress[0]; TSectionStress(1) = sectionResultantStress[1]; TSectionStress(2) = sectionResultantStress[2];

	TSectionTangent(0, 0) = sectionTangent[0][0]; TSectionTangent(0, 1) = sectionTangent[0][1]; TSectionTangent(0, 2) = sectionTangent[0][2];
	TSectionTangent(1, 0) = sectionTangent[1][0]; TSectionTangent(1, 1) = sectionTangent[1][1]; TSectionTangent(1, 2) = sectionTangent[1][2];
	TSectionTangent(2, 0) = sectionTangent[2][0]; TSectionTangent(2, 1) = sectionTangent[2][1]; TSectionTangent(2, 2) = sectionTangent[2][2];

	return 0;
}

const Vector& ReinforcedConcreteLayeredMembraneSection::getSectionDeformation(void)
{
	return TSectionStrain;
}

const Vector& ReinforcedConcreteLayeredMembraneSection::getStressResultant(void)
{
	// newName = TSectionStress/espesor_total
	return TSectionStress;
}

const Matrix& ReinforcedConcreteLayeredMembraneSection::getSectionTangent(void)
{
	return TSectionTangent;
}

const Vector& ReinforcedConcreteLayeredMembraneSection::getCommittedStrain(void)
{
	return CSectionStrain;
}

const Vector& ReinforcedConcreteLayeredMembraneSection::getCommittedStress(void)
{
	return CSectionStress;
}

double ReinforcedConcreteLayeredMembraneSection::getRho(void)
{
	double rhoH = 0.0;

	for (int ic = 0; ic < numberConcreteLayers; ic++) {
		rhoH += (TheConcrete2DMaterial[ic]->getRho()) * t[ic];
	}

	return rhoH;
}

double ReinforcedConcreteLayeredMembraneSection::getEcAvg(void)
{
	DummyStream theDummyStream;

	char aa[80] = "getParameters";
	const char* argv[1];
	argv[0] = aa;

	double EcAvg = 0.0;

	for (int ic = 0; ic < numberConcreteLayers; ic++) {
		
		Response* theResponse = TheConcrete2DMaterial[ic]->setResponse(argv, 1, theDummyStream);

		if (theResponse == 0) {
			opserr << "ReinforcedConcreteLayeredMembraneSection::ReinforcedConcreteLayeredMembraneSection - failed to get concrete parameters for ReinforcedConcreteLayeredMembraneSection with tag: " << this->getTag() << "\n";
			exit(-1);
		}

		// Get OrthotropicRotatingAngleConcreteT2DMaterial01 concrete parameters
		theResponse->getResponse();
		Information& theInfoPar = theResponse->getInformation();
		const Vector& ParNDConc = theInfoPar.getData();

		// Get concrete young's modulus and multiplied by the thickness of each layer
		EcAvg += (ParNDConc[1]) * t[ic];

		delete theResponse;
	}

	EcAvg = EcAvg / h;

	return EcAvg;
}

const Matrix& ReinforcedConcreteLayeredMembraneSection::getInitialTangent(void)
{
	// Get the initial membrane section stiffness
	double initialSectionTangent[3][3];			// membrane section stiffness: [Dm] = h*[MaterialTangent] is the membrane stiffness Tangent = Sum(Tangent_ic(zTop_ic - zBottom_ic)) + Sum(Tangent_is(t_is))
	initialSectionTangent[0][0] = 0.0; initialSectionTangent[0][1] = 0.0;	initialSectionTangent[0][2] = 0.0;
	initialSectionTangent[1][0] = 0.0; initialSectionTangent[1][1] = 0.0;	initialSectionTangent[1][2] = 0.0;
	initialSectionTangent[2][0] = 0.0; initialSectionTangent[2][1] = 0.0;	initialSectionTangent[2][2] = 0.0;

	InitialTangent.Zero();

	// Add the section stiffness contribution from each concrete layer
	for (int ic = 0; ic < numberConcreteLayers; ic++) {
		// Get concrete initial tangent stiffness
		const Matrix& layerConcreteTangent = TheConcrete2DMaterial[ic]->getInitialTangent();

		double Eco00 = layerConcreteTangent(0, 0); double Eco01 = layerConcreteTangent(0, 1); double Eco02 = layerConcreteTangent(0, 2);
		double Eco10 = layerConcreteTangent(1, 0); double Eco11 = layerConcreteTangent(1, 1); double Eco12 = layerConcreteTangent(1, 2);
		double Eco20 = layerConcreteTangent(2, 0); double Eco21 = layerConcreteTangent(2, 1); double Eco22 = layerConcreteTangent(2, 2);

		// Calculate the concrete contribution to section tangent
		initialSectionTangent[0][0] += t[ic] * Eco00; initialSectionTangent[0][1] += t[ic] * Eco01;	initialSectionTangent[0][2] += t[ic] * Eco02;
		initialSectionTangent[1][0] += t[ic] * Eco10; initialSectionTangent[1][1] += t[ic] * Eco11;	initialSectionTangent[1][2] += t[ic] * Eco12;
		initialSectionTangent[2][0] += t[ic] * Eco20; initialSectionTangent[2][1] += t[ic] * Eco21;	initialSectionTangent[2][2] += t[ic] * Eco22;
	}

	// Add the section stiffness contribution from each reinforced steel layer
	for (int iRs = 0; iRs < numberReinforcedSteelLayers; iRs++) {
		// Get reinforced steel initial tangent stiffness
		const Matrix& layerReinforcedSteelTangent = TheReinforcedSteel2DMaterial[iRs]->getInitialTangent();

		double Eso00 = layerReinforcedSteelTangent(0, 0); double Eso01 = layerReinforcedSteelTangent(0, 1); double Eso02 = layerReinforcedSteelTangent(0, 2);
		double Eso10 = layerReinforcedSteelTangent(1, 0); double Eso11 = layerReinforcedSteelTangent(1, 1); double Eso12 = layerReinforcedSteelTangent(1, 2);
		double Eso20 = layerReinforcedSteelTangent(2, 0); double Eso21 = layerReinforcedSteelTangent(2, 1); double Eso22 = layerReinforcedSteelTangent(2, 2);
		// Calculate the reinforced steel contribution to section tangent
		initialSectionTangent[0][0] += h * Eso00; initialSectionTangent[0][1] += h * Eso01; initialSectionTangent[0][2] += h * Eso02;
		initialSectionTangent[1][0] += h * Eso10; initialSectionTangent[1][1] += h * Eso11; initialSectionTangent[1][2] += h * Eso12;
		initialSectionTangent[2][0] += h * Eso20; initialSectionTangent[2][1] += h * Eso21; initialSectionTangent[2][2] += h * Eso22;
	}
	
	InitialTangent(0, 0) = initialSectionTangent[0][0]; InitialTangent(0, 1) = initialSectionTangent[0][1]; InitialTangent(0, 2) = initialSectionTangent[0][2];
	InitialTangent(1, 0) = initialSectionTangent[1][0]; InitialTangent(1, 1) = initialSectionTangent[1][1]; InitialTangent(1, 2) = initialSectionTangent[1][2];
	InitialTangent(2, 0) = initialSectionTangent[2][0]; InitialTangent(2, 1) = initialSectionTangent[2][1]; InitialTangent(2, 2) = initialSectionTangent[2][2];

	return InitialTangent;
}

SectionForceDeformation* ReinforcedConcreteLayeredMembraneSection::getCopy(void)
{
	ReinforcedConcreteLayeredMembraneSection* theCopy = new ReinforcedConcreteLayeredMembraneSection(this->getTag(),
		numberReinforcedSteelLayers,
		numberConcreteLayers,
		TheReinforcedSteel2DMaterial,
		TheConcrete2DMaterial,
		t);

	return theCopy;
}

const ID& ReinforcedConcreteLayeredMembraneSection::getType(void)
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

int ReinforcedConcreteLayeredMembraneSection::getOrder(void) const
{
	return 3;
}

int ReinforcedConcreteLayeredMembraneSection::commitState(void)
{
	int success = 0;

	// Commit the state for each concrete layer 
	for (int iC = 0; iC < numberConcreteLayers; iC++)
		success += TheConcrete2DMaterial[iC]->commitState();

	// Commit the state for each reinforced steel layer
	for (int iRS = 0; iRS < numberReinforcedSteelLayers; iRS++)
		success += TheReinforcedSteel2DMaterial[iRS]->commitState();

	// Commit the history variables
	CSectionStrain = TSectionStrain;
	CSectionStress = TSectionStress;
	CSectionTangent = TSectionTangent;

	// Set the crack pattern for the section
	this->setCrackPattern();

	return success;
}

int ReinforcedConcreteLayeredMembraneSection::revertToLastCommit(void)
{
	int success = 0;

	// Revert to last commit the state of each concrete layer
	for (int iC = 0; iC < numberConcreteLayers; iC++)
		success += TheConcrete2DMaterial[iC]->revertToLastCommit();

	// Revert to last commit the state of each reinforced steel layer
	for (int iRS = 0; iRS < numberReinforcedSteelLayers; iRS++)
		success += TheReinforcedSteel2DMaterial[iRS]->revertToLastCommit();

	// Revert the history variables to last commit
	TSectionStrain = CSectionStrain;
	TSectionStress = CSectionStress;
	TSectionTangent = CSectionTangent;

	return success;
}

int ReinforcedConcreteLayeredMembraneSection::revertToStart(void)
{
	int success = 0;

	// Revert to start the state of each concrete layer
	for (int iC = 0; iC < numberConcreteLayers; iC++)
		success += TheConcrete2DMaterial[iC]->revertToStart();

	// Revert to start the state of each reinforced steel layer
	for (int iRS = 0; iRS < numberReinforcedSteelLayers; iRS++)
		success += TheReinforcedSteel2DMaterial[iRS]->revertToStart();

	// Revert the state variables to start
	crackPattern.Zero();
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

int ReinforcedConcreteLayeredMembraneSection::sendSelf(int commitTag, Channel& theChannel)
{
	int res = 0;
	int dataTag = this->getDbTag();

	static ID iData(3);
	iData(0) = this->getTag();
	iData(1) = numberReinforcedSteelLayers;
	iData(2) = numberConcreteLayers;

	res += theChannel.sendID(dataTag, commitTag, iData);
	if (res < 0) {
		opserr << "WARNING ReinforcedConcreteLayeredMembraneSection::sendSelf() - " << this->getTag() << " failed to send ID data" << endln;
		return res;
	}

	if (numberReinforcedSteelLayers > 0)
	{
		Vector vecDataRS(2 * numberReinforcedSteelLayers);
		int counter = 0;
		int i;
		int matDbTag;
		for (i = 0; i < numberReinforcedSteelLayers; i++)
		{
			vecDataRS(counter++) = (double)TheReinforcedSteel2DMaterial[i]->getClassTag();
			matDbTag = TheReinforcedSteel2DMaterial[i]->getDbTag();
			if (matDbTag == 0) {
				matDbTag = theChannel.getDbTag();
				if (matDbTag != 0)
					TheReinforcedSteel2DMaterial[i]->setDbTag(matDbTag);
			}
			vecDataRS(counter++) = (double)matDbTag;
		}

		res += theChannel.sendVector(dataTag, commitTag, vecDataRS);
		if (res < 0) {
			opserr << "WARNING ReinforcedConcreteLayeredMembraneSection::sendSelf() - " << this->getTag() << " failed to send Vector data" << endln;
			return res;
		}

		for (i = 0; i < numberReinforcedSteelLayers; i++) {
			res += TheReinforcedSteel2DMaterial[i]->sendSelf(commitTag, theChannel);
			if (res < 0) {
				opserr << "WARNING ReinforcedConcreteLayeredMembraneSection::sendSelf() - " << this->getTag() << " failed to send its Material" << endln;
				return res;
			}
		}
	}

	if (numberConcreteLayers > 0)
	{
		Vector vecDataC(3 * numberConcreteLayers);
		int counter = 0;
		int i;

		for (i = 0; i < numberConcreteLayers; i++) {
			vecDataC(counter++) = t[i];
		}

		int matDbTag;
		for (i = 0; i < numberConcreteLayers; i++)
		{
			vecDataC(counter++) = (double)TheConcrete2DMaterial[i]->getClassTag();
			matDbTag = TheConcrete2DMaterial[i]->getDbTag();
			if (matDbTag == 0) {
				matDbTag = theChannel.getDbTag();
				if (matDbTag != 0)
					TheConcrete2DMaterial[i]->setDbTag(matDbTag);
			}
			vecDataC(counter++) = (double)matDbTag;
		}

		res += theChannel.sendVector(dataTag, commitTag, vecDataC);
		if (res < 0) {
			opserr << "WARNING ReinforcedConcreteLayeredMembraneSection::sendSelf() - " << this->getTag() << " failed to send Vector data" << endln;
			return res;
		}

		for (i = 0; i < numberConcreteLayers; i++) {
			res += TheConcrete2DMaterial[i]->sendSelf(commitTag, theChannel);
			if (res < 0) {
				opserr << "WARNING ReinforcedConcreteLayeredMembraneSection::sendSelf() - " << this->getTag() << " failed to send its Material" << endln;
				return res;
			}
		}
	}

	return res;
}

int ReinforcedConcreteLayeredMembraneSection::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	int res = 0;

	int dataTag = this->getDbTag();

	static ID iData(3);
	res += theChannel.recvID(dataTag, commitTag, iData);

	if (res < 0) {
		opserr << "WARNING ReinforcedConcreteLayeredMembraneSection::recvSelf() - " << this->getTag() << " failed to receive ID data" << endln;
		return res;
	}

	this->setTag(iData(0));

	int i;
	if (numberReinforcedSteelLayers != iData(1))
	{
		numberReinforcedSteelLayers = iData(1);
		if (TheReinforcedSteel2DMaterial != 0)
		{
			for (i = 0; i < numberReinforcedSteelLayers; i++)
			{
				if (TheReinforcedSteel2DMaterial[i] != 0) delete TheReinforcedSteel2DMaterial[i];
			}
			delete[] TheReinforcedSteel2DMaterial;
		}
		TheReinforcedSteel2DMaterial = new NDMaterial * [numberReinforcedSteelLayers];
		for (i = 0; i < numberReinforcedSteelLayers; i++)
			TheReinforcedSteel2DMaterial[i] = nullptr;
	}

	if (numberConcreteLayers != iData(2))
	{
		numberConcreteLayers = iData(2);
		if (t != 0) delete[] t;
		t = new double[numberConcreteLayers];
		if (TheConcrete2DMaterial != 0)
		{
			for (i = 0; i < numberConcreteLayers; i++)
			{
				if (TheConcrete2DMaterial[i] != 0) delete TheConcrete2DMaterial[i];
			}
			delete[] TheConcrete2DMaterial;
		}
		TheConcrete2DMaterial = new NDMaterial * [numberConcreteLayers];
		for (i = 0; i < numberConcreteLayers; i++)
			TheConcrete2DMaterial[i] = nullptr;
	}

	if (numberReinforcedSteelLayers > 0)
	{
		Vector vecDataRS(2 * numberReinforcedSteelLayers);
		res += theChannel.recvVector(dataTag, commitTag, vecDataRS);
		if (res < 0) {
			opserr << "WARNING ReinforcedConcreteLayeredMembraneSection::recvSelf() - " << this->getTag() << " failed to receive Vector data" << endln;
			return res;
		}
		int counter = 0;
		for (i = 0; i < numberReinforcedSteelLayers; i++) {
			int matClassTag = (int)vecDataRS(counter++);
			if (TheReinforcedSteel2DMaterial[i] == nullptr || TheReinforcedSteel2DMaterial[i]->getClassTag() != matClassTag) {
				if (TheReinforcedSteel2DMaterial[i]) delete TheReinforcedSteel2DMaterial[i];
				TheReinforcedSteel2DMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
				if (TheReinforcedSteel2DMaterial[i] == nullptr) {
					opserr << "ReinforcedConcreteLayeredMembraneSection::recvSelf() - Broker could no create NDMaterial of class type " << matClassTag << endln;
					return -1;
				}
			}
			TheReinforcedSteel2DMaterial[i]->setDbTag((int)vecDataRS(counter++));
			res += TheReinforcedSteel2DMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
			if (res < 0) {
				opserr << "ReinforcedConcreteLayeredMembraneSection::recvSelf() - material " << i << " failed to recv itself" << endln;
				return res;
			}
		}
	}

	if (numberConcreteLayers > 0)
	{
		Vector vecDataC(3 * numberConcreteLayers);
		res += theChannel.recvVector(dataTag, commitTag, vecDataC);
		if (res < 0) {
			opserr << "WARNING ReinforcedConcreteLayeredMembraneSection::recvSelf() - " << this->getTag() << " failed to receive Vector data" << endln;
			return res;
		}
		int counter = 0;
		for (i = 0; i < numberConcreteLayers; i++)
		{
			t[i] = vecDataC[counter++];
		}
		for (i = 0; i < numberConcreteLayers; i++) {
			int matClassTag = (int)vecDataC(counter++);
			if (TheConcrete2DMaterial[i] == nullptr || TheConcrete2DMaterial[i]->getClassTag() != matClassTag) {
				if (TheConcrete2DMaterial[i]) delete TheConcrete2DMaterial[i];
				TheConcrete2DMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
				if (TheConcrete2DMaterial[i] == nullptr) {
					opserr << "ReinforcedConcreteLayeredMembraneSection::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
					return -1;
				}
			}
			TheConcrete2DMaterial[i]->setDbTag((int)vecDataC(counter++));
			res += TheConcrete2DMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
			if (res < 0) {
				opserr << "ReinforcedConcreteLayeredMembraneSection::recvSelf() - material " << i << " failed to recv itself" << endln;
				return res;
			}
		}
	}

	return res;
}

void ReinforcedConcreteLayeredMembraneSection::Print(OPS_Stream& s, int flag)
{
	s << "ReinforcedConcreteLayeredMembraneSection tag: " << this->getTag() << endln;
	s << "Total thickness h = " << h << endln;

	for (int iC = 0; iC < numberConcreteLayers; iC++) {
		s << "Concrete layer " << iC + 1 << ", thickness t = " << t[iC] << endln;
		TheConcrete2DMaterial[iC]->Print(s, flag);
		s << endln;
	}

	for (int iRS = 0; iRS < numberReinforcedSteelLayers; iRS++) {
		s << "Reinforced Steel layer " << iRS + 1 << endln;
		TheReinforcedSteel2DMaterial[iRS]->Print(s, flag);
		s << endln;
	}

}

Response* ReinforcedConcreteLayeredMembraneSection::setResponse(const char** argv, int argc, OPS_Stream& s)
{
	Response* theResponse = 0;

	if (strcmp(argv[0], "panel_avg_stress") == 0 || strcmp(argv[0], "Panel_AvgStress") == 0) {
		s.tag("SectionOutput");
		s.attr("secType", "ReinforcedConcreteLayeredMembraneSection");
		s.attr("secTag", this->getTag());
		s.tag("ResponseType", "sigma11");
		s.tag("ResponseType", "sigma22");
		s.tag("ResponseType", "sigma12");
		s.endTag();

		Vector data1(3);
		data1.Zero();

		theResponse = new MaterialResponse(this, 1, data1);

	}
	else if (strcmp(argv[0], "cracking_pattern") == 0 || strcmp(argv[0], "Cracking_pattern") == 0) {
		s.tag("SectionOutput");
		s.attr("secType", "ReinforcedConcreteLayeredMembraneSection");
		s.attr("secTag", this->getTag());
		s.tag("ResponseType", "CP1");
		s.tag("ResponseType", "CP2");
		s.tag("ResponseType", "CP3");
		s.tag("ResponseType", "CP4");
		s.endTag();

		Vector data2(4);
		data2.Zero();

		theResponse = new MaterialResponse(this, 2, data2);

	}
	else if (strcmp(argv[0], "thetaPD_angle") == 0 || strcmp(argv[0], "ThetaPD") == 0 || strcmp(argv[0],"thetaPD") == 0) {
		s.tag("SectionOutput");
		s.attr("secType", "ReinforcedConcreteLayeredMembraneSection");
		s.attr("secTag", this->getTag());
		s.tag("ResponseType", "thetaPD");
		s.endTag();

		theResponse = new MaterialResponse(this, 3, 0.0);

	}
	else if (strcmp(argv[0], "getBendingParameters") == 0) {
		s.tag("SectionOutput");
		s.attr("secType", "ReinforcedConcreteLayeredMembraneSection");
		s.attr("secTag", this->getTag());
		s.tag("ResponseType", "Eave");
		s.tag("ResponseType", "Tave");
		s.endTag();

		Vector data3(2);
		data3.Zero();

		theResponse = new MaterialResponse(this, 4, data3);
	}
	else if (strcmp(argv[0], "panel_force") == 0 || strcmp(argv[0], "Panel_Force") == 0 || strcmp(argv[0], "Panel_force") == 0) {
		s.tag("SectionOutput");
		s.attr("secType", "ReinforcedConcreteLayeredMembraneSection");
		s.attr("secTag", this->getTag());
		s.tag("ResponseType", "Nxx");
		s.tag("ResponseType", "Nyy");
		s.tag("ResponseType", "Nxy");
		s.endTag();

		Vector data4(3);
		data4.Zero();

		theResponse = new MaterialResponse(this, 5, data4);
	}
	else if (strcmp(argv[0], "panel_strain") == 0 || strcmp(argv[0], "Panel_Strain") == 0 || strcmp(argv[0], "Panel_strain") == 0) {
		s.tag("SectionOutput");
		s.attr("secType", "ReinforcedConcreteLayeredMembraneSection");
		s.attr("secTag", this->getTag());
		s.tag("ResponseType", "eps11");
		s.tag("ResponseType", "eps22");
		s.tag("ResponseType", "eps12");
		s.endTag();

		Vector data5(3);
		data5.Zero();

		theResponse = new MaterialResponse(this, 6, data5);
	}
	else if (strcmp(argv[0], "CLayer") == 0 || strcmp(argv[0], "ConcLayer") == 0) {
		if (argc != 3) {
			opserr << "ReinforcedConcreteLayeredMembraneSection::setResponse() - need to specify more data\n";
			return 0;
		}
		int pointNum = atoi(argv[1]);
		if (pointNum > 0 && pointNum <= numberConcreteLayers) {
			s.tag("LayerOutput");
			s.attr("number", pointNum);
			s.attr("thickness", t[pointNum]);

			theResponse = TheConcrete2DMaterial[pointNum - 1]->setResponse(&argv[2], argc - 2, s);

		}
	}
	else if (strcmp(argv[0], "RSLayer") == 0 || strcmp(argv[0], "ReinfSteelLayer") == 0) {
		if (argc != 3) {
			opserr << "ReinforcedConcreteLayeredMembraneSection::setResponse() - need to specify more data\n";
			return 0;
		}
		int pointNum = atoi(argv[1]);
		if (pointNum > 0 && pointNum <= numberReinforcedSteelLayers) {
			s.tag("LayerOutput");
			s.attr("number", pointNum);

			theResponse = TheReinforcedSteel2DMaterial[pointNum - 1]->setResponse(&argv[2], argc - 2, s);

		}
	}
	else {

		return this->SectionForceDeformation::setResponse(argv, argc, s);
	}


	return theResponse;
}


int ReinforcedConcreteLayeredMembraneSection::getResponse(int responseID, Information& info)
{
	if (responseID == 1) {
		return info.setVector(this->getSectionStressAvg());
	}
	else if (responseID == 2) {
		return info.setVector(this->getCrackPattern());
	}
	else if (responseID == 3) {
		return info.setDouble(this->getThetaPDAngle());
	}
	else if (responseID == 4) {
		return info.setVector(this->getBendingParameters());
	}
	else if (responseID == 5) {
		return info.setVector(this->getCommittedStress());
	}
	else if (responseID == 6) {
		return info.setVector(this->getCommittedStrain());
	}
	else {
		return 0;

	}
}

// Function that returns bending parameters - added for MEFI3D by Maria Jose Nunez, UChile
Vector ReinforcedConcreteLayeredMembraneSection::getBendingParameters(void)
{
	Vector input_par(2);

	input_par.Zero();

	input_par(0) = this->getEcAvg();
	input_par(1) = h;

	return input_par;
}

void ReinforcedConcreteLayeredMembraneSection::setCrackPattern(void)
{
	// Get strain concrete propierties (ec, ecr)
	double* strainAtFt = new double[numberConcreteLayers];
	double* strainAtFc = new double[numberConcreteLayers];

	DummyStream theDummyStream;
	char aa[80] = "getParameters";
	const char* argv[1];
	argv[0] = aa;

	for (int i = 0; i < numberConcreteLayers; i++)
	{
		Response* theResponse = TheConcrete2DMaterial[i]->setResponse(argv, 1, theDummyStream);

		if (theResponse == 0) {
			opserr << "ReinforcedConcreteLayeredMembraneSection::ReinforcedConcreteLayeredMembraneSection - failed to get strain concrete parameters for RCLMS section with tag: " << this->getTag() << "\n";
			exit(-1);
		}

		theResponse->getResponse();
		Information& theInfoPar = theResponse->getInformation();
		const Vector& ParNDConc = theInfoPar.getData();

		strainAtFc[i] = ParNDConc[2];
		strainAtFt[i] = ParNDConc[3];

		delete theResponse;
	}
	 
	double maxVal = numeric_limits<double>::lowest();
	double minVal = numeric_limits<double>::max();

	for (int i = 0; i < numberConcreteLayers; i++) {
		if (strainAtFc[i] > maxVal) {
			maxVal = strainAtFc[i];
		}
		if (strainAtFt[i] < minVal) {
			minVal = strainAtFt[i];
		}
	}

	ec = maxVal;
	ecr = minVal;


	delete[] strainAtFc;
	delete[] strainAtFt;

	// Calculate the Principal Direction of Strain
	this->calculateStrainPrincipalDirections01();

	// Calculate the Poisson Ratios
	if (crackPattern(0) == 1 || crackPattern(2) == 1) {
		isConcreteCracked = true;
	}
	else {
		isConcreteCracked = false;
	}

	this->calculatePoissonRatios(strainPrincipalDirection(0), strainPrincipalDirection(1));

	// Define the Poisson Ratio Matrix V = [ 1/(1-nu12*nu21) nu12/(1-nu12*nu21) 0 ; nu21/(1-nu12*nu21) 1/(1-nu12*nu21) 0 ; 0 0 1 ]
	double nu[2];
	nu[0] = poissonRatios(0); nu[1] = poissonRatios(1);
	double oneOverOneMinusNu12Nu21 = 1.0 / (1.0 - nu[0] * nu[1]);

	double V[3][3] = { {oneOverOneMinusNu12Nu21,nu[0] * oneOverOneMinusNu12Nu21,0.0},{nu[1] * oneOverOneMinusNu12Nu21,oneOverOneMinusNu12Nu21,0.0},{0.0,0.0,1.0} };

	// Calculate the newUniaxialStrain in the orientation of the Principal Direction of Strain
	double newUniaxialStrainPD[3] = { 0.0,0.0,0.0 };

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			newUniaxialStrainPD[i] = newUniaxialStrainPD[i] + (strainPrincipalDirection(j) * V[i][j]);
		}
	}

	// Set the CrackPattern Propierties
	// Check if the Concrete has crushing
	// Select the minimal strain in the principal direction
	double minStrain = min(newUniaxialStrainPD[0], newUniaxialStrainPD[1]);
	if (minStrain <= ec || (crackPattern(2) == 1 && minStrain < 0)) {
		crackPattern(2) = 1;
		crackPattern(3) = thetaPrincipalDirection;
	}

	// Check if the concrete has been cracked in the direction 1
	double maxStrain = max(newUniaxialStrainPD[0], newUniaxialStrainPD[1]);
	if (maxStrain >= ecr || (crackPattern(0) == 1 && maxStrain >= 0)) {
		crackPattern(0) = 1;
		crackPattern(1) = thetaPrincipalDirection;
	}

}

const Vector& ReinforcedConcreteLayeredMembraneSection::getCrackPattern(void)
{
	return crackPattern;
}

double ReinforcedConcreteLayeredMembraneSection::getThetaPDAngle(void)
{
	return thetaPrincipalDirection;
}

Vector ReinforcedConcreteLayeredMembraneSection::getSectionStressAvg(void)
{
	Vector SectionStressAvg(3);

	SectionStressAvg.Zero();

	SectionStressAvg(0) = TSectionStress(0) / h;
	SectionStressAvg(1) = TSectionStress(1) / h;
	SectionStressAvg(2) = TSectionStress(2) / h;

	return SectionStressAvg;
}

void ReinforcedConcreteLayeredMembraneSection::calculateStrainPrincipalDirections01(void)
{
	double strain_vec[3];		//ex, ey, gamma
	double doubleThetaPD, cos2Theta, sin2Theta;

	// Get strain values from strain of element
	strain_vec[0] = TSectionStrain(0);		// ex
	strain_vec[1] = TSectionStrain(1);		// ey
	strain_vec[2] = TSectionStrain(2);		// gxy

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

void ReinforcedConcreteLayeredMembraneSection::calculateAngle01(double cosTheta, double sinTheta, double& theta)
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

void ReinforcedConcreteLayeredMembraneSection::calculatePoissonRatios(double e1, double e2)
{
	double nu[2];

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