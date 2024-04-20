// Code written/implemented by: Fabian Rojas B.
//								Maria Jose Nunez
//
// Created: 07/2023
// 
// Description: This file contains the LayeredMembraneSection class definition
// A LayeredMembraneSection is a subclass of the sectionForceDeformation class and corresponds to the abstract representation
// for the stress-strain behavior for a Reinforced Concrete Layer Membrane Element in the Finite Element Method or Structural Analysis. 
//
// Reference:
// 1. Rojas, F., Anderson, J. C., Massone, L. M. (2016). A nonlinear quadrilateral layered membrane element with drilling degrees of freedom for 
// the modeling of reinforced concrete walls. Engineering Structures, 124, 521-538.
//
// Source: \OpenSees\SRC\material\section\LayeredMembraneSection
//
// Rev: 1.0

#include <LayeredMembraneSection.h>
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
void* OPS_LayeredMembraneSection()
{
	int numArgs = OPS_GetNumRemainingInputArgs();

	// Parse the script for material parameters
	if (numArgs < 7) {
		opserr << "Want: LayeredMembraneSection $secTag $total_thickness $nLayers -mat{$Material_tags} -thick{$Thicknesses} <-Eout $YoungModulus>" << endln;
		return 0;
	}

	int tag, nLayers;
	double totalThickness;
	int* matTags;								
	double* thickness;
	const char* str;
	NDMaterial** theMats;
	double Eaverage = 0.0;

	// section tag
	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) {
		opserr << "WARNING invalid section LayeredMembraneSection tag" << endln;
		return 0;
	}

	// section total thickness
	if (OPS_GetDouble(&numdata, &totalThickness) != 0) {
		opserr << "WARNING invalid tota_thickness" << endln;
		opserr << "LayeredMembraneSection section: " << tag << endln;
		return 0;
	}

	// layers number
	if (OPS_GetIntInput(&numdata, &nLayers) < 0) {
		opserr << "WARNING invalid nLayers" << endln;
		opserr << "LayeredMembraneSection section: " << tag << endln;
		return 0;
	}

	theMats = new NDMaterial * [nLayers];
	thickness = new double[nLayers];
	matTags = new int[nLayers];

	numArgs = OPS_GetNumRemainingInputArgs();

	while (numArgs > 0) {
		str = OPS_GetString();
		if (strcmp(str, "-mat") == 0) {
			numdata = nLayers;
			if (OPS_GetIntInput(&numdata, matTags) != 0) {
				opserr << "Invalid nDMaterial tag for LayeredMembraneSection " << tag << endln;
				return 0;
			}
			for (int i = 0; i < nLayers; i++) {
				theMats[i] = 0;
				theMats[i] = OPS_getNDMaterial(matTags[i]);
				if (theMats[i] == 0) {
					opserr << "Invalid nDMaterial tag " << matTags[i] << " for LayeredMembraneSection " << tag << endln;
					return 0;
				}
			}
		}
		else if (strcmp(str, "-thick") == 0) {
			numdata = nLayers;
			if (OPS_GetDoubleInput(&numdata, thickness) != 0) {
				opserr << "Invalid thickness parameter for LayeredMembraneSection " << tag << endln;
				return 0;
			}
		}
		else if (strcmp(str, "-Eout") == 0) {
			numdata = 1;
			if (OPS_GetDouble(&numdata, &Eaverage) != 0) {
				opserr << "Invalid Eout parameter for LayeredMembraneSection " << tag << endln;
				return 0;
			}
		}
		else {
			opserr << "WARNING: Invalid option " << str << " in LayeredMembraneSection " << tag << endln;
			return 0;
		}

		numArgs = OPS_GetNumRemainingInputArgs();
	}

	SectionForceDeformation* theSection = new LayeredMembraneSection(tag, totalThickness, nLayers, theMats,thickness,Eaverage);

	// Clean uo dynamic memory
	if (theMats != 0)
		delete[] theMats;
	if (matTags != 0)
		delete[] matTags;
	if (thickness != 0)
		delete[] thickness;

	return theSection;
}

//static vector and matrices
ID LayeredMembraneSection::array(3);

// Full constructor
LayeredMembraneSection::LayeredMembraneSection(int tag,		                                    // section tag
	double totalThickness,																		// section total thickness
	int nLayers,																				// number of layers
	NDMaterial** MaterialObjects,																// array of nDMaterial tags
	double* Thickness,																			// array of layers thicknesses
	double Eaverage)																			// modulus of elasticity																		

	:SectionForceDeformation(tag, SEC_TAG_LayeredMembraneSection),
	t_total(totalThickness), numberLayers(nLayers), Eave(Eaverage),
	CSectionStrain(3), CSectionStress(3), CSectionTangent(3, 3),
	TSectionStrain(3), TSectionStress(3), TSectionTangent(3,3), InitialTangent(3,3)
{
	// Set initial values
	for (int i = 0; i < 3; i++) {
		CSectionStrain(i) = 0.0;
		CSectionStress(i) = 0.0;
		TSectionStrain(i) = 0.0;
		TSectionStress(i) = 0.0;
	}

	// Check materials input
	if (MaterialObjects == 0) {
		opserr << "LayeredMembraneSection::LayeredMembraneSection() - Null nDMaterial array passed.\n";
		exit(-1);
	}

	// Allocate memory for layers thicknesses
	t = new double[numberLayers];

	for (int i = 0; i < numberLayers; i++) {
		t[i] = Thickness[i];
	}

	// Allocate memory for the NDMaterials
	The2DMaterials = new NDMaterial * [numberLayers];

	if (The2DMaterials == 0) {
		opserr << "LayeredMembraneSection::LayeredMembraneSection() - Failed to allocate pointers for ND materials.\n";
		exit(-1);
	}

	// Get copies of the ND materials
	for (int i = 0; i < numberLayers; i++) {
		if (MaterialObjects[i] == 0) {
			opserr << "LayeredMembraneSection::LayeredMembraneSection() - Null ND material pointer passed.\n";
			exit(-1);
		}

		The2DMaterials[i] = MaterialObjects[i]->getCopy("PlaneStress2D");

		if (The2DMaterials[i] == 0) {
			opserr << "LayeredMembraneSection::LayeredMembraneSection() - Failed to copy ND material.\n";
			exit(-1);
		}
	}

}

// Blank constructor (constructor for blank object that recvSelf needs to be invoked upon) (constructor which should be invoked by an FEM_ObjectBroker only)
LayeredMembraneSection::LayeredMembraneSection() :
	SectionForceDeformation(0, SEC_TAG_LayeredMembraneSection),
	t_total(0), numberLayers(0), Eave(0),
	CSectionStrain(3), CSectionStress(3), CSectionTangent(3, 3),
	TSectionStrain(3), TSectionStress(3), TSectionTangent(3, 3), InitialTangent(3, 3)
{
	The2DMaterials = 0;

	this->revertToStart();
}

// Destructor (clean up memory the LayeredMembraneSection objects allocates)
LayeredMembraneSection::~LayeredMembraneSection()
{
	if (t != 0) delete[]t;

	if (The2DMaterials != 0) {
		for (int i = 0; i < numberLayers; i++)
			if (The2DMaterials[i] != 0)
				delete The2DMaterials[i];
		delete[]The2DMaterials;
	}
}

int LayeredMembraneSection::setTrialSectionDeformation(const Vector& newTrialSectionStrain)
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


	// Set the strain at each layer
	for (int i = 0; i < numberLayers; i++) {
		// Get nDMaterial stress and tangent stiffness
		The2DMaterials[i]->setTrialStrain(newTrialSectionStrain);
		const Vector& stressNDM = The2DMaterials[i]->getStress();
		const Matrix& tangentNDM = The2DMaterials[i]->getTangent();

		double sig11 = stressNDM(0);
		double sig22 = stressNDM(1);
		double sig12 = stressNDM(2);
		// Calculate the layer contribution to section resultant stress
		sectionResultantStress[0] += t[i] * sig11;
		sectionResultantStress[1] += t[i] * sig22;
		sectionResultantStress[2] += t[i] * sig12;

		double E00 = tangentNDM(0, 0); double E01 = tangentNDM(0, 1); double E02 = tangentNDM(0, 2);
		double E10 = tangentNDM(1, 0); double E11 = tangentNDM(1, 1); double E12 = tangentNDM(1, 2);
		double E20 = tangentNDM(2, 0); double E21 = tangentNDM(2, 1); double E22 = tangentNDM(2, 2);
		// Calculate the layer contribution to section tangent
		sectionTangent[0][0] += t[i] * E00; sectionTangent[0][1] += t[i] * E01;	sectionTangent[0][2] += t[i] * E02;
		sectionTangent[1][0] += t[i] * E10; sectionTangent[1][1] += t[i] * E11;	sectionTangent[1][2] += t[i] * E12;
		sectionTangent[2][0] += t[i] * E20; sectionTangent[2][1] += t[i] * E21;	sectionTangent[2][2] += t[i] * E22;
	}
	
	// Set the stress and tangent
	TSectionStress(0) = sectionResultantStress[0]; TSectionStress(1) = sectionResultantStress[1]; TSectionStress(2) = sectionResultantStress[2];

	TSectionTangent(0, 0) = sectionTangent[0][0]; TSectionTangent(0, 1) = sectionTangent[0][1]; TSectionTangent(0, 2) = sectionTangent[0][2];
	TSectionTangent(1, 0) = sectionTangent[1][0]; TSectionTangent(1, 1) = sectionTangent[1][1]; TSectionTangent(1, 2) = sectionTangent[1][2];
	TSectionTangent(2, 0) = sectionTangent[2][0]; TSectionTangent(2, 1) = sectionTangent[2][1]; TSectionTangent(2, 2) = sectionTangent[2][2];

	return 0;
}

const Vector& LayeredMembraneSection::getSectionDeformation(void)
{
	return TSectionStrain;
}

const Vector& LayeredMembraneSection::getStressResultant(void)
{
	return TSectionStress;
}

const Matrix& LayeredMembraneSection::getSectionTangent(void)
{
	return TSectionTangent;
}

const Vector& LayeredMembraneSection::getCommittedStrain(void)
{
	return CSectionStrain;
}

const Vector& LayeredMembraneSection::getCommittedStress(void)
{
	return CSectionStress;
}

double LayeredMembraneSection::getRho(void)
{
	double rho = 0.0;

	for (int i = 0; i < numberLayers; i++) {
		rho += (The2DMaterials[i]->getRho()) * t[i];
	}

	return rho;
}


const Matrix& LayeredMembraneSection::getInitialTangent(void)
{
	// Get the initial membrane section stiffness
	double initialSectionTangent[3][3];			// membrane section stiffness: [Dm] = h*[MaterialTangent] is the membrane stiffness Tangent = Sum(Tangent_ic(zTop_ic - zBottom_ic)) + Sum(Tangent_is(t_is))
	initialSectionTangent[0][0] = 0.0; initialSectionTangent[0][1] = 0.0;	initialSectionTangent[0][2] = 0.0;
	initialSectionTangent[1][0] = 0.0; initialSectionTangent[1][1] = 0.0;	initialSectionTangent[1][2] = 0.0;
	initialSectionTangent[2][0] = 0.0; initialSectionTangent[2][1] = 0.0;	initialSectionTangent[2][2] = 0.0;

	InitialTangent.Zero();

	// Add the section stiffness contribution from each layer
	for (int i = 0; i < numberLayers; i++) {
		// Get material initial tangent stiffness
		const Matrix& layerMaterialTangent = The2DMaterials[i]->getInitialTangent();
		
		double Eo00 = layerMaterialTangent(0, 0); double Eo01 = layerMaterialTangent(0, 1); double Eo02 = layerMaterialTangent(0, 2);
		double Eo10 = layerMaterialTangent(1, 0); double Eo11 = layerMaterialTangent(1, 1); double Eo12 = layerMaterialTangent(1, 2);
		double Eo20 = layerMaterialTangent(2, 0); double Eo21 = layerMaterialTangent(2, 1); double Eo22 = layerMaterialTangent(2, 2);

		// Calculate the material contribution to section tangent
		initialSectionTangent[0][0] += t[i] * Eo00; initialSectionTangent[0][1] += t[i] * Eo01;	initialSectionTangent[0][2] += t[i] * Eo02;
		initialSectionTangent[1][0] += t[i] * Eo10; initialSectionTangent[1][1] += t[i] * Eo11;	initialSectionTangent[1][2] += t[i] * Eo12;
		initialSectionTangent[2][0] += t[i] * Eo20; initialSectionTangent[2][1] += t[i] * Eo21;	initialSectionTangent[2][2] += t[i] * Eo22;
	}

	InitialTangent(0, 0) = initialSectionTangent[0][0]; InitialTangent(0, 1) = initialSectionTangent[0][1]; InitialTangent(0, 2) = initialSectionTangent[0][2];
	InitialTangent(1, 0) = initialSectionTangent[1][0]; InitialTangent(1, 1) = initialSectionTangent[1][1]; InitialTangent(1, 2) = initialSectionTangent[1][2];
	InitialTangent(2, 0) = initialSectionTangent[2][0]; InitialTangent(2, 1) = initialSectionTangent[2][1]; InitialTangent(2, 2) = initialSectionTangent[2][2];

	return InitialTangent;
}

SectionForceDeformation* LayeredMembraneSection::getCopy(void)
{
	LayeredMembraneSection* theCopy = new LayeredMembraneSection(this->getTag(),
		t_total,
		numberLayers,
		The2DMaterials,
		t,
		Eave);

	return theCopy;
}

const ID& LayeredMembraneSection::getType(void)
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

int LayeredMembraneSection::getOrder(void) const
{
	return 3;
}

int LayeredMembraneSection::commitState(void)
{
	int success = 0;

	// Commit the state for each layer
	for (int i = 0; i < numberLayers; i++) {
		success += The2DMaterials[i]->commitState();
	}

	// Commit the history variables
	CSectionStrain = TSectionStrain;
	CSectionStress = TSectionStress;
	CSectionTangent = TSectionTangent;

	return success;
}

int LayeredMembraneSection::revertToLastCommit(void)
{
	int success = 0;

	// Revert to last commit the state of each layer
	for (int i = 0; i < numberLayers; i++) {
		success += The2DMaterials[i]->revertToLastCommit();
	}

	// Revert the history variables to last commit
	TSectionStrain = CSectionStrain;
	TSectionStress = CSectionStress;
	TSectionTangent = TSectionTangent;

	return success;
}

int LayeredMembraneSection::revertToStart(void)
{
	int success = 0;

	// Revert to start the state of each layer
	for (int i = 0; i < numberLayers; i++)
		success += The2DMaterials[i]->revertToStart();

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

int LayeredMembraneSection::sendSelf(int commitTag, Channel& theChannel)
{
	int res = 0;
	int dataTag = this->getDbTag();

	//static ID iData(4);
	static Vector iData(4);
	iData(0) = this->getTag();
	iData(1) = t_total;
	iData(2) = numberLayers;
	iData(3) = Eave;

	//res += theChannel.sendID(dataTag, commitTag, iData);
	res += theChannel.sendVector(dataTag, commitTag, iData);
	if (res < 0) {
		opserr << "WARNING LayeredMembraneSection:: sendSelf() - " << this->getTag() << " failed to send vector data" << endln;
		return res;
	}

	if (numberLayers > 0)
	{
		Vector vecData(3 * numberLayers);
		int counter = 0;
		int i;

		for (i = 0; i < numberLayers; i++) {
			vecData(counter++) = t[i];
		}

		int matDbTag;
		for (i = 0; i < numberLayers; i++) {
			vecData(counter++) = (double)The2DMaterials[i]->getClassTag();
			matDbTag = The2DMaterials[i]->getDbTag();
			if (matDbTag == 0) {
				matDbTag = theChannel.getDbTag();
				if (matDbTag != 0)
					The2DMaterials[i]->setDbTag(matDbTag);
			}
			vecData(counter++) = (double)matDbTag;
		}
		res += theChannel.sendVector(dataTag, commitTag, vecData);
		if (res < 0) {
			opserr << "WARNING LayeredMembraneSection::sendSelf() - " << this->getTag() << " failed to send Vector data" << endln;
			return res;
		}

		for (i = 0; i < numberLayers; i++) {
			res += The2DMaterials[i]->sendSelf(commitTag, theChannel);
			if (res < 0) {
				opserr << "WARNING LayeredMembraneSection::sendSelf() - " << this->getTag() << " failed to send its Material" << endln;
				return res;
			}
		}
	}

	return res;
}

int LayeredMembraneSection::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	int res = 0;

	int dataTag = this->getDbTag();

	static Vector iData(4);
	res += theChannel.recvVector(dataTag, commitTag, iData);

	if (res < 0) {
		opserr << "WARNING LayeredMembraneSection::recvSelf() - " << this->getTag() << " failed to receive vector data" << endln;
		return res;
	}

	this->setTag((int)iData(0));
	t_total = iData(1);

	int i;
	if (numberLayers != (int)iData(2))
	{
		numberLayers = (int)iData(2);
		if (t != 0) delete[]t;
		t = new double[numberLayers];
		if (The2DMaterials != 0)
		{
			for (i = 0; i < numberLayers; i++)
			{
				if (The2DMaterials[i] != 0)delete The2DMaterials[i];
			}
			delete[]The2DMaterials;
		}
		The2DMaterials = new NDMaterial * [numberLayers];
		for (i = 0; i < numberLayers; i++)
			The2DMaterials[i] = nullptr;
	}

	Eave	= iData(3);

	if (numberLayers > 0)
	{
		Vector vecData(3 * numberLayers);
		res += theChannel.recvVector(dataTag, commitTag, vecData);
		if (res < 0) {
			opserr << "WARNING LayeredMembraneSection::recvSelf() - " << this->getTag() << " failed to receive Vector data" << endln;
			return res;
		}
		int counter = 0;
		for (i = 0; i < numberLayers; i++)
		{
			t[i] = vecData[counter++];
		}
		for (i = 0; i < numberLayers; i++) {
			int matClassTag = (int)vecData(counter++);
			if (The2DMaterials[i] == nullptr || The2DMaterials[i]->getClassTag() != matClassTag) {
				if (The2DMaterials[i])delete The2DMaterials[i];
				The2DMaterials[i] = theBroker.getNewNDMaterial(matClassTag);
				if (The2DMaterials[i] == nullptr) {
					opserr << "LayeredMembraneSection::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
					return -1;
				}
			}
			The2DMaterials[i]->setDbTag((int)vecData(counter++));
			res += The2DMaterials[i]->recvSelf(commitTag, theChannel, theBroker);
			if (res < 0) {
				opserr << "LayeredMembraneSection::recvSelf() - material " << i << " failed to recv itself" << endln;
				return res;
			}
		}
	}

	return res;
}

void LayeredMembraneSection::Print(OPS_Stream& s, int flag)
{
	s << "LayeredMembraneSection tag: " << this->getTag() << endln;
	s << "Total thickness h = " << t_total << endln;
	
	for (int i = 0; i < numberLayers; i++) {
		s << "Layer " << i + 1 << ", thickness t = " << t[i] << endln;
		The2DMaterials[i]->Print(s, flag);
		s << endln;
	}

}

Response* LayeredMembraneSection::setResponse(const char** argv, int argc, OPS_Stream& s)
{
	Response* theResponse = 0;
	if (strcmp(argv[0], "panel_strain") == 0 || strcmp(argv[0], "Panel_Strain") == 0 || strcmp(argv[0], "Panel_strain") == 0) {
		s.tag("SectionOutput");
		s.attr("secType", "LayeredMembraneSection");
		s.attr("secTag", this->getTag());
		s.tag("ResponseType", "eps11");
		s.tag("ResponseType", "eps22");
		s.tag("ResponseType", "eps12");
		s.endTag();

		Vector data1(3);
		data1.Zero();

		theResponse = new MaterialResponse(this, 1,data1);

	}
	else if (strcmp(argv[0], "panel_force") == 0 || strcmp(argv[0], "Panel_Force") == 0 || strcmp(argv[0], "Panel_force") == 0) {
		s.tag("SectionOutput");
		s.attr("secType", "LayeredMembraneSection");
		s.attr("secTag", this->getTag());
		s.tag("ResponseType", "Nxx");
		s.tag("ResponseType", "Nyy");
		s.tag("ResponseType", "Nxy");
		s.endTag();

		Vector data2(3);
		data2.Zero();

		theResponse = new MaterialResponse(this,2,data2);

	}
	else if (strcmp(argv[0], "getBendingParameters") == 0) {
		s.tag("SectionOutput");
		s.attr("secType", "LayeredMembraneSection");
		s.attr("secTag", this->getTag());
		s.tag("ResponseType", "Eave");
		s.tag("ResponseType", "Tave");
		s.endTag();

		Vector data3(2);
		data3.Zero();

		theResponse = new MaterialResponse(this, 3, data3);
	}
	else if (strcmp(argv[0], "panel_avg_stress") == 0 || strcmp(argv[0], "Panel_AvgStress") == 0) {
		s.tag("SectionOutput");
		s.attr("secType", "LayeredMembraneSection");
		s.attr("secTag", this->getTag());
		s.tag("ResponseType", "sigma11");
		s.tag("ResponseType", "sigma22");
		s.tag("ResponseType", "sigma12");
		s.endTag();

		Vector data4(3);
		data4.Zero();

		theResponse = new MaterialResponse(this, 4, data4);
	}
	else if (strcmp(argv[0], "Layer") == 0 || strcmp(argv[0], "layer") == 0 || strcmp(argv[0], "matLayer") == 0 || strcmp(argv[0], "MatLayer") == 0) {
		if (argc != 3) {
			opserr << "LayeredMembraneSection::setResponse() - need to specify more data\n";
			return 0;
		}
		int pointNum = atoi(argv[1]);
		if (pointNum > 0 && pointNum <= numberLayers) {
			s.tag("LayerOutput");
			s.attr("number", pointNum);
			s.attr("thickness", t[pointNum]);

			theResponse = The2DMaterials[pointNum - 1]->setResponse(&argv[2], argc - 2, s);
		}
	}
	else {

		return this->SectionForceDeformation::setResponse(argv, argc, s);
	}

	return theResponse;
}


int LayeredMembraneSection::getResponse(int responseID, Information& info)
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
	else if (responseID == 4) {
		return info.setVector(this->getSectionStressAvg());
	}
	else {
		return 0;
	}
}

// Function that returns bending parameters - added for MEFI3D by Maria Jose Nunez, UChile
Vector LayeredMembraneSection::getBendingParameters(void)
{
	Vector input_par(2);

	input_par.Zero();

	input_par(0) = Eave;
	input_par(1) = t_total;				

	return input_par;
}

// Function below return values for recorder
Vector LayeredMembraneSection::getSectionStressAvg(void)
{
	Vector SectionStressAvg(3);
	SectionStressAvg.Zero();

	SectionStressAvg(0) = TSectionStress(0) / t_total;
	SectionStressAvg(1) = TSectionStress(1) / t_total;
	SectionStressAvg(2) = TSectionStress(2) / t_total;

	return SectionStressAvg;
}