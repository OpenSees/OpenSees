#include "stdafx.h"
#include "BrickWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Elements;
using namespace OpenSees::Elements::BeamIntegrations;
using namespace OpenSees::Materials::Sections;
using namespace OpenSees::Materials::Uniaxials;

BbarBrickWrapper::BbarBrickWrapper(int tag,
	int node1,
	int node2,
	int node3,
	int node4,
	int node5,
	int node6,
	int node7,
	int node8,
	NDMaterialWrapper^ theMaterial,
	double b1, double b2, double b3) {
	_Element = new BbarBrick(tag,
		node1,
		node2,
		node3,
		node4,
		node5,
		node6,
		node7,
		node8,
		*theMaterial->_NDMaterial,
		b1, b2, b3);
}

BbarBrickWithSensitivityWrapper::BbarBrickWithSensitivityWrapper(int tag,
	int node1,
	int node2,
	int node3,
	int node4,
	int node5,
	int node6,
	int node7,
	int node8,
	NDMaterialWrapper^ theMaterial,
	double b1, double b2, double b3) {
	_Element = new BbarBrickWithSensitivity(tag,
		node1,
		node2,
		node3,
		node4,
		node5,
		node6,
		node7,
		node8,
		*theMaterial->_NDMaterial,
		b1, b2, b3);
}

BrickWrapper::BrickWrapper(int tag,
	int node1,
	int node2,
	int node3,
	int node4,
	int node5,
	int node6,
	int node7,
	int node8,
	NDMaterialWrapper^ theMaterial,
	double b1, double b2, double b3) {
	_Element = new Brick(tag,
		node1,
		node2,
		node3,
		node4,
		node5,
		node6,
		node7,
		node8,
		*theMaterial->_NDMaterial,
		b1, b2, b3);
}

Twenty_Node_BrickWrapper::Twenty_Node_BrickWrapper(int tag,
	int node1,
	int node2,
	int node3,
	int node4,
	int node5,
	int node6,
	int node7,
	int node8,
	int node9,
	int node10,
	int node11,
	int node12,
	int node13,
	int node14,
	int node15,
	int node16,
	int node17,
	int node18,
	int node19,
	int node20,
	NDMaterialWrapper^ theMaterial,
	double b1, double b2, double b3) {
	_Element = new Twenty_Node_Brick(tag,
		node1,
		node2,
		node3,
		node4,
		node5,
		node6,
		node7,
		node8,
		node9,
		node10,
		node11,
		node12,
		node13,
		node14,
		node15,
		node16,
		node17,
		node18,
		node19,
		node20,
		*theMaterial->_NDMaterial,
		b1, b2, b3);
}

