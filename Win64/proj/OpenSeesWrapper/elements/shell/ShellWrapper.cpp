#include "stdafx.h"
#include "ShellWrapper.h"

using namespace System::Runtime::InteropServices;
using namespace OpenSees;
using namespace OpenSees::Elements;
using namespace OpenSees::Materials::NDMaterials;



ShellANDeSWrapper::ShellANDeSWrapper(int element_number,
	int node_numb_1, int node_numb_2, int node_numb_3, double t,
	double E, double nu, double rho)
{
	_Element = new ShellANDeS(element_number,
		node_numb_1, node_numb_2, node_numb_3, t,
		E, nu, rho);
}

ShellDKGQWrapper::ShellDKGQWrapper(int tag,
	int node1,
	int node2,
	int node3,
	int node4,
	SectionForceDeformationWrapper^ theMaterial)
{
	_Element = new ShellDKGQ(tag,
		node1,
		node2,
		node3,
		node4,
		*theMaterial->_SectionForceDeformation);
}

ShellDKGTWrapper::ShellDKGTWrapper(int tag,
	int node1,
	int node2,
	int node3,
	SectionForceDeformationWrapper^ theMaterial)
{
	_Element = new ShellDKGT(tag,
		node1,
		node2,
		node3,
		*theMaterial->_SectionForceDeformation);
}

ShellMITC4Wrapper::ShellMITC4Wrapper(int tag,
	int node1,
	int node2,
	int node3,
	int node4,
	SectionForceDeformationWrapper^ theMaterial, bool updateBasis)
{
	_Element = new ShellMITC4(tag,
		node1,
		node2,
		node3,
		node4,
		*theMaterial->_SectionForceDeformation, updateBasis);
}

ShellMITC4ThermalWrapper::ShellMITC4ThermalWrapper(int tag,
	int node1,
	int node2,
	int node3,
	int node4,
	SectionForceDeformationWrapper^ theMaterial)
{
	_Element = new ShellMITC4Thermal(tag,
		node1,
		node2,
		node3,
		node4,
		*theMaterial->_SectionForceDeformation);
}

ShellMITC9Wrapper::ShellMITC9Wrapper(int tag,
	int node1,
	int node2,
	int node3,
	int node4,
	int node5,
	int node6,
	int node7,
	int node8,
	int node9,
	SectionForceDeformationWrapper^ theMaterial)
{
	_Element = new ShellMITC9(tag,
		node1,
		node2,
		node3,
		node4,
		node5,
		node6,
		node7,
		node8,
		node9,
		*theMaterial->_SectionForceDeformation);
}

ShellNLDKGQWrapper::ShellNLDKGQWrapper(int tag,
	int node1,
	int node2,
	int node3,
	int node4,
	SectionForceDeformationWrapper^ theMaterial)
{
	_Element = new ShellNLDKGQ(tag,
		node1,
		node2,
		node3,
		node4,
		*theMaterial->_SectionForceDeformation);
}

ShellNLDKGQThermalWrapper::ShellNLDKGQThermalWrapper(int tag,
	int node1,
	int node2,
	int node3,
	int node4,
	SectionForceDeformationWrapper^ theMaterial)
{
	_Element = new ShellNLDKGQThermal(tag,
		node1,
		node2,
		node3,
		node4,
		*theMaterial->_SectionForceDeformation);
}

ShellNLDKGTWrapper::ShellNLDKGTWrapper(int tag,
	int node1,
	int node2,
	int node3,
	SectionForceDeformationWrapper^ theMaterial)
{
	_Element = new ShellNLDKGT(tag,
		node1,
		node2,
		node3,
		*theMaterial->_SectionForceDeformation);
}

