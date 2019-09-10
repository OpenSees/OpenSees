#include "stdafx.h"
#include "QuadWrapper.h"

using namespace System::Runtime::InteropServices;
using namespace OpenSees;
using namespace OpenSees::Elements;
using namespace OpenSees::Materials::NDMaterials;


FourNodeQuad3dWrapper::FourNodeQuad3dWrapper(int tag, int Nd1, int Nd2, int Nd3, int Nd4,
	NDMaterialWrapper^ nDMaterial, FourNodeQuadType^ type, double t)
{
	_Element = new FourNodeQuad3d(tag, Nd1, Nd2, Nd3, Nd4, *nDMaterial->_NDMaterial, (char*)(void*)Marshal::StringToHGlobalAnsi(type->ToString()), t);
}

ConstantPressureVolumeQuadWrapper::ConstantPressureVolumeQuadWrapper(int tag,
	int Nd1, int Nd2, int Nd3, int Nd4,
	NDMaterialWrapper^ nDMaterial,
	double t)
{
	_Element = new ConstantPressureVolumeQuad(tag,
		Nd1, Nd2, Nd3, Nd4,
		*nDMaterial->_NDMaterial,
		t);
}

EnhancedQuadWrapper::EnhancedQuadWrapper(int tag,
	int Nd1, int Nd2, int Nd3, int Nd4,
	NDMaterialWrapper^ nDMaterial,
	FourNodeQuadType^ type, double t)
{
	_Element = new EnhancedQuad(tag,
		Nd1, Nd2, Nd3, Nd4,
		*nDMaterial->_NDMaterial,
		(char*)(void*)Marshal::StringToHGlobalAnsi(type->ToString()), t);
}

FourNodeQuadWrapper::FourNodeQuadWrapper(int tag,
	int Nd1, int Nd2, int Nd3, int Nd4,
	NDMaterialWrapper^ nDMaterial,
	FourNodeQuadType^ type, double t)
{
	_Element = new FourNodeQuad(tag,
		Nd1, Nd2, Nd3, Nd4,
		*nDMaterial->_NDMaterial,
		(char*)(void*)Marshal::StringToHGlobalAnsi(type->ToString()), t);
}

FourNodeQuadWrapper::FourNodeQuadWrapper(int tag,
	int Nd1, int Nd2, int Nd3, int Nd4,
	NDMaterialWrapper^ nDMaterial,
	FourNodeQuadType^ type, double t, double pressure,
	double rho,
	double b1, double b2)
{
	_Element = new FourNodeQuad(tag,
		Nd1, Nd2, Nd3, Nd4,
		*nDMaterial->_NDMaterial,
		(char*)(void*)Marshal::StringToHGlobalAnsi(type->ToString()), t, pressure,
		rho,
		b1, b2);
}

FourNodeQuadWithSensitivityWrapper::FourNodeQuadWithSensitivityWrapper(int tag,
	int Nd1, int Nd2, int Nd3, int Nd4,
	NDMaterialWrapper^ nDMaterial,
	FourNodeQuadType^ type, double t)
{
	_Element = new FourNodeQuadWithSensitivity(tag,
		Nd1, Nd2, Nd3, Nd4,
		*nDMaterial->_NDMaterial,
		(char*)(void*)Marshal::StringToHGlobalAnsi(type->ToString()), t);
}

FourNodeQuadWithSensitivityWrapper::FourNodeQuadWithSensitivityWrapper(int tag,
	int Nd1, int Nd2, int Nd3, int Nd4,
	NDMaterialWrapper^ nDMaterial,
	FourNodeQuadType^ type, double t, double pressure,
	double rho,
	double b1, double b2)
{
	_Element = new FourNodeQuadWithSensitivity(tag,
		Nd1, Nd2, Nd3, Nd4,
		*nDMaterial->_NDMaterial,
		(char*)(void*)Marshal::StringToHGlobalAnsi(type->ToString()), t, pressure,
		rho,
		b1, b2);
}

NineNodeMixedQuadWrapper::NineNodeMixedQuadWrapper(int tag,
	int node1,
	int node2,
	int node3,
	int node4,
	int node5,
	int node6,
	int node7,
	int node8,
	int node9,
	NDMaterialWrapper^ nDMaterial)
{
	_Element = new NineNodeMixedQuad(tag,
		node1,
		node2,
		node3,
		node4,
		node5,
		node6,
		node7,
		node8,
		node9,
		*nDMaterial->_NDMaterial);
}


