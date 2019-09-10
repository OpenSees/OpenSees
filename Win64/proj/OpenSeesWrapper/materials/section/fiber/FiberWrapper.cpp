#include "stdafx.h"
#include "FiberWrapper.h"

using namespace OpenSees::Materials::Sections::Repres;

FiberWrapper::FiberWrapper()
{

}

UniaxialFiber2dWrapper::UniaxialFiber2dWrapper(int tag, UniaxialMaterialWrapper^ theMat, double area, double position)
{
	_Fiber = new UniaxialFiber2d(tag, *theMat->_UniaxialMaterial, area, position);
}

UniaxialFiber3dWrapper::UniaxialFiber3dWrapper(int tag, UniaxialMaterialWrapper^ theMat, double area,
	VectorWrapper^ position)
{
	_Fiber = new UniaxialFiber3d(tag, *theMat->_UniaxialMaterial, area, *position->_Vector);
}

NDFiber2dWrapper::NDFiber2dWrapper(int tag, NDMaterialWrapper^ theMat, double Area, double position) {
	_Fiber = new NDFiber2d(tag, *theMat->_NDMaterial, Area, position);
}

NDFiber3dWrapper::NDFiber3dWrapper(int tag, NDMaterialWrapper^ theMat, double Area, double y, double z) {
	_Fiber = new NDFiber3d(tag, *theMat->_NDMaterial, Area, y, z);
}

