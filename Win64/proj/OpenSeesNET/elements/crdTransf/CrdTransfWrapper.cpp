#include "stdafx.h"
#include "CrdTransfWrapper.h"

using namespace OpenSees::Elements::CrdTransfs;
CrdTransfWrapper::CrdTransfWrapper()
{
	 
}

// 2d CrdTransf

LinearCrdTransf2dWrapper::LinearCrdTransf2dWrapper()
{
	_CrdTransf = new LinearCrdTransf2d();
}

CorotCrdTransf2dWrapper::CorotCrdTransf2dWrapper()
{
	_CrdTransf = new CorotCrdTransf2d();
}

PDeltaCrdTransf2dWrapper::PDeltaCrdTransf2dWrapper()
{
	_CrdTransf = new PDeltaCrdTransf2d();
}

// 3d CrdTransf
LinearCrdTransf3dWrapper::LinearCrdTransf3dWrapper(VectorWrapper^ vectorXZPlane)
{
	_CrdTransf = new LinearCrdTransf3d(0, *vectorXZPlane->_Vector);
}

LinearCrdTransf3dWrapper::LinearCrdTransf3dWrapper(VectorWrapper^ vectorXZPlane, VectorWrapper^ rigJntOffsetI, VectorWrapper^ rigJntOffsetJ)
{
	_CrdTransf = new LinearCrdTransf3d(0, *vectorXZPlane->_Vector, *rigJntOffsetI->_Vector, *rigJntOffsetJ->_Vector);
}

CorotCrdTransf3dWrapper::CorotCrdTransf3dWrapper(VectorWrapper^ vectorXZPlane, VectorWrapper^ rigJntOffsetI, VectorWrapper^ rigJntOffsetJ)
{
	_CrdTransf = new CorotCrdTransf3d(0, *vectorXZPlane->_Vector, *rigJntOffsetI->_Vector, *rigJntOffsetJ->_Vector);
}

PDeltaCrdTransf3dWrapper::PDeltaCrdTransf3dWrapper(VectorWrapper^ vectorXZPlane)
{
	_CrdTransf = new PDeltaCrdTransf3d(0, *vectorXZPlane->_Vector);
}

PDeltaCrdTransf3dWrapper::PDeltaCrdTransf3dWrapper(VectorWrapper^ vectorXZPlane, VectorWrapper^ rigJntOffsetI, VectorWrapper^ rigJntOffsetJ)
{
	_CrdTransf = new PDeltaCrdTransf3d(0, *vectorXZPlane->_Vector, *rigJntOffsetI->_Vector, *rigJntOffsetJ->_Vector);
}

CorotCrdTransfWarping2dWrapper::CorotCrdTransfWarping2dWrapper()
{
	_CrdTransf = new CorotCrdTransfWarping2d();
}

CorotCrdTransfWarping2dWrapper::CorotCrdTransfWarping2dWrapper(VectorWrapper^ rigJntOffsetI, VectorWrapper^ rigJntOffsetJ)
{
	_CrdTransf = new CorotCrdTransfWarping2d(0, *rigJntOffsetI->_Vector, *rigJntOffsetI->_Vector);
}

LinearCrdTransf2dIntWrapper::LinearCrdTransf2dIntWrapper()
{
	_CrdTransf = new LinearCrdTransf2dInt(0);
}

LinearCrdTransf2dIntWrapper::LinearCrdTransf2dIntWrapper(VectorWrapper^ rigJntOffsetI, VectorWrapper^ rigJntOffsetJ)
{
	_CrdTransf = new LinearCrdTransf2dInt(0, *rigJntOffsetI->_Vector, *rigJntOffsetI->_Vector);
}

