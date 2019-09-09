#include "stdafx.h"
#include "LoadWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Components::Loads;

LoadWrapper::LoadWrapper(int tag)
	:DomainComponentWrapper(tag)
{

}

LoadWrapper::~LoadWrapper()
{
	if (_Load != 0)
		delete _Load;
}

NodalLoadWrapper::NodalLoadWrapper(int tag, int node, VectorWrapper^ values, bool isLoadConstant)
	:LoadWrapper(0)
{
	Vector* _vector = values->_Vector;
	_NodalLoad = new NodalLoad(tag, node, *_vector, isLoadConstant);
}

Beam2dPartialUniformLoadWrapper::Beam2dPartialUniformLoadWrapper(int tag, double wTrans, double wAxial, double aL, double bL, int eleTag)
	:ElementalLoadWrapper(tag)
{
	_ElementalLoad = new Beam2dPartialUniformLoad(tag, wTrans, wAxial, aL, bL, eleTag);
}

Beam2dPointLoadWrapper::Beam2dPointLoadWrapper(int tag, double Pt, double x, int eleTag, double Pa)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new Beam2dPointLoad(tag, Pt, x, eleTag, Pa);
}

Beam2dTempLoadWrapper::Beam2dTempLoadWrapper(int tag, double Ttop1, double Tbot1, double Ttop2,
	double Tbot2,
	int theElementTag)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new Beam2dTempLoad(tag, Ttop1, Tbot1, Ttop2, Tbot2, theElementTag);
}

Beam2dTempLoadWrapper::Beam2dTempLoadWrapper(int tag, double Tuniform,
	int theElementTag)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new Beam2dTempLoad(tag, Tuniform, theElementTag);
}

Beam2dTempLoadWrapper::Beam2dTempLoadWrapper(int tag, double Ttop, double Tbot,
	int theElementTag)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new Beam2dTempLoad(tag, Ttop, Tbot, theElementTag);
}

Beam2dTempLoadWrapper::Beam2dTempLoadWrapper(int tag, int theElementTag)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new Beam2dTempLoad(tag, theElementTag);
}

Beam2dThermalActionWrapper::Beam2dThermalActionWrapper(int tag, double t1, double locY1, double t2, double locY2,
	double t3, double locY3, double t4, double locY4,
	double t5, double locY5, double t6, double locY6,
	double t7, double locY7, double t8, double locY8,
	double t9, double locY9,
	int theElementTag)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new Beam2dThermalAction(tag, t1, locY1, t2, locY2,
		t3, locY3, t4, locY4,
		t5, locY5, t6, locY6,
		t7, locY7, t8, locY8,
		t9, locY9,
		theElementTag);
}

Beam2dThermalActionWrapper::Beam2dThermalActionWrapper(int tag, double locY1, double locY2,
	TimeSeriesWrapper^ theSeries, int theElementTag)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new Beam2dThermalAction(tag, locY1, locY2,
		theSeries->_TimeSeries, theElementTag);
}

Beam2dThermalActionWrapper::Beam2dThermalActionWrapper(int tag, VectorWrapper^ locs,
	TimeSeriesWrapper^ theSeries, int theElementTag)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new Beam2dThermalAction(tag, *locs->_Vector,
		theSeries->_TimeSeries, theElementTag);
}

Beam2dThermalActionWrapper::Beam2dThermalActionWrapper(int tag, int theElementTag)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new Beam2dThermalAction(tag, theElementTag);
}

Beam2dUniformLoadWrapper::Beam2dUniformLoadWrapper(int tag, double wTrans, double wAxial, int eleTag)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new Beam2dUniformLoad(tag, wTrans, wAxial, eleTag);
}

Beam3dPointLoadWrapper::Beam3dPointLoadWrapper(int tag, double Py, double Pz, double x, int eleTag, double Pa)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new Beam3dPointLoad(tag, Py, Pz, x, eleTag, Pa);
}

Beam3dThermalActionWrapper::Beam3dThermalActionWrapper(int tag, double t1, double locY1, double t2, double locY2,
	double t3, double locY3, double t4, double locY4,
	double t5, double locY5, double t6, double t7, double locZ1,
	double t8, double t9, double locZ2, double t10, double t11, double locZ3,
	double t12, double t13, double locZ4, double t14, double t15, double locZ5,
	int theElementTag)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new Beam3dThermalAction(tag, t1, locY1, t2, locY2,
		t3, locY3, t4, locY4,
		t5, locY5, t6, t7, locZ1,
		t8, t9, locZ2, t10, t11, locZ3,
		t12, t13, locZ4, t14, t15, locZ5,
		theElementTag);
}


Beam3dThermalActionWrapper::Beam3dThermalActionWrapper(int tag, double t1, double locY1, double t2, double locY2,
	double t3, double locY3, double t4, double locY4,
	double t5, double locY5, double t6, double locY6,
	double t7, double locY7, double t8, double locY8,
	double t9, double locY9,
	int theElementTag)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new Beam3dThermalAction(tag, t1, locY1, t2, locY2,
		t3, locY3, t4, locY4,
		t5, locY5, t6, locY6,
		t7, locY7, t8, locY8,
		t9, locY9,
		theElementTag);
}

Beam3dThermalActionWrapper::Beam3dThermalActionWrapper(int tag, double locY1, double locY2, double locZ1, double Z2,
	TimeSeriesWrapper^ theSeries,
	int theElementTag)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new Beam3dThermalAction(tag, locY1, locY2, locZ1, Z2,
		theSeries->_TimeSeries,
		theElementTag);
}


Beam3dThermalActionWrapper::Beam3dThermalActionWrapper(int tag, VectorWrapper^ locs, TimeSeriesWrapper^ theSeries,
	int theElementTag)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new Beam3dThermalAction(tag, *locs->_Vector, theSeries->_TimeSeries,
		theElementTag);
}

Beam3dThermalActionWrapper::Beam3dThermalActionWrapper(int tag, int theElementTag)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new Beam3dThermalAction(tag, theElementTag);
}

Beam3dUniformLoadWrapper::Beam3dUniformLoadWrapper(int tag, double wy, double wz, double wx, int eleTag)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new Beam3dUniformLoad(tag, wy, wz, wx, eleTag);
}

BrickSelfWeightWrapper::BrickSelfWeightWrapper(int tag, int eleTag)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new BrickSelfWeight(tag, eleTag);
}

ShellThermalActionWrapper::ShellThermalActionWrapper(int tag, double t1, double locY1, double t2, double locY2,
	double t3, double locY3, double t4, double locY4,
	double t5, double locY5, double t6, double locY6,
	double t7, double locY7, double t8, double locY8,
	double t9, double locY9,
	int theElementTag)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new ShellThermalAction(tag, t1, locY1, t2, locY2,
		t3, locY3, t4, locY4,
		t5, locY5, t6, locY6,
		t7, locY7, t8, locY8,
		t9, locY9,
		theElementTag);
}

ShellThermalActionWrapper::ShellThermalActionWrapper(int tag, double t1, double locY1, double t2, double locY2,
	double t3, double locY3, double t4, double locY4,
	double t5, double locY5, int theElementTag)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new ShellThermalAction(tag, t1, locY1, t2, locY2,
		t3, locY3, t4, locY4,
		t5, locY5, theElementTag);
}

ShellThermalActionWrapper::ShellThermalActionWrapper(int tag, double t1, double locY1, double t2, double locY2,
	int theElementTag)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new ShellThermalAction(tag, t1, locY1, t2, locY2,
		theElementTag);
}

ShellThermalActionWrapper::ShellThermalActionWrapper(int tag, double locY1, double locY2,
	TimeSeriesWrapper^ theSeries, int theElementTag)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new ShellThermalAction(tag, locY1, locY2,
		theSeries->_TimeSeries, theElementTag);
}

ShellThermalActionWrapper::ShellThermalActionWrapper(int tag, int theElementTag)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new ShellThermalAction(tag, theElementTag);
}


NodalThermalActionWrapper::NodalThermalActionWrapper(int tag, int theNodeTag,
	VectorWrapper^ locy, TimeSeriesWrapper^ theSeries, VectorWrapper^ crds)
	: NodalLoadWrapper(tag)
{
	_NodalLoad = new NodalThermalAction(tag, theNodeTag,
		*locy->_Vector, theSeries->_TimeSeries, crds->_Vector);
}

NodalThermalActionWrapper::NodalThermalActionWrapper(int tag, int theNodeTag,
	double locY1, double locY2, double locZ1, double locZ2,
	TimeSeriesWrapper^ theSeries, VectorWrapper^ crds)
	: NodalLoadWrapper(tag)
{
	_NodalLoad = new NodalThermalAction(tag, theNodeTag,
		locY1, locY2, locZ1, locZ2,
		theSeries->_TimeSeries, crds->_Vector);
}

NodalThermalActionWrapper::NodalThermalActionWrapper(int tag, int theNodeTag,
	double t1, double locY1, double t2, double locY2, VectorWrapper^ crds)
	: NodalLoadWrapper(tag)
{
	_NodalLoad = new NodalThermalAction(tag, theNodeTag,
		t1, locY1, t2, locY2, crds->_Vector);
}

SelfWeightWrapper::SelfWeightWrapper(int tag, double xFact, double yFact, double zFact, int eleTag)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new SelfWeight(tag, xFact, yFact, zFact, eleTag);
}

SurfaceLoaderWrapper::SurfaceLoaderWrapper(int tag, int eleTag)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new SurfaceLoader(tag, eleTag);
}

ThermalActionWrapperWrapper::ThermalActionWrapperWrapper(int tag, int EleTag,
	NodalThermalActionWrapper^ theNodalTA1, NodalThermalActionWrapper^ theNodalTA2)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new ThermalActionWrapper(tag, EleTag,
		(NodalThermalAction*)theNodalTA1->_NodalLoad, (NodalThermalAction*)theNodalTA2->_NodalLoad);
}

ThermalActionWrapperWrapper::ThermalActionWrapperWrapper(int tag, int EleTag,
	NodalThermalActionWrapper^ theNodalTA1, NodalThermalActionWrapper^ theNodalTA2,
	NodalThermalActionWrapper^ theNodalTA3)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new ThermalActionWrapper(tag, EleTag,
		(NodalThermalAction*)theNodalTA1->_NodalLoad, (NodalThermalAction*)theNodalTA2->_NodalLoad,
		(NodalThermalAction*)theNodalTA3->_NodalLoad);
}


ThermalActionWrapperWrapper::ThermalActionWrapperWrapper(int tag, int EleTag,
	NodalThermalActionWrapper^ theNodalTA1, NodalThermalActionWrapper^ theNodalTA2,
	NodalThermalActionWrapper^ theNodalTA3, NodalThermalActionWrapper^ theNodalTA4)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new ThermalActionWrapper(tag, EleTag,
		(NodalThermalAction*)theNodalTA1->_NodalLoad, (NodalThermalAction*)theNodalTA2->_NodalLoad,
		(NodalThermalAction*)theNodalTA3->_NodalLoad, (NodalThermalAction*)theNodalTA4->_NodalLoad);
}

ThermalActionWrapperWrapper::ThermalActionWrapperWrapper(int tag, int EleTag,
	NodalThermalActionWrapper^ theNodalTA1, NodalThermalActionWrapper^ theNodalTA2,
	NodalThermalActionWrapper^ theNodalTA3, NodalThermalActionWrapper^ theNodalTA4,
	NodalThermalActionWrapper^ theNodalTA5)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new ThermalActionWrapper(tag, EleTag,
		(NodalThermalAction*)theNodalTA1->_NodalLoad, (NodalThermalAction*)theNodalTA2->_NodalLoad,
		(NodalThermalAction*)theNodalTA3->_NodalLoad, (NodalThermalAction*)theNodalTA4->_NodalLoad
		, (NodalThermalAction*)theNodalTA5->_NodalLoad);
}

ThermalActionWrapperWrapper::ThermalActionWrapperWrapper(int tag, int EleTag,
	NodalThermalActionWrapper^ theNodalTA1, NodalThermalActionWrapper^ theNodalTA2,
	NodalThermalActionWrapper^ theNodalTA3, NodalThermalActionWrapper^ theNodalTA4,
	NodalThermalActionWrapper^ theNodalTA5, NodalThermalActionWrapper^ theNodalTA6)
	: ElementalLoadWrapper(tag)
{
	_ElementalLoad = new ThermalActionWrapper(tag, EleTag,
		(NodalThermalAction*)theNodalTA1->_NodalLoad, (NodalThermalAction*)theNodalTA2->_NodalLoad,
		(NodalThermalAction*)theNodalTA3->_NodalLoad, (NodalThermalAction*)theNodalTA4->_NodalLoad,
		(NodalThermalAction*)theNodalTA5->_NodalLoad, (NodalThermalAction*)theNodalTA6->_NodalLoad);
}



