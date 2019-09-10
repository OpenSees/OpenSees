#include "stdafx.h"

#include <GroundMotion.h>
#include <TimeSeries.h>
#include "GroudMotionWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Components::GroundMotions;
using namespace System::Runtime::InteropServices;

GroundMotionWrapper::GroundMotionWrapper()
{

}

GroundMotionWrapper::GroundMotionWrapper(TimeSeriesWrapper^ dispSeries,
	TimeSeriesWrapper^ velSeries,
	TimeSeriesWrapper^ accelSeries,
	TimeSeriesIntegratorWrapper^ theIntegrator,
	double dTintegration, double fact)
{
	TimeSeries* _dispSeries = 0;
	TimeSeries* _velSeries = 0;
	TimeSeries* _accelSeries = 0;
	TimeSeriesIntegrator* _theIntegrator = 0;

	if (dispSeries != nullptr)
		_dispSeries = dispSeries->_TimeSeries;

	if (velSeries != nullptr)
		_velSeries = velSeries->_TimeSeries;

	if (accelSeries != nullptr)
		_accelSeries = accelSeries->_TimeSeries;

	if (theIntegrator != nullptr)
		_theIntegrator = theIntegrator->_TimeSeriesIntegrator;

	_GroundMotion = new GroundMotion(_dispSeries, _velSeries, _accelSeries, _theIntegrator, dTintegration, fact);
}

GroundMotionWrapper::~GroundMotionWrapper()
{
	if (_GroundMotion != 0)
		delete _GroundMotion;

}

GroundMotionRecordWrapper::GroundMotionRecordWrapper()
{

}

GroundMotionRecordWrapper::GroundMotionRecordWrapper(String^ fileNameAccel, double timeStep, double fact, double dTintegration)
{
	char* accl = (char*)(void*)Marshal::StringToHGlobalAnsi(fileNameAccel);
	_GroundMotion = new GroundMotionRecord(accl, timeStep, fact, dTintegration);

}

GroundMotionRecordWrapper::GroundMotionRecordWrapper(String^ fileNameAccel, String^ fileNameTime, double timeStep, double fact, double dTintegration)
{
	char* accl = (char*)(void*)Marshal::StringToHGlobalAnsi(fileNameAccel);
	char* time = (char*)(void*)Marshal::StringToHGlobalAnsi(fileNameTime);
	_GroundMotion = new GroundMotionRecord(accl, time, fact, dTintegration);
}

GroundMotionRecordWrapper::~GroundMotionRecordWrapper()
{
	if (_GroundMotion != 0)
		delete _GroundMotion;
}


InterpolatedGroundMotionWrapper::InterpolatedGroundMotionWrapper()
{

}

InterpolatedGroundMotionWrapper::InterpolatedGroundMotionWrapper(array<GroundMotionWrapper^>^ groundMotions, VectorWrapper^  fact, bool destroyMotions, double deltaPeak)
{
	GroundMotion** gms = new GroundMotion*[groundMotions->Length];
	for (int i = 0; i < groundMotions->Length; i++)
	{
		gms[i] = groundMotions[i]->_GroundMotion;
	}
	_GroundMotion = new InterpolatedGroundMotion(gms, *fact->_Vector, destroyMotions, deltaPeak);
}

InterpolatedGroundMotionWrapper::~InterpolatedGroundMotionWrapper()
{
	if (_GroundMotion != 0)
		delete _GroundMotion;
}

