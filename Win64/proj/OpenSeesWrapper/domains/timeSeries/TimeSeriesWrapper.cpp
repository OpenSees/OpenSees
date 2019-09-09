#include "stdafx.h"
#include "TimeSeriesWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Components::Timeseries;
using namespace System;
using namespace System::Runtime::InteropServices;

TimeSeriesWrapper::TimeSeriesWrapper()
{
	
}

LinearSeriesWrapper::LinearSeriesWrapper()
{
	this->_TimeSeries = new  LinearSeries();
}

PathSeriesWrapper::PathSeriesWrapper(VectorWrapper^ theLoadPath,
	double theTimeIncr,
	double theFactor,
	bool useLast,
	bool prependZero,
	double tStart)
{
	this->_TimeSeries = new  PathSeries(0, *theLoadPath->_Vector, theTimeIncr, theFactor, useLast, prependZero, tStart);
}

PathTimeSeriesWrapper::PathTimeSeriesWrapper(VectorWrapper^ thePath,
	VectorWrapper^ theTime,
	double cfactor,
	bool useLast)
{
	this->_TimeSeries = new  PathTimeSeries(0, *thePath->_Vector, *theTime->_Vector, cfactor, useLast);
}

ConstantSeriesWrapper::ConstantSeriesWrapper(double cFactor)
{
	this->_TimeSeries = new  ConstantSeries(0, cFactor);
}

PathTimeSeriesThermalWrapper::PathTimeSeriesThermalWrapper(String^ fileName, int DataNum, bool tempOut, double cfactor)
{
	char* fname = (char*)(void*)Marshal::StringToHGlobalAnsi(fileName);
	this->_TimeSeries = new  PathTimeSeriesThermal(0, fname, DataNum, tempOut, cfactor);
}

PathTimeSeriesThermalWrapper::PathTimeSeriesThermalWrapper(int DataNum, bool tempOut, double cfactor)
{
	this->_TimeSeries = new  PathTimeSeriesThermal(0,  DataNum, tempOut, cfactor);
}

PeerMotionWrapper::PeerMotionWrapper(String^ earthquake,
	String^ station,
	String^ responseType,
	double cfactor)
{
	char* cearthquake = (char*)(void*)Marshal::StringToHGlobalAnsi(earthquake);
	char* cstation = (char*)(void*)Marshal::StringToHGlobalAnsi(station);
	char* cresponseType = (char*)(void*)Marshal::StringToHGlobalAnsi(responseType);

	this->_TimeSeries = new  PeerMotion(0, cearthquake, cstation, cresponseType, cfactor);
}

PulseSeriesWrapper::PulseSeriesWrapper(double tStart,
	double tFinish,
	double period,
	double pWidth,
	double phaseShift,
	double cFactor,
	double zeroShift)
{
	this->_TimeSeries = new  PulseSeries(0,  tStart,
		 tFinish,
		 period,
		 pWidth,
		 phaseShift,
		 cFactor,
		 zeroShift);
}

RectangularSeriesWrapper::RectangularSeriesWrapper(double tStart, double tFinal, double cFactor)
{
	this->_TimeSeries = new  RectangularSeries(0, tStart, tFinal, cFactor);
}

TriangleSeriesWrapper::TriangleSeriesWrapper(double tStart,
	double tFinish,
	double period,
	double phaseShift,
	double cFactor,
	double zeroShift)
{
	this->_TimeSeries = new  TriangleSeries(0,  tStart,
		 tFinish,
		 period,
		 phaseShift,
		 cFactor,
		 zeroShift);
}

TrigSeriesWrapper::TrigSeriesWrapper(double tStart,
	double tFinish,
	double period,
	double phaseShift,
	double cFactor,
	double zeroShift)
{
	this->_TimeSeries = new  TrigSeries(0, tStart,
		tFinish,
		period,
		phaseShift,
		cFactor,
		zeroShift);
}

// Integrators
TimeSeriesIntegratorWrapper::TimeSeriesIntegratorWrapper()
{

}

TrapezoidalTimeSeriesIntegratorWrapper::TrapezoidalTimeSeriesIntegratorWrapper()
{
	_TimeSeriesIntegrator = new TrapezoidalTimeSeriesIntegrator();
}

SimpsonTimeSeriesIntegratorWrapper::SimpsonTimeSeriesIntegratorWrapper()
{
	_TimeSeriesIntegrator = new SimpsonTimeSeriesIntegrator();
}

