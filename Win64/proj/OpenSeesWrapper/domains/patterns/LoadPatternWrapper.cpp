#include "stdafx.h"
#include <LoadPattern.h>
#include "LoadPatternWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Components;
using namespace OpenSees::Components::Timeseries;
using namespace OpenSees::Components::GroundMotions;
using namespace OpenSees::Components::LoadPatterns;
LoadPatternWrapper::LoadPatternWrapper(int tag)
{
	_LoadPattern = new LoadPattern(tag);
}

LoadPatternWrapper::LoadPatternWrapper()
{
	
}

void
LoadPatternWrapper::SetTimeSeries(TimeSeriesWrapper^ timeSeries)
{
	TimeSeries* _timeSeries = timeSeries->_TimeSeries;
	_LoadPattern->setTimeSeries(_timeSeries);
}

EarthquakePatternWrapper::EarthquakePatternWrapper(int tag)
	:LoadPatternWrapper(tag)
{

}

UniformExcitationWrapper::UniformExcitationWrapper(int tag
	, GroundMotionWrapper^ groudMotion, int dof, double vel0, double fact)
	:EarthquakePatternWrapper(tag)
{
	_LoadPattern = new UniformExcitation(*groudMotion->_GroundMotion, dof, tag, vel0, fact);
}

FireLoadPatternWrapper::FireLoadPatternWrapper(int tag)
	:LoadPatternWrapper(tag)
{
	_LoadPattern = new FireLoadPattern(tag);
}

void
FireLoadPatternWrapper::SetFireTimeSeries(TimeSeriesWrapper^ theSeries1, TimeSeriesWrapper^ theSeries2,
	TimeSeriesWrapper^ theSeries3, TimeSeriesWrapper^ theSeries4, TimeSeriesWrapper^ theSeries5,
	TimeSeriesWrapper^ theSeries6, TimeSeriesWrapper^ theSeries7, TimeSeriesWrapper^ theSeries8, TimeSeriesWrapper^ theSeries9) {

	((FireLoadPattern*)_LoadPattern)->setFireTimeSeries(theSeries1->_TimeSeries, theSeries2->_TimeSeries,
		theSeries3->_TimeSeries, theSeries4->_TimeSeries, theSeries5->_TimeSeries,
		theSeries6->_TimeSeries, theSeries7->_TimeSeries, theSeries8->_TimeSeries, theSeries9->_TimeSeries);
}

MultiSupportPatternWrapper::MultiSupportPatternWrapper(int tag)
	:LoadPatternWrapper(tag)
{
	_LoadPattern = new MultiSupportPattern(tag);
}