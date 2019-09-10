#pragma once

#include <GroundMotion.h>
#include <GroundMotionRecord.h>
#include <InterpolatedGroundMotion.h>
#include "../timeSeries/TimeSeriesWrapper.h"
#include "../../matrix/VectorWrapper.h"
#include "../../actors/IMovableObjectWrapper.h"
using namespace System;
using namespace OpenSees::Components::Timeseries;
namespace OpenSees {
	namespace Components {
		namespace GroundMotions {
		public ref class GroundMotionWrapper : IMovableObjectWrapper
		{
		public:
			GroundMotionWrapper();
			GroundMotionWrapper(TimeSeriesWrapper^ dispSeries,
				TimeSeriesWrapper^ velSeries,
				TimeSeriesWrapper^ accelSeries,
				TimeSeriesIntegratorWrapper^ theIntegrator,
				double dTintegration, double fact);

			double GetDuration() {
				return _GroundMotion->getDuration();
			}

			double GetPeakAccel() {
				return _GroundMotion->getPeakAccel();
			}

			double GetPeakVel() {
				return _GroundMotion->getPeakVel();
			}

			double GetPeakDisp() {
				return _GroundMotion->getPeakDisp();
			}

			double GetAccel(double time) {
				return _GroundMotion->getAccel(time);
			}

			double GetVel(double time) {
				return _GroundMotion->getVel(time);
			}

			double GetDisp(double time) {
				return _GroundMotion->getDisp(time);
			}

			void SetIntegrator(TimeSeriesIntegratorWrapper^ timeSeriesIntegrator) {
				_GroundMotion->setIntegrator(timeSeriesIntegrator->_TimeSeriesIntegrator);
			}

			TimeSeriesWrapper^ Integrate(TimeSeriesWrapper^ timeSeries, double delta) {
				TimeSeries* ret = _GroundMotion->integrate(timeSeries->_TimeSeries, delta);
				TimeSeriesWrapper^ wret = gcnew TimeSeriesWrapper();
				wret->_TimeSeries = ret;
				return wret;
			}

			~GroundMotionWrapper();

		internal:
			GroundMotion * _GroundMotion;
		};

		public ref class GroundMotionRecordWrapper : GroundMotionWrapper
		{
		public:
			GroundMotionRecordWrapper();
			GroundMotionRecordWrapper(String^ fileNameAccel, double timeStep,double fact, double dTintegration);
			GroundMotionRecordWrapper(String^ fileNameAccel,String^ fileNameTime, double timeStep, double fact, double dTintegration);
			~GroundMotionRecordWrapper();

		internal:
		};

		public ref class InterpolatedGroundMotionWrapper : GroundMotionWrapper
		{
		public:
			InterpolatedGroundMotionWrapper();
			InterpolatedGroundMotionWrapper(array<GroundMotionWrapper^>^ groundMotions, VectorWrapper^  fact, bool destroyMotions, double deltaPeak);
			~InterpolatedGroundMotionWrapper();

		internal:
		};
		}
	}
}
