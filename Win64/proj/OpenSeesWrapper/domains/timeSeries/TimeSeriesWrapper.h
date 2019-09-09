#pragma once
#include <TimeSeries.h>
#include <TimeSeriesIntegrator.h>
#include <LinearSeries.h>
#include <PathSeries.h>
#include <PathTimeSeries.h>
#include <TimeSeriesIntegrator.h>
#include <TrapezoidalTimeSeriesIntegrator.h>
#include <ConstantSeries.h>
//#include <DiscretizedRandomProcessSeries.h>
#include <PathTimeSeriesThermal.h>
#include <PeerMotion.h>
#include <PulseSeries.h>
#include <RectangularSeries.h>
#include <SimpsonTimeSeriesIntegrator.h>
//#include <SimulatedRandomProcessSeries.h>
#include <TriangleSeries.h>
#include <TrigSeries.h>


#include "../../taggeds/TaggedObjectWrapper.h"
#include "../../actors/IMovableObjectWrapper.h"
#include "../../matrix/VectorWrapper.h"

using namespace System;
namespace OpenSees {
	namespace Components {
		namespace Timeseries {
			public ref class TimeSeriesWrapper : TaggedObjectWrapper, IMovableObjectWrapper
			{
			public:
				TimeSeriesWrapper();
				double GetFactor(double time) {
					return _TimeSeries->getFactor(time);
				}
				~TimeSeriesWrapper() {
					if (_TimeSeries != 0)
						delete _TimeSeries;
				};

			internal:
				TimeSeries * _TimeSeries;
			};

			public ref class PathTimeSeriesWrapper : TimeSeriesWrapper
			{
			public:
				PathTimeSeriesWrapper(VectorWrapper^ thePath,
					VectorWrapper^ theTime,
					double cfactor,
					bool useLast);
				~PathTimeSeriesWrapper() {
					if (_TimeSeries != 0)
						delete _TimeSeries;
				};

			internal:

			};

			public ref class PathSeriesWrapper : TimeSeriesWrapper
			{
			public:
				PathSeriesWrapper(VectorWrapper^ theLoadPath,
					double theTimeIncr,
					double theFactor,
					bool useLast,
					bool prependZero,
					double tStart);
				~PathSeriesWrapper() {
					if (_TimeSeries != 0)
						delete _TimeSeries;
				};

			internal:

			};

			public ref class LinearSeriesWrapper : TimeSeriesWrapper
			{
			public:
				LinearSeriesWrapper();
				~LinearSeriesWrapper() {
					if (_TimeSeries != 0)
						delete _TimeSeries;
				};

			internal:

			};

			public ref class ConstantSeriesWrapper : TimeSeriesWrapper
			{
			public:
				ConstantSeriesWrapper(double cFactor);
				~ConstantSeriesWrapper() {
					if (_TimeSeries != 0)
						delete _TimeSeries;
				};

			internal:

			};

			public ref class PathTimeSeriesThermalWrapper : TimeSeriesWrapper
			{
			public:
				PathTimeSeriesThermalWrapper(String^ fileName, int DataNum, bool tempOut, double cfactor);
				PathTimeSeriesThermalWrapper(int DataNum, bool tempOut, double cfactor);
				~PathTimeSeriesThermalWrapper() {
					if (_TimeSeries != 0)
						delete _TimeSeries;
				};

			internal:

			};

			public ref class PeerMotionWrapper : TimeSeriesWrapper
			{
			public:
				PeerMotionWrapper(String^ earthquake,
					String^ station,
					String^ responseType,
					double cfactor);
				PeerMotionWrapper() {};
				~PeerMotionWrapper() {
					if (_TimeSeries != 0)
						delete _TimeSeries;
				};

			internal:

			};

			public ref class PulseSeriesWrapper : TimeSeriesWrapper
			{
			public:
				PulseSeriesWrapper(double tStart,
					double tFinish,
					double period,
					double pWidth,
					double phaseShift,
					double cFactor,
					double zeroShift);
				PulseSeriesWrapper() {};
				~PulseSeriesWrapper() {
					if (_TimeSeries != 0)
						delete _TimeSeries;
				};

			internal:

			};

			public ref class RectangularSeriesWrapper : TimeSeriesWrapper
			{
			public:
				RectangularSeriesWrapper(double tStart, double tFinal, double cFactor);
				RectangularSeriesWrapper() {};
				~RectangularSeriesWrapper() {
					if (_TimeSeries != 0)
						delete _TimeSeries;
				};

			internal:

			};

			public ref class TriangleSeriesWrapper : TimeSeriesWrapper
			{
			public:
				TriangleSeriesWrapper(double tStart,
					double tFinish,
					double period,
					double phaseShift,
					double cFactor,
					double zeroShift);
				TriangleSeriesWrapper() {};
				~TriangleSeriesWrapper() {
					if (_TimeSeries != 0)
						delete _TimeSeries;
				};

			internal:

			};

			public ref class TrigSeriesWrapper : TimeSeriesWrapper
			{
			public:
				TrigSeriesWrapper(double tStart,
					double tFinish,
					double period,
					double phaseShift,
					double cFactor,
					double zeroShift);
				TrigSeriesWrapper() {};
				~TrigSeriesWrapper() {
					if (_TimeSeries != 0)
						delete _TimeSeries;
				};

			internal:

			};

			// Integrators
			public ref class TimeSeriesIntegratorWrapper abstract : IMovableObjectWrapper
			{
			public:
				TimeSeriesIntegratorWrapper();
				~TimeSeriesIntegratorWrapper() {
					if (_TimeSeriesIntegrator != 0)
						delete _TimeSeriesIntegrator;
				};

			internal:
				TimeSeriesIntegrator * _TimeSeriesIntegrator;
			};

			public ref class TrapezoidalTimeSeriesIntegratorWrapper : TimeSeriesIntegratorWrapper
			{
			public:
				TrapezoidalTimeSeriesIntegratorWrapper();
				~TrapezoidalTimeSeriesIntegratorWrapper() {
					if (_TimeSeriesIntegrator != 0)
						delete _TimeSeriesIntegrator;
				};

			internal:

			};

			public ref class SimpsonTimeSeriesIntegratorWrapper : TimeSeriesIntegratorWrapper
			{
			public:
				SimpsonTimeSeriesIntegratorWrapper();
				~SimpsonTimeSeriesIntegratorWrapper() {
					if (_TimeSeriesIntegrator != 0)
						delete _TimeSeriesIntegrator;
				};

			internal:

			};
		}
	}
}
