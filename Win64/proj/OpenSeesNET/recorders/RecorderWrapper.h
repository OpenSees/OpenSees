#pragma once
#include <response\Response.h>
#include <recorder\DamageRecorder.h>
#include <recorder\NodeRecorder.h>
#include <recorder\DriftRecorder.h>
#include <recorder\ElementRecorder.h>
#include <recorder\EnvelopeDriftRecorder.h>
#include <recorder\EnvelopeElementRecorder.h>
#include <recorder\EnvelopeNodeRecorder.h>
#include <recorder\MaxNodeDispRecorder.h>
#include <recorder\EnvelopeDriftRecorder.h>
#include <recorder\NormElementRecorder.h>
#include <recorder\NormEnvelopeElementRecorder.h>
#include <recorder\PatternRecorder.h>


#include "../taggeds/TaggedObjectWrapper.h"
#include "../actors/IMovableObjectWrapper.h"
#include "../matrix/IDWrapper.h"
#include "../domains/domain/BaseDomainWrapper.h"
#include "../domains/timeSeries/TimeSeriesWrapper.h"
#include "../damage/DamageModelWrapper.h"
#include "../handlers/HandlerWrapper.h"


using namespace System;
using namespace OpenSees;
using namespace OpenSees::DamageModels;
using namespace OpenSees::Components;
using namespace OpenSees::Components::Timeseries;
using namespace OpenSees::Handlers;

namespace OpenSees {
	namespace Recorders {
		

		public ref class RecorderWrapper abstract : TaggedObjectWrapper, IMovableObjectWrapper
		{
		public:
			RecorderWrapper() {};
			~RecorderWrapper() {};
		internal:
			Recorder * _Recorder;
		private:

		};

		public ref class DamageRecorderWrapper : RecorderWrapper
		{
		public:
			DamageRecorderWrapper(int elemid, IDWrapper^ secIDs, int dofid, DamageModelWrapper^ dmgPtr,
				BaseDomainWrapper^ theDomainPtr,
				bool echotimeflag, double deltat, OPS_StreamWrapper^ theOutputStream);
			~DamageRecorderWrapper() {};
		};

		public ref class NodeRecorderWrapper : RecorderWrapper
		{
		public:
			NodeRecorderWrapper(IDWrapper^ theDof,
				IDWrapper^ theNodes,
				int sensitivity,
				String^ dataToStore,
				BaseDomainWrapper^ theDomain,
				OPS_StreamWrapper^ theOutputHandler);
			~NodeRecorderWrapper() {};

		};

		public ref class DriftRecorderWrapper : RecorderWrapper
		{
		public:
			DriftRecorderWrapper(int ndI, int ndJ, int dof, int perpDirn,
				BaseDomainWrapper^ theDomain,
				OPS_StreamWrapper^ theHandler,
				bool echoTime,
				double dT);

			DriftRecorderWrapper(int ndI, int ndJ, int dof, int perpDirn,
				BaseDomainWrapper^ theDomain,
				OPS_StreamWrapper^ theHandler);

			DriftRecorderWrapper(IDWrapper^ ndI, IDWrapper^ ndJ, int dof, int perpDirn,
				BaseDomainWrapper^ theDomain,
				OPS_StreamWrapper^ theHandler,
				bool echoTime,
				double dT);

			DriftRecorderWrapper(IDWrapper^ ndI, IDWrapper^ ndJ, int dof, int perpDirn,
				BaseDomainWrapper^ theDomain,
				OPS_StreamWrapper^ theHandler);

			~DriftRecorderWrapper() {};

			

		};

		public ref class ElementRecorderWrapper : RecorderWrapper
		{
		public:
			ElementRecorderWrapper(IDWrapper^ eleID,
				array<String^>^ argv,
				bool echoTime,
				BaseDomainWrapper^ theDomain,
				OPS_StreamWrapper^ theOutputHandler,
				double deltaT,
				IDWrapper^ dof);
			ElementRecorderWrapper(IDWrapper^ eleID,
				array<String^>^ argv,
				bool echoTime,
				BaseDomainWrapper^ theDomain,
				OPS_StreamWrapper^ theOutputHandler);
			~ElementRecorderWrapper() {};

		};

		public ref class EnvelopeDriftRecorderWrapper : RecorderWrapper
		{
		public:
			EnvelopeDriftRecorderWrapper(int ndI, int ndJ, int dof, int perpDirn,
				BaseDomainWrapper^ theDomain,
				OPS_StreamWrapper^ theHandler,
				bool echoTime);
			EnvelopeDriftRecorderWrapper(IDWrapper^ ndI, IDWrapper^ ndJ, int dof, int perpDirn,
				BaseDomainWrapper^ theDomain,
				OPS_StreamWrapper^ theHandler,
				bool echoTime);
			~EnvelopeDriftRecorderWrapper() {};
		};

		public ref class EnvelopeElementRecorderWrapper : RecorderWrapper
		{
		public:
			EnvelopeElementRecorderWrapper(IDWrapper^ eleID,
				array<String^>^ argv,
				BaseDomainWrapper^ theDomain,
				OPS_StreamWrapper^ theHandler,
				double deltaT,
				bool echoTimeFlag,
				IDWrapper^ dof);

			EnvelopeElementRecorderWrapper(IDWrapper^ eleID,
				array<String^>^ argv,
				BaseDomainWrapper^ theDomain,
				OPS_StreamWrapper^ theHandler);

			~EnvelopeElementRecorderWrapper() {};
		};

		public ref class EnvelopeNodeRecorderWrapper : RecorderWrapper
		{
		public:
			EnvelopeNodeRecorderWrapper(IDWrapper^ theDof,
				IDWrapper^ theNodes,
				String^ dataToStore,
				BaseDomainWrapper^ theDomain,
				OPS_StreamWrapper^ theHandler,
				double deltaT,
				bool echoTimeFlag,
				array<TimeSeriesWrapper^>^ theTimeSeries);
			EnvelopeNodeRecorderWrapper(IDWrapper^ theDof,
				IDWrapper^ theNodes,
				String^ dataToStore,
				BaseDomainWrapper^ theDomain,
				OPS_StreamWrapper^ theHandler);
			~EnvelopeNodeRecorderWrapper() {};
		};

		public ref class MaxNodeDispRecorderWrapper : RecorderWrapper
		{
		public:
			MaxNodeDispRecorderWrapper(int dof, IDWrapper^ theNodes, BaseDomainWrapper^ theDomain);
			~MaxNodeDispRecorderWrapper() {};
		};

		public ref class NormElementRecorderWrapper : RecorderWrapper
		{
		public:
			NormElementRecorderWrapper(IDWrapper^ eleID,
				array<String^>^ argv,
				bool echoTimeFlag,
				BaseDomainWrapper^ theDomain,
				OPS_StreamWrapper^ theHandler,
				double deltaT,
				IDWrapper^ dof);

			NormElementRecorderWrapper(IDWrapper^ eleID,
				array<String^>^ argv,
				bool echoTimeFlag,
				BaseDomainWrapper^ theDomain,
				OPS_StreamWrapper^ theHandler);

			~NormElementRecorderWrapper() {};
		};

		public ref class NormEnvelopeElementRecorderWrapper : RecorderWrapper
		{
		public:
			NormEnvelopeElementRecorderWrapper(IDWrapper^ eleID,
				array<String^>^ argv,
				BaseDomainWrapper^ theDomain,
				OPS_StreamWrapper^ theHandler,
				double deltaT,
				bool echoTimeFlag,
				IDWrapper^ dof);

			NormEnvelopeElementRecorderWrapper(IDWrapper^ eleID,
				array<String^>^ argv,
				BaseDomainWrapper^ theDomain,
				OPS_StreamWrapper^ theHandler);

			~NormEnvelopeElementRecorderWrapper() {};
		};

		public ref class PatternRecorderWrapper : RecorderWrapper
		{
		public:
			PatternRecorderWrapper(int thePattern,
				BaseDomainWrapper^ theDomain,
				String^ argv,
				double deltaT,
				int startFlag);

			~PatternRecorderWrapper() {};
		};
	}
}