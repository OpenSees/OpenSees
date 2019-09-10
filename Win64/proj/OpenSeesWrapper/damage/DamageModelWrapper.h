#pragma once
#include <DamageModel.h>
#include <HystereticEnergy.h>
#include <Kratzig.h>
#include <Mehanny.h>
#include <NormalizedPeak.h>
#include <ParkAng.h>
#include "../actors/IMovableObjectWrapper.h"
#include "../taggeds/TaggedObjectWrapper.h"
namespace OpenSees {
	namespace DamageModels {
		public ref class DamageModelWrapper :TaggedObjectWrapper, IMovableObjectWrapper
		{
		public:
			DamageModelWrapper();
			~DamageModelWrapper();
		internal:
			DamageModel * _DamageModel;
		private:

		};

		public ref class HystereticEnergyWrapper :DamageModelWrapper
		{
		public:
			HystereticEnergyWrapper(int tag, double Etot, double Cpow);
			~HystereticEnergyWrapper();
		internal:

		private:

		};

		public ref class KratzigWrapper :DamageModelWrapper
		{
		public:
			KratzigWrapper(int tag, double ultimatePosVal, double ultimateNegVal);
			~KratzigWrapper();
		internal:

		private:

		};

		public ref class MehannyWrapper :DamageModelWrapper
		{
		public:
			MehannyWrapper(int tag, double alpha, double beta, double gamma,
				double ultimatePosValue, double ultimateNegValue, double abstol, double reltol, double posmodifier, double negmodifier);
			~MehannyWrapper();
		internal:

		private:

		};

		public enum class DamageTypeWrapper {
			Force,
			Deformation,
			PlasticDefo,
			TotalEnergy,
			PlasticEnergy,
		};

		public ref class NormalizedPeakWrapper :DamageModelWrapper
		{
		public:
			NormalizedPeakWrapper(int tag, double maxVal, double minVal, DamageTypeWrapper^ type);
			~NormalizedPeakWrapper();
		internal:

		private:

		};


		public ref class ParkAngWrapper :DamageModelWrapper
		{
		public:
			ParkAngWrapper(int tag, double deltaU, double beta, double sigmaY);
			~ParkAngWrapper();
		internal:

		private:

		};
	}
}
