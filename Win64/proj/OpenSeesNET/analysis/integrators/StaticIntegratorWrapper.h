#pragma once

#include <StaticIntegrator.h>
#include <LoadControl.h>
#include <DisplacementControl.h>
#include <ArcLength.h>
#include <ArcLength1.h>
#include <HSConstraint.h>
#include <LoadPath.h>
#include <MinUnbalDispNorm.h>
#include <EQPath.h>

#include "IncrementalIntegratorWrapper.h"
#include "../../domains/domain/DomainWrapper.h"
#include "../../matrix/VectorWrapper.h"

namespace OpenSees {
	namespace Integrators {
		namespace Static {
			public ref class StaticIntegratorWrapper abstract : IncrementalIntegratorWrapper
			{
			public:
				StaticIntegratorWrapper();
				~StaticIntegratorWrapper();
			internal:
				StaticIntegrator * _StaticIntegrator;
			private:

			};

			public ref class LoadControlWrapper : StaticIntegratorWrapper
			{
			public:
				LoadControlWrapper(double deltaLambda, int numIncr, double minLambda, double maxlambda);
				~LoadControlWrapper();
			internal:

			private:

			};

			public ref class DisplacementControlWrapper : StaticIntegratorWrapper
			{
			public:
				DisplacementControlWrapper(int node, int dof,
					double increment,
					DomainWrapper^ domain,
					int numIncr,
					double min, double max);
				~DisplacementControlWrapper();
			internal:

			private:

			};

			public ref class ArcLengthWrapper : StaticIntegratorWrapper
			{
			public:
				ArcLengthWrapper(double arcLength, double alpha);
				~ArcLengthWrapper();
			internal:

			private:

			};

			public ref class ArcLength1Wrapper : StaticIntegratorWrapper
			{
			public:
				ArcLength1Wrapper(double arcLength, double alpha);
				~ArcLength1Wrapper();
			internal:

			private:

			};

			public ref class HSConstraintWrapper : StaticIntegratorWrapper
			{
			public:
				HSConstraintWrapper(double arcLength, double psi_u, double psi_f, double u_ref);
				~HSConstraintWrapper();
			internal:

			private:

			};

			public ref class LoadPathWrapper : StaticIntegratorWrapper
			{
			public:
				LoadPathWrapper();
				LoadPathWrapper(VectorWrapper^ theLoadPath);
				~LoadPathWrapper();
			internal:

			private:

			};

			public ref class MinUnbalDispNormWrapper : StaticIntegratorWrapper
			{
			public:
				MinUnbalDispNormWrapper(double lambda1, int specNumIterStep, double dlambda1min, double dlambda1max, int signFirstStepMethod);
				~MinUnbalDispNormWrapper();
			internal:

			private:

			};

			public ref class QuadraticMethodWrapper : StaticIntegratorWrapper
			{
			public:
				QuadraticMethodWrapper(double arcLength, int type);
				~QuadraticMethodWrapper();
			internal:

			private:

			};
		}
	}
}
