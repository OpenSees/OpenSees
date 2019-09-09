#pragma once

#include <LineSearch.h>
#include "../../actors/IMovableObjectWrapper.h"


namespace OpenSees {
	namespace Algorithms {
		namespace LineSearchs {
			public ref class LineSearchWrapper abstract : IMovableObjectWrapper
			{
			public:
				LineSearchWrapper();
				~LineSearchWrapper();
			internal:
				LineSearch* _LineSearch;
			};



			public ref class BisectionLineSearchWrapper :LineSearchWrapper
			{
			public:
				BisectionLineSearchWrapper();
				BisectionLineSearchWrapper(double tolerance,
					int maxIter,
					double minEta,
					double maxEta,
					int flag);
				~BisectionLineSearchWrapper();
			internal:

			};

			public ref class InitialInterpolatedLineSearchWrapper :LineSearchWrapper
			{
			public:
				InitialInterpolatedLineSearchWrapper();
				InitialInterpolatedLineSearchWrapper(double tolerance,
					int maxIter,
					double minEta,
					double maxEta,
					int printFlag);
				~InitialInterpolatedLineSearchWrapper();
			internal:

			};

			public ref class SecantLineSearchWrapper :LineSearchWrapper
			{
			public:
				SecantLineSearchWrapper();
				SecantLineSearchWrapper(double tolerance,
					int maxIter,
					double minEta,
					double maxEta,
					int printFlag);
				~SecantLineSearchWrapper();
			internal:

			};

			public ref class RegulaFalsiLineSearchWrapper :LineSearchWrapper
			{
			public:
				RegulaFalsiLineSearchWrapper();
				RegulaFalsiLineSearchWrapper(double tolerance,
					int maxIter,
					double minEta,
					double maxEta,
					int printFlag);
				~RegulaFalsiLineSearchWrapper();
			internal:

			};
		}
		
	}
}