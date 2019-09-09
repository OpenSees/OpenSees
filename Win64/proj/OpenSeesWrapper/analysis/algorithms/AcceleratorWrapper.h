#pragma once

#include <Accelerator.h>
#include "../../actors/IMovableObjectWrapper.h"


namespace OpenSees {
	namespace Algorithms {
		namespace Accelerators {
			public ref class AcceleratorWrapper abstract : IMovableObjectWrapper
			{
			public:
				AcceleratorWrapper();
				~AcceleratorWrapper();
			internal:
				Accelerator* _Accelerator;
			};



			/*public ref class DifferenceAcceleratorWrapper :AcceleratorWrapper
			{
			public:
				DifferenceAcceleratorWrapper(int maxDim, int tangent);
				~DifferenceAcceleratorWrapper();
			internal:

			};*/

			public ref class DifferenceAccelerator2Wrapper :AcceleratorWrapper
			{
			public:
				DifferenceAccelerator2Wrapper(int maxDim, int tangent);
				~DifferenceAccelerator2Wrapper();
			internal:

			};

			public ref class KrylovAcceleratorWrapper :AcceleratorWrapper
			{
			public:
				KrylovAcceleratorWrapper(int maxDim, int tangent);
				~KrylovAcceleratorWrapper();
			internal:

			};

			public ref class KrylovAccelerator2Wrapper :AcceleratorWrapper
			{
			public:
				KrylovAccelerator2Wrapper(int maxDim, int tangent);
				~KrylovAccelerator2Wrapper();
			internal:

			};

			public ref class MillerAcceleratorWrapper :AcceleratorWrapper
			{
			public:
				MillerAcceleratorWrapper(int maxDim , double tol ,
					int tangent);
				~MillerAcceleratorWrapper();
			internal:

			};

			public ref class PeriodicAcceleratorWrapper :AcceleratorWrapper
			{
			public:
				PeriodicAcceleratorWrapper(int iter, int tangent);
				~PeriodicAcceleratorWrapper();
			internal:

			};

			public ref class RaphsonAcceleratorWrapper :AcceleratorWrapper
			{
			public:
				RaphsonAcceleratorWrapper(int tangent);
				~RaphsonAcceleratorWrapper();
			internal:

			};

			public ref class SecantAccelerator1Wrapper :AcceleratorWrapper
			{
			public:
				SecantAccelerator1Wrapper(int maxIter, int tangent);
				~SecantAccelerator1Wrapper();
			internal:

			};

			public ref class SecantAccelerator2Wrapper :AcceleratorWrapper
			{
			public:
				SecantAccelerator2Wrapper(int maxIter, int tangent);
				~SecantAccelerator2Wrapper();
			internal:

			};

			public ref class SecantAccelerator3Wrapper :AcceleratorWrapper
			{
			public:
				SecantAccelerator3Wrapper(int maxIter, int tangent);
				~SecantAccelerator3Wrapper();
			internal:

			};
		}
		
	}
}