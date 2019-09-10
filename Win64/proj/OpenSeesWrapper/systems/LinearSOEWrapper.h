#pragma once
#include <LinearSOE.h>
#include "SystemOfEqnWrapper.h"

namespace OpenSees {
	namespace Systems {
		namespace Linears {
			public ref class LinearSOEWrapper abstract : SystemOfEqnWrapper
			{
			public:
				LinearSOEWrapper();
				~LinearSOEWrapper();
			internal:
				LinearSOE * _LinearSOE;

			};
		}
	}
}
