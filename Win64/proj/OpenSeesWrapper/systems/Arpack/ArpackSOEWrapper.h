#pragma once
#include "../EigenSOEWrapper.h"
#include <EigenSOE.h>

namespace OpenSees {
	namespace Systems {
		namespace Eigens {
			public ref class ArpackSOEWrapper : EigenSOEWrapper
			{
			public:
				ArpackSOEWrapper(double shift);
				ArpackSOEWrapper();
				~ArpackSOEWrapper();
			internal:

			private:

			};
		}
	}
}
