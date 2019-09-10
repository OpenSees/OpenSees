#pragma once
#include <GraphNumberer.h>
#include <RCM.h>
#include <PlainNumberer.h>
#include "../actors/IMovableObjectWrapper.h"

namespace OpenSees {
	namespace GraphNumberers {
		public ref class GraphNumbererWrapper : IMovableObjectWrapper
		{
		public:
			GraphNumbererWrapper();
			~GraphNumbererWrapper();
		internal:
			GraphNumberer * _GraphNumberer;
		private:

		};

		public ref class RCMWrapper : GraphNumbererWrapper
		{
		public:
			RCMWrapper(bool GPS);
			~RCMWrapper();
		internal:

		private:

		};
	}
}
