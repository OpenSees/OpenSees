#pragma once
#include <DOF_Numberer.h>
#include "../../graphs/GraphNumbererWrapper.h"
#include "../../actors/IMovableObjectWrapper.h"
using namespace OpenSees::GraphNumberers;
namespace OpenSees {
	namespace Numberers {
		public ref class DOF_NumbererWrapper : IMovableObjectWrapper
		{
		public:
			DOF_NumbererWrapper(GraphNumbererWrapper^ graphNumberer);
			~DOF_NumbererWrapper();
		internal:
			DOF_Numberer * _DOFNumberer;
		private:

		};

		public ref class PlainNumbererWrapper : DOF_NumbererWrapper
		{
		public:
			PlainNumbererWrapper();
			~PlainNumbererWrapper();
		internal:

		private:

		};
	}
}
