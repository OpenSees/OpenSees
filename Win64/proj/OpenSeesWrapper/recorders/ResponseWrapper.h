#pragma once
#include <response.h>


using namespace System;

namespace OpenSees {
	namespace Recorders {
		public ref class ResponseWrapper 
		{
		public:
			ResponseWrapper();
			~ResponseWrapper() {
				if (_Response != 0)
					delete _Response;
			};
		internal:
			Response *_Response;
		};
	}
}