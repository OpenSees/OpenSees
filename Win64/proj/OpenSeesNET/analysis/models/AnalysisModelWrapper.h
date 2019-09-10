#pragma once
#include <AnalysisModel.h>
#include "../../actors/IMovableObjectWrapper.h"


namespace OpenSees {

		public ref class AnalysisModelWrapper : IMovableObjectWrapper
		{
		public:
			AnalysisModelWrapper();
			int GetNumEqn() {
				return _AnalysisModel->getNumEqn();
			}
			~AnalysisModelWrapper();
		internal:
			AnalysisModel * _AnalysisModel;
		private:

		};
	}
