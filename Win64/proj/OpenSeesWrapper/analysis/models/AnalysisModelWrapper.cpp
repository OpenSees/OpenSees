#include "stdafx.h"
#include <AnalysisModel.h>
#include "AnalysisModelWrapper.h"

using namespace OpenSees;
AnalysisModelWrapper::AnalysisModelWrapper()
{
	_AnalysisModel = new AnalysisModel();
}

AnalysisModelWrapper::~AnalysisModelWrapper()
{
	if (_AnalysisModel != 0)
		delete _AnalysisModel;
}

