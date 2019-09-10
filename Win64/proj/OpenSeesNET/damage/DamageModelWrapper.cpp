#include "stdafx.h"
#include "DamageModelWrapper.h"
using namespace OpenSees;
using namespace OpenSees::DamageModels;
using namespace System::Runtime::InteropServices;
DamageModelWrapper::DamageModelWrapper()
{
	
}

DamageModelWrapper::~DamageModelWrapper()
{
	if (_DamageModel != 0)
		delete _DamageModel;
}


HystereticEnergyWrapper::HystereticEnergyWrapper(int tag, double Etot, double Cpow)
{
	_DamageModel = new HystereticEnergy(tag, Etot, Cpow);
}

HystereticEnergyWrapper::~HystereticEnergyWrapper()
{
	if (_DamageModel != 0)
		delete _DamageModel;
}

KratzigWrapper::KratzigWrapper(int tag, double ultimatePosVal, double ultimateNegVal)
{
	_DamageModel = new Kratzig(tag, ultimatePosVal, ultimateNegVal);
}

KratzigWrapper::~KratzigWrapper()
{
	if (_DamageModel != 0)
		delete _DamageModel;
}


MehannyWrapper::MehannyWrapper(int tag, double alpha, double beta, double gamma,
	double ultimatePosValue, double ultimateNegValue, double abstol, double reltol, double posmodifier, double negmodifier)
{

	_DamageModel = new Mehanny(tag, alpha, beta, gamma, ultimatePosValue, ultimateNegValue, abstol, reltol, posmodifier, negmodifier);
}

MehannyWrapper::~MehannyWrapper()
{
	if (_DamageModel != 0)
		delete _DamageModel;
}


NormalizedPeakWrapper::NormalizedPeakWrapper(int tag, double maxVal, double minVal, DamageTypeWrapper^ type)
{
	System::String^ strType = type->ToString();
	char* str2 = (char*)(void*)Marshal::StringToHGlobalAnsi(strType);
	_DamageModel = new NormalizedPeak(tag, maxVal, minVal, str2);
}

NormalizedPeakWrapper::~NormalizedPeakWrapper()
{
	if (_DamageModel != 0)
		delete _DamageModel;
}


ParkAngWrapper::ParkAngWrapper(int tag, double deltaU, double beta, double sigmaY)
{
	_DamageModel = new ParkAng(tag, deltaU, beta, sigmaY);
}

ParkAngWrapper::~ParkAngWrapper()
{
	if (_DamageModel != 0)
		delete _DamageModel;
}


