#include "stdafx.h"
#include "BeamIntegrationWrapper.h"


using namespace OpenSees::Elements::BeamIntegrations;

BeamIntegrationWrapper::BeamIntegrationWrapper(BeamIntegrationType^ type)
{
	_type = *type;
	switch (*type)
	{
	case BeamIntegrationType::Legendre: {
		_BeamIntegration = new LegendreBeamIntegration();

		break;
	}
	case BeamIntegrationType::Lobatto: {
		_BeamIntegration = new LobattoBeamIntegration();
		break;
	}
	case BeamIntegrationType::Radau: {
		_BeamIntegration = new RadauBeamIntegration();
		break;
	}

	case BeamIntegrationType::CompositeSimpson: {
		_BeamIntegration = new CompositeSimpsonBeamIntegration();
		break;
	}
	case BeamIntegrationType::NewtonCotes: {
		_BeamIntegration = new NewtonCotesBeamIntegration();
		break;
	}
	case BeamIntegrationType::Trapezoidal: {
		_BeamIntegration = new TrapezoidalBeamIntegration();
		break;
	}
	default:
		_BeamIntegration = new LegendreBeamIntegration();
		break;
	}
}


DistHingeIntegrationWrapper::DistHingeIntegrationWrapper(double lpI, double lpJ, BeamIntegrationWrapper^ bi)
{
	_BeamIntegration = new DistHingeIntegration(lpI, lpJ, *bi->_BeamIntegration);
}

FixedLocationBeamIntegrationWrapper::FixedLocationBeamIntegrationWrapper(int nIP, VectorWrapper^ pt)
{
	_BeamIntegration = new FixedLocationBeamIntegration(nIP, *pt->_Vector);
}

HingeEndpointBeamIntegrationWrapper::HingeEndpointBeamIntegrationWrapper(double lpI, double lpJ)
{
	_BeamIntegration = new HingeEndpointBeamIntegration(lpI, lpJ);
}

HingeMidpointBeamIntegrationWrapper::HingeMidpointBeamIntegrationWrapper(double lpI, double lpJ)
{
	_BeamIntegration = new HingeMidpointBeamIntegration(lpI, lpJ);
}

HingeRadauBeamIntegrationWrapper::HingeRadauBeamIntegrationWrapper(double lpI, double lpJ)
{
	_BeamIntegration = new HingeRadauBeamIntegration(lpI, lpJ);
}

HingeRadauTwoBeamIntegrationWrapper::HingeRadauTwoBeamIntegrationWrapper(double lpI, double lpJ)
{
	_BeamIntegration = new HingeRadauTwoBeamIntegration(lpI, lpJ);
}

LowOrderBeamIntegrationWrapper::LowOrderBeamIntegrationWrapper(int nIP, VectorWrapper^ pt, int nc, VectorWrapper^ wt)
{
	_BeamIntegration = new LowOrderBeamIntegration(nIP, *pt->_Vector, nc, *wt->_Vector);
}

MidDistanceBeamIntegrationWrapper::MidDistanceBeamIntegrationWrapper(int nIP, VectorWrapper^ pt)
{
	_BeamIntegration = new MidDistanceBeamIntegration(nIP, *pt->_Vector);
}

RegularizedHingeIntegrationWrapper::RegularizedHingeIntegrationWrapper(BeamIntegrationWrapper^ bi,
	double lpI, double lpJ,
	double epsI, double epsJ)
{
	_BeamIntegration = new RegularizedHingeIntegration(*bi->_BeamIntegration,lpI, lpJ, epsI, epsJ);
}

UserDefinedBeamIntegrationWrapper::UserDefinedBeamIntegrationWrapper(int nIP, VectorWrapper^ pt, VectorWrapper^ wt)
{
	_BeamIntegration = new UserDefinedBeamIntegration(nIP, *pt->_Vector, *wt->_Vector);
}

UserDefinedHingeIntegrationWrapper::UserDefinedHingeIntegrationWrapper(int npL, VectorWrapper^ ptL, VectorWrapper^ wtL,
	int npR, VectorWrapper^ ptR, VectorWrapper^ wtR)
{
	_BeamIntegration = new UserDefinedHingeIntegration(npL, *ptL->_Vector, *wtL->_Vector,
		npR, *ptR->_Vector, *wtR->_Vector);
}

