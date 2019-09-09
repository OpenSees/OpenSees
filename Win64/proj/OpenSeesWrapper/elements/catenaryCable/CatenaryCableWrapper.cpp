#include "stdafx.h"
#include "CatenaryCableWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Elements;
using namespace OpenSees::Elements::BeamIntegrations;
using namespace OpenSees::Materials::Sections;
using namespace OpenSees::Materials::Uniaxials;

CatenaryCableWrapper::CatenaryCableWrapper(int tag, int node1, int node2, double weight, double E, double A, double L0, double alpha, double temperature_change, double rho, double error_tol, int Nsubsteps, int massType) {
	_Element = new CatenaryCable( tag,  node1,  node2,  weight,  E,  A,  L0,  alpha,  temperature_change,  rho,  error_tol,  Nsubsteps,  massType);
}


