#include "stdafx.h"
#include "AdapterWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Elements;
using namespace OpenSees::Elements::BeamIntegrations;
using namespace OpenSees::Materials::Sections;
using namespace OpenSees::Materials::Uniaxials;

AdapterWrapper::AdapterWrapper(int tag, IDWrapper^ nodes, IDWrapper^ dof,
	MatrixWrapper^ stif, int ipPort,
	int addRayleigh, MatrixWrapper^ mass) {
	_Element = new Adapter(tag, *nodes->_ID, dof->_ID,
		*stif->_Matrix, ipPort,
		addRayleigh, mass != nullptr ? mass->_Matrix : 0);
}

ActuatorCorotWrapper::ActuatorCorotWrapper(int tag, int dim, int Nd1, int Nd2,
	double EA, int ipPort, int addRayleigh, double rho) {
	_Element = new ActuatorCorot(tag,  dim,  Nd1,  Nd2,
		 EA,  ipPort,  addRayleigh,  rho);
}

ActuatorWrapper::ActuatorWrapper(int tag, int dim, int Nd1, int Nd2,
	double EA, int ipPort, int addRayleigh, double rho) {
	_Element = new Actuator(tag, dim, Nd1, Nd2,
		EA, ipPort, addRayleigh, rho);
}