#include "stdafx.h"
#include "Quad4FiberOverlayWrapper.h"


using namespace OpenSees;
using namespace OpenSees::Elements;
using namespace OpenSees::Elements::BeamIntegrations;
using namespace OpenSees::Materials::Sections;
using namespace OpenSees::Materials::Uniaxials;


Quad4FiberOverlayWrapper::Quad4FiberOverlayWrapper(int tag, int nd1, int nd2, int nd3, int nd4,
	UniaxialMaterialWrapper^ m, double Af, double beta1, double beta2) {

	_Element = new Quad4FiberOverlay(tag, nd1, nd2, nd3, nd4,
		*m->_UniaxialMaterial, Af, beta1, beta2);
}





