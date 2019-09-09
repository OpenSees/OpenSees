#include "stdafx.h"
#include "SurfaceLoadWrapper.h"


using namespace OpenSees;
using namespace OpenSees::Elements;
using namespace OpenSees::Elements::BeamIntegrations;
using namespace OpenSees::Materials::Sections;
using namespace OpenSees::Materials::Uniaxials;


SurfaceLoadWrapper::SurfaceLoadWrapper(int tag, int Nd1, int Nd2, int Nd3, int Nd4, double pressure) {

	_Element = new SurfaceLoad(tag, Nd1, Nd2, Nd3, Nd4, pressure);
}

TriSurfaceLoadWrapper::TriSurfaceLoadWrapper(int tag, int Nd1, int Nd2, int Nd3, double pressure) {

	_Element = new TriSurfaceLoad(tag, Nd1, Nd2, Nd3, pressure);
}





