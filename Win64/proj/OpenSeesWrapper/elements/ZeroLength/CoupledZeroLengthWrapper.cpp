#include "stdafx.h"
#include "CoupledZeroLengthWrapper.h"

using namespace System::Runtime::InteropServices;
using namespace OpenSees;
using namespace OpenSees::Elements;
using namespace OpenSees::Elements::BeamIntegrations;
using namespace OpenSees::Materials::Sections;
using namespace OpenSees::Materials::Uniaxials;



CoupledZeroLengthWrapper::CoupledZeroLengthWrapper(int tag,
	int Nd1, int Nd2,
	UniaxialMaterialWrapper^ theMaterial,
	int direction1, int direction2,
	int doRayleighDamping) {

	_Element = new CoupledZeroLength( tag,
		 Nd1,  Nd2,
		*theMaterial->_UniaxialMaterial,
		 direction1,  direction2,
		 doRayleighDamping);
}





