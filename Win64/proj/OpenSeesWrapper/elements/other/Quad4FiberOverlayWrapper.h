#pragma once



#include <UWelements/Quad4FiberOverlay.h>



#include "../../taggeds/TaggedObjectWrapper.h"
#include "../../actors/IMovableObjectWrapper.h"
#include "../ElementWrapper.h"
#include "../crdTransf/CrdTransfWrapper.h"
#include "../../materials/section/SectionForceDeformationWrapper.h"
#include "../../materials/uniaxial/UniaxialMaterialWrapper.h"
#include "../../materials/ndmaterial/NDMaterialWrapper.h"
#include "../beamColumn/beamIntegration/BeamIntegrationWrapper.h"
#include "../../matrix/IDWrapper.h"
#include "../../matrix/MatrixWrapper.h"
#include "../../matrix/VectorWrapper.h"
#include "../../domains/domain/DomainWrapper.h"
#include "../../damage/DamageModelWrapper.h"


using namespace System;
using namespace OpenSees;
using namespace OpenSees::DamageModels;
using namespace OpenSees::Components;
using namespace OpenSees::Elements::BeamIntegrations;
using namespace OpenSees::Elements::CrdTransfs;
using namespace OpenSees::Materials::Sections;
using namespace OpenSees::Materials::Uniaxials;
using namespace OpenSees::Materials::NDMaterials;

namespace OpenSees {
	namespace Elements {
		

		public ref class Quad4FiberOverlayWrapper : ElementWrapper
		{

		public:
			Quad4FiberOverlayWrapper(int tag, int nd1, int nd2, int nd3, int nd4,
				UniaxialMaterialWrapper^ m, double Af, double beta1, double beta2);

			Quad4FiberOverlayWrapper() {};
			~Quad4FiberOverlayWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

	}

}
