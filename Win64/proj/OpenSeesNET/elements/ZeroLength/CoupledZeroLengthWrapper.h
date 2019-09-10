#pragma once



#include <zeroLength\CoupledZeroLength.h>




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
		public enum class  PlaneElementType { PlaneStrain, PlaneStress };

#pragma region zeroLength
		public ref class CoupledZeroLengthWrapper : ElementWrapper
		{

		public:
			CoupledZeroLengthWrapper(int tag,
				int Nd1, int Nd2,
				UniaxialMaterialWrapper^ theMaterial,
				int direction1, int direction2,
				int doRayleighDamping);

			CoupledZeroLengthWrapper() {};
			~CoupledZeroLengthWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};
#pragma endregion


	}

}
