#pragma once



#include <surfaceLoad\SurfaceLoad.h>
#include <surfaceLoad\TriSurfaceLoad.h>




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
		
#pragma region surfaceLoad
		public ref class SurfaceLoadWrapper : ElementWrapper
		{

		public:
			SurfaceLoadWrapper(int tag, int Nd1, int Nd2, int Nd3, int Nd4, double pressure);

			SurfaceLoadWrapper() {};
			~SurfaceLoadWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class TriSurfaceLoadWrapper : ElementWrapper
		{

		public:
			TriSurfaceLoadWrapper(int tag, int Nd1, int Nd2, int Nd3, double pressure);

			TriSurfaceLoadWrapper() {};
			~TriSurfaceLoadWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};
#pragma endregion

	}

}
