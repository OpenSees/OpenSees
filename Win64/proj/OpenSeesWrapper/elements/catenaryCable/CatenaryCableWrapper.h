#pragma once

#include <catenaryCable/CatenaryCable.h>


#include "../ElementWrapper.h"
#include "../../materials/section/SectionForceDeformationWrapper.h"
#include "../../materials/uniaxial/UniaxialMaterialWrapper.h"
#include "../../materials/ndmaterial/NDMaterialWrapper.h"
#include "../beamColumn/beamIntegration/BeamIntegrationWrapper.h"
#include "../../matrix/IDWrapper.h"
#include "../../matrix/MatrixWrapper.h"
#include "../../matrix/VectorWrapper.h"


using namespace OpenSees;
using namespace OpenSees::Elements::BeamIntegrations;
using namespace OpenSees::Materials::Sections;
using namespace OpenSees::Materials::Uniaxials;
using namespace OpenSees::Materials::NDMaterials;

namespace OpenSees {
	namespace Elements {
		public ref class CatenaryCableWrapper : ElementWrapper
		{
		public:
			CatenaryCableWrapper(int tag, int node1, int node2, double weight, double E, double A, double L0, double alpha, double temperature_change, double rho, double error_tol, int Nsubsteps, int massType);
			CatenaryCableWrapper() {};
			~CatenaryCableWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};


	}

}
