#pragma once

#include <adapter/Adapter.h>
#include <adapter/ActuatorCorot.h>
#include <adapter/Actuator.h>

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
		public ref class AdapterWrapper : ElementWrapper
		{
		public:
			AdapterWrapper(int tag, IDWrapper^ nodes, IDWrapper^ dof,
				MatrixWrapper^ stif, int ipPort,
				int addRayleigh, MatrixWrapper^ mass);
			AdapterWrapper() {};
			~AdapterWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class ActuatorCorotWrapper : ElementWrapper
		{
		public:
			ActuatorCorotWrapper(int tag, int dim, int Nd1, int Nd2,
				double EA, int ipPort, int addRayleigh, double rho);
			ActuatorCorotWrapper() {};
			~ActuatorCorotWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class ActuatorWrapper : ElementWrapper
		{
		public:
			ActuatorWrapper(int tag, int dim, int Nd1, int Nd2,
				double EA, int ipPort, int addRayleigh, double rho);
			ActuatorWrapper() {};
			~ActuatorWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

	}

}
