#pragma once

#include <brick/BbarBrick.h>
#include <brick/BbarBrickWithSensitivity.h>
#include <brick/Brick.h>
#include <brick/Twenty_Node_Brick.h>


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
		public ref class BbarBrickWrapper : ElementWrapper
		{
		public:
			BbarBrickWrapper(int tag,
				int node1,
				int node2,
				int node3,
				int node4,
				int node5,
				int node6,
				int node7,
				int node8,
				NDMaterialWrapper^ theMaterial,
				double b1, double b2, double b3);
			BbarBrickWrapper() {};
			~BbarBrickWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class BbarBrickWithSensitivityWrapper : ElementWrapper
		{
		public:
			BbarBrickWithSensitivityWrapper(int tag,
				int node1,
				int node2,
				int node3,
				int node4,
				int node5,
				int node6,
				int node7,
				int node8,
				NDMaterialWrapper^ theMaterial,
				double b1, double b2, double b3);
			BbarBrickWithSensitivityWrapper() {};
			~BbarBrickWithSensitivityWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class BrickWrapper : ElementWrapper
		{
		public:
			BrickWrapper(int tag,
				int node1,
				int node2,
				int node3,
				int node4,
				int node5,
				int node6,
				int node7,
				int node8,
				NDMaterialWrapper^ theMaterial,
				double b1, double b2, double b3);
			BrickWrapper() {};
			~BrickWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class Twenty_Node_BrickWrapper : ElementWrapper
		{
		public:
			Twenty_Node_BrickWrapper(int tag,
				int node1,
				int node2,
				int node3,
				int node4,
				int node5,
				int node6,
				int node7,
				int node8,
				int node9,
				int node10,
				int node11,
				int node12,
				int node13,
				int node14,
				int node15,
				int node16,
				int node17,
				int node18,
				int node19,
				int node20,
				NDMaterialWrapper^ theMaterial,
				double b1, double b2, double b3);
			Twenty_Node_BrickWrapper() {};
			~Twenty_Node_BrickWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

	}

}
