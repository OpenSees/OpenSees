#pragma once

#include <ConstantPressureVolumeQuad.h>
#include <EnhancedQuad.h>
#include <FourNodeQuad.h>
#include <FourNodeQuad3d.h>
#include <FourNodeQuadWithSensitivity.h>
#include <NineNodeMixedQuad.h>

#include "../ElementWrapper.h"
#include "../../materials/ndmaterial/NDMaterialWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Materials::NDMaterials;
namespace OpenSees {
	namespace Elements {
		public enum class  FourNodeQuadType { PlaneStrain, PlaneStress };

		public ref class ConstantPressureVolumeQuadWrapper : ElementWrapper
		{
		public:
			ConstantPressureVolumeQuadWrapper(int tag,
				int Nd1, int Nd2, int Nd3, int Nd4,
				NDMaterialWrapper^ nDMaterial,
				double t);
			~ConstantPressureVolumeQuadWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class EnhancedQuadWrapper : ElementWrapper
		{
		public:
			EnhancedQuadWrapper(int tag,
				int Nd1, int Nd2, int Nd3, int Nd4,
				NDMaterialWrapper^ nDMaterial,
				FourNodeQuadType^ type, double t);
			~EnhancedQuadWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class FourNodeQuadWrapper : ElementWrapper
		{
		public:
			FourNodeQuadWrapper(int tag,
				int Nd1, int Nd2, int Nd3, int Nd4,
				NDMaterialWrapper^ nDMaterial,
				FourNodeQuadType^ type, double t);
			FourNodeQuadWrapper(int tag,
				int Nd1, int Nd2, int Nd3, int Nd4,
				NDMaterialWrapper^ nDMaterial,
				FourNodeQuadType^ type, double t, double pressure,
				double rho,
				double b1, double b2);
			~FourNodeQuadWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class FourNodeQuadWithSensitivityWrapper : ElementWrapper
		{
		public:
			FourNodeQuadWithSensitivityWrapper(int tag,
				int Nd1, int Nd2, int Nd3, int Nd4,
				NDMaterialWrapper^ nDMaterial,
				FourNodeQuadType^ type, double t);
			FourNodeQuadWithSensitivityWrapper(int tag,
				int Nd1, int Nd2, int Nd3, int Nd4,
				NDMaterialWrapper^ nDMaterial,
				FourNodeQuadType^ type, double t, double pressure,
				double rho,
				double b1, double b2);
			~FourNodeQuadWithSensitivityWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class FourNodeQuad3dWrapper : ElementWrapper
		{
		public:
			FourNodeQuad3dWrapper(int tag, 
				int Nd1, int Nd2, int Nd3, int Nd4,
				NDMaterialWrapper^ nDMaterial,
				FourNodeQuadType^ type,
				double t);
			~FourNodeQuad3dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class NineNodeMixedQuadWrapper : ElementWrapper
		{
		public:
			NineNodeMixedQuadWrapper(int tag,
				int node1,
				int node2,
				int node3,
				int node4,
				int node5,
				int node6,
				int node7,
				int node8,
				int node9,
				NDMaterialWrapper^ nDMaterial);
			~NineNodeMixedQuadWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};
	}
	
}
