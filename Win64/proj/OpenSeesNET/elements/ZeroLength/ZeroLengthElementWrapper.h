#pragma once



//#include <zeroLength\CoupledZeroLength.h>
#include <zeroLength\ZeroLength.h>
#include <zeroLength\ZeroLengthContact2D.h>
#include <zeroLength\ZeroLengthContact3D.h>
#include <zeroLength\ZeroLengthContactNTS2D.h>
#include <zeroLength\ZeroLengthImpact3D.h>
#include <zeroLength\ZeroLengthInterface2D.h>
#include <zeroLength\ZeroLengthND.h>
#include <zeroLength\ZeroLengthRocking.h>
#include <zeroLength\ZeroLengthSection.h>



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
		/*public ref class CoupledZeroLengthWrapper : ElementWrapper
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
		};*/

		public ref class ZeroLengthWrapper : ElementWrapper
		{

		public:
			ZeroLengthWrapper(int tag,
				int dimension,
				int Nd1, int Nd2,
				VectorWrapper^ x,
				VectorWrapper^ yprime,
				UniaxialMaterialWrapper^ theMaterial,
				int direction,
				int doRayleighDamping);

			ZeroLengthWrapper(int tag,
				int dimension,
				int Nd1, int Nd2,
				VectorWrapper^ x,
				VectorWrapper^ yprime,
				UniaxialMaterialWrapper^ theMaterial,
				UniaxialMaterialWrapper^ theDampingMaterial,
				int direction);

			ZeroLengthWrapper(int tag,
				int dimension,
				int Nd1, int Nd2,
				VectorWrapper^ x,
				VectorWrapper^ yprime,
				int n1dMat,
				array<UniaxialMaterialWrapper^>^ theMaterials,
				IDWrapper^ direction,
				int doRaylieghDamping);

			ZeroLengthWrapper(int tag,
				int dimension,
				int Nd1, int Nd2,
				VectorWrapper^ x,
				VectorWrapper^ yprime,
				int n1dMat,
				array<UniaxialMaterialWrapper^>^ theMaterials,
				array<UniaxialMaterialWrapper^>^ theDampMaterial,
				IDWrapper^ direction,
				int doRaylieghDamping);

			ZeroLengthWrapper() {};
			~ZeroLengthWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class ZeroLengthContact2DWrapper : ElementWrapper
		{

		public:
			ZeroLengthContact2DWrapper(int tag, int Nd1, int Nd2,
				double Kn, double Kt, double fRatio,
				VectorWrapper^ normal);

			ZeroLengthContact2DWrapper() {};
			~ZeroLengthContact2DWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class ZeroLengthContact3DWrapper : ElementWrapper
		{

		public:
			ZeroLengthContact3DWrapper(int tag,
				int Nd1, int Nd2,
				int direction, double Kn, double Kt, double fRatio, double c,
				double originX, double originY);

			ZeroLengthContact3DWrapper() {};
			~ZeroLengthContact3DWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class ZeroLengthContactNTS2DWrapper : ElementWrapper
		{

		public:
			ZeroLengthContactNTS2DWrapper(int tag, int sNdNum, int mNdNum, IDWrapper^ Nodes,
				double Kn, double Kt, double fRatio);

			ZeroLengthContactNTS2DWrapper() {};
			~ZeroLengthContactNTS2DWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class ZeroLengthImpact3DWrapper : ElementWrapper
		{

		public:
			ZeroLengthImpact3DWrapper(int tag,
				int Nd1, int Nd2,
				int direction,
				double initGapInput, double fRatio, double Kt,
				double Kn, double Kn2Input, double Delta_yInput,
				double c);

			ZeroLengthImpact3DWrapper() {};
			~ZeroLengthImpact3DWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class ZeroLengthInterface2DWrapper : ElementWrapper
		{

		public:
			ZeroLengthInterface2DWrapper(int tag, int sNdNum, int mNdNum, int sDof, int mDof, IDWrapper^ Nodes,
				double Kn, double Kt, double fRatio);

			ZeroLengthInterface2DWrapper() {};
			~ZeroLengthInterface2DWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class ZeroLengthNDWrapper : ElementWrapper
		{

		public:
			ZeroLengthNDWrapper(int tag,
				int dimension,
				int Nd1, int Nd2,
				VectorWrapper^ x,
				VectorWrapper^ yprime,
				NDMaterialWrapper^ theNDMaterial);

			ZeroLengthNDWrapper(int tag,
				int dimension,
				int Nd1, int Nd2,
				VectorWrapper^ x,
				VectorWrapper^ yprime,
				NDMaterialWrapper^ theNDMaterial, UniaxialMaterialWrapper^ the1DMaterial);

			ZeroLengthNDWrapper() {};
			~ZeroLengthNDWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class ZeroLengthRockingWrapper : ElementWrapper
		{

		public:
			ZeroLengthRockingWrapper(int tag,
				int dimension,
				int Nd1, int Nd2,
				VectorWrapper^ x,
				VectorWrapper^ yprime,
				double kr, double R,
				double theta, double kappa,
				double xi, double dispTol, double velTol);

			ZeroLengthRockingWrapper() {};
			~ZeroLengthRockingWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

		public ref class ZeroLengthSectionWrapper : ElementWrapper
		{

		public:
			ZeroLengthSectionWrapper(int tag,
				int dimension,
				int Nd1, int Nd2,
				VectorWrapper^ x,
				VectorWrapper^ yprime,
				SectionForceDeformationWrapper^ theSection,
				int doRayleighDamping);

			ZeroLengthSectionWrapper() {};
			~ZeroLengthSectionWrapper() {
				if (_Element != 0)
					delete _Element;
			};
		};

#pragma endregion


	}

}
