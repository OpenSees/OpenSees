#pragma once

#include <integration/SectionIntegration.h>
#include <integration/RCSectionIntegration.h>
#include <integration/RCTBeamSectionIntegration.h>
#include <integration/WideFlangeSectionIntegration.h>

#include <ID.h>
#include <Fiber.h>
#include <repres\cell\Cell.h>
#include <repres\reinfBar\ReinfBar.h>

#include <UniaxialFiber2d.h>
#include <UniaxialFiber3d.h>
#include <UniaxialMaterial.h>

#include <ElasticSection3d.h>
#include <FiberSection2d.h>
#include <FiberSection3d.h>
#include <dispBeamColumnInt/FiberSection2dInt.h>

#include "../../../actors/IMovableObjectWrapper.h"
#include "../../uniaxial/UniaxialMaterialWrapper.h"
#include "../fiber/FiberWrapper.h"
#include "../repres/sect/SectionRepresWrapper.h"

using namespace System;
using namespace System::Collections::Generic;
using namespace OpenSees::Materials;
using namespace OpenSees::Materials::Uniaxials;
using namespace OpenSees::Materials::Sections::Repres;
using namespace OpenSees::Materials::Sections;

namespace OpenSees {
	namespace Materials {
		namespace Sections {
			namespace Integrations {
				public ref class SectionIntegrationWrapper abstract : IMovableObjectWrapper
				{
				public:
					SectionIntegrationWrapper() {};
					~SectionIntegrationWrapper() {
						if (_SectionIntegration != 0)
							delete _SectionIntegration;
					};
				internal:
					SectionIntegration * _SectionIntegration;
				};

				public ref class RCSectionIntegrationWrapper : SectionIntegrationWrapper
				{
				public:
					RCSectionIntegrationWrapper(double d, double b, double Atop, double Abottom,
						double Aside, double cover,
						int Nfcore, int Nfcover, int Nfs);
					~RCSectionIntegrationWrapper() {
						if (_SectionIntegration != 0)
							delete _SectionIntegration;
					};
				internal:
					
				};

				public ref class RCTBeamSectionIntegrationWrapper : SectionIntegrationWrapper
				{
				public:
					RCTBeamSectionIntegrationWrapper(double d, double bw, double beff, double hf, double Atop, double Abottom, double flcov,
						double wcov, int Nflcover, int Nwcover, int Nflcore, int Nwcore, int NsteelTop, int NsteelBottom);
					~RCTBeamSectionIntegrationWrapper() {
						if (_SectionIntegration != 0)
							delete _SectionIntegration;
					};
				internal:

				};

				public ref class WideFlangeSectionIntegrationWrapper : SectionIntegrationWrapper
				{
				public:
					WideFlangeSectionIntegrationWrapper(double d, double tw, double bf, double tf,
						int Nfdw, int Nftf);
					~WideFlangeSectionIntegrationWrapper() {
						if (_SectionIntegration != 0)
							delete _SectionIntegration;
					};
				internal:

				};
			}
		}
	}
}
