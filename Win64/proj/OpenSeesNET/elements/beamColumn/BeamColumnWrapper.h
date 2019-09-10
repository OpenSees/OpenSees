#pragma once

#include <dispBeamColumn/AxEqDispBeamColumn2d.h>
#include <dispBeamColumn/DispBeamColumn2d.h>
#include <dispBeamColumn/DispBeamColumn2dThermal.h>
#include <dispBeamColumn/DispBeamColumn2dWithSensitivity.h>
#include <dispBeamColumn/DispBeamColumn3d.h>
#include <dispBeamColumn/DispBeamColumn3dThermal.h>
#include <dispBeamColumn/DispBeamColumn3dWithSensitivity.h>

#include <dispBeamColumnInt/DispBeamColumn2dInt.h>

#include <ElasticBeam2d.h>
#include <ElasticBeam3d.h>
#include <ElasticTimoshenkoBeam2d.h>
#include <ElasticTimoshenkoBeam3d.h>
#include <ModElasticBeam2d.h>

#include <ElasticForceBeamColumn2d.h>
#include <ElasticForceBeamColumn3d.h>
#include <ElasticForceBeamColumnWarping2d.h>
#include <ForceBeamColumn2d.h>
#include <ForceBeamColumn2dThermal.h>
#include <ForceBeamColumn3d.h>
#include <ForceBeamColumnCBDI2d.h>
#include <ForceBeamColumnWarping2d.h>

#include <mvlem/MVLEM.h>
#include <mvlem/SFI_MVLEM.h>

#include "../ElementWrapper.h"
#include "../../materials/section/SectionForceDeformationWrapper.h"
#include "../../materials/uniaxial/UniaxialMaterialWrapper.h"
#include "../../materials/ndmaterial/NDMaterialWrapper.h"
#include "beamIntegration/BeamIntegrationWrapper.h"
#include "../crdTransf/CrdTransfWrapper.h"


using namespace OpenSees;
using namespace OpenSees::Elements::CrdTransfs;
using namespace OpenSees::Elements::BeamIntegrations;
using namespace OpenSees::Materials::Sections;
using namespace OpenSees::Materials::Uniaxials;
using namespace OpenSees::Materials::NDMaterials;

namespace OpenSees {
	namespace Elements {

		public interface class ISectionForceDeformationBeamColumnWrapper
		{
		public:
			array<SectionForceDeformationWrapper^>^ GetSections();
			BeamIntegrationWrapper^ GetBeamIntegration();
			CrdTransfWrapper^ GetCrdTransf();
		};

		public ref class AxEqDispBeamColumn2dWrapper : ElementWrapper
		{
		public:
			AxEqDispBeamColumn2dWrapper(int tag, int nd1, int nd2,
				array<SectionForceDeformationWrapper^>^ s,
				BeamIntegrationWrapper^ bi, CrdTransfWrapper^ coordTransf, double tol,
				double rho, int cMass, int maxIters);
			AxEqDispBeamColumn2dWrapper() {};
			~AxEqDispBeamColumn2dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class DispBeamColumn2dWrapper : ElementWrapper, ISectionForceDeformationBeamColumnWrapper
		{
		public:
			DispBeamColumn2dWrapper(int tag, int nd1, int nd2,
				array<SectionForceDeformationWrapper^>^ sections,
				BeamIntegrationWrapper^ beamIntegr,
				CrdTransfWrapper^ coordTransf, double rho, int cMass);
			DispBeamColumn2dWrapper() {};

			virtual array<SectionForceDeformationWrapper^>^ __clrcall GetSections() sealed { 
				DispBeamColumn2d* ele = (DispBeamColumn2d*)_Element;
				SectionForceDeformation** sections = ele->theSections;
				int numSections = ele->numSections;
				array<SectionForceDeformationWrapper^>^ _sections = gcnew array<SectionForceDeformationWrapper^>(numSections);
				for (int i = 0; i < numSections; i++)
				{
					SectionForceDeformationWrapper^  sec = gcnew SectionForceDeformationWrapper();
					sec->_SectionForceDeformation = sections[i];
					sec->_TaggedObject = sections[i];
					_sections[i] = sec;
				}

				return _sections;
			};

			virtual BeamIntegrationWrapper^ __clrcall GetBeamIntegration() sealed {
				DispBeamColumn2d* ele = (DispBeamColumn2d*)_Element;
				BeamIntegrationWrapper^ beamInt = gcnew BeamIntegrationWrapper();
				beamInt->_BeamIntegration = ele->beamInt;
				return beamInt;
			};

			virtual CrdTransfWrapper^ __clrcall GetCrdTransf() sealed {
				DispBeamColumn2d* ele = (DispBeamColumn2d*)_Element;
				CrdTransfWrapper^ crdTransf = gcnew CrdTransfWrapper();
				crdTransf->_CrdTransf = ele->crdTransf;
				return crdTransf;
			};

			~DispBeamColumn2dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class DispBeamColumn2dThermalWrapper : ElementWrapper
		{
		public:
			DispBeamColumn2dThermalWrapper(int tag, int nd1, int nd2,
				array<SectionForceDeformationWrapper^>^ sections,
				BeamIntegrationWrapper^ beamIntegr,
				CrdTransfWrapper^ coordTransf, double rho);
			DispBeamColumn2dThermalWrapper() {};

			~DispBeamColumn2dThermalWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class DispBeamColumn2dWithSensitivityWrapper : ElementWrapper
		{
		public:
			DispBeamColumn2dWithSensitivityWrapper(int tag, int nd1, int nd2,
				array<SectionForceDeformationWrapper^>^ sections,
				BeamIntegrationWrapper^ beamIntegr,
				CrdTransfWrapper^ coordTransf, double rho);
			DispBeamColumn2dWithSensitivityWrapper() {};

			~DispBeamColumn2dWithSensitivityWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class DispBeamColumn3dWrapper : ElementWrapper, ISectionForceDeformationBeamColumnWrapper
		{
		public:
			DispBeamColumn3dWrapper(int tag, int nodeI, int nodeJ,
				int numSections, array<SectionForceDeformationWrapper^>^ sections,
				BeamIntegrationWrapper^ beamIntegr,
				CrdTransfWrapper^ coordTransf);
			DispBeamColumn3dWrapper(int tag, int nodeI, int nodeJ,
				int numSections, SectionForceDeformationWrapper^ section,
				BeamIntegrationWrapper^ beamIntegr,
				CrdTransfWrapper^ coordTransf);
			virtual array<SectionForceDeformationWrapper^>^ __clrcall GetSections() sealed {
				DispBeamColumn3d* ele = (DispBeamColumn3d*)_Element;
				SectionForceDeformation** sections = ele->theSections;
				int numSections = ele->numSections;
				array<SectionForceDeformationWrapper^>^ _sections = gcnew array<SectionForceDeformationWrapper^>(numSections);
				for (int i = 0; i < numSections; i++)
				{
					SectionForceDeformationWrapper^  sec = gcnew SectionForceDeformationWrapper();
					sec->_SectionForceDeformation = sections[i];
					sec->_TaggedObject = sections[i];
					_sections[i] = sec;
				}
				return _sections;
			};

			virtual BeamIntegrationWrapper^ __clrcall GetBeamIntegration() sealed {
				DispBeamColumn3d* ele = (DispBeamColumn3d*)_Element;
				BeamIntegrationWrapper^ beamInt = gcnew BeamIntegrationWrapper();
				beamInt->_BeamIntegration = ele->beamInt;
				return beamInt;
			};

			virtual CrdTransfWrapper^ __clrcall GetCrdTransf() sealed {
				DispBeamColumn3d* ele = (DispBeamColumn3d*)_Element;
				CrdTransfWrapper^ crdTransf = gcnew CrdTransfWrapper();
				crdTransf->_CrdTransf = ele->crdTransf;
				return crdTransf;
			};
			~DispBeamColumn3dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class DispBeamColumn3dThermalWrapper : ElementWrapper
		{
		public:
			DispBeamColumn3dThermalWrapper(int tag, int nd1, int nd2,
				array<SectionForceDeformationWrapper^>^ sections,
				BeamIntegrationWrapper^ beamIntegr,
				CrdTransfWrapper^ coordTransf);
			~DispBeamColumn3dThermalWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class DispBeamColumn3dWithSensitivityWrapper : ElementWrapper
		{
		public:
			DispBeamColumn3dWithSensitivityWrapper(int tag, int nd1, int nd2,
				array<SectionForceDeformationWrapper^>^ sections,
				BeamIntegrationWrapper^ beamIntegr,
				CrdTransfWrapper^ coordTransf,double r);
			~DispBeamColumn3dWithSensitivityWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class DispBeamColumn2dIntWrapper : ElementWrapper
		{
		public:
			DispBeamColumn2dIntWrapper(int tag, int nodeI, int nodeJ,
				array<SectionForceDeformationWrapper^>^ sections,
				CrdTransfWrapper^ coordTransf, double C, double rho);

			~DispBeamColumn2dIntWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class ElasticBeam2dWrapper : ElementWrapper
		{
		public:
			ElasticBeam2dWrapper(int tag, double A, double E, double I,
				int Nd1, int Nd2, CrdTransfWrapper^ theTransf,
				double alpha, double d,
				double rho, int cMass);
			ElasticBeam2dWrapper() {};

			~ElasticBeam2dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class ElasticBeam3dWrapper : ElementWrapper
		{
		public:
			ElasticBeam3dWrapper(int tag, double A, double E, double G,
				double Jx, double Iy, double Iz,
				int Nd1, int Nd2, CrdTransfWrapper^ theTransf,
				double rho, int cMass,
				int sectionTag);

			ElasticBeam3dWrapper(int tag, int Nd1, int Nd2, SectionForceDeformationWrapper^ section,
				CrdTransfWrapper^ theTransf, double rho, int cMass);

			ElasticBeam3dWrapper() {};

			~ElasticBeam3dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class ElasticTimoshenkoBeam2dWrapper : ElementWrapper
		{
		public:
			ElasticTimoshenkoBeam2dWrapper(int tag, int Nd1, int Nd2, double E, double G,
				double A, double Iz, double Avy, CrdTransfWrapper^ theTransf,
				double rho, int cMass);

			ElasticTimoshenkoBeam2dWrapper() {};

			~ElasticTimoshenkoBeam2dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class ElasticTimoshenkoBeam3dWrapper : ElementWrapper
		{
		public:
			ElasticTimoshenkoBeam3dWrapper(int tag, int Nd1, int Nd2, double E, double G,
				double A, double Jx, double Iy, double Iz, double Avy, double Avz,
				CrdTransfWrapper^ theTransf, double rho, int cMass);

			ElasticTimoshenkoBeam3dWrapper() {};

			~ElasticTimoshenkoBeam3dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class ModElasticBeam2dWrapper : ElementWrapper
		{
		public:
			ModElasticBeam2dWrapper(int tag, double A, double E, double I,
				int Nd1, int Nd2, double K11, double K33,
				double K44, CrdTransfWrapper^ theTransf,
				double alpha, double d ,
				double rho , int cMass);

			ModElasticBeam2dWrapper() {};

			~ModElasticBeam2dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class ElasticForceBeamColumn2dWrapper : ElementWrapper
		{
		public:
			ElasticForceBeamColumn2dWrapper(int tag, int nodeI, int nodeJ,
				array<SectionForceDeformationWrapper^>^ sec,
				BeamIntegrationWrapper^ beamIntegr,
				CrdTransfWrapper^ coordTransf, double rho);

			ElasticForceBeamColumn2dWrapper() {};

			~ElasticForceBeamColumn2dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class ElasticForceBeamColumn3dWrapper : ElementWrapper
		{
		public:
			ElasticForceBeamColumn3dWrapper(int tag, int nodeI, int nodeJ,
				array<SectionForceDeformationWrapper^>^ sec,
				BeamIntegrationWrapper^ beamIntegr,
				CrdTransfWrapper^ coordTransf, double rho);

			ElasticForceBeamColumn3dWrapper() {};

			~ElasticForceBeamColumn3dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class ElasticForceBeamColumnWarping2dWrapper : ElementWrapper
		{
		public:
			ElasticForceBeamColumnWarping2dWrapper(int tag, int nodeI, int nodeJ,
				array<SectionForceDeformationWrapper^>^ sec,
				BeamIntegrationWrapper^ beamIntegr,
				CrdTransfWrapper^ coordTransf, double rho);

			ElasticForceBeamColumnWarping2dWrapper() {};

			~ElasticForceBeamColumnWarping2dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class ForceBeamColumn2dWrapper : ElementWrapper, ISectionForceDeformationBeamColumnWrapper
		{
		public:
			ForceBeamColumn2dWrapper(int tag, int nodeI, int nodeJ,
				int numSections, array<SectionForceDeformationWrapper^>^ sections,
				BeamIntegrationWrapper^ beamIntegr,
				CrdTransfWrapper^ coordTransf);
			ForceBeamColumn2dWrapper(int tag, int nodeI, int nodeJ,
				int numSections, SectionForceDeformationWrapper^ section,
				BeamIntegrationWrapper^ beamIntegr,
				CrdTransfWrapper^ coordTransf);
			ForceBeamColumn2dWrapper(int tag, int nodeI, int nodeJ,
				int numSections, array<SectionForceDeformationWrapper^>^ sections,
				BeamIntegrationWrapper^ beamIntegr,
				CrdTransfWrapper^ coordTransf, double rho,
				int maxNumIters, double tolerance);

			virtual array<SectionForceDeformationWrapper^>^ __clrcall GetSections() sealed {
				ForceBeamColumn2d* ele = (ForceBeamColumn2d*)_Element;
				if (ele == 0) return nullptr;
				SectionForceDeformation** sections = ele->sections;
				int numSections = ele->numSections;
				array<SectionForceDeformationWrapper^>^ _sections = gcnew array<SectionForceDeformationWrapper^>(numSections);
				for (int i = 0; i < numSections; i++)
				{
					SectionForceDeformationWrapper^  sec = gcnew SectionForceDeformationWrapper();
					sec->_SectionForceDeformation = sections[i];
					sec->_TaggedObject = sections[i];
					_sections[i] = sec;
					
				}
				return _sections;
			};

			virtual BeamIntegrationWrapper^ __clrcall GetBeamIntegration() sealed {
				ForceBeamColumn2d* ele = (ForceBeamColumn2d*)_Element;
				BeamIntegrationWrapper^ beamInt = gcnew BeamIntegrationWrapper();
				beamInt->_BeamIntegration = ele->beamIntegr;
				return beamInt;
			};

			virtual CrdTransfWrapper^ __clrcall GetCrdTransf() sealed {
				ForceBeamColumn2d* ele = (ForceBeamColumn2d*)_Element;
				CrdTransfWrapper^ crdTransf = gcnew CrdTransfWrapper();
				crdTransf->_CrdTransf = ele->crdTransf;
				return crdTransf;
			};

			~ForceBeamColumn2dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class ForceBeamColumn2dThermalWrapper : ElementWrapper
		{
		public:
			ForceBeamColumn2dThermalWrapper() {};
			ForceBeamColumn2dThermalWrapper(int tag, int nodeI, int nodeJ,
				int numSections, array<SectionForceDeformationWrapper^>^ sections,
				BeamIntegrationWrapper^ beamIntegr,
				CrdTransfWrapper^ coordTransf, double rho,
				int maxNumIters, double tolerance);

			~ForceBeamColumn2dThermalWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class ForceBeamColumn3dWrapper : ElementWrapper, ISectionForceDeformationBeamColumnWrapper
		{
		public:
			ForceBeamColumn3dWrapper(int tag, int nodeI, int nodeJ,
				int numSections, array<SectionForceDeformationWrapper^>^ sections,
				BeamIntegrationWrapper^ beamIntegr,
				CrdTransfWrapper^ coordTransf);
			ForceBeamColumn3dWrapper(int tag, int nodeI, int nodeJ,
				int numSections, SectionForceDeformationWrapper^ section,
				BeamIntegrationWrapper^ beamIntegr,
				CrdTransfWrapper^ coordTransf);
			ForceBeamColumn3dWrapper(int tag, int nodeI, int nodeJ,
				int numSections, array<SectionForceDeformationWrapper^>^ sections,
				BeamIntegrationWrapper^ beamIntegr,
				CrdTransfWrapper^ coordTransf, double rho,
				int maxNumIters, double tolerance);
			ForceBeamColumn3dWrapper(int tag, int nodeI, int nodeJ,
				int numSections, SectionForceDeformationWrapper^ section,
				BeamIntegrationWrapper^ beamIntegr,
				CrdTransfWrapper^ coordTransf, double rho,
				int maxNumIters, double tolerance);

			virtual array<SectionForceDeformationWrapper^>^ __clrcall GetSections() sealed {
				ForceBeamColumn3d* ele = (ForceBeamColumn3d*)_Element;
				if (ele == 0) return nullptr;
				SectionForceDeformation** sections = ele->sections;
				int numSections = ele->numSections;
				array<SectionForceDeformationWrapper^>^ _sections = gcnew array<SectionForceDeformationWrapper^>(numSections);
				for (int i = 0; i < numSections; i++)
				{
					SectionForceDeformationWrapper^  sec = gcnew SectionForceDeformationWrapper();
					sec->_SectionForceDeformation = sections[i];
					sec->_TaggedObject = sections[i];
					_sections[i] = sec;

				}
				return _sections;
			};

			virtual BeamIntegrationWrapper^ __clrcall GetBeamIntegration() sealed {
				ForceBeamColumn3d* ele = (ForceBeamColumn3d*)_Element;
				BeamIntegrationWrapper^ beamInt = gcnew BeamIntegrationWrapper();
				beamInt->_BeamIntegration = ele->beamIntegr;
				return beamInt;
			};

			virtual CrdTransfWrapper^ __clrcall GetCrdTransf() sealed {
				ForceBeamColumn3d* ele = (ForceBeamColumn3d*)_Element;
				CrdTransfWrapper^ crdTransf = gcnew CrdTransfWrapper();
				crdTransf->_CrdTransf = ele->crdTransf;
				return crdTransf;
			};

			~ForceBeamColumn3dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class ForceBeamColumnWarping2dWrapper : ElementWrapper
		{
		public:
			ForceBeamColumnWarping2dWrapper() {};
			ForceBeamColumnWarping2dWrapper(int tag, int nodeI, int nodeJ,
				int numSections, array<SectionForceDeformationWrapper^>^ sections,
				BeamIntegrationWrapper^ beamIntegr,
				CrdTransfWrapper^ coordTransf, double rho,
				int maxNumIters, double tolerance);
			~ForceBeamColumnWarping2dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class ForceBeamColumnCBDI2dWrapper : ElementWrapper
		{
		public:
			ForceBeamColumnCBDI2dWrapper() {};
			ForceBeamColumnCBDI2dWrapper(int tag, int nodeI, int nodeJ,
				int numSections, array<SectionForceDeformationWrapper^>^ sections,
				BeamIntegrationWrapper^ beamIntegr,
				CrdTransfWrapper^ coordTransf, double rho, bool includeShear,
				int maxNumIters, double tolerance);
			~ForceBeamColumnCBDI2dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class MVLEMWrapper : ElementWrapper
		{
		public:
			MVLEMWrapper(int tag,
				double Dens,
				int Nd1, int Nd2,
				array<UniaxialMaterialWrapper^>^ materialsConcrete,
				array<UniaxialMaterialWrapper^>^ materialsSteel,
				array<UniaxialMaterialWrapper^>^ materialsShear,
				array<double>^ Rho,
				array<double>^ thickness,
				array<double>^ width,
				int mm,
				double cc);
			MVLEMWrapper() {};
			~MVLEMWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class SFI_MVLEMWrapper : ElementWrapper
		{
		public:
			SFI_MVLEMWrapper(int tag,					// element tag
				int Nd1, int Nd2,				// end node tags          
				array<NDMaterialWrapper^>^ materials,			// array of material tags
				array<double>^ Thickness,				// array of macro-fiber thickness
				array<double>^ Width,					// array of macro-fiber widths
				int mm,							// number of macro-fibers (RC panels)
				double cc);
			SFI_MVLEMWrapper() {};
			~SFI_MVLEMWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};
	}

}
