#include "stdafx.h"
#include "BeamColumnWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Elements;
using namespace OpenSees::Elements::CrdTransfs;
using namespace OpenSees::Elements::BeamIntegrations;
using namespace OpenSees::Materials::Sections;
using namespace OpenSees::Materials::Uniaxials;

ForceBeamColumn3dWrapper::ForceBeamColumn3dWrapper(int tag, int nodeI, int nodeJ,
	int numSections, array<SectionForceDeformationWrapper^>^ sections,
	BeamIntegrationWrapper^ beamIntegr,
	CrdTransfWrapper^ coordTransf)
{
	SectionForceDeformation** _sections = new SectionForceDeformation*[numSections];
	for (int i = 0; i < numSections; i++)
		_sections[i] = sections[i]->_SectionForceDeformation->getCopy();
	_Element = new ForceBeamColumn3d(tag, nodeI, nodeJ, numSections, _sections, *beamIntegr->_BeamIntegration, *coordTransf->_CrdTransf);
}

ForceBeamColumn3dWrapper::ForceBeamColumn3dWrapper(int tag, int nodeI, int nodeJ,
	int numSections, SectionForceDeformationWrapper^ section,
	BeamIntegrationWrapper^ beamIntegr,
	CrdTransfWrapper^ coordTransf)
{
	SectionForceDeformation** _sections = new SectionForceDeformation*[numSections];
	for (int i = 0; i < numSections; i++)
		_sections[i] = section->_SectionForceDeformation->getCopy();
	_Element = new ForceBeamColumn3d(tag, nodeI, nodeJ, numSections, _sections, *beamIntegr->_BeamIntegration, *coordTransf->_CrdTransf);
}

ForceBeamColumn3dWrapper::ForceBeamColumn3dWrapper(int tag, int nodeI, int nodeJ,
	int numSections, array<SectionForceDeformationWrapper^>^ sections,
	BeamIntegrationWrapper^ beamIntegr,
	CrdTransfWrapper^ coordTransf, double rho,
	int maxNumIters, double tolerance) {

	SectionForceDeformation** _sections = new SectionForceDeformation*[numSections];
	for (int i = 0; i < numSections; i++)
		_sections[i] = sections[i]->_SectionForceDeformation->getCopy();
	_Element = new ForceBeamColumn3d(tag, nodeI, nodeJ, numSections, _sections, *beamIntegr->_BeamIntegration,
		*coordTransf->_CrdTransf,rho,maxNumIters,tolerance);

}

ForceBeamColumn3dWrapper::ForceBeamColumn3dWrapper(int tag, int nodeI, int nodeJ,
	int numSections, SectionForceDeformationWrapper^ section,
	BeamIntegrationWrapper^ beamIntegr,
	CrdTransfWrapper^ coordTransf, double rho,
	int maxNumIters, double tolerance) {

	SectionForceDeformation** _sections = new SectionForceDeformation*[numSections];
	for (int i = 0; i < numSections; i++)
		_sections[i] = section->_SectionForceDeformation->getCopy();
	_Element = new ForceBeamColumn3d(tag, nodeI, nodeJ, numSections, _sections, *beamIntegr->_BeamIntegration,
		*coordTransf->_CrdTransf, rho, maxNumIters, tolerance);
}




ForceBeamColumn2dWrapper::ForceBeamColumn2dWrapper(int tag, int nodeI, int nodeJ,
	int numSections, array<SectionForceDeformationWrapper^>^ sections,
	BeamIntegrationWrapper^ beamIntegr,
	CrdTransfWrapper^ coordTransf)
{
	SectionForceDeformation** _sections = new SectionForceDeformation*[numSections];
	for (int i = 0; i < numSections; i++)
		_sections[i] = sections[i]->_SectionForceDeformation;
	_Element = new ForceBeamColumn2d(tag, nodeI, nodeJ, numSections, _sections, *beamIntegr->_BeamIntegration, *coordTransf->_CrdTransf);
}

ForceBeamColumn2dWrapper::ForceBeamColumn2dWrapper(int tag, int nodeI, int nodeJ,
	int numSections, SectionForceDeformationWrapper^ section,
	BeamIntegrationWrapper^ beamIntegr,
	CrdTransfWrapper^ coordTransf)
{
	SectionForceDeformation** _sections = new SectionForceDeformation*[numSections];
	for (int i = 0; i < numSections; i++)
		_sections[i] = section->_SectionForceDeformation->getCopy();
	_Element = new ForceBeamColumn2d(tag, nodeI, nodeJ, numSections, _sections, *beamIntegr->_BeamIntegration, *coordTransf->_CrdTransf);
}

ForceBeamColumn2dWrapper::ForceBeamColumn2dWrapper(int tag, int nodeI, int nodeJ,
	int numSections, array<SectionForceDeformationWrapper^>^ sections,
	BeamIntegrationWrapper^ beamIntegr,
	CrdTransfWrapper^ coordTransf, double rho,
	int maxNumIters, double tolerance) {

	SectionForceDeformation** _sections = new SectionForceDeformation*[numSections];
	for (int i = 0; i < numSections; i++)
		_sections[i] = sections[i]->_SectionForceDeformation;
	_Element = new ForceBeamColumn2d(tag, nodeI, nodeJ, numSections, _sections, *beamIntegr->_BeamIntegration,
		*coordTransf->_CrdTransf, rho, maxNumIters, tolerance);

}


ElasticBeam2dWrapper::ElasticBeam2dWrapper(int tag, double A, double E, double I, int Nd1, int Nd2,
	CrdTransfWrapper ^ theTransf, double alpha, double d,
	double rho, int cMass)
{
	_Element = new ElasticBeam2d(tag, A, E, I, Nd1, Nd2, *theTransf->_CrdTransf, alpha, d, rho, cMass);
}


DispBeamColumn3dWrapper::DispBeamColumn3dWrapper(int tag, int nodeI, int nodeJ,
	int numSections, array<SectionForceDeformationWrapper^>^ sections,
	BeamIntegrationWrapper^ beamIntegr,
	CrdTransfWrapper^ coordTransf)
{
	SectionForceDeformation** _sections = new SectionForceDeformation*[numSections];
	for (int i = 0; i < numSections; i++)
		_sections[i] = sections[i]->_SectionForceDeformation;
	_Element = new DispBeamColumn3d(tag, nodeI, nodeJ, numSections, _sections, *beamIntegr->_BeamIntegration, *coordTransf->_CrdTransf);
}

DispBeamColumn3dWrapper::DispBeamColumn3dWrapper(int tag, int nodeI, int nodeJ,
	int numSections, SectionForceDeformationWrapper^ section,
	BeamIntegrationWrapper^ beamIntegr,
	CrdTransfWrapper^ coordTransf)
{
	SectionForceDeformation** _sections = new SectionForceDeformation*[numSections];
	for (int i = 0; i < numSections; i++)
		_sections[i] = section->_SectionForceDeformation->getCopy();
	_Element = new DispBeamColumn3d(tag, nodeI, nodeJ, numSections, _sections, *beamIntegr->_BeamIntegration, *coordTransf->_CrdTransf);
}


AxEqDispBeamColumn2dWrapper::AxEqDispBeamColumn2dWrapper(int tag, int nd1, int nd2,
	array<SectionForceDeformationWrapper^>^ s,
	BeamIntegrationWrapper^ bi, CrdTransfWrapper^ coordTransf, double tol,
	double rho, int cMass, int maxIters)
{
	SectionForceDeformation** _sections = new SectionForceDeformation*[s->Length];
	for (int i = 0; i < s->Length; i++)
		_sections[i] = s[i]->_SectionForceDeformation->getCopy();
	_Element = new AxEqDispBeamColumn2d(tag, nd1, nd2, s->Length, _sections, *bi->_BeamIntegration, *coordTransf->_CrdTransf, tol,rho,cMass,maxIters);
}

DispBeamColumn2dWrapper::DispBeamColumn2dWrapper(int tag, int nd1, int nd2,
	array<SectionForceDeformationWrapper^>^ s,
	BeamIntegrationWrapper^ bi, CrdTransfWrapper^ coordTransf,
	double rho, int cMass)
{
	SectionForceDeformation** _sections = new SectionForceDeformation*[s->Length];
	for (int i = 0; i < s->Length; i++)
		_sections[i] = s[i]->_SectionForceDeformation->getCopy();
	_Element = new DispBeamColumn2d(tag, nd1, nd2, s->Length, _sections, *bi->_BeamIntegration, *coordTransf->_CrdTransf, rho, cMass);
}

DispBeamColumn2dThermalWrapper::DispBeamColumn2dThermalWrapper(int tag, int nd1, int nd2,
	array<SectionForceDeformationWrapper^>^ s,
	BeamIntegrationWrapper^ bi, CrdTransfWrapper^ coordTransf,
	double rho)
{
	SectionForceDeformation** _sections = new SectionForceDeformation*[s->Length];
	for (int i = 0; i < s->Length; i++)
		_sections[i] = s[i]->_SectionForceDeformation->getCopy();
	_Element = new DispBeamColumn2dThermal(tag, nd1, nd2, s->Length, _sections, *bi->_BeamIntegration, *coordTransf->_CrdTransf, rho);
}

DispBeamColumn2dWithSensitivityWrapper::DispBeamColumn2dWithSensitivityWrapper(int tag, int nd1, int nd2,
	array<SectionForceDeformationWrapper^>^ s,
	BeamIntegrationWrapper^ bi, CrdTransfWrapper^ coordTransf,
	double rho)
{
	SectionForceDeformation** _sections = new SectionForceDeformation*[s->Length];
	for (int i = 0; i < s->Length; i++)
		_sections[i] = s[i]->_SectionForceDeformation->getCopy();
	_Element = new DispBeamColumn2dWithSensitivity(tag, nd1, nd2, s->Length, _sections, *bi->_BeamIntegration, *coordTransf->_CrdTransf, rho);
}


DispBeamColumn3dThermalWrapper::DispBeamColumn3dThermalWrapper(int tag, int nd1, int nd2,
	array<SectionForceDeformationWrapper^>^ s,
	BeamIntegrationWrapper^ bi, CrdTransfWrapper^ coordTransf)
{
	SectionForceDeformation** _sections = new SectionForceDeformation*[s->Length];
	for (int i = 0; i < s->Length; i++)
		_sections[i] = s[i]->_SectionForceDeformation->getCopy();
	_Element = new DispBeamColumn3dThermal(tag, nd1, nd2, s->Length, _sections, *bi->_BeamIntegration, *coordTransf->_CrdTransf);
}

DispBeamColumn3dWithSensitivityWrapper::DispBeamColumn3dWithSensitivityWrapper(int tag, int nd1, int nd2,
	array<SectionForceDeformationWrapper^>^ s,
	BeamIntegrationWrapper^ bi, CrdTransfWrapper^ coordTransf, double r)
{
	SectionForceDeformation** _sections = new SectionForceDeformation*[s->Length];
	for (int i = 0; i < s->Length; i++)
		_sections[i] = s[i]->_SectionForceDeformation->getCopy();
	_Element = new DispBeamColumn3dWithSensitivity(tag, nd1, nd2, s->Length, _sections, *bi->_BeamIntegration, *coordTransf->_CrdTransf, r);
}

DispBeamColumn2dIntWrapper::DispBeamColumn2dIntWrapper(int tag, int nd1, int nd2,
	array<SectionForceDeformationWrapper^>^ s,
	CrdTransfWrapper^ coordTransf, double c, double rho)
{
	SectionForceDeformation** _sections = new SectionForceDeformation*[s->Length];
	for (int i = 0; i < s->Length; i++)
		_sections[i] = s[i]->_SectionForceDeformation->getCopy();
	_Element = new DispBeamColumn2dInt(tag, nd1, nd2, s->Length, _sections, *coordTransf->_CrdTransf, c,rho);
}

ElasticBeam3dWrapper::ElasticBeam3dWrapper(int tag, double A, double E, double G,
	double Jx, double Iy, double Iz,
	int Nd1, int Nd2, CrdTransfWrapper^ theTransf,
	double rho, int cMass,
	int sectionTag)
{
	_Element = new ElasticBeam3d(tag, A, E, G,
		Jx, Iy, Iz,
		Nd1, Nd2, *theTransf->_CrdTransf,
		rho, cMass,
		sectionTag);
}

ElasticBeam3dWrapper::ElasticBeam3dWrapper(int tag, int Nd1, int Nd2, SectionForceDeformationWrapper^ section,
	CrdTransfWrapper^ theTransf, double rho, int cMass)
{
	_Element = new ElasticBeam3d(tag, Nd1, Nd2, section->_SectionForceDeformation, *theTransf->_CrdTransf, rho, cMass);
}

ElasticTimoshenkoBeam2dWrapper::ElasticTimoshenkoBeam2dWrapper(int tag, int Nd1, int Nd2, double E, double G,
	double A, double Iz, double Avy, CrdTransfWrapper^ theTransf,
	double rho, int cMass)
{
	_Element = new ElasticTimoshenkoBeam2d(tag, Nd1, Nd2, E, G,
		A, Iz, Avy, *theTransf->_CrdTransf,
		rho, cMass);
}


ElasticTimoshenkoBeam3dWrapper::ElasticTimoshenkoBeam3dWrapper(int tag, int Nd1, int Nd2, double E, double G,
	double A, double Jx, double Iy, double Iz, double Avy, double Avz,
	CrdTransfWrapper^ theTransf, double rho, int cMass)
{
	_Element = new ElasticTimoshenkoBeam3d( tag,  Nd1,  Nd2,  E,  G,
		 A,  Jx,  Iy,  Iz,  Avy,  Avz,
		 *theTransf->_CrdTransf,  rho,  cMass);
}


ModElasticBeam2dWrapper::ModElasticBeam2dWrapper(int tag, double A, double E, double I,
	int Nd1, int Nd2, double K11, double K33,
	double K44, CrdTransfWrapper^ theTransf,
	double alpha, double d,
	double rho, int cMass)
{
	_Element = new ModElasticBeam2d( tag,  A,  E,  I,
		 Nd1,  Nd2,  K11,  K33,
		 K44, *theTransf->_CrdTransf,
		 alpha,  d,
		 rho,  cMass);
}

ElasticForceBeamColumn2dWrapper::ElasticForceBeamColumn2dWrapper(int tag, int nodeI, int nodeJ,
	array<SectionForceDeformationWrapper^>^ s,
	BeamIntegrationWrapper^ beamIntegr,
	CrdTransfWrapper^ coordTransf, double rho)
{
	SectionForceDeformation** _sections = new SectionForceDeformation*[s->Length];
	for (int i = 0; i < s->Length; i++)
		_sections[i] = s[i]->_SectionForceDeformation->getCopy();
	_Element = new ElasticForceBeamColumn2d(tag, nodeI, nodeJ, s->Length, _sections, *beamIntegr->_BeamIntegration, *coordTransf->_CrdTransf, rho);
}


ElasticForceBeamColumn3dWrapper::ElasticForceBeamColumn3dWrapper(int tag, int nodeI, int nodeJ,
	array<SectionForceDeformationWrapper^>^ s,
	BeamIntegrationWrapper^ beamIntegr,
	CrdTransfWrapper^ coordTransf, double rho)
{
	SectionForceDeformation** _sections = new SectionForceDeformation*[s->Length];
	for (int i = 0; i < s->Length; i++)
		_sections[i] = s[i]->_SectionForceDeformation->getCopy();
	_Element = new ElasticForceBeamColumn3d(tag, nodeI, nodeJ, s->Length, _sections, *beamIntegr->_BeamIntegration, *coordTransf->_CrdTransf, rho);
}

ElasticForceBeamColumnWarping2dWrapper::ElasticForceBeamColumnWarping2dWrapper(int tag, int nodeI, int nodeJ,
	array<SectionForceDeformationWrapper^>^ s,
	BeamIntegrationWrapper^ beamIntegr,
	CrdTransfWrapper^ coordTransf, double rho)
{
	SectionForceDeformation** _sections = new SectionForceDeformation*[s->Length];
	for (int i = 0; i < s->Length; i++)
		_sections[i] = s[i]->_SectionForceDeformation->getCopy();
	_Element = new ElasticForceBeamColumnWarping2d(tag, nodeI, nodeJ, s->Length, _sections, *beamIntegr->_BeamIntegration, *coordTransf->_CrdTransf, rho);
}


ForceBeamColumn2dThermalWrapper::ForceBeamColumn2dThermalWrapper(int tag, int nodeI, int nodeJ,
	int numSections, array<SectionForceDeformationWrapper^>^ sections,
	BeamIntegrationWrapper^ beamIntegr,
	CrdTransfWrapper^ coordTransf, double rho,
	int maxNumIters, double tolerance)
{
	SectionForceDeformation** _sections = new SectionForceDeformation*[sections->Length];
	for (int i = 0; i < sections->Length; i++)
		_sections[i] = sections[i]->_SectionForceDeformation->getCopy();
	_Element = new ForceBeamColumn2dThermal(tag, nodeI, nodeJ, sections->Length, _sections, *beamIntegr->_BeamIntegration, *coordTransf->_CrdTransf, rho, maxNumIters, tolerance);
}


ForceBeamColumnWarping2dWrapper::ForceBeamColumnWarping2dWrapper(int tag, int nodeI, int nodeJ,
	int numSections, array<SectionForceDeformationWrapper^>^ sections,
	BeamIntegrationWrapper^ beamIntegr,
	CrdTransfWrapper^ coordTransf, double rho,
	int maxNumIters, double tolerance)
{
	SectionForceDeformation** _sections = new SectionForceDeformation*[sections->Length];
	for (int i = 0; i < sections->Length; i++)
		_sections[i] = sections[i]->_SectionForceDeformation->getCopy();
	_Element = new ForceBeamColumnWarping2d(tag, nodeI, nodeJ, sections->Length, _sections, *beamIntegr->_BeamIntegration, *coordTransf->_CrdTransf, rho, maxNumIters, tolerance);
}

ForceBeamColumnCBDI2dWrapper::ForceBeamColumnCBDI2dWrapper(int tag, int nodeI, int nodeJ,
	int numSections, array<SectionForceDeformationWrapper^>^ sections,
	BeamIntegrationWrapper^ beamIntegr,
	CrdTransfWrapper^ coordTransf, double rho, bool includeShear,
	int maxNumIters, double tolerance)
{
	SectionForceDeformation** _sections = new SectionForceDeformation*[sections->Length];
	for (int i = 0; i < sections->Length; i++)
		_sections[i] = sections[i]->_SectionForceDeformation->getCopy();
	_Element = new ForceBeamColumnCBDI2d(tag, nodeI, nodeJ, sections->Length, _sections, *beamIntegr->_BeamIntegration, *coordTransf->_CrdTransf, rho, includeShear, maxNumIters, tolerance);
}


MVLEMWrapper::MVLEMWrapper(int tag,
	double Dens,
	int Nd1, int Nd2,
	array<UniaxialMaterialWrapper^>^ materialsConcrete,
	array<UniaxialMaterialWrapper^>^ materialsSteel,
	array<UniaxialMaterialWrapper^>^ materialsShear,
	array<double>^ Rho,
	array<double>^ thickness,
	array<double>^ width,
	int mm,
	double cc)
{
	UniaxialMaterial** concs = new UniaxialMaterial*[materialsConcrete->Length];
	for (int i = 0; i < materialsConcrete->Length; i++)
	{
		concs[i] = materialsConcrete[i]->_UniaxialMaterial;
	}

	UniaxialMaterial** stls = new UniaxialMaterial*[materialsSteel->Length];
	for (int i = 0; i < materialsSteel->Length; i++)
	{
		stls[i] = materialsSteel[i]->_UniaxialMaterial;
	}

	UniaxialMaterial** shrs = new UniaxialMaterial*[materialsShear->Length];
	for (int i = 0; i < materialsShear->Length; i++)
	{
		shrs[i] = materialsShear[i]->_UniaxialMaterial;
	}

	double* _rho = new double[Rho->Length];
	for (int i = 0; i < Rho->Length; i++)
	{
		_rho[i] = Rho[i];
	}

	double* _thickness = new double[thickness->Length];
	for (int i = 0; i < thickness->Length; i++)
	{
		_thickness[i] = thickness[i];
	}

	double* _width = new double[width->Length];
	for (int i = 0; i < width->Length; i++)
	{
		_width[i] = width[i];
	}
	_Element = new MVLEM(tag, Dens, Nd1, Nd2, concs, stls, shrs, _rho, _thickness, _width, mm, cc);
}

SFI_MVLEMWrapper::SFI_MVLEMWrapper(int tag,					// element tag
	int Nd1, int Nd2,				// end node tags          
	array<NDMaterialWrapper^>^ materials,			// array of material tags
	array<double>^ Thickness,				// array of macro-fiber thickness
	array<double>^ Width,					// array of macro-fiber widths
	int mm,							// number of macro-fibers (RC panels)
	double cc)
{
	NDMaterial** mats = new NDMaterial*[materials->Length];
	for (int i = 0; i < materials->Length; i++)
	{
		mats[i] = materials[i]->_NDMaterial;
	}

	double* _Thickness = new double[Thickness->Length];
	for (int i = 0; i < Thickness->Length; i++)
	{
		_Thickness[i] = Thickness[i];
	}

	double* _Width = new double[Width->Length];
	for (int i = 0; i < Width->Length; i++)
	{
		_Width[i] = Width[i];
	}

	
	_Element = new SFI_MVLEM(tag, Nd1, Nd2, mats, _Thickness, _Width, mm, cc);
}
