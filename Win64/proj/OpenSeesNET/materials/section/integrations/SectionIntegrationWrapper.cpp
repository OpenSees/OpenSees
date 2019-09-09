#include "stdafx.h"
#include "SectionIntegrationWrapper.h"

using namespace System;
using namespace System::Collections::Generic;
using namespace OpenSees;
using namespace OpenSees::Materials;
using namespace OpenSees::Materials::Uniaxials;
using namespace OpenSees::Materials::Sections;
using namespace OpenSees::Materials::Sections::Integrations;
using namespace OpenSees::Materials::Sections::Repres;


RCSectionIntegrationWrapper::RCSectionIntegrationWrapper(double d, double b, double Atop, double Abottom,
	double Aside, double cover,
	int Nfcore, int Nfcover, int Nfs) {
	_SectionIntegration = new RCSectionIntegration(d, b, Atop, Abottom, Aside, cover, Nfcore, Nfcover, Nfs);
}

RCTBeamSectionIntegrationWrapper::RCTBeamSectionIntegrationWrapper(double d, double bw, double beff, double hf, double Atop, double Abottom, double flcov,
	double wcov, int Nflcover, int Nwcover, int Nflcore, int Nwcore, int NsteelTop, int NsteelBottom) {
	_SectionIntegration = new RCTBeamSectionIntegration(d, bw, beff, hf, Atop, Abottom, flcov, wcov, Nflcover, Nwcover, Nflcore, Nwcore, NsteelTop, NsteelBottom);
}

WideFlangeSectionIntegrationWrapper::WideFlangeSectionIntegrationWrapper(double d, double tw, double bf, double tf,
	int Nfdw, int Nftf) {
	_SectionIntegration = new WideFlangeSectionIntegration(d, tw, bf, tf, Nfdw, Nftf);
}



