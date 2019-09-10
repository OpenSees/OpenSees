#include "stdafx.h"
#include "BearingWrapper.h"

using namespace System;
using namespace OpenSees;
using namespace OpenSees::Elements;
using namespace OpenSees::Elements::BeamIntegrations;
using namespace OpenSees::Materials::Sections;
using namespace OpenSees::Materials::Uniaxials;
using namespace OpenSees::Materials::NDMaterials;

UniaxialMaterial** array2pointer(array<UniaxialMaterialWrapper^>^ theMaterials) {
	UniaxialMaterial** mats = new UniaxialMaterial*[theMaterials->Length];
	for (int i = 0; i < theMaterials->Length; i++)
	{
		mats[i] = theMaterials[i]->_UniaxialMaterial;
	}

	return mats;
}

double* array2pointer1(array<double>^ arr) {
	double* rets = new double[arr->Length];
	for (int i = 0; i < arr->Length; i++)
	{
		rets[i] = arr[i];
	}

	return rets;
}

ElastomericBearingBoucWen2dWrapper::ElastomericBearingBoucWen2dWrapper(int tag, int Nd1, int Nd2,
	double kInit, double qd, double alpha1,
	array<UniaxialMaterialWrapper^>^ theMaterials)
{
	_Element = new ElastomericBearingBoucWen2d(tag, Nd1, Nd2,
		kInit, qd, alpha1,
		array2pointer(theMaterials));
}

ElastomericBearingBoucWen2dWrapper::ElastomericBearingBoucWen2dWrapper(int tag, int Nd1, int Nd2,
	double kInit, double qd, double alpha1,
	array<UniaxialMaterialWrapper^>^ theMaterials,
	VectorWrapper^ y, VectorWrapper^ x,
	double alpha2, double mu,
	double eta, double beta,
	double gamma, double shearDistI,
	int addRayleigh, double mass,
	int maxIter, double tol) {
	_Element = new ElastomericBearingBoucWen2d(tag, Nd1, Nd2,
		kInit, qd, alpha1,
		array2pointer(theMaterials),
		*y->_Vector, *y->_Vector,
		alpha2, mu,
		eta, beta,
		gamma, shearDistI,
		addRayleigh, mass,
		maxIter, tol);

}


ElastomericBearingBoucWen3dWrapper::ElastomericBearingBoucWen3dWrapper(int tag, int Nd1, int Nd2,
	double kInit, double qd, double alpha1,
	array<UniaxialMaterialWrapper^>^ theMaterials, VectorWrapper^ y)
{
	_Element = new ElastomericBearingBoucWen3d(tag, Nd1, Nd2,
		kInit, qd, alpha1,
		array2pointer(theMaterials), *y->_Vector);
}

ElastomericBearingBoucWen3dWrapper::ElastomericBearingBoucWen3dWrapper(int tag, int Nd1, int Nd2,
	double kInit, double qd, double alpha1,
	array<UniaxialMaterialWrapper^>^ theMaterials,
	VectorWrapper^ y, VectorWrapper^ x,
	double alpha2, double mu,
	double eta, double beta,
	double gamma, double shearDistI,
	int addRayleigh, double mass,
	int maxIter, double tol) {
	_Element = new ElastomericBearingBoucWen3d(tag, Nd1, Nd2,
		kInit, qd, alpha1,
		array2pointer(theMaterials),
		*y->_Vector, *y->_Vector,
		alpha2, mu,
		eta, beta,
		gamma, shearDistI,
		addRayleigh, mass,
		maxIter, tol);

}


ElastomericBearingBoucWenMod3dWrapper::ElastomericBearingBoucWenMod3dWrapper(int tag, int Nd1, int Nd2,
	double kInit, double fy,
	double Gr, double Kbulk, double D1, double D2,
	double ts, double tr, int n, double alpha1)
{
	_Element = new ElastomericBearingBoucWenMod3d(tag, Nd1, Nd2,
		kInit, fy,
		Gr, Kbulk, D1, D2,
		ts, tr, n, alpha1);
}


ElastomericBearingBoucWenMod3dWrapper::ElastomericBearingBoucWenMod3dWrapper(int tag, int Nd1, int Nd2,
	double kInit, double fy,
	double Gr, double Kbulk, double D1, double D2,
	double ts, double tr, int n, double alpha1,
	double alpha2, double mu, double eta,
	double beta, double gamma,
	double a1, double a2, double T,
	double b1, double b2, double b3, double b4,
	VectorWrapper^ y, VectorWrapper^ x,
	double shearDistI,
	int addRayleigh, double mass,
	int maxIter, double tol) {
	_Element = new ElastomericBearingBoucWenMod3d(tag, Nd1, Nd2,
		kInit, fy,
		Gr, Kbulk, D1, D2,
		ts, tr, n, alpha1,
		alpha2, mu, eta,
		beta, gamma,
		a1, a2, T,
		b1, b2, b3, b4,
		*y->_Vector, *x->_Vector,
		shearDistI,
		addRayleigh, mass,
		maxIter, tol);

}

ElastomericBearingPlasticity2dWrapper::ElastomericBearingPlasticity2dWrapper(int tag, int Nd1, int Nd2,
	double kInit, double qd, double alpha1,
	array<UniaxialMaterialWrapper^>^ theMaterials) {
	_Element = new ElastomericBearingPlasticity2d(tag, Nd1, Nd2,
		kInit, qd, alpha1,
		array2pointer(theMaterials));
}

ElastomericBearingPlasticity2dWrapper::ElastomericBearingPlasticity2dWrapper(int tag, int Nd1, int Nd2,
	double kInit, double qd, double alpha1,
	array<UniaxialMaterialWrapper^>^ theMaterials,
	VectorWrapper^ y, VectorWrapper^ x,
	double alpha2, double mu,
	double shearDistI,
	int addRayleigh, double mass) {
	_Element = new ElastomericBearingPlasticity2d(tag, Nd1, Nd2,
		kInit, qd, alpha1,
		array2pointer(theMaterials), *y->_Vector, *x->_Vector, alpha2, mu,
		shearDistI,
		addRayleigh, mass);
}

//ElastomericBearingPlasticity2dWrapper::ElastomericBearingPlasticity2dWrapper(int tag, int Nd1, int Nd2,
//	double kInit, double qd, double alpha1,
//	array<UniaxialMaterialWrapper^>^ theMaterials) {
//	_Element = new ElastomericBearingPlasticity2d(tag, Nd1, Nd2,
//		kInit, qd, alpha1,
//		array2pointer(theMaterials));
//}
//
//ElastomericBearingPlasticity2dWrapper::ElastomericBearingPlasticity2dWrapper(int tag, int Nd1, int Nd2,
//	double kInit, double qd, double alpha1,
//	array<UniaxialMaterialWrapper^>^ theMaterials,
//	VectorWrapper^ y, VectorWrapper^ x,
//	double alpha2, double mu,
//	double shearDistI,
//	int addRayleigh, double mass){
//	_Element = new ElastomericBearingPlasticity2d(tag, Nd1, Nd2,
//		kInit, qd, alpha1,
//		array2pointer(theMaterials),
//		*y->_Vector, *x->_Vector,
//		alpha2, mu,
//		shearDistI,
//		addRayleigh, mass);
//}

ElastomericBearingPlasticity3dWrapper::ElastomericBearingPlasticity3dWrapper(int tag, int Nd1, int Nd2,
	double kInit, double qd, double alpha1,
	array<UniaxialMaterialWrapper^>^ theMaterials, VectorWrapper^ y) {
	_Element = new ElastomericBearingPlasticity3d(tag, Nd1, Nd2,
		kInit, qd, alpha1,
		array2pointer(theMaterials), *y->_Vector);
}

ElastomericBearingPlasticity3dWrapper::ElastomericBearingPlasticity3dWrapper(int tag, int Nd1, int Nd2,
	double kInit, double qd, double alpha1,
	array<UniaxialMaterialWrapper^>^ theMaterials,
	VectorWrapper^ y, VectorWrapper^ x,
	double alpha2, double mu,
	double shearDistI,
	int addRayleigh, double mass) {
	_Element = new ElastomericBearingPlasticity3d(tag, Nd1, Nd2,
		kInit, qd, alpha1,
		array2pointer(theMaterials),
		*y->_Vector, *x->_Vector,
		alpha2, mu,
		shearDistI,
		addRayleigh, mass);
}

ElastomericBearingUFRP2dWrapper::ElastomericBearingUFRP2dWrapper(int tag, int Nd1, int Nd2, double uy,
	double a1, double a2, double a3, double a4, double a5,
	double b, double c, array<UniaxialMaterialWrapper^>^ theMaterials) {
	_Element = new ElastomericBearingUFRP2d(tag, Nd1, Nd2, uy,
		a1, a2, a3, a4, a5,
		b, c,
		array2pointer(theMaterials));
}

ElastomericBearingUFRP2dWrapper::ElastomericBearingUFRP2dWrapper(int tag, int Nd1, int Nd2, double uy,
	double a1, double a2, double a3, double a4, double a5,
	double b, double c, array<UniaxialMaterialWrapper^>^ theMaterials,
	VectorWrapper^ y, VectorWrapper^ x,
	double eta, double beta, double gamma,
	double shearDistI, int addRayleigh, double mass,
	int maxIter, double tol) {
	_Element = new ElastomericBearingUFRP2d(tag, Nd1, Nd2, uy,
		a1, a2, a3, a4, a5,
		b, c,
		array2pointer(theMaterials), *y->_Vector, *x->_Vector,
		eta, beta, gamma,
		shearDistI, addRayleigh, mass,
		maxIter, tol);
}


ElastomericXWrapper::ElastomericXWrapper(int eleTag, int Nd1, int Nd2, double qd, double alpha, double Gr, double Kbulk,
	double D1, double D2, double ts, double tr, double n, VectorWrapper^ y) {
	_Element = new ElastomericX(eleTag, Nd1, Nd2, qd, alpha, Gr, Kbulk,
		D1, D2, ts, tr, n, *y->_Vector);
}

ElastomericXWrapper::ElastomericXWrapper(int eleTag, int Nd1, int Nd2, double qd, double alpha, double Gr, double Kbulk,
	double D1, double D2, double ts, double tr, double n, VectorWrapper^ y, VectorWrapper^ x,
	double kc, double PhiM, double ac, double sDratio, double m,
	double cd, double tc, int tag1, int tag2, int tag3, int tag4) {
	_Element = new ElastomericX( eleTag,  Nd1,  Nd2,  qd,  alpha,  Gr,  Kbulk,
		 D1,  D2,  ts,  tr,  n, *y->_Vector, *x->_Vector,
		 kc,  PhiM,  ac,  sDratio,  m,
		 cd,  tc,  tag1,  tag2,  tag3,  tag4);
}


HDRWrapper::HDRWrapper(int tag, int Nd1, int Nd2, double Gr, double Kbulk, double D1, double D2, double ts,
	double tr, int n, double a1, double a2, double a3, double b1, double b2, double b3,
	double c1, double c2, double c3, double c4, VectorWrapper^ y) {
	_Element = new HDR(tag, Nd1, Nd2, Gr, Kbulk, D1, D2, ts,
		tr, n, a1, a2, a3, b1, b2, b3,
		c1, c2, c3, c4, *y->_Vector);
}

HDRWrapper::HDRWrapper(int tag, int Nd1, int Nd2, double Gr, double Kbulk, double D1, double D2, double ts,
	double tr, int n, double a1, double a2, double a3, double b1, double b2, double b3,
	double c1, double c2, double c3, double c4, VectorWrapper^ y, VectorWrapper^ x,
	double kc, double PhiM, double ac, double sDratio, double m,
	double tc) {
	_Element = new HDR(tag, Nd1, Nd2, Gr, Kbulk, D1, D2, ts,
		tr, n, a1, a2, a3, b1, b2, b3,
		c1, c2, c3, c4, *y->_Vector, *x->_Vector,
		kc, PhiM, ac, sDratio, m,
		tc);
}

LeadRubberXWrapper::LeadRubberXWrapper(int eleTag, int Nd1, int Nd2, double qd, double alpha, double Gr, double Kbulk,
	double D1, double D2, double ts, double tr, double n, VectorWrapper^ y) {

	_Element = new LeadRubberX(eleTag, Nd1, Nd2, qd, alpha, Gr, Kbulk,
		D1, D2, ts, tr, n, *y->_Vector);
}

LeadRubberXWrapper::LeadRubberXWrapper(int eleTag, int Nd1, int Nd2, double qd, double alpha, double Gr, double Kbulk,
	double D1, double D2, double ts, double tr, double n, VectorWrapper^ y, VectorWrapper^ x,
	double kc, double PhiM, double ac, double sDratio, double m,
	double cd, double tc, double qL, double cL, double kS,
	double aS, int tag1, int tag2, int tag3, int tag4, int tag5) {

	_Element = new LeadRubberX(eleTag, Nd1, Nd2, qd, alpha, Gr, Kbulk,
		D1, D2, ts, tr, n, *y->_Vector, *x->_Vector,
		kc, PhiM, ac, sDratio, m,
		cd, tc, qL, cL, kS,
		aS, tag1, tag2, tag3, tag4, tag5);
}

CoulombWrapper::CoulombWrapper(int tag, double mu) {

	_FrictionModel = new Coulomb( tag,  mu);
}

VelDependentWrapper::VelDependentWrapper(int tag, double muSlow, double muFast, double transRate) {

	_FrictionModel = new VelDependent( tag,  muSlow,  muFast,  transRate);
}

VelDepMultiLinearWrapper::VelDepMultiLinearWrapper(int tag,
	VectorWrapper^ velocityPoints,
	VectorWrapper^ frictionPoints) {

	_FrictionModel = new VelDepMultiLinear(tag, *velocityPoints->_Vector, *frictionPoints->_Vector);
}

VelNormalFrcDepWrapper::VelNormalFrcDepWrapper(int tag,
	double aSlow, double nSlow, double aFast, double nFast,
	double alpha0, double alpha1, double alpha2, double maxMuFact) {

	_FrictionModel = new VelNormalFrcDep(tag,
		aSlow, nSlow, aFast, nFast,
		alpha0, alpha1, alpha2, maxMuFact);
}

VelPressureDepWrapper::VelPressureDepWrapper(int tag, double muSlow, double muFast0, double A,
	double deltaMu, double alpha, double transRate) {

	_FrictionModel = new VelPressureDep(tag, muSlow, muFast0, A,
		deltaMu, alpha, transRate);
}

FlatSliderSimple2dWrapper::FlatSliderSimple2dWrapper(int tag, int Nd1, int Nd2,
	FrictionModelWrapper^ theFrnMdl, double kInit,
	array<UniaxialMaterialWrapper^>^ theMaterials) {

	_Element = new FlatSliderSimple2d(tag, Nd1, Nd2,
		*theFrnMdl->_FrictionModel, kInit,
		array2pointer(theMaterials));
}

FlatSliderSimple2dWrapper::FlatSliderSimple2dWrapper(int tag, int Nd1, int Nd2,
	FrictionModelWrapper^ theFrnMdl, double kInit,
	array<UniaxialMaterialWrapper^>^ theMaterials,
	VectorWrapper^ y, VectorWrapper^ x,
	double shearDistI,
	int addRayleigh, double mass,
	int maxIter, double tol) {

	_Element = new FlatSliderSimple2d(tag, Nd1, Nd2,
		*theFrnMdl->_FrictionModel, kInit,
		array2pointer(theMaterials),
		*y->_Vector, *x->_Vector,
		shearDistI,
		addRayleigh, mass,
		maxIter, tol);
}

FlatSliderSimple3dWrapper::FlatSliderSimple3dWrapper(int tag, int Nd1, int Nd2,
	FrictionModelWrapper^ theFrnMdl, double kInit,
	array<UniaxialMaterialWrapper^>^ theMaterials) {

	_Element = new FlatSliderSimple3d(tag, Nd1, Nd2,
		*theFrnMdl->_FrictionModel, kInit,
		array2pointer(theMaterials));
}

FlatSliderSimple3dWrapper::FlatSliderSimple3dWrapper(int tag, int Nd1, int Nd2,
	FrictionModelWrapper^ theFrnMdl, double kInit,
	array<UniaxialMaterialWrapper^>^ theMaterials,
	VectorWrapper^ y, VectorWrapper^ x,
	double shearDistI,
	int addRayleigh, double mass,
	int maxIter, double tol) {

	_Element = new FlatSliderSimple3d(tag, Nd1, Nd2,
		*theFrnMdl->_FrictionModel, kInit,
		array2pointer(theMaterials),
		*y->_Vector, *x->_Vector,
		shearDistI,
		addRayleigh, mass,
		maxIter, tol);
}

FPBearingPTVWrapper::FPBearingPTVWrapper(int tag, int Nd1, int Nd2, double MuReference,
	int IsPDependent, double refP, int IsTDependent, double Diffusivity_Steel,
	double Conductivity_Steel, int IsVDependent, double rate_v_mu, double Reff,
	double r_Contact, double kInit,
	UniaxialMaterialWrapper^ theMatA, UniaxialMaterialWrapper^ theMatB,
	UniaxialMaterialWrapper^ theMatC, UniaxialMaterialWrapper^ theMatD) {

	_Element = new FPBearingPTV( tag,  Nd1,  Nd2,  MuReference,
		 IsPDependent,  refP,  IsTDependent,  Diffusivity_Steel,
		 Conductivity_Steel,  IsVDependent,  rate_v_mu,  Reff,
		 r_Contact,  kInit,
		*theMatA->_UniaxialMaterial, *theMatB->_UniaxialMaterial,
		*theMatC->_UniaxialMaterial, *theMatD->_UniaxialMaterial);
}

FPBearingPTVWrapper::FPBearingPTVWrapper(int tag, int Nd1, int Nd2, double MuReference,
	int IsPDependent, double refP, int IsTDependent, double Diffusivity_Steel,
	double Conductivity_Steel, int IsVDependent, double rate_v_mu, double Reff,
	double r_Contact, double kInit,
	UniaxialMaterialWrapper^ theMatA, UniaxialMaterialWrapper^ theMatB,
	UniaxialMaterialWrapper^ theMatC, UniaxialMaterialWrapper^ theMatD,
	VectorWrapper^ x, VectorWrapper^  y,
	double shearDistI,
	int addRayleigh, double mass,
	int maxIter, double tol, int unit) {

	_Element = new FPBearingPTV(tag, Nd1, Nd2, MuReference,
		IsPDependent, refP, IsTDependent, Diffusivity_Steel,
		Conductivity_Steel, IsVDependent, rate_v_mu, Reff,
		r_Contact, kInit,
		*theMatA->_UniaxialMaterial, *theMatB->_UniaxialMaterial,
		*theMatC->_UniaxialMaterial, *theMatD->_UniaxialMaterial,
		*x->_Vector, *y->_Vector,
		shearDistI,
		addRayleigh, mass,
		maxIter, tol, unit);
}

MultiFP2dWrapper::MultiFP2dWrapper(int tag,
	int Nd1, int Nd2,
	UniaxialMaterialWrapper^ theFrictionModel,
	UniaxialMaterialWrapper^ theVerticalModel,
	double w0, int axialCase) {

	_Element = new MultiFP2d(tag,
		Nd1, Nd2,
		theFrictionModel->_UniaxialMaterial,
		theVerticalModel->_UniaxialMaterial,
		w0, axialCase);
}

MultiFP2dWrapper::MultiFP2dWrapper(int tag,
	int Nd1, int Nd2,
	int type,
	VectorWrapper^ R,
	VectorWrapper^ h,
	VectorWrapper^ D,
	VectorWrapper^ d,
	VectorWrapper^ mu,
	double Kvert,
	double w0, int axialCase) {

	_Element = new MultiFP2d( tag,
		 Nd1,  Nd2,
		 type,
		*R->_Vector,
		*h->_Vector,
		*D->_Vector,
		*d->_Vector,
		*mu->_Vector,
		 Kvert,
		 w0,  axialCase);
}

RJWatsonEQS2dWrapper::RJWatsonEQS2dWrapper(int tag, int Nd1, int Nd2,
	FrictionModelWrapper^ theFrnMdl, double kInit,
	array<UniaxialMaterialWrapper^>^ theMaterials) {

	_Element = new RJWatsonEQS2d(tag, Nd1, Nd2,
		*theFrnMdl->_FrictionModel, kInit,
		array2pointer(theMaterials));
}

RJWatsonEQS2dWrapper::RJWatsonEQS2dWrapper(int tag, int Nd1, int Nd2,
	FrictionModelWrapper^ theFrnMdl, double kInit,
	array<UniaxialMaterialWrapper^>^ theMaterials,
	VectorWrapper^  y, VectorWrapper^  x,
	double shearDistI,
	int addRayleigh, double mass,
	int maxIter, double tol,
	double kFactUplift) {

	_Element = new RJWatsonEQS2d(tag, Nd1, Nd2,
		*theFrnMdl->_FrictionModel, kInit,
		array2pointer(theMaterials), *y->_Vector, *x->_Vector,
		shearDistI,
		addRayleigh, mass,
		maxIter, tol,
		kFactUplift);
}

RJWatsonEQS3dWrapper::RJWatsonEQS3dWrapper(int tag, int Nd1, int Nd2,
	FrictionModelWrapper^ theFrnMdl, double kInit,
	array<UniaxialMaterialWrapper^>^ theMaterials) {

	_Element = new RJWatsonEQS3d(tag, Nd1, Nd2,
		*theFrnMdl->_FrictionModel, kInit,
		array2pointer(theMaterials));
}

RJWatsonEQS3dWrapper::RJWatsonEQS3dWrapper(int tag, int Nd1, int Nd2,
	FrictionModelWrapper^ theFrnMdl, double kInit,
	array<UniaxialMaterialWrapper^>^ theMaterials,
	VectorWrapper^  y, VectorWrapper^  x,
	double shearDistI,
	int addRayleigh, double mass,
	int maxIter, double tol,
	double kFactUplift) {

	_Element = new RJWatsonEQS3d(tag, Nd1, Nd2,
		*theFrnMdl->_FrictionModel, kInit,
		array2pointer(theMaterials), *y->_Vector, *x->_Vector,
		shearDistI,
		addRayleigh, mass,
		maxIter, tol,
		kFactUplift);
}


SingleFPSimple2dWrapper::SingleFPSimple2dWrapper(int tag, int Nd1, int Nd2,
	FrictionModelWrapper^ theFrnMdl, double Reff, double kInit,
	array<UniaxialMaterialWrapper^>^ theMaterials) {

	_Element = new SingleFPSimple2d(tag, Nd1, Nd2,
		*theFrnMdl->_FrictionModel, Reff, kInit,
		array2pointer(theMaterials));
}

SingleFPSimple2dWrapper::SingleFPSimple2dWrapper(int tag, int Nd1, int Nd2,
	FrictionModelWrapper^ theFrnMdl, double Reff, double kInit,
	array<UniaxialMaterialWrapper^>^ theMaterials,
	VectorWrapper^  y, VectorWrapper^  x,
	double shearDistI, int addRayleigh,
	int inclVertDisp , double mass,
	int maxIter, double tol,
	double kFactUplift) {

	_Element = new SingleFPSimple2d(tag, Nd1, Nd2,
		*theFrnMdl->_FrictionModel, Reff, kInit,
		array2pointer(theMaterials), *y->_Vector, *x->_Vector,
		shearDistI, addRayleigh,
		inclVertDisp, mass,
		maxIter, tol,
		kFactUplift);
}

SingleFPSimple3dWrapper::SingleFPSimple3dWrapper(int tag, int Nd1, int Nd2,
	FrictionModelWrapper^ theFrnMdl, double Reff, double kInit,
	array<UniaxialMaterialWrapper^>^ theMaterials) {

	_Element = new SingleFPSimple3d(tag, Nd1, Nd2,
		*theFrnMdl->_FrictionModel, Reff, kInit,
		array2pointer(theMaterials));
}

SingleFPSimple3dWrapper::SingleFPSimple3dWrapper(int tag, int Nd1, int Nd2,
	FrictionModelWrapper^ theFrnMdl, double Reff, double kInit,
	array<UniaxialMaterialWrapper^>^ theMaterials,
	VectorWrapper^  y, VectorWrapper^  x,
	double shearDistI,
	int addRayleigh, int inclVertDisp, double mass,
	int maxIter, double tol,
	double kFactUplift) {

	_Element = new SingleFPSimple3d(tag, Nd1, Nd2,
		*theFrnMdl->_FrictionModel, Reff, kInit,
		array2pointer(theMaterials), *y->_Vector, *x->_Vector,
		shearDistI,
		addRayleigh, inclVertDisp, mass,
		maxIter, tol,
		kFactUplift);
}



TFP_BearingWrapper::TFP_BearingWrapper(int tag,
	int Nd1, int Nd2,
	array<double>^ r,
	array<double>^ dio,
	array<double>^ di,
	array<double>^ mu,
	array<double>^ h,
	double H0,
	double a,
	double K,
	double vYield) {

	_Element = new TFP_Bearing(tag,
		Nd1, Nd2,
		array2pointer1(r),
		array2pointer1(dio),
		array2pointer1(di),
		array2pointer1(mu),
		array2pointer1(h),
		H0,
		a,
		K,
		vYield);
}

TFP_Bearing2dWrapper::TFP_Bearing2dWrapper(int tag,
	int Nd1, int Nd2,
	array<double>^ r,
	array<double>^ dio,
	array<double>^ di,
	array<double>^ mu,
	array<double>^ h,
	double H0,
	double a,
	double K,
	double vYield) {

	_Element = new TFP_Bearing2d(tag,
		Nd1, Nd2,
		array2pointer1(r),
		array2pointer1(dio),
		array2pointer1(di),
		array2pointer1(mu),
		array2pointer1(h),
		H0,
		a,
		K,
		vYield);
}

TPB1DWrapper::TPB1DWrapper(int tag,
	int Nd1,
	int Nd2,
	int dir,
	array<double>^ mu,
	array<double>^ R,
	array<double>^ h,
	array<double>^ D,
	array<double>^ d,
	double W) {

	_Element = new TPB1D(tag,
		Nd1,
		Nd2,
		dir,
		array2pointer1(mu),
		array2pointer1(R),
		array2pointer1(h),
		array2pointer1(D),
		array2pointer1(d),
		W);
}

TripleFrictionPendulumWrapper::TripleFrictionPendulumWrapper(int tag,
	int Nd1, int Nd2,
	array<FrictionModelWrapper^>^ theFrnMdls,
	array<UniaxialMaterialWrapper^>^ theMaterials,
	double L1,
	double L2,
	double L3,
	double Ubar1,
	double Ubar2,
	double Ubar3,
	double W,
	double Uy,
	double Kvt,
	double minFv,
	double tol) {

	FrictionModel** rets = new FrictionModel*[theFrnMdls->Length];
	for (int i = 0; i < theFrnMdls->Length; i++)
	{
		rets[i] = theFrnMdls[i]->_FrictionModel;
	}

	_Element = new TripleFrictionPendulum(tag,
		Nd1, Nd2,
		rets,
		array2pointer(theMaterials),
		L1,
		L2,
		L3,
		Ubar1,
		Ubar2,
		Ubar3,
		W,
		Uy,
		Kvt,
		minFv,
		tol);
}

KikuchiBearingWrapper::KikuchiBearingWrapper(int Tag, int Nd1, int Nd2,
	int Shape, double Size, double TotalRubber, double TotalHeight,
	int NMSS, UniaxialMaterialWrapper^ MatMSS, double LimDisp,
	int NMNS, UniaxialMaterialWrapper^ MatMNS, double Lambda,
	VectorWrapper^ OriYp, VectorWrapper^ OriX, double Mass,
	bool IfPDInput, bool IfTilt,
	double AdjCi, double AdjCj,
	bool IfBalance, double LimFo, double LimFi, int NIter) {

	_Element = new KikuchiBearing(Tag, Nd1, Nd2,
		Shape, Size, TotalRubber, TotalHeight,
		NMSS, MatMSS->_UniaxialMaterial, LimDisp,
		NMNS, MatMNS->_UniaxialMaterial, Lambda,
		*OriYp->_Vector, *OriX->_Vector, Mass,
		IfPDInput, IfTilt,
		AdjCi, AdjCj,
		IfBalance, LimFo, LimFi, NIter);
}

MultipleNormalSpringWrapper::MultipleNormalSpringWrapper(int Tag, int Nd1, int Nd2,
	int NDivide,
	UniaxialMaterialWrapper^ Material,
	int Shape,
	double Size,
	double Lambda,
	VectorWrapper^  OriYp, VectorWrapper^  OriX,
	double Mass) {

	_Element = new MultipleNormalSpring(Tag, Nd1, Nd2,
		NDivide,
		Material->_UniaxialMaterial,
		Shape,
		Size,
		Lambda,
		*OriYp->_Vector, *OriX->_Vector,
		Mass);
}

MultipleShearSpringWrapper::MultipleShearSpringWrapper(int Tag, int Nd1, int Nd2,
	int NSpring,
	UniaxialMaterialWrapper^ Material,
	double LimDisp,
	VectorWrapper^  OriYp, VectorWrapper^  OriX,
	double Mass) {

	_Element = new MultipleShearSpring(Tag, Nd1, Nd2,
		NSpring,
		Material->_UniaxialMaterial,
		LimDisp,
		*OriYp->_Vector, *OriX->_Vector,
		Mass);
}

MultipleShearSpringWrapper::MultipleShearSpringWrapper(int Tag, int Nd1, int Nd2,
	array<UniaxialMaterialWrapper^>^ Materials,
	int NSpring,
	double LimDisp,
	VectorWrapper^  OriYp, VectorWrapper^  OriX,
	double Mass) {

	_Element = new MultipleShearSpring(Tag, Nd1, Nd2,
		array2pointer(Materials),
		NSpring,
		LimDisp,
		*OriYp->_Vector, *OriX->_Vector,
		Mass);
}

YamamotoBiaxialHDRWrapper::YamamotoBiaxialHDRWrapper(int Tag, int Nd1, int Nd2, int Tp, double DDo, double DDi, double Hr,
	double Cr, double Cs,
	VectorWrapper^  OriYp, VectorWrapper^  OriX,
	double Mass) {

	_Element = new YamamotoBiaxialHDR(Tag, Nd1, Nd2, Tp, DDo, DDi, Hr,
		Cr, Cs,
		*OriYp->_Vector, *OriX->_Vector,
		Mass);
}