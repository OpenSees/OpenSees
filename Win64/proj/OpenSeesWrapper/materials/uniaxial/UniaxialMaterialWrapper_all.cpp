#include "stdafx.h"
#include "UniaxialMaterialWrapper.h"
#include "UniaxialMaterialWrapper_all.h"

using namespace OpenSees;
using namespace OpenSees::Materials;
using namespace OpenSees::Materials::Uniaxials;
using namespace System::Runtime::InteropServices;

double* array2pointer5(array<double>^ arr) {
	double* rets = new double[arr->Length];
	for (int i = 0; i < arr->Length; i++)
	{
		rets[i] = arr[i];
	}

	return rets;
}

UniaxialMaterial** array2pointer5(array<UniaxialMaterialWrapper^>^ theMaterials) {
	UniaxialMaterial** mats = new UniaxialMaterial*[theMaterials->Length];
	for (int i = 0; i < theMaterials->Length; i++)
	{
		mats[i] = theMaterials[i]->_UniaxialMaterial;
	}

	return mats;
}

ConcreteL01Wrapper::ConcreteL01Wrapper(int tag, double fpc, double eco) {
	_UniaxialMaterial = new ConcreteL01(tag, fpc, eco);
}

ConcreteZ01Wrapper::ConcreteZ01Wrapper(int tag, double fpc, double eco) {
	_UniaxialMaterial = new ConcreteZ01(tag, fpc, eco);
}

SteelZ01Wrapper::SteelZ01Wrapper(int tag, double fy, double E0, double fpc, double rou,
	double ac,
	double rc) {
	_UniaxialMaterial = new SteelZ01(tag, fy, E0, fpc, rou, ac, rc);
}

TendonL01Wrapper::TendonL01Wrapper(int tag, double fy, double E0, double fpu, double rou, double epsp,
	double ac,
	double rc) {
	_UniaxialMaterial = new TendonL01(tag, fy, E0, fpu, rou, epsp, ac, rc);
}

//NDMaterial** array2pointer5(array<NDMaterialWrapper^>^ theMaterials) {
//	NDMaterial** mats = new NDMaterial*[theMaterials->Length];
//	for (int i = 0; i < theMaterials->Length; i++)
//	{
//		mats[i] = theMaterials[i]->_NDMaterial;
//	}
//
//	return mats;
//}

//UniaxialMaterialWrapper::UniaxialMaterialWrapper()
//{
//
//}
//
//UniaxialMaterialWrapper::UniaxialMaterialWrapper(int tag)
//{
//
//}

ArctangentBackboneWrapper::ArctangentBackboneWrapper(int tag, double K1, double gammaY, double alpha) {
	_HystereticBackbone = new ArctangentBackbone(tag, K1, gammaY, alpha);
}

ManderBackboneWrapper::ManderBackboneWrapper(int tag, double fc, double epsc, double Ec) {
	_HystereticBackbone = new ManderBackbone(tag, fc, epsc, Ec);
}

RaynorBackboneWrapper::RaynorBackboneWrapper(int tag, double es, double f1, double f2, double epsh, double epsm, double c1, double ey) {
	_HystereticBackbone = new 				RaynorBackbone(tag, es, f1, f2, epsh, epsm, c1, ey);
}

ReeseSandBackboneWrapper::ReeseSandBackboneWrapper(int tag, double kx, double ym, double pm,
	double yu, double pu) {
	_HystereticBackbone = new ReeseSandBackbone(tag, kx, ym, pm, yu, pu);
}

ReeseSoftClayBackboneWrapper::ReeseSoftClayBackboneWrapper(int tag, double pu, double y50, double n) {
	_HystereticBackbone = new ReeseSoftClayBackbone(tag, pu, y50, n);
}

ReeseStiffClayBelowWSWrapper::ReeseStiffClayBelowWSWrapper(int tag, double esi, double y, double as, double pc) {
	_HystereticBackbone = new ReeseStiffClayBelowWS(tag, esi, y, as, pc);
}

TrilinearBackboneWrapper::TrilinearBackboneWrapper(int tag, double e1, double s1,
	double e2, double s2, double e3, double s3) {
	_HystereticBackbone = new TrilinearBackbone(tag, e1, s1, e2, s2, e3, s3);
}

TrilinearBackboneWrapper::TrilinearBackboneWrapper(int tag, double e1, double s1,
	double e2, double s2) {
	_HystereticBackbone = new TrilinearBackbone(tag, e1, s1, e2, s2);
}

BackboneMaterialWrapper::BackboneMaterialWrapper(int tag, HystereticBackboneWrapper^ backbone) {
	_UniaxialMaterial = new BackboneMaterial(tag, *backbone->_HystereticBackbone);
}

DrainBilinearMaterialWrapper::DrainBilinearMaterialWrapper(int tag,
	double E, double fyp, double fyn, double alpha,
	double ecaps, double ecapk, double ecapa, double ecapd,
	double cs, double ck, double ca, double cd,
	double capSlope, double capDispP, double capDispN, double res, double beto) {
	_UniaxialMaterial = new DrainBilinearMaterial(tag, E, fyp, fyn, alpha, ecaps, ecapk, ecapa, ecapd, cs, ck, ca, cd, capSlope, capDispP, capDispN, res, beto);
}

DrainBilinearMaterialWrapper::DrainBilinearMaterialWrapper(int tag, VectorWrapper^ input, double beto) {
	_UniaxialMaterial = new DrainBilinearMaterial(tag, *input->_Vector, beto);
}

DrainClough1MaterialWrapper::DrainClough1MaterialWrapper(int tag,
	double E, double fyp, double fyn, double alpha,
	double ecaps, double ecapk, double ecapa, double ecapd,
	double cs, double ck, double ca, double cd,
	double capSlope, double capDispP, double capDispN, double res, double beto) {
	_UniaxialMaterial = new DrainClough1Material(tag, E, fyp, fyn, alpha, ecaps, ecapk, ecapa, ecapd, cs, ck, ca, cd, capSlope, capDispP, capDispN, res, beto);
}

DrainClough1MaterialWrapper::DrainClough1MaterialWrapper(int tag, VectorWrapper^ input, double beto) {
	_UniaxialMaterial = new DrainClough1Material(tag, *input->_Vector, beto);
}

DrainClough2MaterialWrapper::DrainClough2MaterialWrapper(int tag,
	double E, double fyp, double fyn, double alpha,
	double ecaps, double ecapk, double ecapa, double ecapd,
	double cs, double ck, double ca, double cd,
	double capSlope, double capDispP, double capDispN, double res, double beto) {
	_UniaxialMaterial = new DrainClough2Material(tag, E, fyp, fyn, alpha, ecaps, ecapk, ecapa, ecapd, cs, ck, ca, cd, capSlope, capDispP, capDispN, res, beto);
}

DrainClough2MaterialWrapper::DrainClough2MaterialWrapper(int tag, VectorWrapper^ input, double beto) {
	_UniaxialMaterial = new DrainClough2Material(tag, *input->_Vector, beto);
}

DrainHardeningMaterialWrapper::DrainHardeningMaterialWrapper(int tag,
	double E, double sigY, double Hiso, double Hkin, double beto) {
	_UniaxialMaterial = new DrainHardeningMaterial(tag, E, sigY, Hiso, Hkin, beto);
}

DrainPinch1MaterialWrapper::DrainPinch1MaterialWrapper(int tag,
	double E, double fyp, double fyn, double alpha,
	double ecaps, double ecapk, double ecapa, double ecapd,
	double cs, double ck, double ca, double cd,
	double capSlope, double capDispP, double capDispN,
	double fpp, double fpn, double pinch, double res, double beto) {
	_UniaxialMaterial = new DrainPinch1Material(tag, E, fyp, fyn, alpha, ecaps, ecapk, ecapa, ecapd, cs, ck, ca, cd, capSlope, capDispP, capDispN, fpp, fpn, pinch, res, beto);
}

FedeasBond1MaterialWrapper::FedeasBond1MaterialWrapper(int tag,
	double u1p, double q1p, double u2p, double u3p, double q3p,
	double u1n, double q1n, double u2n, double u3n, double q3n,
	double s0, double bb) {
	_UniaxialMaterial = new FedeasBond1Material(tag, u1p, q1p, u2p, u3p, q3p, u1n, q1n, u2n, u3n, q3n, s0, bb);
}

FedeasBond1MaterialWrapper::FedeasBond1MaterialWrapper(int tag, VectorWrapper^ input) {
	_UniaxialMaterial = new FedeasBond1Material(tag, *input->_Vector);
}

FedeasBond2MaterialWrapper::FedeasBond2MaterialWrapper(int tag,
	double u1p, double q1p, double u2p, double u3p, double q3p,
	double u1n, double q1n, double u2n, double u3n, double q3n,
	double s0, double bb, double alp, double aln) {
	_UniaxialMaterial = new FedeasBond2Material(tag, u1p, q1p, u2p, u3p, q3p, u1n, q1n, u2n, u3n, q3n, s0, bb, alp, aln);
}

FedeasBond2MaterialWrapper::FedeasBond2MaterialWrapper(int tag, VectorWrapper^ data) {
	_UniaxialMaterial = new FedeasBond2Material(tag, *data->_Vector);
}

FedeasConcr1MaterialWrapper::FedeasConcr1MaterialWrapper(int tag,
	double fc, double ec, double fu, double eu) {
	_UniaxialMaterial = new FedeasConcr1Material(tag, fc, ec, fu, eu);
}

FedeasConcr1MaterialWrapper::FedeasConcr1MaterialWrapper(int tag, VectorWrapper^ data) {
	_UniaxialMaterial = new FedeasConcr1Material(tag, *data->_Vector);
}

FedeasConcr2MaterialWrapper::FedeasConcr2MaterialWrapper(int tag,
	double fc, double ec, double fu, double eu,
	double ratio, double ft, double Ets) {
	_UniaxialMaterial = new FedeasConcr2Material(tag, fc, ec, fu, eu, ratio, ft, Ets);
}

FedeasConcr2MaterialWrapper::FedeasConcr2MaterialWrapper(int tag, VectorWrapper^ data) {
	_UniaxialMaterial = new FedeasConcr2Material(tag, *data->_Vector);
}

FedeasConcr3MaterialWrapper::FedeasConcr3MaterialWrapper(int tag, double fc, double ec, double fu, double eu,
	double rat, double ft, double epst0, double ft0, double beta, double epstu) {
	_UniaxialMaterial = new FedeasConcr3Material(tag, fc, ec, fu, eu, rat, ft, epst0, ft0, beta, epstu);
}

FedeasConcr3MaterialWrapper::FedeasConcr3MaterialWrapper(int tag, VectorWrapper^ data) {
	_UniaxialMaterial = new FedeasConcr3Material(tag, *data->_Vector);
}

FedeasHardeningMaterialWrapper::FedeasHardeningMaterialWrapper(int tag,
	double E, double sigmaY, double Hiso, double Hkin) {
	_UniaxialMaterial = new FedeasHardeningMaterial(tag, E, sigmaY, Hiso, Hkin);
}

FedeasHardeningMaterialWrapper::FedeasHardeningMaterialWrapper(int tag, VectorWrapper^ data) {
	_UniaxialMaterial = new FedeasHardeningMaterial(tag, *data->_Vector);
}

FedeasHyster1MaterialWrapper::FedeasHyster1MaterialWrapper(int tag,
	double mom1p, double rot1p, double mom2p, double rot2p,
	double mom1n, double rot1n, double mom2n, double rot2n,
	double pinchX, double pinchY, double damfc1, double damfc2) {
	_UniaxialMaterial = new FedeasHyster1Material(tag, mom1p, rot1p, mom2p, rot2p, mom1n, rot1n, mom2n, rot2n, pinchX, pinchY, damfc1, damfc2);
}

FedeasHyster1MaterialWrapper::FedeasHyster1MaterialWrapper(int tag, VectorWrapper^ data) {
	_UniaxialMaterial = new FedeasHyster1Material(tag, *data->_Vector);
}

FedeasHyster2MaterialWrapper::FedeasHyster2MaterialWrapper(int tag,
	double mom1p, double rot1p, double mom2p, double rot2p,
	double mom3p, double rot3p, double mom1n, double rot1n,
	double mom2n, double rot2n, double mom3n, double rot3n,
	double pinchX, double pinchY, double damfc1, double damfc2) {
	_UniaxialMaterial = new FedeasHyster2Material(tag, mom1p, rot1p, mom2p, rot2p, mom3p, rot3p, mom1n, rot1n, mom2n, rot2n, mom3n, rot3n, pinchX, pinchY, damfc1, damfc2);
}

FedeasHyster2MaterialWrapper::FedeasHyster2MaterialWrapper(int tag,
	double mom1p, double rot1p, double mom2p, double rot2p,
	double mom1n, double rot1n, double mom2n, double rot2n,
	double pinchX, double pinchY, double damfc1, double damfc2) {
	_UniaxialMaterial = new FedeasHyster2Material(tag, mom1p, rot1p, mom2p, rot2p, mom1n, rot1n, mom2n, rot2n, pinchX, pinchY, damfc1, damfc2);
}

FedeasHyster2MaterialWrapper::FedeasHyster2MaterialWrapper(int tag, VectorWrapper^ data) {
	_UniaxialMaterial = new FedeasHyster2Material(tag, *data->_Vector);
}


FedeasSteel1MaterialWrapper::FedeasSteel1MaterialWrapper(int tag,
	double fy, double E0, double b) {
	_UniaxialMaterial = new FedeasSteel1Material(tag, fy, E0, b);
}


FedeasSteel1MaterialWrapper::
FedeasSteel1MaterialWrapper(int tag,
	double fy, double E0, double b,
	double a1, double a2, double a3, double a4) {
	_UniaxialMaterial = new
		FedeasSteel1Material(tag, fy, E0, b, a1, a2, a3, a4);
}

FedeasSteel1MaterialWrapper::FedeasSteel1MaterialWrapper(int tag, VectorWrapper^ data) {
	_UniaxialMaterial = new FedeasSteel1Material(tag, *data->_Vector);
}

FedeasSteel2MaterialWrapper::FedeasSteel2MaterialWrapper(int tag,
	double fy, double E0, double b,
	double R0, double cR1, double cR2,
	double a1, double a2, double a3, double a4) {
	_UniaxialMaterial = new FedeasSteel2Material(tag, fy, E0, b, R0, cR1, cR2, a1, a2, a3, a4);
}

FedeasSteel2MaterialWrapper::FedeasSteel2MaterialWrapper(int tag,
	double fy, double E0, double b,
	double R0, double cR1, double cR2) {
	_UniaxialMaterial = new FedeasSteel2Material(tag, fy, E0, b, R0, cR1, cR2);
}

FedeasSteel2MaterialWrapper::FedeasSteel2MaterialWrapper(int tag,
	double fy, double E0, double b) {
	_UniaxialMaterial = new FedeasSteel2Material(tag, fy, E0, b);
}

FedeasSteel2MaterialWrapper::FedeasSteel2MaterialWrapper(int tag, VectorWrapper^ data) {
	_UniaxialMaterial = new FedeasSteel2Material(tag, *data->_Vector);
}

PlasticDamageMaterialWrapper::PlasticDamageMaterialWrapper(int tag, double E, double Ft, double Fc, double ft_max,
	double fcy, double fc_max, double kt_crit, double Relax) {
	_UniaxialMaterial = new PlasticDamageMaterial(tag, E, Ft, Fc, ft_max, fcy, fc_max, kt_crit, Relax);
}

PlasticDamageMaterialWrapper::PlasticDamageMaterialWrapper(int tag, VectorWrapper^ data) {
	_UniaxialMaterial = new PlasticDamageMaterial(tag, *data->_Vector);
}


LimitStateMaterialWrapper::LimitStateMaterialWrapper(int tag,
	double mom1p, double rot1p, double mom2p, double rot2p,
	double mom3p, double rot3p,
	double mom1n, double rot1n, double mom2n, double rot2n,
	double mom3n, double rot3n,
	double pinchX, double pinchY,
	double damfc1, double damfc2,
	double beta) {
	_UniaxialMaterial = new LimitStateMaterial(tag, mom1p, rot1p, mom2p, rot2p, mom3p, rot3p, mom1n, rot1n, mom2n, rot2n, mom3n, rot3n, pinchX, pinchY, damfc1, damfc2, beta);
}

LimitStateMaterialWrapper::LimitStateMaterialWrapper(int tag,
	double mom1p, double rot1p, double mom2p, double rot2p,
	double mom1n, double rot1n, double mom2n, double rot2n,
	double pinchX, double pinchY,
	double damfc1, double damfc2,
	double beta) {
	_UniaxialMaterial = new
		LimitStateMaterial(tag, mom1p, rot1p, mom2p, rot2p, mom1n, rot1n, mom2n, rot2n, pinchX, pinchY, damfc1, damfc2, beta);
}

LimitStateMaterialWrapper::LimitStateMaterialWrapper(int tag,
	double mom1p, double rot1p, double mom2p, double rot2p,
	double mom3p, double rot3p,
	double mom1n, double rot1n, double mom2n, double rot2n,
	double mom3n, double rot3n,
	double pinchX, double pinchY,
	double damfc1, double damfc2,
	double beta, LimitCurveWrapper^ theCurve,
	int curveType, int degrade) {
	_UniaxialMaterial = new
		LimitStateMaterial(tag, mom1p, rot1p, mom2p, rot2p, mom3p, rot3p, mom1n, rot1n, mom2n, rot2n, mom3n, rot3n, pinchX, pinchY, damfc1, damfc2, beta, *theCurve->_LimitCurve, curveType, degrade);
}

PinchingLimitStateMaterialWrapper::PinchingLimitStateMaterialWrapper(int matTag,
	int nodeT, int nodeB, int drftAx, double Kelas, int crvTyp, int crvTag,
	double YpinchUPN, double YpinchRPN, double XpinchRPN,
	double YpinchUNP, double YpinchRNP, double XpinchRNP,
	double dmgStrsLimE, double dmgDispMax,
	double dmgE1, double dmgE2, double dmgE3, double dmgE4, double dmgELim,
	double dmgU1, double dmgU2, double dmgU3, double dmgU4, double dmgULim,
	double dmgR1, double dmgR2, double dmgR3, double dmgR4, double dmgRLim, double dmgRCyc,
	double dmgS1, double dmgS2, double dmgS3, double dmgS4, double dmgSLim, double dmgSCyc,
	int eleTag, double b, double d, double h, double a, double st, double As, double Acc,
	double ld, double db, double rhot, double fc, double fy, double fyt,
	BaseDomainWrapper^ theDomain, NodeWrapper^ theNodeT, NodeWrapper^ theNodeB, LimitCurveWrapper^ theCurve, BaseElementWrapper^ theElement) {
	_UniaxialMaterial = new PinchingLimitStateMaterial(matTag, nodeT, nodeB, drftAx, Kelas, crvTyp, crvTag, YpinchUPN, YpinchRPN, XpinchRPN, YpinchUNP, YpinchRNP, XpinchRNP, dmgStrsLimE, dmgDispMax, dmgE1, dmgE2, dmgE3, dmgE4, dmgELim, dmgU1, dmgU2, dmgU3, dmgU4, dmgULim, dmgR1, dmgR2, dmgR3, dmgR4, dmgRLim, dmgRCyc, dmgS1, dmgS2, dmgS3, dmgS4, dmgSLim, dmgSCyc, eleTag, b, d, h, a, st, As, Acc, ld, db, rhot, fc, fy, fyt, theDomain->_Domain, theNodeT->_Node, theNodeB->_Node, *theCurve->_LimitCurve, theElement->_Element);
}

PySimple1Wrapper::PySimple1Wrapper(int tag, int classtag, int soilType, double pult, double y50,
	double drag, double dashpot) {
	_UniaxialMaterial = new PySimple1(tag, classtag, soilType, pult, y50, drag, dashpot);
}

PySimple2Wrapper::PySimple2Wrapper(int tag, int classtag, int soilType, double pult, double y50,
	double drag, double dashpot) {
	_UniaxialMaterial = new PySimple2(tag, classtag, soilType, pult, y50, drag, dashpot);
}

PySimple3Wrapper::PySimple3Wrapper(int tag, int classtag, double pult, double pyield, double kmax, double Hmod, double dashpot) {
	throw gcnew System::NotImplementedException();
	//_UniaxialMaterial = new PySimple3(tag, classtag, pult, pyield, kmax, Hmod, dashpot);
}

QzSimple1Wrapper::QzSimple1Wrapper(int tag, int qzType, double Qult, double z50, double suction,
	double dashpot) {
	_UniaxialMaterial = new QzSimple1(tag, qzType, Qult, z50, suction, dashpot);
}

QzSimple2Wrapper::QzSimple2Wrapper(int tag, int qzType, double Qult, double z50, double suction,
	double dashpot) {
	_UniaxialMaterial = new QzSimple2(tag, qzType, Qult, z50, suction, dashpot);
}

TzSimple1Wrapper::TzSimple1Wrapper(int tag, int classtag, int tzType, double tult, double z50, double dashpot) {
	_UniaxialMaterial = new TzSimple1(tag, classtag, tzType, tult, z50, dashpot);
}

TzSimple2Wrapper::TzSimple2Wrapper(int tag, int classtag, int tzType, double tult, double z50, double dashpot) {
	_UniaxialMaterial = new TzSimple2(tag, classtag, tzType, tult, z50, dashpot);
}

TzLiq1Wrapper::TzLiq1Wrapper(int tag, int classtag, int tzType, double tult, double z50,
	double dashpot, int solidElem1, int solidElem2, BaseDomainWrapper^ theDomain) {
	_UniaxialMaterial = new TzLiq1(tag, classtag, tzType, tult, z50, dashpot, solidElem1, solidElem2, theDomain->_Domain);
}

TzLiq1Wrapper::TzLiq1Wrapper(int tag, int classtag, int tzType, double tult, double z50,
	double dashpot, BaseDomainWrapper^ theDomain, TimeSeriesWrapper^ theSeries) {
	_UniaxialMaterial = new 		TzLiq1(tag, classtag, tzType, tult, z50, dashpot, theDomain->_Domain, theSeries->_TimeSeries);
}

PyLiq1Wrapper::PyLiq1Wrapper(int tag, int classtag, int soilType, double pult, double y50, double drag,
	double dashpot, double pRes, int solidElem1, int solidElem2, BaseDomainWrapper^ theDomain) {
	_UniaxialMaterial = new PyLiq1(tag, classtag, soilType, pult, y50, drag, dashpot, pRes, solidElem1, solidElem2, theDomain->_Domain);
}

PyLiq1Wrapper::
PyLiq1Wrapper(int tag, int classtag, int soilType, double pult, double y50, double drag,
	double dashpot, double pRes, BaseDomainWrapper^ theDomain, TimeSeriesWrapper^ theSeries) {
	_UniaxialMaterial = new
		PyLiq1(tag, classtag, soilType, pult, y50, drag, dashpot, pRes, theDomain->_Domain, theSeries->_TimeSeries);
}


BilinWrapper::BilinWrapper(int tag,
	double Ke0, double As, double AsNeg, double My_pos, double My_neg,
	double LamdaS, double LamdaD, double LamdaA, double LamdaK, double Cs,
	double Cd, double Ca, double Ck, double Thetap_pos, double Thetap_neg,
	double Thetapc_pos, double Thetapc_neg, double K, double KNeg, double Thetau_pos,
	double Thetau_neg, double PDPlus, double PDNeg, double nFactor) {
	_UniaxialMaterial = new Bilin(tag, Ke0, As, AsNeg, My_pos, My_neg, LamdaS, LamdaD, LamdaA, LamdaK, Cs, Cd, Ca, Ck, Thetap_pos, Thetap_neg, Thetapc_pos, Thetapc_neg, K, KNeg, Thetau_pos, Thetau_neg, PDPlus, PDNeg, nFactor);
}

BilinWrapper::BilinWrapper(int tag,
	double Ke0, double As, double AsNeg, double My_pos, double My_neg,
	double LamdaS, double LamdaD, double LamdaA, double LamdaK, double Cs,
	double Cd, double Ca, double Ck, double Thetap_pos, double Thetap_neg,
	double Thetapc_pos, double Thetapc_neg, double K, double KNeg, double Thetau_pos,
	double Thetau_neg, double PDPlus, double PDNeg) {
	_UniaxialMaterial = new
		Bilin(tag, Ke0, As, AsNeg, My_pos, My_neg, LamdaS, LamdaD, LamdaA, LamdaK, Cs, Cd, Ca, Ck, Thetap_pos, Thetap_neg, Thetapc_pos, Thetapc_neg, K, KNeg, Thetau_pos, Thetau_neg, PDPlus, PDNeg);
}

Bilin02Wrapper::Bilin02Wrapper(int tag,
	double Ke0, double As, double AsNeg, double My_pos, double My_neg,
	double LamdaS, double LamdaD, double LamdaA, double LamdaK, double Cs,
	double Cd, double Ca, double Ck, double Thetap_pos, double Thetap_neg,
	double Thetapc_pos, double Thetapc_neg, double K, double KNeg, double Thetau_pos,
	double Thetau_neg, double PDPlus, double PDNeg, double nFactor) {
	_UniaxialMaterial = new Bilin02(tag, Ke0, As, AsNeg, My_pos, My_neg, LamdaS, LamdaD, LamdaA, LamdaK, Cs, Cd, Ca, Ck, Thetap_pos, Thetap_neg, Thetapc_pos, Thetapc_neg, K, KNeg, Thetau_pos, Thetau_neg, PDPlus, PDNeg, nFactor);
}

Bilin02Wrapper::Bilin02Wrapper(int tag,
	double Ke0, double As, double AsNeg, double My_pos, double My_neg,
	double LamdaS, double LamdaD, double LamdaA, double LamdaK, double Cs,
	double Cd, double Ca, double Ck, double Thetap_pos, double Thetap_neg,
	double Thetapc_pos, double Thetapc_neg, double K, double KNeg, double Thetau_pos,
	double Thetau_neg, double PDPlus, double PDNeg) {
	_UniaxialMaterial = new
		Bilin02(tag, Ke0, As, AsNeg, My_pos, My_neg, LamdaS, LamdaD, LamdaA, LamdaK, Cs, Cd, Ca, Ck, Thetap_pos, Thetap_neg, Thetapc_pos, Thetapc_neg, K, KNeg, Thetau_pos, Thetau_neg, PDPlus, PDNeg);
}

BilinearWrapper::BilinearWrapper(int tag, VectorWrapper^ inputParam, DamageModelWrapper^ strength, DamageModelWrapper^ stiffness, DamageModelWrapper^ capping) {
	_UniaxialMaterial = new Bilinear(tag, *inputParam->_Vector, strength->_DamageModel, stiffness->_DamageModel, capping->_DamageModel);
}

CloughWrapper::CloughWrapper(int tag, VectorWrapper^ inputParam) {
	_UniaxialMaterial = new Clough(tag, *inputParam->_Vector);
}

CloughDamageWrapper::CloughDamageWrapper(int tag, VectorWrapper^ inputParam, DamageModelWrapper^ strength, DamageModelWrapper^ stiffness, DamageModelWrapper^  accelerated, DamageModelWrapper^ capping) {
	_UniaxialMaterial = new CloughDamage(tag, *inputParam->_Vector, strength->_DamageModel, stiffness->_DamageModel, accelerated->_DamageModel, capping->_DamageModel);
}

CloughHenryWrapper::CloughHenryWrapper(int tag, VectorWrapper^ inputParam) {
	_UniaxialMaterial = new CloughHenry(tag, *inputParam->_Vector);
}

PinchingWrapper::PinchingWrapper(int tag, VectorWrapper^ inputParam) {
	_UniaxialMaterial = new Pinching(tag, *inputParam->_Vector);
}

PinchingDamageWrapper::PinchingDamageWrapper(int tag, VectorWrapper^ inputParam, DamageModelWrapper^ strength, DamageModelWrapper^ stiffness, DamageModelWrapper^  accelerated, DamageModelWrapper^ capping) {
	_UniaxialMaterial = new PinchingDamage(tag, *inputParam->_Vector, strength->_DamageModel, stiffness->_DamageModel, accelerated->_DamageModel, capping->_DamageModel);
}

AxialSpWrapper::AxialSpWrapper(int tag, double sce, double fty, double fcy,
	double bte, double bty, double bcy, double fcr) {
	_UniaxialMaterial = new AxialSp(tag, sce, fty, fcy, bte, bty, bcy, fcr);
}

AxialSpHDWrapper::AxialSpHDWrapper(int tag, double sce, double fty, double fcy, double bte,
	double bty, double bth, double bcy, double fcr, double ath) {
	_UniaxialMaterial = new AxialSpHD(tag, sce, fty, fcy, bte, bty, bth, bcy, fcr, ath);
}

BarSlipMaterialWrapper::BarSlipMaterialWrapper(int tag,
	double fc, double fy, double Es, double fu,
	double Eh, double db, double ld, int nbars, double width, double depth,
	int bsflag, int type) {
	_UniaxialMaterial = new BarSlipMaterial(tag, fc, fy, Es, fu, Eh, db, ld, nbars, width, depth, bsflag, type);
}

BarSlipMaterialWrapper::BarSlipMaterialWrapper(int tag,
	double fc, double fy, double Es, double fu,
	double Eh, double db, double ld, int nbars, double width, double depth,
	int bsflag, int type, int damage, int unit) {
	_UniaxialMaterial = new
		BarSlipMaterial(tag, fc, fy, Es, fu, Eh, db, ld, nbars, width, depth, bsflag, type, damage, unit);
}


BilinearOilDamperWrapper::BilinearOilDamperWrapper(int tag, double K, double C, double Fr, double p, double LGap, double NM, double RelTol, double AbsTol, double MaxHalf) {
	_UniaxialMaterial = new BilinearOilDamper(tag, K, C, Fr, p, LGap, NM, RelTol, AbsTol, MaxHalf);
}

Bond_SP01Wrapper::Bond_SP01Wrapper(int tag, double fy, double sy, double fu, double su, double Kz, double R, double Cd, double db, double fc, double la) {
	_UniaxialMaterial = new Bond_SP01(tag, fy, sy, fu, su, Kz, R, Cd, db, fc, la);
}

Bond_SP01Wrapper::Bond_SP01Wrapper(int tag, double fy, double sy, double fu, double su, double Kz, double R) {
	_UniaxialMaterial = new
		Bond_SP01(tag, fy, sy, fu, su, Kz, R);
}

BoucWenMaterialWrapper::BoucWenMaterialWrapper(int tag,
	double alpha,
	double ko,
	double n,
	double gamma,
	double beta,
	double Ao,
	double deltaA,
	double deltaNu,
	double deltaEta,
	double tolerance,
	int maxNumIter) {
	_UniaxialMaterial = new BoucWenMaterial(tag, alpha, ko, n, gamma, beta, Ao, deltaA, deltaNu, deltaEta, tolerance, maxNumIter);
}

BoucWenOriginalWrapper::BoucWenOriginalWrapper(int tag,
	double Ei,
	double fy,
	double alphaL) {
	_UniaxialMaterial = new BoucWenOriginal(tag, Ei, fy, alphaL);
}


BoucWenOriginalWrapper::

BoucWenOriginalWrapper(int tag,
	double Ei,
	double fy,
	double alphaL,
	double alphaNL,
	double mu,
	double eta,
	double beta,
	double gamma,
	double tol,
	int maxIter) {
	_UniaxialMaterial = new

		BoucWenOriginal(tag, Ei, fy, alphaL, alphaNL, mu, eta, beta, gamma, tol, maxIter);
}

CableMaterialWrapper::CableMaterialWrapper(int tag, double Prestress, double E, double unitWeightEff, double L_Element) {
	_UniaxialMaterial = new CableMaterial(tag, Prestress, E, unitWeightEff, L_Element);
}

CastWrapper::CastWrapper(int tag, double nLegs, double bo,
	double h, double fy, double eo, double l, double b,
	double R0, double cR1, double cR2,
	double a1, double a2, double a3, double a4) {
	_UniaxialMaterial = new Cast(tag, nLegs, bo, h, fy, eo, l, b, R0, cR1, cR2, a1, a2, a3, a4);
}

CastWrapper::CastWrapper(int tag, double nLegs, double bo, double h,
	double fy, double eo, double l, double b,
	double R0, double cR1, double cR2) {
	_UniaxialMaterial = new
		Cast(tag, nLegs, bo, h, fy, eo, l, b, R0, cR1, cR2);
}

CastWrapper::CastWrapper(int tag, double nLegs, double bo, double h,
	double fy, double eo, double l, double b) {
	_UniaxialMaterial = new
		Cast(tag, nLegs, bo, h, fy, eo, l, b);
}

CFSSSWPWrapper::CFSSSWPWrapper(int tag,
	double hight, int width, double fuf, double fyf,
	double tf, double Af, double fus, double fys, double ts,
	double np, double ds, double Vs, double screw_Spacing, double A, double L) {
	_UniaxialMaterial = new CFSSSWP(tag, hight, width, fuf, fyf, tf, Af, fus, fys, ts, np, ds, Vs, screw_Spacing, A, L);
}

CFSWSWPWrapper::CFSWSWPWrapper(int tag,
	double hight, int width, double fuf,
	double tf, double Ife,
	double Ifi, double ts,
	double np, double ds, double Vs,
	double screw_Spacing, double nc, double type, double A, double L) {
	_UniaxialMaterial = new CFSWSWP(tag, hight, width, fuf, tf, Ife, Ifi, ts, np, ds, Vs, screw_Spacing, nc, type, A, L);
}

Concrete01Wrapper::Concrete01Wrapper(int tag, double fpc, double eco, double fpcu, double ecu)
{
	_UniaxialMaterial = new Concrete01(tag, fpc, eco, fpcu, ecu);
}

Concrete01WithSITCWrapper::Concrete01WithSITCWrapper(int tag, double fpc, double eco, double fpcu, double ecu, double endStrainSITC) {
	_UniaxialMaterial = new Concrete01WithSITC(tag, fpc, eco, fpcu, ecu, endStrainSITC);
}



Concrete02Wrapper::Concrete02Wrapper(int tag, double fpc, double eco, double fpcu, double ecu, double rat, double ft, double Ets)
{
	_UniaxialMaterial = new Concrete02(tag, fpc, eco, fpcu, ecu, rat, ft, Ets);
}

Concrete02ThermalWrapper::Concrete02ThermalWrapper(int tag, double _fc, double _epsc0, double _fcu,
	double _epscu, double _rat, double _ft, double _Ets) {
	_UniaxialMaterial = new Concrete02Thermal(tag, _fc, _epsc0, _fcu, _epscu, _rat, _ft, _Ets);
}

Concrete04Wrapper::Concrete04Wrapper(int tag, double fpc, double eco, double ecu, double Ec0, double fct, double etu) {
	_UniaxialMaterial = new Concrete04(tag, fpc, eco, ecu, Ec0, fct, etu);
}

Concrete04Wrapper::Concrete04Wrapper(int tag, double fpc, double eco, double ecu, double Ec0, double fct, double etu, double beta) {
	_UniaxialMaterial = new
		Concrete04(tag, fpc, eco, ecu, Ec0, fct, etu, beta);
}

Concrete04Wrapper::Concrete04Wrapper(int tag, double fpc, double eco, double ecu, double Ec0) {
	_UniaxialMaterial = new
		Concrete04(tag, fpc, eco, ecu, Ec0);
}


Concrete06Wrapper::Concrete06Wrapper(int tag, double fc, double eo, double r, double k, double alphaC, double fcr, double ecr, double b, double alphaT) {
	_UniaxialMaterial = new Concrete06(tag, fc, eo, r, k, alphaC, fcr, ecr, b, alphaT);
}

Concrete07Wrapper::Concrete07Wrapper(int tag, double FPC, double EPSC0, double EC, double FPT, double ESPT0, double XCRP, double XCRN, double R) {
	_UniaxialMaterial = new Concrete07(tag, FPC, EPSC0, EC, FPT, ESPT0, XCRP, XCRN, R);
}

ConcreteCMWrapper::ConcreteCMWrapper(int tag, double fpcc, double epcc, double Ec, double rc, double xcrn, double ft, double et, double rt, double xcrp) {
	_UniaxialMaterial = new ConcreteCM(tag, fpcc, epcc, Ec, rc, xcrn, ft, et, rt, xcrp);
}

ConcreteCMWrapper::ConcreteCMWrapper(int tag, double fpcc, double epcc, double Ec, double rc, double xcrn, double ft, double et, double rt, double xcrp, int mon) {
	_UniaxialMaterial = new
		ConcreteCM(tag, fpcc, epcc, Ec, rc, xcrn, ft, et, rt, xcrp, mon);
}

ConcreteCMWrapper::ConcreteCMWrapper(int tag, double fpcc, double epcc, double Ec, double rc, double xcrn, double ft, double et, double rt, double xcrp, int Gap, int dummy) {
	_UniaxialMaterial = new
		ConcreteCM(tag, fpcc, epcc, Ec, rc, xcrn, ft, et, rt, xcrp, Gap, dummy);
}

ConcreteDWrapper::ConcreteDWrapper(int tag, double fc0, double ec0, double ft0, double eptt0, double Ec0,
	double alphac0, double alphat0, double cesp0, double etap0) {
	_UniaxialMaterial = new ConcreteD(tag, fc0, ec0, ft0, eptt0, Ec0, alphac0, alphat0, cesp0, etap0);
}

ConcreteDWrapper::ConcreteDWrapper(int tag, double fc0, double ec0, double ft0, double eptt0, double Ec0,
	double alphac0, double alphat0) {
	_UniaxialMaterial = new ConcreteD(tag, fc0, ec0, ft0, eptt0, Ec0, alphac0, alphat0);
}

ConcreteECThermalWrapper::ConcreteECThermalWrapper(int tag, double _fc, double _epsc0, double _fcu,
	double _epscu, double _rat, double _ft, double _Ets) {
	_UniaxialMaterial = new ConcreteECThermal(tag, _fc, _epsc0, _fcu, _epscu, _rat, _ft, _Ets);
}

ConcreteSakaiKawashimaWrapper::ConcreteSakaiKawashimaWrapper(int tag, double YMx, double Sigcc, double EPScc) {
	_UniaxialMaterial = new ConcreteSakaiKawashima(tag, YMx, Sigcc, EPScc);
}

ConcretewBetaWrapper::ConcretewBetaWrapper(int tag, double fpc, double ec0, double fcint, double ecint, double fcres, double ecres, double fct, double ftint, double etint, double ftres, double etres, double lambda, double alpha, double bint, double etbint, double bres, double etbres, double M, double E0, double fcc, double ecc) {
	_UniaxialMaterial = new ConcretewBeta(tag, fpc, ec0, fcint, ecint, fcres, ecres, fct, ftint, etint, ftres, etres, lambda, alpha, bint, etbint, bres, etbres, M, E0, fcc, ecc);
}

ConfinedConcrete01Wrapper::ConfinedConcrete01Wrapper(int tag, VectorWrapper^ eps, VectorWrapper^ sigmac) {
	throw gcnew System::NotImplementedException();
	//_UniaxialMaterial = new ConfinedConcrete01(tag, eps->_Vector, sigmac->_Vector);
}

ConfinedConcrete01Wrapper::ConfinedConcrete01Wrapper(int tag, int secType, int dim, VectorWrapper^ semiLength,
	VectorWrapper^ phis, VectorWrapper^ S,
	VectorWrapper^ fyh, VectorWrapper^ Es0,
	VectorWrapper^ haRatio, VectorWrapper^ mueps,
	VectorWrapper^ As, VectorWrapper^ Is,
	double rhos, double fpc, double stRatio, double Ec,
	int epscuOption, double epscu, double epscuLimit,
	int nuOption, double nuc, double phiLon, int concrType,
	int aggrType, double tol, int maxNumIter) {
	throw gcnew System::NotImplementedException();

	//_UniaxialMaterial = new		ConfinedConcrete01(tag, secType, dim, *semiLength->_Vector, *phis->_Vector, *S->_Vector, *fyh->_Vector, *Es0->_Vector, *haRatio->_Vector, *mueps->_Vector, *As->_Vector, *Is->_Vector, rhos, fpc, stRatio, Ec, epscuOption, epscu, epscuLimit, nuOption, nuc, phiLon, concrType, aggrType, tol, maxNumIter);
}

DamperMaterialWrapper::DamperMaterialWrapper(int tag,
	UniaxialMaterialWrapper^ theMaterial) {
	_UniaxialMaterial = new DamperMaterial(tag, theMaterial->_UniaxialMaterial);
}

Dodd_RestrepoWrapper::Dodd_RestrepoWrapper(int tag,
	double Fy,
	double Fsu,
	double ESH,
	double ESU,
	double Youngs,
	double ESHI,
	double FSHI,
	double OmegaFac,
	double Conv) {
	_UniaxialMaterial = new Dodd_Restrepo(tag, Fy, Fsu, ESH, ESU, Youngs, ESHI, FSHI, OmegaFac, Conv);
}

ECC01Wrapper::ECC01Wrapper(int tag, double SIGT0, double EPST0, double SIGT1, double EPST1, double EPST2, double SIGC0,
	double EPSC0, double EPSC1, double ALPHAT1, double ALPHAT2, double ALPHAC, double ALPHACU, double BETAT, double BETAC) {
	_UniaxialMaterial = new ECC01(tag, SIGT0, EPST0, SIGT1, EPST1, EPST2, SIGC0, EPSC0, EPSC1, ALPHAT1, ALPHAT2, ALPHAC, ALPHACU, BETAT, BETAC);
}

Elastic2MaterialWrapper::Elastic2MaterialWrapper(int tag, double E, double eta) {
	_UniaxialMaterial = new Elastic2Material(tag, E, eta);
}

ElasticBilinWrapper::ElasticBilinWrapper(int tag, double E1, double E2, double eps2) {
	_UniaxialMaterial = new ElasticBilin(tag, E1, E2, eps2);
}

ElasticBilinWrapper::ElasticBilinWrapper(int tag, double E1P, double E2P, double epsP, double E1N, double E2N, double eps2N) {
	_UniaxialMaterial = new
		ElasticBilin(tag, E1P, E2P, epsP, E1N, E2N, eps2N);
}

ElasticMaterialWrapper::ElasticMaterialWrapper(int tag, double E, double eta)
{
	_UniaxialMaterial = new ElasticMaterial(tag, E, eta);
}

//ElasticMaterialThermalWrapper::ElasticMaterialThermalWrapper(int tag, double Epos, double alpha, double et, double Eneg, int softindex) {
//	_UniaxialMaterial = new ElasticMaterialThermal(tag, Epos, alpha, et, Eneg, softindex);
//}
//
//ElasticMaterialThermalWrapper::ElasticMaterialThermalWrapper() {
//	_UniaxialMaterial = new ElasticMaterialThermal();
//}

ElasticMaterialWrapper::ElasticMaterialWrapper(int tag, double Epos, double eta, double Eneg) {
	_UniaxialMaterial = new ElasticMaterial(tag, Epos, eta, Eneg);
}



ElasticMultiLinearWrapper::ElasticMultiLinearWrapper(int tag,
	VectorWrapper^ strainPoints,
	VectorWrapper^ stressPoints,
	double eta) {
	_UniaxialMaterial = new ElasticMultiLinear(tag, *strainPoints->_Vector, *stressPoints->_Vector, eta);
}

ElasticPPMaterialWrapper::ElasticPPMaterialWrapper(int tag, double E, double eyp) {
	_UniaxialMaterial = new ElasticPPMaterial(tag, E, eyp);
}

ElasticPPMaterialWrapper::ElasticPPMaterialWrapper(int tag, double E, double eyp, double eyn, double ezero) {
	_UniaxialMaterial = new
		ElasticPPMaterial(tag, E, eyp, eyn, ezero);
}

ENTMaterialWrapper::ENTMaterialWrapper(int tag, double E) {
	_UniaxialMaterial = new ENTMaterial(tag, E);
}

ENTMaterialWrapper::ENTMaterialWrapper(int tag, double E, double a, double b) {
	_UniaxialMaterial = new
		ENTMaterial(tag, E, a, b);
}

EPPGapMaterialWrapper::EPPGapMaterialWrapper(int tag, double E, double fy, double gap, double eta, int damage) {
	_UniaxialMaterial = new EPPGapMaterial(tag, E, fy, gap, eta, damage);
}

FatigueMaterialWrapper::FatigueMaterialWrapper(int tag, UniaxialMaterialWrapper^ material) {
	_UniaxialMaterial = new FatigueMaterial(tag, *material->_UniaxialMaterial);
}

FatigueMaterialWrapper::
FatigueMaterialWrapper(int tag, UniaxialMaterialWrapper^ material, double Dmax,
	double E0,
	double m,
	double minStrain,
	double maxStrain) {
	_UniaxialMaterial = new		FatigueMaterial(tag, *material->_UniaxialMaterial, Dmax, E0, m, minStrain, maxStrain);
}

FRPConfinedConcreteWrapper::FRPConfinedConcreteWrapper(int tag, double fpc1, double fpc2, double epsc0, double D, double c, double Ej, double Sj, double tj, double eju, double S, double fyl, double fyh, double dlong, double dtrans, double Es, double v0, double k, double useBuck) {
	_UniaxialMaterial = new FRPConfinedConcrete(tag, fpc1, fpc2, epsc0, D, c, Ej, Sj, tj, eju, S, fyl, fyh, dlong, dtrans, Es, v0, k, useBuck);
}

FRPConfinedConcrete02Wrapper::FRPConfinedConcrete02Wrapper(int tag, double fc0, double Ec, double ec0, double t, double Efrp, double eps_h_rup, double R, double ft, double Ets, int Unit) {
	_UniaxialMaterial = new FRPConfinedConcrete02(tag, fc0, Ec, ec0, t, Efrp, eps_h_rup, R, ft, Ets, Unit);
}

FRPConfinedConcrete02Wrapper::
FRPConfinedConcrete02Wrapper(int tag, double fc0, double Ec, double ec0, double fcc, double ecu, double ft, double Ets, int Unit) {
	_UniaxialMaterial = new		FRPConfinedConcrete02(tag, fc0, Ec, ec0, fcc, ecu, ft, Ets, Unit);
}

FRPConfinedConcrete02Wrapper::
FRPConfinedConcrete02Wrapper(int tag, double fc0, double Ec, double ec0, double ft, double Ets, int Unit) {
	_UniaxialMaterial = new		FRPConfinedConcrete02(tag, fc0, Ec, ec0, ft, Ets, Unit);
}

HardeningMaterialWrapper::HardeningMaterialWrapper(int tag, double E, double sigmaY,
	double K, double H, double eta) {
	_UniaxialMaterial = new HardeningMaterial(tag, E, sigmaY, K, H, eta);
}

HookGapWrapper::HookGapWrapper(int tag, double E, double gap) {
	_UniaxialMaterial = new HookGap(tag, E, gap);
}

HookGapWrapper::HookGapWrapper(int tag, double E, double gapN, double gapP) {
	_UniaxialMaterial = new		HookGap(tag, E, gapN, gapP);
}


HyperbolicGapMaterialWrapper::HyperbolicGapMaterialWrapper(int tag, double Kmax, double Kur, double Rf, double Fult, double gap) {
	_UniaxialMaterial = new HyperbolicGapMaterial(tag, Kmax, Kur, Rf, Fult, gap);
}

HystereticMaterialWrapper::HystereticMaterialWrapper(int tag,
	double mom1p, double rot1p, double mom2p, double rot2p,
	double mom3p, double rot3p,
	double mom1n, double rot1n, double mom2n, double rot2n,
	double mom3n, double rot3n,
	double pinchX, double pinchY,
	double damfc1, double damfc2,
	double beta) {
	_UniaxialMaterial = new HystereticMaterial(tag, mom1p, rot1p, mom2p, rot2p, mom3p, rot3p, mom1n, rot1n, mom2n, rot2n, mom3n, rot3n, pinchX, pinchY, damfc1, damfc2, beta);
}

HystereticMaterialWrapper::HystereticMaterialWrapper(int tag,
	double mom1p, double rot1p, double mom2p, double rot2p,
	double mom1n, double rot1n, double mom2n, double rot2n,
	double pinchX, double pinchY,
	double damfc1, double damfc2,
	double beta) {
	_UniaxialMaterial = new		HystereticMaterial(tag, mom1p, rot1p, mom2p, rot2p, mom1n, rot1n, mom2n, rot2n, pinchX, pinchY, damfc1, damfc2, beta);
}

ImpactMaterialWrapper::ImpactMaterialWrapper(int tag, double K1, double K2, double Delta_y, double gap) {
	_UniaxialMaterial = new ImpactMaterial(tag, K1, K2, Delta_y, gap);
}

InitStrainMaterialWrapper::InitStrainMaterialWrapper(int tag, UniaxialMaterialWrapper^ material, double epsInit) {
	_UniaxialMaterial = new InitStrainMaterial(tag, *material->_UniaxialMaterial, epsInit);
}

InitStressMaterialWrapper::InitStressMaterialWrapper(int tag, UniaxialMaterialWrapper^ material, double sigInit) {
	_UniaxialMaterial = new InitStressMaterial(tag, *material->_UniaxialMaterial, sigInit);
}

KikuchiAikenHDRWrapper::KikuchiAikenHDRWrapper(int tag, int tp, double ar, double hr,
	double cg, double ch, double cu, double rs, double rf) {
	_UniaxialMaterial = new KikuchiAikenHDR(tag, tp, ar, hr, cg, ch, cu, rs, rf);
}

KikuchiAikenLRBWrapper::KikuchiAikenLRBWrapper(int tag, int type, double ar, double hr, double gr, double ap, double tp,
	double alph, double beta, double temp, double rk, double rq, double rs, double rf) {
	_UniaxialMaterial = new KikuchiAikenLRB(tag, type, ar, hr, gr, ap, tp, alph, beta, temp, rk, rq, rs, rf);
}

MaxwellWrapper::MaxwellWrapper(int tag, double K, double C, double Alpha, double L, int returnD) {
	_UniaxialMaterial = new Maxwell(tag, K, C, Alpha, L, returnD);
}

MinMaxMaterialWrapper::MinMaxMaterialWrapper(int tag, UniaxialMaterialWrapper^ theMaterial, double min, double max) {
	_UniaxialMaterial = new MinMaxMaterial(tag, *theMaterial->_UniaxialMaterial, min, max);
}

ModIMKPeakOrientedWrapper::ModIMKPeakOrientedWrapper(int tag, double Ke0, double AlfanPos, double AlfanNeg, double My_pos, double My_neg,			// Updated: Filipe Ribeiro and Andre Barbosa
	double Ls, double Ld, double La, double Lk, double Cs, double Cd, double Ca, double Ck,
	double ThetaPpos, double ThetaPneg, double ThetaPCpos, double ThetaPCneg,
	double ResfacPos, double ResfacNeg, double FracDispPos, double FracDispNeg,
	double DPos, double DNeg, double nFactor) {
	_UniaxialMaterial = new ModIMKPeakOriented(tag, Ke0, AlfanPos, AlfanNeg, My_pos, My_neg, Ls, Ld, La, Lk, Cs, Cd, Ca, Ck, ThetaPpos, ThetaPneg, ThetaPCpos, ThetaPCneg, ResfacPos, ResfacNeg, FracDispPos, FracDispNeg, DPos, DNeg, nFactor);
}

ModIMKPeakOrientedWrapper::ModIMKPeakOrientedWrapper(int tag, double Ke0, double AlfanPos, double AlfanNeg, double My_pos, double My_neg,			// Updated: Filipe Ribeiro and Andre Barbosa 
	double Ls, double Ld, double La, double Lk, double Cs, double Cd, double Ca, double Ck,
	double ThetaPpos, double ThetaPneg, double ThetaPCpos, double ThetaPCneg,
	double ResfacPos, double ResfacNeg, double FracDispPos, double FracDispNeg,
	double DPos, double DNeg) {
	_UniaxialMaterial = new		ModIMKPeakOriented(tag, Ke0, AlfanPos, AlfanNeg, My_pos, My_neg, Ls, Ld, La, Lk, Cs, Cd, Ca, Ck, ThetaPpos, ThetaPneg, ThetaPCpos, ThetaPCneg, ResfacPos, ResfacNeg, FracDispPos, FracDispNeg, DPos, DNeg);
}

ModIMKPeakOriented02Wrapper::ModIMKPeakOriented02Wrapper(int tag, double Ke0, double AlfanPos, double AlfanNeg, double My_pos, double My_neg,
	double Ls, double Ld, double La, double Lk, double Cs, double Cd, double Ca, double Ck,
	double ThetaPpos, double ThetaPneg, double ThetaPCpos, double ThetaPCneg,
	double ResfacPos, double ResfacNeg, double FracDispPos, double FracDispNeg,
	double DPos, double DNeg, double C_Fp, double C_Fn, double nFactor) {
	_UniaxialMaterial = new ModIMKPeakOriented02(tag, Ke0, AlfanPos, AlfanNeg, My_pos, My_neg, Ls, Ld, La, Lk, Cs, Cd, Ca, Ck, ThetaPpos, ThetaPneg, ThetaPCpos, ThetaPCneg, ResfacPos, ResfacNeg, FracDispPos, FracDispNeg, DPos, DNeg, C_Fp, C_Fn, nFactor);
}

ModIMKPeakOriented02Wrapper::ModIMKPeakOriented02Wrapper(int tag, double Ke0, double AlfanPos, double AlfanNeg, double My_pos, double My_neg,
	double Ls, double Ld, double La, double Lk, double Cs, double Cd, double Ca, double Ck,
	double ThetaPpos, double ThetaPneg, double ThetaPCpos, double ThetaPCneg,
	double ResfacPos, double ResfacNeg, double FracDispPos, double FracDispNeg,
	double DPos, double DNeg, double C_Fp, double C_Fn) {
	_UniaxialMaterial = new		ModIMKPeakOriented02(tag, Ke0, AlfanPos, AlfanNeg, My_pos, My_neg, Ls, Ld, La, Lk, Cs, Cd, Ca, Ck, ThetaPpos, ThetaPneg, ThetaPCpos, ThetaPCneg, ResfacPos, ResfacNeg, FracDispPos, FracDispNeg, DPos, DNeg, C_Fp, C_Fn);
}

ModIMKPeakOriented02Wrapper::ModIMKPeakOriented02Wrapper(int tag, double Ke0, double AlfanPos, double AlfanNeg, double My_pos, double My_neg,
	double Ls, double Ld, double La, double Lk, double Cs, double Cd, double Ca, double Ck,
	double ThetaPpos, double ThetaPneg, double ThetaPCpos, double ThetaPCneg,
	double ResfacPos, double ResfacNeg, double FracDispPos, double FracDispNeg,
	double DPos, double DNeg) {
	_UniaxialMaterial = new		ModIMKPeakOriented02(tag, Ke0, AlfanPos, AlfanNeg, My_pos, My_neg, Ls, Ld, La, Lk, Cs, Cd, Ca, Ck, ThetaPpos, ThetaPneg, ThetaPCpos, ThetaPCneg, ResfacPos, ResfacNeg, FracDispPos, FracDispNeg, DPos, DNeg);
}

ModIMKPeakOriented02Wrapper::ModIMKPeakOriented02Wrapper(int tag, double Ke0, double AlfanPos, double AlfanNeg, double My_pos, double My_neg,
	double Ls, double Ld, double La, double Lk, double Cs, double Cd, double Ca, double Ck,
	double ThetaPpos, double ThetaPneg, double ThetaPCpos, double ThetaPCneg,
	double ResfacPos, double ResfacNeg, double FracDispPos, double FracDispNeg,
	double DPos, double DNeg, double nFactor) {
	_UniaxialMaterial = new		ModIMKPeakOriented02(tag, Ke0, AlfanPos, AlfanNeg, My_pos, My_neg, Ls, Ld, La, Lk, Cs, Cd, Ca, Ck, ThetaPpos, ThetaPneg, ThetaPCpos, ThetaPCneg, ResfacPos, ResfacNeg, FracDispPos, FracDispNeg, DPos, DNeg, nFactor);
}

ModIMKPinchingWrapper::ModIMKPinchingWrapper(int tag, double Ke0, double AlfanPos, double AlfanNeg, double My_pos, double My_neg, double FprPos, double FprNeg, double A_pinch,		// Updated: Filipe Ribeiro and Andre Barbosa
	double Ls, double Ld, double La, double Lk, double Cs, double Cd, double Ca, double Ck,
	double ThetaPpos, double ThetaPneg, double ThetaPCpos, double ThetaPCneg,
	double ResfacPos, double ResfacNeg, double FracDispPos, double FracDispNeg,
	double DPos, double DNeg, double nFactor) {
	_UniaxialMaterial = new ModIMKPinching(tag, Ke0, AlfanPos, AlfanNeg, My_pos, My_neg, FprPos, FprNeg, A_pinch, Ls, Ld, La, Lk, Cs, Cd, Ca, Ck, ThetaPpos, ThetaPneg, ThetaPCpos, ThetaPCneg, ResfacPos, ResfacNeg, FracDispPos, FracDispNeg, DPos, DNeg, nFactor);
}

ModIMKPinchingWrapper::ModIMKPinchingWrapper(int tag, double Ke0, double AlfanPos, double AlfanNeg, double My_pos, double My_neg, double FprPos, double FprNeg, double A_pinch,		// Updated: Filipe Ribeiro and Andre Barbosa
	double Ls, double Ld, double La, double Lk, double Cs, double Cd, double Ca, double Ck,
	double ThetaPpos, double ThetaPneg, double ThetaPCpos, double ThetaPCneg,
	double ResfacPos, double ResfacNeg, double FracDispPos, double FracDispNeg,
	double DPos, double DNeg) {
	_UniaxialMaterial = new		ModIMKPinching(tag, Ke0, AlfanPos, AlfanNeg, My_pos, My_neg, FprPos, FprNeg, A_pinch, Ls, Ld, La, Lk, Cs, Cd, Ca, Ck, ThetaPpos, ThetaPneg, ThetaPCpos, ThetaPCneg, ResfacPos, ResfacNeg, FracDispPos, FracDispNeg, DPos, DNeg);
}

ModIMKPinching02Wrapper::ModIMKPinching02Wrapper(int tag, double Ke0, double AlfanPos, double AlfanNeg, double My_pos, double My_neg, double FprPos, double FprNeg, double A_pinch,		// Updated: Filipe Ribeiro and Andre Barbosa
	double Ls, double Ld, double La, double Lk, double Cs, double Cd, double Ca, double Ck,
	double ThetaPpos, double ThetaPneg, double ThetaPCpos, double ThetaPCneg,
	double ResfacPos, double ResfacNeg, double FracDispPos, double FracDispNeg,
	double DPos, double DNeg, double nFactor) {
	_UniaxialMaterial = new ModIMKPinching02(tag, Ke0, AlfanPos, AlfanNeg, My_pos, My_neg, FprPos, FprNeg, A_pinch, Ls, Ld, La, Lk, Cs, Cd, Ca, Ck, ThetaPpos, ThetaPneg, ThetaPCpos, ThetaPCneg, ResfacPos, ResfacNeg, FracDispPos, FracDispNeg, DPos, DNeg, nFactor);
}

ModIMKPinching02Wrapper::ModIMKPinching02Wrapper(int tag, double Ke0, double AlfanPos, double AlfanNeg, double My_pos, double My_neg, double FprPos, double FprNeg, double A_pinch,		// Updated: Filipe Ribeiro and Andre Barbosa
	double Ls, double Ld, double La, double Lk, double Cs, double Cd, double Ca, double Ck,
	double ThetaPpos, double ThetaPneg, double ThetaPCpos, double ThetaPCneg,
	double ResfacPos, double ResfacNeg, double FracDispPos, double FracDispNeg,
	double DPos, double DNeg) {
	_UniaxialMaterial = new		ModIMKPinching02(tag, Ke0, AlfanPos, AlfanNeg, My_pos, My_neg, FprPos, FprNeg, A_pinch, Ls, Ld, La, Lk, Cs, Cd, Ca, Ck, ThetaPpos, ThetaPneg, ThetaPCpos, ThetaPCneg, ResfacPos, ResfacNeg, FracDispPos, FracDispNeg, DPos, DNeg);
}

MultiLinearWrapper::MultiLinearWrapper(int tag, VectorWrapper^ s, VectorWrapper^ e) {
	_UniaxialMaterial = new MultiLinear(tag, *s->_Vector, *e->_Vector);
}

OriginCenteredWrapper::OriginCenteredWrapper(int tag,
	double f1, double e1, double f2,
	double e2, double f3, double e3) {
	_UniaxialMaterial = new OriginCentered(tag, f1, e1, f2, e2, f3, e3);
}

ParallelMaterialWrapper::ParallelMaterialWrapper(int tag,
	int numMaterial,
	array<UniaxialMaterialWrapper^>^ theMaterials,
	VectorWrapper^ theFactors) {
	_UniaxialMaterial = new ParallelMaterial(tag, numMaterial, array2pointer5(theMaterials), theFactors->_Vector);
}

PathIndependentMaterialWrapper::PathIndependentMaterialWrapper(int tag,
	UniaxialMaterialWrapper^ theMaterial) {
	_UniaxialMaterial = new PathIndependentMaterial(tag, *theMaterial->_UniaxialMaterial);
}

Pinching4MaterialWrapper::Pinching4MaterialWrapper(int tag,
	double stress1p, double strain1p, double stress2p, double strain2p,
	double stress3p, double strain3p, double stress4p, double strain4p,
	double stress1n, double strain1n, double stress2n, double strain2n,
	double stress3n, double strain3n, double stress4n, double strain4n,
	double rDispP, double rForceP, double uForceP,
	double rDispN, double rForceN, double uForceN,
	double gammaK1, double gammaK2, double gammaK3,
	double gammaK4, double gammaKLimit,
	double gammaD1, double gammaD2, double gammaD3,
	double gammaD4, double gammaDLimit,
	double gammaF1, double gammaF2, double gammaF3,
	double gammaF4, double gammaFLimit, double gammaE, int DmgCyc) {
	_UniaxialMaterial = new Pinching4Material(tag, stress1p, strain1p, stress2p, strain2p, stress3p, strain3p, stress4p, strain4p, stress1n, strain1n, stress2n, strain2n, stress3n, strain3n, stress4n, strain4n, rDispP, rForceP, uForceP, rDispN, rForceN, uForceN, gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit, gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit, gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit, gammaE, DmgCyc);
}


Pinching4MaterialWrapper::

Pinching4MaterialWrapper(int tag,
	double stress1p, double strain1p, double stress2p, double strain2p,
	double stress3p, double strain3p, double stress4p, double strain4p,
	double rDispP, double rForceP, double uForceP,
	double gammaK1, double gammaK2, double gammaK3,
	double gammaK4, double gammaKLimit,
	double gammaD1, double gammaD2, double gammaD3,
	double gammaD4, double gammaDLimit,
	double gammaF1, double gammaF2, double gammaF3,
	double gammaF4, double gammaFLimit, double gammaE, int DmgCyc) {
	_UniaxialMaterial = new

		Pinching4Material(tag, stress1p, strain1p, stress2p, strain2p, stress3p, strain3p, stress4p, strain4p, rDispP, rForceP, uForceP, gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit, gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit, gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit, gammaE, DmgCyc);
}




pyUCLAWrapper::pyUCLAWrapper(int tag, int soilType, double pult, double y50, double Cd) {
	_UniaxialMaterial = new pyUCLA(tag, soilType, pult, y50, Cd);
}

RambergOsgoodSteelWrapper::RambergOsgoodSteelWrapper(int tag, double fy, double E0, double rezaA, double rezaN) {
	_UniaxialMaterial = new RambergOsgoodSteel(tag, fy, E0, rezaA, rezaN);
}

RambergOsgoodSteelWrapper::
RambergOsgoodSteelWrapper(int tag, double fy, double E0) {
	_UniaxialMaterial = new
		RambergOsgoodSteel(tag, fy, E0);
}

ReinforcingSteelWrapper::ReinforcingSteelWrapper(int tag, double fyield, double fultimate, double youngs, double youngs_hard,
	double estrainhard, double eultimate, int buckModel, double slenderness, double alpha, double r,
	double gama, double Fatigue1, double Fatigue2, double Degrade1,
	double RC1, double RC2, double RC3, double A1, double HardLim) {
	_UniaxialMaterial = new ReinforcingSteel(tag, fyield, fultimate, youngs, youngs_hard, estrainhard, eultimate, buckModel, slenderness, alpha, r, gama, Fatigue1, Fatigue2, Degrade1, RC1, RC2, RC3, A1, HardLim);
}

ResilienceLowWrapper::ResilienceLowWrapper(int tag, double PY, double DPmax, double Pmax, double Ke, double Kd) {
	_UniaxialMaterial = new ResilienceLow(tag, PY, DPmax, Pmax, Ke, Kd);
}

ResilienceMaterialHRWrapper::ResilienceMaterialHRWrapper(int tag, double  DY, double PY, double DPmax, double Pmax, double Ke, double Kd, double k) {
	_UniaxialMaterial = new ResilienceMaterialHR(tag, DY, PY, DPmax, Pmax, Ke, Kd, k);
}

SAWSMaterialWrapper::SAWSMaterialWrapper(int tag,
	double F0, double FI, double DU, double S0,
	double R1, double R2, double R3, double R4,
	double ALPHA, double BETA) {
	_UniaxialMaterial = new SAWSMaterial(tag, F0, FI, DU, S0, R1, R2, R3, R4, ALPHA, BETA);
}


SelfCenteringMaterialWrapper::SelfCenteringMaterialWrapper(int tag, double k1, double k2,
	double ActF, double beta, double SlipDef,
	double BearDef, double rBear) {
	_UniaxialMaterial = new SelfCenteringMaterial(tag, k1, k2, ActF, beta, SlipDef, BearDef, rBear);
}

SeriesMaterialWrapper::SeriesMaterialWrapper(int tag,
	int numMaterial,
	array<UniaxialMaterialWrapper^>^ theMaterials,
	int maxIter, double tol) {
	_UniaxialMaterial = new SeriesMaterial(tag, numMaterial, array2pointer5(theMaterials), maxIter, tol);
}

ShearPanelMaterialWrapper::ShearPanelMaterialWrapper(int tag,
	double stress1p, double strain1p, double stress2p, double strain2p,
	double stress3p, double strain3p, double stress4p, double strain4p,
	double stress1n, double strain1n, double stress2n, double strain2n,
	double stress3n, double strain3n, double stress4n, double strain4n,
	double rDispP, double rForceP, double uForceP,
	double rDispN, double rForceN, double uForceN,
	double gammaK1, double gammaK2, double gammaK3,
	double gammaK4, double gammaKLimit,
	double gammaD1, double gammaD2, double gammaD3,
	double gammaD4, double gammaDLimit,
	double gammaF1, double gammaF2, double gammaF3,
	double gammaF4, double gammaFLimit, double gammaE, double YldStress) {
	_UniaxialMaterial = new ShearPanelMaterial(tag, stress1p, strain1p, stress2p, strain2p, stress3p, strain3p, stress4p, strain4p, stress1n, strain1n, stress2n, strain2n, stress3n, strain3n, stress4n, strain4n, rDispP, rForceP, uForceP, rDispN, rForceN, uForceN, gammaK1, gammaK2, gammaK3, gammaK4, gammaKLimit, gammaD1, gammaD2, gammaD3, gammaD4, gammaDLimit, gammaF1, gammaF2, gammaF3, gammaF4, gammaFLimit, gammaE, YldStress);
}

SimpleFractureMaterialWrapper::SimpleFractureMaterialWrapper(int tag, UniaxialMaterialWrapper^ material, double max) {
	_UniaxialMaterial = new SimpleFractureMaterial(tag, *material->_UniaxialMaterial, max);
}

SmoothPSConcreteWrapper::SmoothPSConcreteWrapper(int tag, double fpc, double fpcu, double Ec, double eco, double ecu, double  eta) {
	_UniaxialMaterial = new SmoothPSConcrete(tag, fpc, fpcu, Ec, eco, ecu, eta);
}

StainlessECThermalWrapper::StainlessECThermalWrapper(int tag, int gradeTag, double fy, double E0, double fu, double sigInit) {
	_UniaxialMaterial = new StainlessECThermal(tag, gradeTag, fy, E0, fu, sigInit);
}

Steel01ThermalWrapper::Steel01ThermalWrapper(int tag, double fy, double e0, double b) {
	_UniaxialMaterial = new Steel01Thermal(tag, fy, e0, b);
}

Steel01ThermalWrapper::Steel01ThermalWrapper(int tag, double fy, double e0, double b,
	double a1, double a2,
	double a3, double a4) {
	_UniaxialMaterial = new
		Steel01Thermal(tag, fy, e0, b, a1, a2, a3, a4);
}

Steel02Wrapper::Steel02Wrapper(int tag,
	double fy, double E0, double b,
	double R0, double cR1, double cR2) {
	_UniaxialMaterial = new Steel02(tag, fy, E0, b, R0, cR1, cR2);
}

Steel02Wrapper::Steel02Wrapper(int tag,
	double fy, double E0, double b,
	double R0, double cR1, double cR2,
	double a1, double a2, double a3, double a4, double sigInit) {
	_UniaxialMaterial = new
		Steel02(tag, fy, E0, b, R0, cR1, cR2, a1, a2, a3, a4, sigInit);
}



Steel01Wrapper::Steel01Wrapper(int tag, double fy, double e0, double b) {
	_UniaxialMaterial = new Steel01(tag, fy, e0, b);
}

Steel01Wrapper::Steel01Wrapper(int tag, double fy, double e0, double b,
	double a1, double a2,
	double a3, double a4) {
	_UniaxialMaterial = new
		Steel01(tag, fy, e0, b, a1, a2, a3, a4);
}

Steel02ThermalWrapper::Steel02ThermalWrapper(int tag,
	double fy, double E0, double b,
	double R0, double cR1, double cR2) {
	_UniaxialMaterial = new Steel02Thermal(tag, fy, E0, b, R0, cR1, cR2);
}

Steel02ThermalWrapper::Steel02ThermalWrapper(int tag,
	double fy, double E0, double b,
	double R0, double cR1, double cR2,
	double a1, double a2, double a3, double a4, double sigInit) {
	_UniaxialMaterial = new
		Steel02Thermal(tag, fy, E0, b, R0, cR1, cR2, a1, a2, a3, a4, sigInit);
}

Steel2Wrapper::Steel2Wrapper(int tag,
	double fy, double E0, double b,
	double R0, double cR1, double cR2,
	double a1, double a2, double a3, double a4, double sigInit) {
	_UniaxialMaterial = new Steel2(tag, fy, E0, b, R0, cR1, cR2, a1, a2, a3, a4, sigInit);
}

Steel2Wrapper::Steel2Wrapper(int tag,
	double fy, double E0, double b,
	double R0, double cR1, double cR2) {
	_UniaxialMaterial = new
		Steel2(tag, fy, E0, b, R0, cR1, cR2);
}

Steel03Wrapper::Steel03Wrapper(int tag, double fy, double E0, double b, double r, double cR1, double cR2,
	double a1, double a2,
	double a3, double a4) {
	_UniaxialMaterial = new Steel03(tag, fy, E0, b, r, cR1, cR2, a1, a2, a3, a4);
}

Steel4Wrapper::Steel4Wrapper(int tag,
	//basics
	double f_y, double E_0,
	//kinematic hardening
	double b_k, double R_0, double r_1, double r_2,
	double b_kc, double R_0c, double r_1c, double r_2c,
	//isotropic hardening
	double b_i, double rho_i, double b_l, double R_i, double l_yp,
	double b_ic, double rho_ic, double b_lc, double R_ic,
	//ultimate strength limit
	double f_u, double R_u, double f_uc, double R_uc,
	//load history memory
	int cycNum,
	//initial stress
	double sig_init) {
	_UniaxialMaterial = new Steel4(tag, f_y, E_0, b_k, R_0, r_1, r_2, b_kc, R_0c, r_1c, r_2c, b_i, rho_i, b_l, R_i, l_yp, b_ic, rho_ic, b_lc, R_ic, f_u, R_u, f_uc, R_uc, cycNum, sig_init);
}

SteelBRBWrapper::SteelBRBWrapper(int tag,
	double E,
	double sigmaY0,
	double sigmaY_T,
	double alpha_T,
	double alpha_C,
	double sigmaY_C,
	double beta_T,
	double beta_C,
	double delta_T,
	double delta_C,
	double Tol) {
	_UniaxialMaterial = new SteelBRB(tag, E, sigmaY0, sigmaY_T, alpha_T, alpha_C, sigmaY_C, beta_T, beta_C, delta_T, delta_C, Tol);
}

SteelECThermalWrapper::SteelECThermalWrapper(int tag, int typeTag, double fy, double E0,
	double a1, double a2,
	double a3, double a4) {
	_UniaxialMaterial = new SteelECThermal(tag, typeTag, fy, E0, a1, a2, a3, a4);
}

SteelMPWrapper::SteelMPWrapper(int tag, double FY, double E, double B, double R, double CR1, double CR2, double A1, double A2) {
	_UniaxialMaterial = new SteelMP(tag, FY, E, B, R, CR1, CR2, A1, A2);
}

SteelMPFWrapper::SteelMPFWrapper(int tag, double sigyieldp, double sigyieldn,
	double E0, double bp, double bn, double R0, double a1, double a2) {
	_UniaxialMaterial = new SteelMPF(tag, sigyieldp, sigyieldn, E0, bp, bn, R0, a1, a2);
}

SteelMPFWrapper::
SteelMPFWrapper(int tag, double sigyieldp, double sigyieldn,
	double E0, double bp, double bn, double R0, double a1, double a2, double a3, double a4, double a5, double a6) {
	_UniaxialMaterial = new
		SteelMPF(tag, sigyieldp, sigyieldn, E0, bp, bn, R0, a1, a2, a3, a4, a5, a6);
}

UniaxialJ2PlasticityWrapper::UniaxialJ2PlasticityWrapper(int tag, double E, double sigmaY,
	double Hkin, double Hiso) {
	_UniaxialMaterial = new UniaxialJ2Plasticity(tag, E, sigmaY, Hkin, Hiso);
}

ViscousDamperWrapper::ViscousDamperWrapper(int tag, double K, double C, double Alpha, double LGap, double NM, double RelTol, double AbsTol, double MaxHalf) {
	_UniaxialMaterial = new ViscousDamper(tag, K, C, Alpha, LGap, NM, RelTol, AbsTol, MaxHalf);
}

ViscousMaterialWrapper::ViscousMaterialWrapper(int tag, double C, double Alpha, double minVel) {
	_UniaxialMaterial = new ViscousMaterial(tag, C, Alpha, minVel);
}

ExternalUniaxialMaterialWrapper::ExternalUniaxialMaterialWrapper(int tag) {
	_UniaxialMaterial = new ExternalUniaxialMaterial(tag);
}


