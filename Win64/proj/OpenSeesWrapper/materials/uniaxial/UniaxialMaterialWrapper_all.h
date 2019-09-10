#pragma once

#include <UniaxialMaterial.h>
#include "UniaxialMaterialWrapper.h"

// external uni material
#include <ExternalUniaxialMaterial.h>

// backbone
#include <backbone\ArctangentBackbone.h>
#include <BackboneMaterial.h>
#include <backbone\HystereticBackbone.h>
#include <backbone\ManderBackbone.h>
#include <backbone\RaynorBackbone.h>
#include <backbone\ReeseSandBackbone.h>
#include <backbone\ReeseSandBackbone.h>
#include <backbone\ReeseSoftClayBackbone.h>
#include <backbone\ReeseStiffClayBelowWS.h>
#include <backbone\TrilinearBackbone.h>

// drain
#include <drain\DrainBilinearMaterial.h>
#include <drain\DrainClough1Material.h>
#include <drain\DrainClough2Material.h>
#include <drain\DrainHardeningMaterial.h>
#include <drain\DrainPinch1Material.h>

// fedeas
#include <fedeas\FedeasBond1Material.h>
#include <fedeas\FedeasBond2Material.h>
#include <fedeas\FedeasConcr1Material.h>
#include <fedeas\FedeasConcr2Material.h>
#include <fedeas\FedeasConcr3Material.h>
#include <fedeas\FedeasHardeningMaterial.h>
#include <fedeas\FedeasHyster1Material.h>
#include <fedeas\FedeasHyster2Material.h>
#include <fedeas\FedeasSteel1Material.h>
#include <fedeas\FedeasSteel2Material.h>

// fedeas
#include <limitstate\LimitStateMaterial.h>
#include <limitstate\PinchingLimitStateMaterial.h>

//py_tz_qz
#include <py/PySimple1.h>
#include <py/PySimple2.h>
//#include <py/PySimple3.h>
#include <py/QzSimple1.h>
#include <py/QzSimple2.h>
#include <py/TzSimple1.h>
#include <py/TzSimple2.h>
#include <py/TzLiq1.h>
#include <py/PyLiq1.h>


//snap
#include <snap/Bilin.h>
#include <snap/Bilin02.h>
#include <snap/Bilinear.h>
#include <snap/Clough.h>
#include <snap/CloughDamage.h>
#include <snap/CloughHenry.h>
#include <snap/Pinching.h>
#include <snap/PinchingDamage.h>

// uniaxial
#include <AxialSp.h>
#include <AxialSpHD.h>
#include <BarSlipMaterial.h>
#include <BilinearOilDamper.h>
#include <Bond_SP01.h>
#include <BoucWenMaterial.h>
#include <BoucWenOriginal.h>
#include <BWBN.h>
#include <CableMaterial.h>
#include <Cast.h>
#include <CFSSSWP.h>
#include <CFSWSWP.h>
#include <Concrete01.h>
#include <Concrete01WithSITC.h>
#include <Concrete02.h>
#include <Concrete02Thermal.h>
#include <Concrete04.h>
#include <Concrete06.h>
#include <Concrete07.h>
#include <ConcreteCM.h>
#include <ConcreteD.h>
#include <ConcreteECThermal.h>
#include <ConcreteSakaiKawashima.h>
#include <ConcretewBeta.h>
#include <ConfinedConcrete01.h>

//#include <vector>

#include <CubicSpline.h>
#include <DamperMaterial.h>
#include <Dodd_Restrepo.h>
#include <DrainMaterial.h>
#include <ECC01.h>
#include <Elastic2Material.h>
#include <ElasticBilin.h>
#include <ElasticMaterial.h>
#include <ElasticMaterialThermal.h>
#include <ElasticMultiLinear.h>
#include <ElasticPPMaterial.h>
#include <ENTMaterial.h>
#include <EPPGapMaterial.h>
#include <FatigueMaterial.h>
#include <FedeasMaterial.h>
#include <FRPConfinedConcrete.h>
#include <FRPConfinedConcrete02.h>
#include <HardeningMaterial.h>
#include <HookGap.h>
#include <HyperbolicGapMaterial.h>
#include <HystereticMaterial.h>
#include <ImpactMaterial.h>
#include <InitStrainMaterial.h>
#include <InitStressMaterial.h>
#include <KikuchiAikenHDR.h>
#include <KikuchiAikenLRB.h>
#include <Maxwell.h>
#include <MinMaxMaterial.h>
#include <ModIMKPeakOriented.h>
#include <ModIMKPeakOriented02.h>
#include <ModIMKPinching.h>
#include <ModIMKPinching02.h>
#include <MultiLinear.h>
#include <NewUniaxialMaterial.h>
#include <OriginCentered.h>
#include <ParallelMaterial.h>
#include <PathIndependentMaterial.h>
#include <Pinching4Material.h>
#include <limitState/PinchingLimitStateMaterial.h>
#include <fedeas/PlasticDamageMaterial.h>
#include <pyUCLA.h>
#include <RambergOsgoodSteel.h>
#include <ReinforcingSteel.h>
#include <ResilienceLow.h>
#include <ResilienceMaterialHR.h>
#include <SAWSMaterial.h>
#include <SelfCenteringMaterial.h>
#include <SeriesMaterial.h>
#include <ShearPanelMaterial.h>
#include <SimpleFractureMaterial.h>
#include <SmoothPSConcrete.h>
#include <StainlessECThermal.h>
#include <Steel01.h>
#include <Steel01Thermal.h>
#include <Steel02.h>
#include <Steel02Thermal.h>
#include <Steel2.h>
#include <Steel03.h>
#include <Steel4.h>
#include <SteelBRB.h>
#include <SteelECThermal.h>
#include <SteelMP.h>
#include <SteelMPF.h>
#include <TriMatrix.h>
#include <UniaxialJ2Plasticity.h>
#include <ViscousDamper.h>
#include <ViscousMaterial.h>
#include <WrapperUniaxialMaterial.h>

#include <reinforcedConcretePlaneStress\ConcreteL01.h>
#include <reinforcedConcretePlaneStress\ConcreteZ01.h>
#include <reinforcedConcretePlaneStress\SteelZ01.h>
#include <reinforcedConcretePlaneStress\TendonL01.h>

#include <WrapperUniaxialMaterial.h>

#include "../MaterialWrapper.h"

#include "LimitCurveWrapper.h"
//#include "../../elements/ZeroLength/CoupledZeroLengthWrapper.h"
#include "../../damage/DamageModelWrapper.h"
#include "../../actors/IMovableObjectWrapper.h"
#include "../../taggeds/TaggedObjectWrapper.h"
#include "../../domains/domain/BaseDomainWrapper.h"
#include "../../domains/timeSeries/TimeSeriesWrapper.h"
#include "../../elements/BaseElementWrapper.h"
#include "../../domains/nodes/NodeWrapper.h"
#include "../../matrix/VectorWrapper.h"
//#include "../../handlers/HandlerWrapper.h"
//#include "../../recorders/ResponseWrapper.h"
//#include "../../recorders/InformationWrapper.h"


using namespace System;
using namespace System::Runtime::InteropServices;

using namespace OpenSees;
//using namespace OpenSees::Recorders;
//using namespace OpenSees::Handlers;
using namespace OpenSees::Components;
using namespace OpenSees::Components::Timeseries;
using namespace OpenSees::Elements;
using namespace OpenSees::DamageModels;
using namespace OpenSees::Materials;
using namespace OpenSees::Materials::Uniaxials;


namespace OpenSees {
	namespace Materials {
		namespace Uniaxials {

			//public ref class UniaxialMaterialWrapper : MaterialWrapper, IMovableObjectWrapper
			//{
			//public:
			//	UniaxialMaterialWrapper();
			//	UniaxialMaterialWrapper(int tag);
			//	~UniaxialMaterialWrapper() {
			//		if (_UniaxialMaterial != 0)
			//			delete _UniaxialMaterial;
			//	};

			//	UniaxialMaterialWrapper^ GetCopy() {
			//		UniaxialMaterial* mat = _UniaxialMaterial->getCopy();
			//		UniaxialMaterialWrapper^ theMat = gcnew UniaxialMaterialWrapper();
			//		theMat->_TaggedObject = mat;
			//		theMat->_UniaxialMaterial = mat;
			//		return theMat;
			//	}

			//	double GetStrain() { return _UniaxialMaterial->getStrain(); };
			//	double GetStrainRate() { return _UniaxialMaterial->getStrainRate(); };
			//	double GetStress() { return _UniaxialMaterial->getStress(); };
			//	double GetTangent() { return _UniaxialMaterial->getTangent(); };
			//	double GetInitialTangent() { return _UniaxialMaterial->getInitialTangent(); };
			//	double GetDampTangent() { return _UniaxialMaterial->getDampTangent(); };
			//	double GetRho() { return _UniaxialMaterial->getRho(); };
			//	int CommitState() {
			//		return _UniaxialMaterial->commitState();
			//	}
			//	int RevertToLastCommit() {
			//		return _UniaxialMaterial->revertToLastCommit();
			//	}
			//	int RevertToStart() {
			//		return _UniaxialMaterial->revertToStart();
			//	}
			//	String^ GetType() {
			//		return gcnew String(_UniaxialMaterial->getClassType());
			//	}
			//	int GetTypeTag() {
			//		return _UniaxialMaterial->getClassTag();
			//	}

			//	/*ResponseWrapper^ SetResponse(array<String^>^ argv, OpenSees::Handlers::OPS_StreamWrapper^ output) {
			//		char ** _argv = new char*[argv->Length];
			//		for (int i = 0; i < argv->Length; i++)
			//		{
			//			_argv[i] = OPS::StringToChar(argv[i]);
			//		}
			//		
			//		Response* response = this->_UniaxialMaterial->setResponse(_argv, argv->Length, *output->_OPS_StreamPtr);
			//		ResponseWrapper^ _response = gcnew ResponseWrapper();
			//		_response->_Response = response;
			//		return _response;
			//	}*/

			//	/*int GetResponse1(int responseID, InformationWrapper^ info) {
			//		return _UniaxialMaterial->getResponse(responseID, *info->_Information);
			//	}*/

			//internal:
			//	UniaxialMaterial * _UniaxialMaterial;
			//};


#pragma region Backbone
			public ref class HystereticBackboneWrapper abstract : TaggedObjectWrapper, IMovableObjectWrapper
			{
			public:
				HystereticBackboneWrapper() {};
				
				~HystereticBackboneWrapper() {
					if (_HystereticBackbone != 0)
						delete _HystereticBackbone;
				};

			internal:
				HystereticBackbone* _HystereticBackbone;
			};

			public ref class ArctangentBackboneWrapper : HystereticBackboneWrapper
			{
			public:
				ArctangentBackboneWrapper(int tag, double K1, double gammaY, double alpha);

				~ArctangentBackboneWrapper() {
					if (_HystereticBackbone != 0)
						delete _HystereticBackbone;
				};
			};


			public ref class ManderBackboneWrapper : HystereticBackboneWrapper
			{
			public:
				ManderBackboneWrapper(int tag, double fc, double epsc, double Ec);

				~ManderBackboneWrapper() {
					if (_HystereticBackbone != 0)
						delete _HystereticBackbone;
				};
			};

			public ref class RaynorBackboneWrapper : HystereticBackboneWrapper
			{
			public:
				RaynorBackboneWrapper(int tag, double es, double f1, double f2, double epsh, double epsm, double c1, double ey);

				~RaynorBackboneWrapper() {
					if (_HystereticBackbone != 0)
						delete _HystereticBackbone;
				};
			};

			public ref class ReeseSandBackboneWrapper : HystereticBackboneWrapper
			{
			public:
				ReeseSandBackboneWrapper(int tag, double kx, double ym, double pm,
					double yu, double pu);

				~ReeseSandBackboneWrapper() {
					if (_HystereticBackbone != 0)
						delete _HystereticBackbone;
				};
			};

			public ref class ReeseSoftClayBackboneWrapper : HystereticBackboneWrapper
			{
			public:
				ReeseSoftClayBackboneWrapper(int tag, double pu, double y50, double n);

				~ReeseSoftClayBackboneWrapper() {
					if (_HystereticBackbone != 0)
						delete _HystereticBackbone;
				};
			};

			public ref class ReeseStiffClayBelowWSWrapper : HystereticBackboneWrapper
			{
			public:
				ReeseStiffClayBelowWSWrapper(int tag, double esi, double y, double as, double pc);

				~ReeseStiffClayBelowWSWrapper() {
					if (_HystereticBackbone != 0)
						delete _HystereticBackbone;
				};
			};

			public ref class TrilinearBackboneWrapper : HystereticBackboneWrapper
			{
			public:
				TrilinearBackboneWrapper(int tag, double e1, double s1,
					double e2, double s2, double e3, double s3);
				TrilinearBackboneWrapper(int tag, double e1, double s1,
					double e2, double s2);

				~TrilinearBackboneWrapper() {
					if (_HystereticBackbone != 0)
						delete _HystereticBackbone;
				};
			};



			public ref class BackboneMaterialWrapper : OpenSees::Materials::Uniaxials::UniaxialMaterialWrapper
			{
			public:
				BackboneMaterialWrapper(int tag, HystereticBackboneWrapper^ backbone);

				~BackboneMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};
			};
#pragma endregion

#pragma region drain
			public ref class DrainMaterialWrapper abstract : UniaxialMaterialWrapper
			{
			public:
				DrainMaterialWrapper() {};
				DrainMaterialWrapper(int tag, int classTag, int numHV, int numData, double beto) {};
				~DrainMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};
			};

			public ref class DrainBilinearMaterialWrapper : DrainMaterialWrapper
			{
			public:
				DrainBilinearMaterialWrapper(int tag,
					double E, double fyp, double fyn, double alpha,
					double ecaps, double ecapk, double ecapa, double ecapd,
					double cs, double ck, double ca, double cd,
					double capSlope, double capDispP, double capDispN, double res, double beto);
				DrainBilinearMaterialWrapper(int tag, VectorWrapper^ input, double beto);
				~DrainBilinearMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};
			};

			public ref class DrainClough1MaterialWrapper : DrainMaterialWrapper
			{
			public:
				DrainClough1MaterialWrapper(int tag,
					double E, double fyp, double fyn, double alpha,
					double ecaps, double ecapk, double ecapa, double ecapd,
					double cs, double ck, double ca, double cd,
					double capSlope, double capDispP, double capDispN, double res, double beto);
				DrainClough1MaterialWrapper(int tag, VectorWrapper^ input, double beto);
				~DrainClough1MaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};
			};

			public ref class DrainClough2MaterialWrapper : DrainMaterialWrapper
			{
			public:
				DrainClough2MaterialWrapper(int tag,
					double E, double fyp, double fyn, double alpha,
					double ecaps, double ecapk, double ecapa, double ecapd,
					double cs, double ck, double ca, double cd,
					double capSlope, double capDispP, double capDispN, double res, double beto);
				DrainClough2MaterialWrapper(int tag, VectorWrapper^ input, double beto);
				~DrainClough2MaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};
			};

			public ref class DrainHardeningMaterialWrapper : DrainMaterialWrapper
			{
			public:
				DrainHardeningMaterialWrapper(int tag,
					double E, double sigY, double Hiso, double Hkin, double beto);
				~DrainHardeningMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};
			};

			public ref class DrainPinch1MaterialWrapper : DrainMaterialWrapper
			{
			public:
				DrainPinch1MaterialWrapper(int tag,
					double E, double fyp, double fyn, double alpha,
					double ecaps, double ecapk, double ecapa, double ecapd,
					double cs, double ck, double ca, double cd,
					double capSlope, double capDispP, double capDispN,
					double fpp, double fpn, double pinch, double res, double beto);
				~DrainPinch1MaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};
			};


#pragma endregion

#pragma region fedeas
			public ref class FedeasMaterialWrapper abstract : UniaxialMaterialWrapper
			{
			public:
				FedeasMaterialWrapper() {};
				FedeasMaterialWrapper(int tag, int classTag, int numHV, int numData) {};
				~FedeasMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};
			};

			public ref class FedeasBond1MaterialWrapper : FedeasMaterialWrapper
			{
			public:
				FedeasBond1MaterialWrapper(int tag,
					double u1p, double q1p, double u2p, double u3p, double q3p,
					double u1n, double q1n, double u2n, double u3n, double q3n,
					double s0, double bb);
				FedeasBond1MaterialWrapper(int tag, VectorWrapper^ input);
				~FedeasBond1MaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};
			};

			public ref class FedeasBond2MaterialWrapper : FedeasMaterialWrapper
			{
			public:
				FedeasBond2MaterialWrapper(int tag,
					double u1p, double q1p, double u2p, double u3p, double q3p,
					double u1n, double q1n, double u2n, double u3n, double q3n,
					double s0, double bb, double alp, double aln);
				FedeasBond2MaterialWrapper(int tag, VectorWrapper^ data);
				~FedeasBond2MaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};
			};

			public ref class FedeasConcr1MaterialWrapper : FedeasMaterialWrapper
			{
			public:
				FedeasConcr1MaterialWrapper(int tag,
					double fc, double ec, double fu, double eu);
				FedeasConcr1MaterialWrapper(int tag, VectorWrapper^ data);
				~FedeasConcr1MaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};
			};

			public ref class FedeasConcr2MaterialWrapper : FedeasMaterialWrapper
			{
			public:
				FedeasConcr2MaterialWrapper(int tag,
					double fc, double ec, double fu, double eu,
					double ratio, double ft, double Ets);
				FedeasConcr2MaterialWrapper(int tag, VectorWrapper^ data);
				~FedeasConcr2MaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};
			};

			public ref class FedeasConcr3MaterialWrapper : FedeasMaterialWrapper
			{
			public:
				FedeasConcr3MaterialWrapper(int tag, double fc, double ec, double fu, double eu,
					double rat, double ft, double epst0, double ft0, double beta, double epstu);
				FedeasConcr3MaterialWrapper(int tag, VectorWrapper^ data);
				~FedeasConcr3MaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};
			};

			public ref class FedeasHardeningMaterialWrapper : FedeasMaterialWrapper
			{
			public:
				FedeasHardeningMaterialWrapper(int tag,
					double E, double sigmaY, double Hiso, double Hkin);
				FedeasHardeningMaterialWrapper(int tag, VectorWrapper^ data);
				~FedeasHardeningMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};
			};

			public ref class FedeasHyster1MaterialWrapper : FedeasMaterialWrapper
			{
			public:
				FedeasHyster1MaterialWrapper(int tag,
					double mom1p, double rot1p, double mom2p, double rot2p,
					double mom1n, double rot1n, double mom2n, double rot2n,
					double pinchX, double pinchY, double damfc1 , double damfc2 );
				FedeasHyster1MaterialWrapper(int tag, VectorWrapper^ data);
				~FedeasHyster1MaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};
			};

			public ref class FedeasHyster2MaterialWrapper : FedeasMaterialWrapper
			{
			public:
				FedeasHyster2MaterialWrapper(int tag,
					double mom1p, double rot1p, double mom2p, double rot2p,
					double mom3p, double rot3p, double mom1n, double rot1n,
					double mom2n, double rot2n, double mom3n, double rot3n,
					double pinchX, double pinchY, double damfc1, double damfc2);
				FedeasHyster2MaterialWrapper(int tag,
					double mom1p, double rot1p, double mom2p, double rot2p,
					double mom1n, double rot1n, double mom2n, double rot2n,
					double pinchX, double pinchY, double damfc1, double damfc2);
				FedeasHyster2MaterialWrapper(int tag, VectorWrapper^ data);
				~FedeasHyster2MaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};
			};

			public ref class FedeasSteel1MaterialWrapper : FedeasMaterialWrapper
			{
			public:
				FedeasSteel1MaterialWrapper(int tag,
					double fy, double E0, double b,
					double a1, double a2, double a3, double a4);
				FedeasSteel1MaterialWrapper(int tag,
					double fy, double E0, double b);
				FedeasSteel1MaterialWrapper(int tag, VectorWrapper^ data);
				~FedeasSteel1MaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};
			};

			public ref class FedeasSteel2MaterialWrapper : FedeasMaterialWrapper
			{
			public:
				FedeasSteel2MaterialWrapper(int tag,
					double fy, double E0, double b,
					double R0, double cR1, double cR2,
					double a1, double a2, double a3, double a4);
				FedeasSteel2MaterialWrapper(int tag,
					double fy, double E0, double b,
					double R0, double cR1, double cR2);
				FedeasSteel2MaterialWrapper(int tag,
					double fy, double E0, double b);
				FedeasSteel2MaterialWrapper(int tag, VectorWrapper^ data);
				~FedeasSteel2MaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};
			};

			public ref class PlasticDamageMaterialWrapper : FedeasMaterialWrapper
			{
			public:
				PlasticDamageMaterialWrapper(int tag, double E, double Ft, double Fc, double ft_max,
					double fcy, double fc_max, double kt_crit, double Relax);
				PlasticDamageMaterialWrapper(int tag, VectorWrapper^ data);
				~PlasticDamageMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};
			};

#pragma endregion

#pragma region LimitState
			
			public ref class LimitStateMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				LimitStateMaterialWrapper(int tag,
					double mom1p, double rot1p, double mom2p, double rot2p,
					double mom3p, double rot3p,
					double mom1n, double rot1n, double mom2n, double rot2n,
					double mom3n, double rot3n,
					double pinchX, double pinchY,
					double damfc1, double damfc2,
					double beta);
				LimitStateMaterialWrapper(int tag,
					double mom1p, double rot1p, double mom2p, double rot2p,
					double mom1n, double rot1n, double mom2n, double rot2n,
					double pinchX, double pinchY,
					double damfc1, double damfc2,
					double beta);
				LimitStateMaterialWrapper(int tag,
					double mom1p, double rot1p, double mom2p, double rot2p,
					double mom3p, double rot3p,
					double mom1n, double rot1n, double mom2n, double rot2n,
					double mom3n, double rot3n,
					double pinchX, double pinchY,
					double damfc1, double damfc2,
					double beta, LimitCurveWrapper^ theCurve,
					int curveType, int degrade);
				~LimitStateMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class PinchingLimitStateMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				PinchingLimitStateMaterialWrapper(int matTag,
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
					BaseDomainWrapper^ theDomain, NodeWrapper^ theNodeT, NodeWrapper^ theNodeB, LimitCurveWrapper^ theCurve, BaseElementWrapper^ theElement);
				PinchingLimitStateMaterialWrapper() {};
				~PinchingLimitStateMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

#pragma endregion

#pragma region py-tz-qz
			public ref class PySimple1Wrapper : UniaxialMaterialWrapper
			{
			public:
				PySimple1Wrapper(int tag, int classtag, int soilType, double pult, double y50,
					double drag, double dashpot);
				PySimple1Wrapper() {};
				~PySimple1Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:
			};

			public ref class PySimple2Wrapper : UniaxialMaterialWrapper
			{
			public:
				PySimple2Wrapper(int tag, int classtag, int soilType, double pult, double y50,
					double drag, double dashpot);
				PySimple2Wrapper() {};
				~PySimple2Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class PySimple3Wrapper : UniaxialMaterialWrapper
			{
			public:
				PySimple3Wrapper(int tag, int classtag, double pult, double pyield, double kmax, double Hmod, double dashpot);
				PySimple3Wrapper() {};
				~PySimple3Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class QzSimple1Wrapper : UniaxialMaterialWrapper
			{
			public:
				QzSimple1Wrapper(int tag, int qzType, double Qult, double z50, double suction,
					double dashpot);
				QzSimple1Wrapper() {};
				~QzSimple1Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class QzSimple2Wrapper : UniaxialMaterialWrapper
			{
			public:
				QzSimple2Wrapper(int tag, int qzType, double Qult, double z50, double suction,
					double dashpot);
				QzSimple2Wrapper() {};
				~QzSimple2Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class TzSimple1Wrapper : UniaxialMaterialWrapper
			{
			public:
				TzSimple1Wrapper(int tag, int classtag, int tzType, double tult, double z50, double dashpot);
				TzSimple1Wrapper() {};
				~TzSimple1Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class TzSimple2Wrapper : UniaxialMaterialWrapper
			{
			public:
				TzSimple2Wrapper(int tag, int classtag, int tzType, double tult, double z50, double dashpot);
				TzSimple2Wrapper() {};
				~TzSimple2Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class TzLiq1Wrapper : TzSimple1Wrapper
			{
			public:
				TzLiq1Wrapper(int tag, int classtag, int tzType, double tult, double z50,
					double dashpot, int solidElem1, int solidElem2, BaseDomainWrapper^ theDomain);
				TzLiq1Wrapper(int tag, int classtag, int tzType, double tult, double z50,
					double dashpot, BaseDomainWrapper^ theDomain, TimeSeriesWrapper^ theSeries);
				TzLiq1Wrapper() {};
				~TzLiq1Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:
			};

			public ref class PyLiq1Wrapper : PySimple1Wrapper
			{
			public:
				PyLiq1Wrapper(int tag, int classtag, int soilType, double pult, double y50, double drag,
					double dashpot, double pRes, int solidElem1, int solidElem2, BaseDomainWrapper^ theDomain);
				PyLiq1Wrapper(int tag, int classtag, int soilType, double pult, double y50, double drag,
					double dashpot, double pRes, BaseDomainWrapper^ theDomain, TimeSeriesWrapper^ theSeries);
				PyLiq1Wrapper() {};
				~PyLiq1Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};


#pragma endregion

#pragma region snap
			public ref class BilinWrapper : UniaxialMaterialWrapper
			{
			public:
				BilinWrapper(int tag,
					double Ke0, double As, double AsNeg, double My_pos, double My_neg,
					double LamdaS, double LamdaD, double LamdaA, double LamdaK, double Cs,
					double Cd, double Ca, double Ck, double Thetap_pos, double Thetap_neg,
					double Thetapc_pos, double Thetapc_neg, double K, double KNeg, double Thetau_pos,
					double Thetau_neg, double PDPlus, double PDNeg, double nFactor);
				BilinWrapper(int tag,
					double Ke0, double As, double AsNeg, double My_pos, double My_neg,
					double LamdaS, double LamdaD, double LamdaA, double LamdaK, double Cs,
					double Cd, double Ca, double Ck, double Thetap_pos, double Thetap_neg,
					double Thetapc_pos, double Thetapc_neg, double K, double KNeg, double Thetau_pos,
					double Thetau_neg, double PDPlus, double PDNeg);
				BilinWrapper() {};
				~BilinWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class Bilin02Wrapper : UniaxialMaterialWrapper
			{
			public:
				Bilin02Wrapper(int tag,
					double Ke0, double As, double AsNeg, double My_pos, double My_neg,
					double LamdaS, double LamdaD, double LamdaA, double LamdaK, double Cs,
					double Cd, double Ca, double Ck, double Thetap_pos, double Thetap_neg,
					double Thetapc_pos, double Thetapc_neg, double K, double KNeg, double Thetau_pos,
					double Thetau_neg, double PDPlus, double PDNeg, double nFactor);
				Bilin02Wrapper(int tag,
					double Ke0, double As, double AsNeg, double My_pos, double My_neg,
					double LamdaS, double LamdaD, double LamdaA, double LamdaK, double Cs,
					double Cd, double Ca, double Ck, double Thetap_pos, double Thetap_neg,
					double Thetapc_pos, double Thetapc_neg, double K, double KNeg, double Thetau_pos,
					double Thetau_neg, double PDPlus, double PDNeg);
				Bilin02Wrapper() {};
				~Bilin02Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class BilinearWrapper : UniaxialMaterialWrapper
			{
			public:
				BilinearWrapper(int tag, VectorWrapper^ inputParam, DamageModelWrapper^ strength, DamageModelWrapper^ stiffness, DamageModelWrapper^ capping);
				BilinearWrapper() {};
				~BilinearWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class CloughWrapper : UniaxialMaterialWrapper
			{
			public:
				CloughWrapper(int tag, VectorWrapper^ inputParam);
				CloughWrapper() {};
				~CloughWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class CloughDamageWrapper : UniaxialMaterialWrapper
			{
			public:
				CloughDamageWrapper(int tag, VectorWrapper^ inputParam, DamageModelWrapper^ strength, DamageModelWrapper^ stiffness, DamageModelWrapper^  accelerated, DamageModelWrapper^ capping);
				CloughDamageWrapper() {};
				~CloughDamageWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class CloughHenryWrapper : UniaxialMaterialWrapper
			{
			public:
				CloughHenryWrapper(int tag, VectorWrapper^ inputParam);
				CloughHenryWrapper() {};
				~CloughHenryWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class PinchingWrapper : UniaxialMaterialWrapper
			{
			public:
				PinchingWrapper(int tag, VectorWrapper^ inputParam);
				PinchingWrapper() {};
				~PinchingWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class PinchingDamageWrapper : UniaxialMaterialWrapper
			{
			public:
				PinchingDamageWrapper(int tag, VectorWrapper^ inputParam, DamageModelWrapper^ strength, DamageModelWrapper^ stiffness, DamageModelWrapper^  accelerated, DamageModelWrapper^ capping);
				PinchingDamageWrapper() {};
				~PinchingDamageWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

#pragma endregion

			public ref class AxialSpWrapper : UniaxialMaterialWrapper
			{
			public:
				AxialSpWrapper(int tag, double sce, double fty, double fcy,
					double bte, double bty, double bcy, double fcr);
				AxialSpWrapper() {};
				~AxialSpWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class AxialSpHDWrapper : UniaxialMaterialWrapper
			{
			public:
				AxialSpHDWrapper(int tag, double sce, double fty, double fcy, double bte,
					double bty, double bth, double bcy, double fcr, double ath);

				AxialSpHDWrapper() {};
				~AxialSpHDWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class BarSlipMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				BarSlipMaterialWrapper(int tag,
					double fc, double fy, double Es, double fu,
					double Eh, double db, double ld, int nbars, double width, double depth,
					int bsflag, int type);
				BarSlipMaterialWrapper(int tag,
					double fc, double fy, double Es, double fu,
					double Eh, double db, double ld, int nbars, double width, double depth,
					int bsflag, int type, int damage, int unit);
				BarSlipMaterialWrapper() {};
				~BarSlipMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};
			
			public ref class BilinearOilDamperWrapper : UniaxialMaterialWrapper
			{
			public:
				BilinearOilDamperWrapper(int tag, double K, double C, double Fr, double p, double LGap, double NM, double RelTol, double AbsTol, double MaxHalf);
				BilinearOilDamperWrapper() {};
				~BilinearOilDamperWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class Bond_SP01Wrapper : UniaxialMaterialWrapper
			{
			public:
				Bond_SP01Wrapper(int tag, double fy, double sy, double fu, double su, double Kz, double R, double Cd, double db, double fc, double la);
				Bond_SP01Wrapper(int tag, double fy, double sy, double fu, double su, double Kz, double R);
				Bond_SP01Wrapper() {};
				~Bond_SP01Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class BoucWenMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				BoucWenMaterialWrapper(int tag,
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
					int maxNumIter);
				BoucWenMaterialWrapper() {};
				~BoucWenMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class BoucWenOriginalWrapper : UniaxialMaterialWrapper
			{
			public:
				BoucWenOriginalWrapper(int tag,
					double Ei,
					double fy,
					double alphaL);

				BoucWenOriginalWrapper(int tag,
					double Ei,
					double fy,
					double alphaL,
					double alphaNL ,
					double mu ,
					double eta ,
					double beta ,
					double gamma ,
					double tol ,
					int maxIter );
				BoucWenOriginalWrapper() {};
				~BoucWenOriginalWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class CableMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				CableMaterialWrapper(int tag, double Prestress, double E, double unitWeightEff, double L_Element);

				CableMaterialWrapper() {};
				~CableMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class CastWrapper : UniaxialMaterialWrapper
			{
			public:
				CastWrapper(int tag, double nLegs, double bo,
					double h, double fy, double eo, double l, double b,
					double R0, double cR1, double cR2,
					double a1, double a2, double a3, double a4);
				CastWrapper(int tag, double nLegs, double bo, double h,
					double fy, double eo, double l, double b,
					double R0, double cR1, double cR2);
				CastWrapper(int tag, double nLegs, double bo, double h,
					double fy, double eo, double l, double b);
				CastWrapper() {};
				~CastWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class CFSSSWPWrapper : UniaxialMaterialWrapper
			{
			public:
				CFSSSWPWrapper(int tag,
					double hight, int width, double fuf, double fyf,
					double tf, double Af, double fus, double fys, double ts,
					double np, double ds, double Vs, double screw_Spacing, double A, double L);
				CFSSSWPWrapper() {};
				~CFSSSWPWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class CFSWSWPWrapper : UniaxialMaterialWrapper
			{
			public:
				CFSWSWPWrapper(int tag,
					double hight, int width, double fuf,
					double tf, double Ife,
					double Ifi, double ts,
					double np, double ds, double Vs,
					double screw_Spacing, double nc, double type, double A, double L);
				CFSWSWPWrapper() {};
				~CFSWSWPWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class Concrete01Wrapper : UniaxialMaterialWrapper
			{
			public:
				Concrete01Wrapper(int tag, double fpc, double eco, double fpcu, double ecu);
				Concrete01Wrapper() {};
				~Concrete01Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class Concrete01WithSITCWrapper : UniaxialMaterialWrapper
			{
			public:
				Concrete01WithSITCWrapper(int tag, double fpc, double eco, double fpcu, double ecu, double endStrainSITC);
				Concrete01WithSITCWrapper() {};
				~Concrete01WithSITCWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class Concrete02Wrapper : UniaxialMaterialWrapper
			{
			public:
				Concrete02Wrapper(int tag, double fpc, double eco, double fpcu, double ecu, double rat, double ft, double Ets);
				Concrete02Wrapper() {};
				~Concrete02Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class Concrete02ThermalWrapper : UniaxialMaterialWrapper
			{
			public:
				Concrete02ThermalWrapper(int tag, double _fc, double _epsc0, double _fcu,
					double _epscu, double _rat, double _ft, double _Ets);
				Concrete02ThermalWrapper() {};
				~Concrete02ThermalWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class Concrete04Wrapper : UniaxialMaterialWrapper
			{
			public:
				Concrete04Wrapper(int tag, double fpc, double eco, double ecu, double Ec0, double fct, double etu);
				Concrete04Wrapper(int tag, double fpc, double eco, double ecu, double Ec0, double fct, double etu, double beta);
				Concrete04Wrapper(int tag, double fpc, double eco, double ecu, double Ec0);
				Concrete04Wrapper() {};
				~Concrete04Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class Concrete06Wrapper : UniaxialMaterialWrapper
			{
			public:
				Concrete06Wrapper(int tag, double fc, double eo, double r, double k, double alphaC, double fcr, double ecr, double b, double alphaT);
				Concrete06Wrapper() {};
				~Concrete06Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class Concrete07Wrapper : UniaxialMaterialWrapper
			{
			public:
				Concrete07Wrapper(int tag, double FPC, double EPSC0, double EC, double FPT, double ESPT0, double XCRP, double XCRN, double R);
				Concrete07Wrapper() {};
				~Concrete07Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class ConcreteCMWrapper : UniaxialMaterialWrapper
			{
			public:
				ConcreteCMWrapper(int tag, double fpcc, double epcc, double Ec, double rc, double xcrn, double ft, double et, double rt, double xcrp);
				ConcreteCMWrapper(int tag, double fpcc, double epcc, double Ec, double rc, double xcrn, double ft, double et, double rt, double xcrp, int mon);
				ConcreteCMWrapper(int tag, double fpcc, double epcc, double Ec, double rc, double xcrn, double ft, double et, double rt, double xcrp, int Gap, int dummy);
				ConcreteCMWrapper() {};
				~ConcreteCMWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class ConcreteDWrapper : UniaxialMaterialWrapper
			{
			public:
				ConcreteDWrapper(int tag, double fc0, double ec0, double ft0, double eptt0, double Ec0,
					double alphac0, double alphat0, double cesp0, double etap0);
				ConcreteDWrapper(int tag, double fc0, double ec0, double ft0, double eptt0, double Ec0,
					double alphac0, double alphat0);
				ConcreteDWrapper() {};
				~ConcreteDWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class ConcreteECThermalWrapper : UniaxialMaterialWrapper
			{
			public:
				ConcreteECThermalWrapper(int tag, double _fc, double _epsc0, double _fcu,
					double _epscu, double _rat, double _ft, double _Ets);
				ConcreteECThermalWrapper() {};
				~ConcreteECThermalWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class ConcreteSakaiKawashimaWrapper : UniaxialMaterialWrapper
			{
			public:
				ConcreteSakaiKawashimaWrapper(int tag, double YMx, double Sigcc, double EPScc);
				ConcreteSakaiKawashimaWrapper() {};
				~ConcreteSakaiKawashimaWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class ConcretewBetaWrapper : UniaxialMaterialWrapper
			{
			public:
				ConcretewBetaWrapper(int tag, double fpc, double ec0, double fcint, double ecint, double fcres, double ecres, double fct, double ftint, double etint, double ftres, double etres, double lambda, double alpha, double bint, double etbint, double bres, double etbres, double M, double E0, double fcc, double ecc);
				ConcretewBetaWrapper() {};
				~ConcretewBetaWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class ConfinedConcrete01Wrapper : UniaxialMaterialWrapper
			{
			public:
				ConfinedConcrete01Wrapper(int tag, VectorWrapper^ eps, VectorWrapper^ sigmac);
				ConfinedConcrete01Wrapper(int tag, int secType, int dim, VectorWrapper^ semiLength,
					VectorWrapper^ phis, VectorWrapper^ S,
					VectorWrapper^ fyh, VectorWrapper^ Es0,
					VectorWrapper^ haRatio, VectorWrapper^ mueps,
					VectorWrapper^ As, VectorWrapper^ Is,
					double rhos, double fpc, double stRatio, double Ec,
					int epscuOption, double epscu, double epscuLimit,
					int nuOption, double nuc, double phiLon, int concrType,
					int aggrType, double tol, int maxNumIter);
				ConfinedConcrete01Wrapper() {};
				~ConfinedConcrete01Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class DamperMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				DamperMaterialWrapper(int tag,
					UniaxialMaterialWrapper^ theMaterial);
				DamperMaterialWrapper() {};
				~DamperMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class Dodd_RestrepoWrapper : UniaxialMaterialWrapper
			{
			public:
				Dodd_RestrepoWrapper(int tag,
					double Fy,
					double Fsu,
					double ESH,
					double ESU,
					double Youngs,
					double ESHI,
					double FSHI,
					double OmegaFac,
					double Conv);
				Dodd_RestrepoWrapper() {};
				~Dodd_RestrepoWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class ECC01Wrapper : UniaxialMaterialWrapper
			{
			public:
				ECC01Wrapper(int tag, double SIGT0, double EPST0, double SIGT1, double EPST1, double EPST2, double SIGC0,
					double EPSC0, double EPSC1, double ALPHAT1, double ALPHAT2, double ALPHAC, double ALPHACU, double BETAT, double BETAC);
				ECC01Wrapper() {};
				~ECC01Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class Elastic2MaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				Elastic2MaterialWrapper(int tag, double E, double eta);
				Elastic2MaterialWrapper() {};
				~Elastic2MaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class ElasticBilinWrapper : UniaxialMaterialWrapper
			{
			public:
				ElasticBilinWrapper(int tag, double E1, double E2, double eps2);
				ElasticBilinWrapper(int tag, double E1P, double E2P, double epsP, double E1N, double E2N, double eps2N);
				ElasticBilinWrapper() {};
				~ElasticBilinWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class ElasticMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				ElasticMaterialWrapper(int tag, double E, double eta);
				ElasticMaterialWrapper(int tag, double Epos, double eta, double Eneg);
				~ElasticMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};
			private:

			};

			/*public ref class ElasticMaterialThermalWrapper : UniaxialMaterialWrapper
			{
			public:
				ElasticMaterialThermalWrapper(int tag, double Epos, double alpha, double et, double Eneg, int softindex);
				ElasticMaterialThermalWrapper();
				~ElasticMaterialThermalWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};*/

			public ref class ElasticMultiLinearWrapper : UniaxialMaterialWrapper
			{
			public:
				ElasticMultiLinearWrapper(int tag,
					VectorWrapper^ strainPoints,
					VectorWrapper^ stressPoints,
					double eta );
				ElasticMultiLinearWrapper() {};
				~ElasticMultiLinearWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class ElasticPPMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				ElasticPPMaterialWrapper(int tag, double E, double eyp);
				ElasticPPMaterialWrapper(int tag, double E, double eyp, double eyn, double ezero);
				ElasticPPMaterialWrapper() {};
				~ElasticPPMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class ENTMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				ENTMaterialWrapper(int tag, double E);
				ENTMaterialWrapper(int tag, double E, double a, double b);
				ENTMaterialWrapper() {};
				~ENTMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class EPPGapMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				EPPGapMaterialWrapper(int tag, double E, double fy, double gap, double eta, int damage);
				EPPGapMaterialWrapper() {};
				~EPPGapMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class FatigueMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				FatigueMaterialWrapper(int tag, UniaxialMaterialWrapper^ material);
				FatigueMaterialWrapper(int tag, UniaxialMaterialWrapper^ material, double Dmax,
					double E0,
					double m,
					double minStrain,
					double maxStrain);
				FatigueMaterialWrapper() {};
				~FatigueMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class FRPConfinedConcreteWrapper : UniaxialMaterialWrapper
			{
			public:
				FRPConfinedConcreteWrapper(int tag, double fpc1, double fpc2, double epsc0, double D, double c, double Ej, double Sj, double tj, double eju, double S, double fyl, double fyh, double dlong, double dtrans, double Es, double v0, double k, double useBuck);
				FRPConfinedConcreteWrapper() {};
				~FRPConfinedConcreteWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class FRPConfinedConcrete02Wrapper : UniaxialMaterialWrapper
			{
			public:
				FRPConfinedConcrete02Wrapper(int tag, double fc0, double Ec, double ec0, double t, double Efrp, double eps_h_rup, double R, double ft, double Ets, int Unit);
				FRPConfinedConcrete02Wrapper(int tag, double fc0, double Ec, double ec0, double fcc, double ecu, double ft, double Ets, int Unit);
				FRPConfinedConcrete02Wrapper(int tag, double fc0, double Ec, double ec0, double ft, double Ets, int Unit);
				FRPConfinedConcrete02Wrapper() {};
				~FRPConfinedConcrete02Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class HardeningMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				HardeningMaterialWrapper(int tag, double E, double sigmaY,
					double K, double H, double eta);
				HardeningMaterialWrapper() {};
				~HardeningMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class HookGapWrapper : UniaxialMaterialWrapper
			{
			public:
				HookGapWrapper(int tag, double E, double gap);
				HookGapWrapper(int tag, double E, double gapN, double gapP);
				HookGapWrapper() {};
				~HookGapWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class HyperbolicGapMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				HyperbolicGapMaterialWrapper(int tag, double Kmax, double Kur, double Rf, double Fult, double gap);
				HyperbolicGapMaterialWrapper() {};
				~HyperbolicGapMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class HystereticMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				HystereticMaterialWrapper(int tag,
					double mom1p, double rot1p, double mom2p, double rot2p,
					double mom3p, double rot3p,
					double mom1n, double rot1n, double mom2n, double rot2n,
					double mom3n, double rot3n,
					double pinchX, double pinchY,
					double damfc1 , double damfc2 ,
					double beta );
				HystereticMaterialWrapper(int tag,
					double mom1p, double rot1p, double mom2p, double rot2p,
					double mom1n, double rot1n, double mom2n, double rot2n,
					double pinchX, double pinchY,
					double damfc1, double damfc2,
					double beta);
				HystereticMaterialWrapper() {};
				~HystereticMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class ImpactMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				ImpactMaterialWrapper(int tag, double K1, double K2, double Delta_y, double gap);
				ImpactMaterialWrapper() {};
				~ImpactMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class InitStrainMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				InitStrainMaterialWrapper(int tag, UniaxialMaterialWrapper^ material, double epsInit);
				InitStrainMaterialWrapper() {};
				~InitStrainMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class InitStressMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				InitStressMaterialWrapper(int tag, UniaxialMaterialWrapper^ material, double sigInit);
				InitStressMaterialWrapper() {};
				~InitStressMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class KikuchiAikenHDRWrapper : UniaxialMaterialWrapper
			{
			public:
				KikuchiAikenHDRWrapper(int tag, int tp, double ar, double hr,
					double cg, double ch, double cu, double rs, double rf);
				KikuchiAikenHDRWrapper() {};
				~KikuchiAikenHDRWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class KikuchiAikenLRBWrapper : UniaxialMaterialWrapper
			{
			public:
				KikuchiAikenLRBWrapper(int tag, int type, double ar, double hr, double gr, double ap, double tp,
					double alph, double beta, double temp, double rk, double rq, double rs, double rf);
				KikuchiAikenLRBWrapper() {};
				~KikuchiAikenLRBWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class MaxwellWrapper : UniaxialMaterialWrapper
			{
			public:
				MaxwellWrapper(int tag, double K, double C, double Alpha, double L, int returnD);
				MaxwellWrapper() {};
				~MaxwellWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class MinMaxMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				MinMaxMaterialWrapper(int tag, UniaxialMaterialWrapper^ theMaterial, double min, double max);
				~MinMaxMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class ModIMKPeakOrientedWrapper : UniaxialMaterialWrapper
			{
			public:
				ModIMKPeakOrientedWrapper(int tag, double Ke0, double AlfanPos, double AlfanNeg, double My_pos, double My_neg,			// Updated: Filipe Ribeiro and Andre Barbosa
					double Ls, double Ld, double La, double Lk, double Cs, double Cd, double Ca, double Ck,
					double ThetaPpos, double ThetaPneg, double ThetaPCpos, double ThetaPCneg,
					double ResfacPos, double ResfacNeg, double FracDispPos, double FracDispNeg,
					double DPos, double DNeg, double nFactor);
				ModIMKPeakOrientedWrapper(int tag, double Ke0, double AlfanPos, double AlfanNeg, double My_pos, double My_neg,			// Updated: Filipe Ribeiro and Andre Barbosa 
					double Ls, double Ld, double La, double Lk, double Cs, double Cd, double Ca, double Ck,
					double ThetaPpos, double ThetaPneg, double ThetaPCpos, double ThetaPCneg,
					double ResfacPos, double ResfacNeg, double FracDispPos, double FracDispNeg,
					double DPos, double DNeg);
				~ModIMKPeakOrientedWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class ModIMKPeakOriented02Wrapper : UniaxialMaterialWrapper
			{
			public:
				ModIMKPeakOriented02Wrapper(int tag, double Ke0, double AlfanPos, double AlfanNeg, double My_pos, double My_neg,
					double Ls, double Ld, double La, double Lk, double Cs, double Cd, double Ca, double Ck,
					double ThetaPpos, double ThetaPneg, double ThetaPCpos, double ThetaPCneg,
					double ResfacPos, double ResfacNeg, double FracDispPos, double FracDispNeg,
					double DPos, double DNeg, double C_Fp, double C_Fn, double nFactor);
				ModIMKPeakOriented02Wrapper(int tag, double Ke0, double AlfanPos, double AlfanNeg, double My_pos, double My_neg,
					double Ls, double Ld, double La, double Lk, double Cs, double Cd, double Ca, double Ck,
					double ThetaPpos, double ThetaPneg, double ThetaPCpos, double ThetaPCneg,
					double ResfacPos, double ResfacNeg, double FracDispPos, double FracDispNeg,
					double DPos, double DNeg, double C_Fp, double C_Fn);
				ModIMKPeakOriented02Wrapper(int tag, double Ke0, double AlfanPos, double AlfanNeg, double My_pos, double My_neg,
					double Ls, double Ld, double La, double Lk, double Cs, double Cd, double Ca, double Ck,
					double ThetaPpos, double ThetaPneg, double ThetaPCpos, double ThetaPCneg,
					double ResfacPos, double ResfacNeg, double FracDispPos, double FracDispNeg,
					double DPos, double DNeg);
				ModIMKPeakOriented02Wrapper(int tag, double Ke0, double AlfanPos, double AlfanNeg, double My_pos, double My_neg,
					double Ls, double Ld, double La, double Lk, double Cs, double Cd, double Ca, double Ck,
					double ThetaPpos, double ThetaPneg, double ThetaPCpos, double ThetaPCneg,
					double ResfacPos, double ResfacNeg, double FracDispPos, double FracDispNeg,
					double DPos, double DNeg, double nFactor);
				~ModIMKPeakOriented02Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class ModIMKPinchingWrapper : UniaxialMaterialWrapper
			{
			public:
				ModIMKPinchingWrapper(int tag, double Ke0, double AlfanPos, double AlfanNeg, double My_pos, double My_neg, double FprPos, double FprNeg, double A_pinch,		// Updated: Filipe Ribeiro and Andre Barbosa
					double Ls, double Ld, double La, double Lk, double Cs, double Cd, double Ca, double Ck,
					double ThetaPpos, double ThetaPneg, double ThetaPCpos, double ThetaPCneg,
					double ResfacPos, double ResfacNeg, double FracDispPos, double FracDispNeg,
					double DPos, double DNeg, double nFactor);
				ModIMKPinchingWrapper(int tag, double Ke0, double AlfanPos, double AlfanNeg, double My_pos, double My_neg, double FprPos, double FprNeg, double A_pinch,		// Updated: Filipe Ribeiro and Andre Barbosa
					double Ls, double Ld, double La, double Lk, double Cs, double Cd, double Ca, double Ck,
					double ThetaPpos, double ThetaPneg, double ThetaPCpos, double ThetaPCneg,
					double ResfacPos, double ResfacNeg, double FracDispPos, double FracDispNeg,
					double DPos, double DNeg);
				
				~ModIMKPinchingWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class ModIMKPinching02Wrapper : UniaxialMaterialWrapper
			{
			public:
				ModIMKPinching02Wrapper(int tag, double Ke0, double AlfanPos, double AlfanNeg, double My_pos, double My_neg, double FprPos, double FprNeg, double A_pinch,		// Updated: Filipe Ribeiro and Andre Barbosa
					double Ls, double Ld, double La, double Lk, double Cs, double Cd, double Ca, double Ck,
					double ThetaPpos, double ThetaPneg, double ThetaPCpos, double ThetaPCneg,
					double ResfacPos, double ResfacNeg, double FracDispPos, double FracDispNeg,
					double DPos, double DNeg, double nFactor);
				ModIMKPinching02Wrapper(int tag, double Ke0, double AlfanPos, double AlfanNeg, double My_pos, double My_neg, double FprPos, double FprNeg, double A_pinch,		// Updated: Filipe Ribeiro and Andre Barbosa
					double Ls, double Ld, double La, double Lk, double Cs, double Cd, double Ca, double Ck,
					double ThetaPpos, double ThetaPneg, double ThetaPCpos, double ThetaPCneg,
					double ResfacPos, double ResfacNeg, double FracDispPos, double FracDispNeg,
					double DPos, double DNeg);

				~ModIMKPinching02Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class MultiLinearWrapper : UniaxialMaterialWrapper
			{
			public:
				MultiLinearWrapper(int tag, VectorWrapper^ s, VectorWrapper^ e);
				
				~MultiLinearWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class NewUniaxialMaterialWrapper abstract : UniaxialMaterialWrapper
			{
			public:
				NewUniaxialMaterialWrapper(int tag) {};
				~NewUniaxialMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class OriginCenteredWrapper : UniaxialMaterialWrapper
			{
			public:
				OriginCenteredWrapper(int tag,
					double f1, double e1, double f2,
					double e2, double f3, double e3);

				~OriginCenteredWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class ParallelMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				ParallelMaterialWrapper(int tag,
					int numMaterial,
					array<UniaxialMaterialWrapper^>^ theMaterials,
					VectorWrapper^ theFactors);

				~ParallelMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class PathIndependentMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				PathIndependentMaterialWrapper(int tag,
					UniaxialMaterialWrapper^ theMaterial);

				~PathIndependentMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class Pinching4MaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				Pinching4MaterialWrapper(int tag,
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
					double gammaF4, double gammaFLimit, double gammaE, int DmgCyc);

				Pinching4MaterialWrapper(int tag,
					double stress1p, double strain1p, double stress2p, double strain2p,
					double stress3p, double strain3p, double stress4p, double strain4p,
					double rDispP, double rForceP, double uForceP,
					double gammaK1, double gammaK2, double gammaK3,
					double gammaK4, double gammaKLimit,
					double gammaD1, double gammaD2, double gammaD3,
					double gammaD4, double gammaDLimit,
					double gammaF1, double gammaF2, double gammaF3,
					double gammaF4, double gammaFLimit, double gammaE, int DmgCyc);

				~Pinching4MaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			

			public ref class pyUCLAWrapper : UniaxialMaterialWrapper
			{
			public:
				pyUCLAWrapper(int tag, int soilType, double pult, double y50, double Cd);
				pyUCLAWrapper() {};
				~pyUCLAWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class RambergOsgoodSteelWrapper : UniaxialMaterialWrapper
			{
			public:
				RambergOsgoodSteelWrapper(int tag, double fy, double E0, double rezaA, double rezaN);
				RambergOsgoodSteelWrapper(int tag, double fy, double E0);
				RambergOsgoodSteelWrapper() {};
				~RambergOsgoodSteelWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class ReinforcingSteelWrapper : UniaxialMaterialWrapper
			{
			public:
				ReinforcingSteelWrapper(int tag, double fyield, double fultimate, double youngs, double youngs_hard,
					double estrainhard, double eultimate, int buckModel, double slenderness, double alpha, double r,
					double gama, double Fatigue1, double Fatigue2, double Degrade1,
					double RC1, double RC2, double RC3, double A1, double HardLim);
				ReinforcingSteelWrapper() {};
				~ReinforcingSteelWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class ResilienceLowWrapper : UniaxialMaterialWrapper
			{
			public:
				ResilienceLowWrapper(int tag, double PY, double DPmax, double Pmax, double Ke, double Kd);
				ResilienceLowWrapper() {};
				~ResilienceLowWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class ResilienceMaterialHRWrapper : UniaxialMaterialWrapper
			{
			public:
				ResilienceMaterialHRWrapper(int tag, double  DY, double PY, double DPmax, double Pmax, double Ke, double Kd, double k);
				ResilienceMaterialHRWrapper() {};
				~ResilienceMaterialHRWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class SAWSMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				SAWSMaterialWrapper(int tag,
					double F0, double FI, double DU, double S0,
					double R1, double R2, double R3, double R4,
					double ALPHA, double BETA);
				SAWSMaterialWrapper() {};
				~SAWSMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class SelfCenteringMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				SelfCenteringMaterialWrapper(int tag, double k1, double k2,
					double ActF, double beta, double SlipDef,
					double BearDef, double rBear);
				SelfCenteringMaterialWrapper() {};
				~SelfCenteringMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class SeriesMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				SeriesMaterialWrapper(int tag,
					int numMaterial,
					array<UniaxialMaterialWrapper^>^ theMaterials,
					int maxIter, double tol);
				SeriesMaterialWrapper() {};
				~SeriesMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class ShearPanelMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				ShearPanelMaterialWrapper(int tag,
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
					double gammaF4, double gammaFLimit, double gammaE, double YldStress);
				ShearPanelMaterialWrapper() {};
				~ShearPanelMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class SimpleFractureMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				SimpleFractureMaterialWrapper(int tag, UniaxialMaterialWrapper^ material, double max);
				SimpleFractureMaterialWrapper() {};
				~SimpleFractureMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class SmoothPSConcreteWrapper : UniaxialMaterialWrapper
			{
			public:
				SmoothPSConcreteWrapper(int tag, double fpc, double fpcu, double Ec, double eco, double ecu, double  eta);
				SmoothPSConcreteWrapper() {};
				~SmoothPSConcreteWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class StainlessECThermalWrapper : UniaxialMaterialWrapper
			{
			public:
				StainlessECThermalWrapper(int tag, int gradeTag, double fy, double E0, double fu, double sigInit);
				StainlessECThermalWrapper() {};
				~StainlessECThermalWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class Steel01Wrapper : UniaxialMaterialWrapper
			{
			public:
				Steel01Wrapper(int tag, double fy, double e0, double b);
				Steel01Wrapper(int tag, double fy, double e0, double b,
					double a1, double a2,
					double a3, double a4);
				~Steel01Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class Steel01ThermalWrapper : UniaxialMaterialWrapper
			{
			public:
				Steel01ThermalWrapper(int tag, double fy, double e0, double b);
				Steel01ThermalWrapper(int tag, double fy, double e0, double b,
					double a1, double a2,
					double a3, double a4);
				~Steel01ThermalWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};
			
			public ref class Steel02Wrapper : UniaxialMaterialWrapper
			{
			public:
				Steel02Wrapper(int tag,
					double fy, double E0, double b,
					double R0, double cR1, double cR2);
				Steel02Wrapper(int tag,
					double fy, double E0, double b,
					double R0, double cR1, double cR2,
					double a1, double a2, double a3, double a4, double sigInit);
				~Steel02Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class Steel02ThermalWrapper : UniaxialMaterialWrapper
			{
			public:
				Steel02ThermalWrapper(int tag,
					double fy, double E0, double b,
					double R0, double cR1, double cR2);
				Steel02ThermalWrapper(int tag,
					double fy, double E0, double b,
					double R0, double cR1, double cR2,
					double a1, double a2, double a3, double a4, double sigInit);
				~Steel02ThermalWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class Steel2Wrapper : UniaxialMaterialWrapper
			{
			public:
				Steel2Wrapper(int tag,
					double fy, double E0, double b,
					double R0, double cR1, double cR2,
					double a1, double a2, double a3, double a4, double sigInit);
				Steel2Wrapper(int tag,
					double fy, double E0, double b,
					double R0, double cR1, double cR2);
				~Steel2Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class Steel03Wrapper : UniaxialMaterialWrapper
			{
			public:
				Steel03Wrapper(int tag, double fy, double E0, double b, double r, double cR1, double cR2,
					double a1 , double a2 ,
					double a3 , double a4 );
				~Steel03Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class Steel4Wrapper : UniaxialMaterialWrapper
			{
			public:
				Steel4Wrapper(int tag,
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
					double sig_init);
				~Steel4Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class SteelBRBWrapper : UniaxialMaterialWrapper
			{
			public:
				SteelBRBWrapper(int tag,
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
					double Tol);
				~SteelBRBWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class SteelECThermalWrapper : UniaxialMaterialWrapper
			{
			public:
				SteelECThermalWrapper(int tag, int typeTag, double fy, double E0,
					double a1 , double a2 ,
					double a3 , double a4);
				~SteelECThermalWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class SteelMPWrapper : UniaxialMaterialWrapper
			{
			public:
				SteelMPWrapper(int tag, double FY, double E, double B, double R, double CR1, double CR2, double A1, double A2);
				~SteelMPWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class SteelMPFWrapper : UniaxialMaterialWrapper
			{
			public:
				SteelMPFWrapper(int tag, double sigyieldp, double sigyieldn,
					double E0, double bp, double bn, double R0, double a1, double a2);
				SteelMPFWrapper(int tag, double sigyieldp, double sigyieldn,
					double E0, double bp, double bn, double R0, double a1, double a2, double a3, double a4, double a5, double a6);
				~SteelMPFWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class UniaxialJ2PlasticityWrapper : UniaxialMaterialWrapper
			{
			public:
				UniaxialJ2PlasticityWrapper(int tag, double E, double sigmaY,
					double Hkin, double Hiso);
				~UniaxialJ2PlasticityWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class ViscousDamperWrapper : UniaxialMaterialWrapper
			{
			public:
				ViscousDamperWrapper(int tag, double K, double C, double Alpha, double LGap, double NM, double RelTol, double AbsTol, double MaxHalf);
				~ViscousDamperWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class ViscousMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				ViscousMaterialWrapper(int tag, double C, double Alpha, double minVel);
				~ViscousMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};

			internal:

			};

			public ref class ConcreteL01Wrapper : UniaxialMaterialWrapper
			{
			public:
				ConcreteL01Wrapper(int tag, double fpc, double eco);
				~ConcreteL01Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};
			};

			public ref class ConcreteZ01Wrapper : UniaxialMaterialWrapper
			{
			public:
				ConcreteZ01Wrapper(int tag, double fpc, double eco);
				~ConcreteZ01Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};
			};

			public ref class SteelZ01Wrapper : UniaxialMaterialWrapper
			{
			public:
				SteelZ01Wrapper(int tag, double fy, double E0, double fpc, double rou,
					double ac,
					double rc);
				~SteelZ01Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};
			};

			public ref class TendonL01Wrapper : UniaxialMaterialWrapper
			{
			public:
				TendonL01Wrapper(int tag, double fy, double E0, double fpu, double rou, double epsp,
					double ac,
					double rc);
				~TendonL01Wrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;
				};
			};

			public delegate int MDSetTrialStrain(double, double);
			delegate int UMDSetTrialStrain(double, double);

			public delegate double MDGetStress();
			delegate double UMDGetStress(void);

			public delegate double MDGetTangent();
			delegate double UMDGetTangent(void);

			public delegate double MDGetInitialTangent();
			delegate double UMDGetInitialTangent(void);

			public delegate double MDGetStrain();
			delegate double UMDGetStrain(void);

			public delegate double MDGetDampTangent();
			delegate double UMDGetDampTangent(void);

			public delegate double MDGetStrainRate();
			delegate double UMDGetStrainRate(void);

			public delegate double MDGetRho();
			delegate double UMDGetRho(void);

			public delegate int MDCommitState();
			delegate int UMDCommitState(void);

			public delegate int MDRevertToLastCommit();
			delegate int UMDRevertToLastCommit(void);

			public delegate int MDRevertToStart();
			delegate int UMDRevertToStart(void);

			public delegate UniaxialMaterialWrapper^ MDGetCopy();
			delegate UniaxialMaterial* UMDGetCopy(void);

			public delegate void MDPrint();
			delegate void UMDPrint(void);

			public ref class ExternalUniaxialMaterialWrapper : UniaxialMaterialWrapper
			{
			public:
				ExternalUniaxialMaterialWrapper(int tag);
				
				~ExternalUniaxialMaterialWrapper() {
					if (_UniaxialMaterial != 0)
						delete _UniaxialMaterial;

					umdSetTrialStrain_GC.Free();
					umdGetStress_GC.Free();
					umdGetTangent_GC.Free();
					umdGetInitialTangent_GC.Free();
					umdGetStrain_GC.Free();
					umdCommitState_GC.Free();
					umdRevertToLastCommit_GC.Free();
					umdRevertToStart_GC.Free();
					umdGetCopy_GC.Free();
					umdPrintSelf_GC.Free();
					umdGetDampTangent_GC.Free();
					umdGetStrainRate_GC.Free();
					umdGetRho_GC.Free();
				};

				void SetLinks(
					MDSetTrialStrain ^ _MDSetTrialStrain,
					MDGetStress ^ _MDGetStress,
					MDGetTangent ^ _MDGetTangent,
					MDGetInitialTangent ^ _MDGetInitialTangent,
					MDGetDampTangent^ _MDGetDampTangent,
					MDGetStrain ^ _MDGetStrain,
					MDGetStrainRate^ __MDGetStrainRate,
					MDGetRho^ _MDGetRho,
					MDCommitState ^ _MDCommitState,
					MDRevertToLastCommit ^ _MDRevertToLastCommit,
					MDRevertToStart ^ _MDRevertToStart,
					MDGetCopy ^ _MDGetCopy,
					MDPrint ^ _MDPrint)
				{
					this->_MDSetTrialStrain = _MDSetTrialStrain;
					UMDSetTrialStrain^ umdSetTrialStrain = gcnew UMDSetTrialStrain(this, &ExternalUniaxialMaterialWrapper::SetTrialStrainMethod);
					umdSetTrialStrain_GC = GCHandle::Alloc(umdSetTrialStrain);
					IntPtr ip = Marshal::GetFunctionPointerForDelegate(umdSetTrialStrain);
					UMSetTrialStrain setTrialStrain = static_cast<UMSetTrialStrain>(ip.ToPointer());

					this->_MDGetStress = _MDGetStress;
					UMDGetStress^ umdGetStress = gcnew UMDGetStress(this, &ExternalUniaxialMaterialWrapper::GetStressMethod);
					umdGetStress_GC = GCHandle::Alloc(umdGetStress);
					ip = Marshal::GetFunctionPointerForDelegate(umdGetStress);
					UMGetStress getStress = static_cast<UMGetStress>(ip.ToPointer());

					this->_MDGetTangent = _MDGetTangent;
					UMDGetTangent^ umdGetTangent = gcnew UMDGetTangent(this, &ExternalUniaxialMaterialWrapper::GetTangentMethod);
					umdGetTangent_GC = GCHandle::Alloc(umdGetTangent);
					ip = Marshal::GetFunctionPointerForDelegate(umdGetTangent);
					UMGetTangent getTangent = static_cast<UMGetTangent>(ip.ToPointer());

					this->_MDGetInitialTangent = _MDGetInitialTangent;
					UMDGetInitialTangent^ umdGetInitialTangent = gcnew UMDGetInitialTangent(this, &ExternalUniaxialMaterialWrapper::GetInitialTangentMethod);
					umdGetInitialTangent_GC = GCHandle::Alloc(umdGetInitialTangent);
					ip = Marshal::GetFunctionPointerForDelegate(umdGetInitialTangent);
					UMGetInitialTangent getInitialTangent = static_cast<UMGetInitialTangent>(ip.ToPointer());

					this->_MDGetStrain = _MDGetStrain;
					UMDGetStrain^ umdGetStrain = gcnew UMDGetStrain(this, &ExternalUniaxialMaterialWrapper::GetStrainMethod);
					umdGetStrain_GC = GCHandle::Alloc(umdGetStrain);
					ip = Marshal::GetFunctionPointerForDelegate(umdGetStrain);
					UMGetStrain getStrain = static_cast<UMGetStrain>(ip.ToPointer());

					this->_MDCommitState = _MDCommitState;
					UMDCommitState^ umdCommitState = gcnew UMDCommitState(this, &ExternalUniaxialMaterialWrapper::CommitStateMethod);
					umdCommitState_GC = GCHandle::Alloc(umdCommitState);
					ip = Marshal::GetFunctionPointerForDelegate(umdCommitState);
					UMCommitState commitState = static_cast<UMCommitState>(ip.ToPointer());

					this->_MDRevertToLastCommit = _MDRevertToLastCommit;
					UMDRevertToLastCommit^ umdRevertToLastCommit = gcnew UMDRevertToLastCommit(this, &ExternalUniaxialMaterialWrapper::RevertToLastCommitMethod);
					umdRevertToLastCommit_GC = GCHandle::Alloc(umdRevertToLastCommit);
					ip = Marshal::GetFunctionPointerForDelegate(umdRevertToLastCommit);
					UMRevertToLastCommit revertToLastCommit = static_cast<UMRevertToLastCommit>(ip.ToPointer());

					this->_MDRevertToStart = _MDRevertToStart;
					UMDRevertToStart^ umdRevertToStart = gcnew UMDRevertToStart(this, &ExternalUniaxialMaterialWrapper::RevertToStartMethod);
					umdRevertToStart_GC = GCHandle::Alloc(umdRevertToStart);
					ip = Marshal::GetFunctionPointerForDelegate(umdRevertToStart);
					UMRevertToStart revertToStart = static_cast<UMRevertToStart>(ip.ToPointer());

					this->_MDGetCopy = _MDGetCopy;
					UMDGetCopy^ umdGetCopy = gcnew UMDGetCopy(this, &ExternalUniaxialMaterialWrapper::GetCopyMethod);
					umdGetCopy_GC = GCHandle::Alloc(umdGetCopy);
					ip = Marshal::GetFunctionPointerForDelegate(umdGetCopy);
					UMGetCopy getCopy = static_cast<UMGetCopy>(ip.ToPointer());

					this->_MDPrint = _MDPrint;
					UMDPrint^ umdPrintSelf = gcnew UMDPrint(this, &ExternalUniaxialMaterialWrapper::PrintMethod);
					umdPrintSelf_GC = GCHandle::Alloc(umdPrintSelf);
					ip = Marshal::GetFunctionPointerForDelegate(umdPrintSelf);
					UMPrint printSelf = static_cast<UMPrint>(ip.ToPointer());

					this->_MDGetDampTangent = _MDGetDampTangent;
					UMDGetDampTangent^ umdGetDampTangent = gcnew UMDGetDampTangent(this, &ExternalUniaxialMaterialWrapper::GetDampTangentMethod);
					umdGetDampTangent_GC = GCHandle::Alloc(umdGetDampTangent);
					ip = Marshal::GetFunctionPointerForDelegate(umdGetDampTangent);
					UMGetDampTangent getDampTangent = static_cast<UMGetDampTangent>(ip.ToPointer());

					this->_MDGetStrainRate = _MDGetStrainRate;
					UMDGetStrainRate^ umdGetStrainRate = gcnew UMDGetStrainRate(this, &ExternalUniaxialMaterialWrapper::GetStrainRateMethod);
					umdGetStrainRate_GC = GCHandle::Alloc(umdGetStrainRate);
					ip = Marshal::GetFunctionPointerForDelegate(umdGetStrainRate);
					UMGetStrainRate getStrainRate = static_cast<UMGetStrainRate>(ip.ToPointer());

					this->_MDGetRho = _MDGetRho;
					UMDGetRho^ umdGetRho = gcnew UMDGetRho(this, &ExternalUniaxialMaterialWrapper::GetRhoMethod);
					umdGetRho_GC = GCHandle::Alloc(umdGetRho);
					ip = Marshal::GetFunctionPointerForDelegate(umdGetRho);
					UMGetRho getRho = static_cast<UMGetRho>(ip.ToPointer());

					((ExternalUniaxialMaterial*)_UniaxialMaterial)->SetLinks(setTrialStrain, getStress, getTangent, getInitialTangent,
						getDampTangent, getStrain, getStrainRate, getRho, commitState,
						revertToLastCommit, revertToStart, getCopy, printSelf);

				};
			private:
				GCHandle umdSetTrialStrain_GC;
				GCHandle umdGetStress_GC;
				GCHandle umdGetTangent_GC;
				GCHandle umdGetInitialTangent_GC;
				GCHandle umdGetStrain_GC;
				GCHandle umdCommitState_GC;
				GCHandle umdRevertToLastCommit_GC;
				GCHandle umdRevertToStart_GC;
				GCHandle umdGetCopy_GC;
				GCHandle umdPrintSelf_GC;

				GCHandle umdGetDampTangent_GC;
				GCHandle umdGetStrainRate_GC;
				GCHandle umdGetRho_GC;

				MDSetTrialStrain^ _MDSetTrialStrain;
				MDGetStress^ _MDGetStress;
				MDGetTangent^ _MDGetTangent;
				MDGetInitialTangent^ _MDGetInitialTangent;
				MDGetStrain^ _MDGetStrain;
				MDCommitState^ _MDCommitState;
				MDRevertToLastCommit^ _MDRevertToLastCommit;
				MDRevertToStart^ _MDRevertToStart;
				MDGetCopy^ _MDGetCopy;
				MDPrint^ _MDPrint;
				MDGetDampTangent^ _MDGetDampTangent;
				MDGetStrainRate^ _MDGetStrainRate;
				MDGetRho^ _MDGetRho;
			//protected:

			internal:

				double GetDampTangentMethod(void) {
					return this->_MDGetDampTangent();
				};
				double GetStrainRateMethod(void) {
					return this->_MDGetStrainRate();
				};
				double GetRhoMethod(void) {
					return this->_MDGetRho();
				};

				int SetTrialStrainMethod(double strain, double strainRate) {
					return this->_MDSetTrialStrain(strain, strainRate);
				}

				double GetStressMethod() {
					return this->_MDGetStress();
				}

				double GetTangentMethod() {
					return this->_MDGetTangent();
				}

				double GetInitialTangentMethod() {
					return this->_MDGetInitialTangent();
				}

				double GetStrainMethod() {
					return this->_MDGetInitialTangent();
				}

				int CommitStateMethod() {
					return this->_MDCommitState();
				}

				int RevertToLastCommitMethod() {
					return this->_MDRevertToLastCommit();
				}

				int RevertToStartMethod() {
					return this->_MDRevertToStart();
				}
				
				UniaxialMaterial* GetCopyMethod() {
					UniaxialMaterialWrapper^ copy = _MDGetCopy();
					return copy->_UniaxialMaterial;
				}

				void PrintMethod() {
					return this->_MDPrint();
				}
			
			};
			

			

			

			
		}
	}
	
}
