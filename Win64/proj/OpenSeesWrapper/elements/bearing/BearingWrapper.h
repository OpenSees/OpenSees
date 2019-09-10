#pragma once

#include <elastomericBearing/ElastomericBearingBoucWen2d.h>
#include <elastomericBearing/ElastomericBearingBoucWen3d.h>
#include <elastomericBearing/ElastomericBearingBoucWenMod3d.h>
#include <elastomericBearing/ElastomericBearingPlasticity2d.h>
#include <elastomericBearing/ElastomericBearingPlasticity3d.h>
#include <elastomericBearing/ElastomericBearingUFRP2d.h>
#include <elastomericBearing/ElastomericX.h>
#include <elastomericBearing/HDR.h>
#include <elastomericBearing/LeadRubberX.h>

#include <frictionBearing\frictionModel\FrictionModel.h>
#include <frictionBearing\frictionModel\Coulomb.h>
#include <frictionBearing\frictionModel\FrictionResponse.h>
#include <frictionBearing\frictionModel\VelDependent.h>
#include <frictionBearing\frictionModel\VelDepMultiLinear.h>
#include <frictionBearing\frictionModel\VelNormalFrcDep.h>
#include <frictionBearing\frictionModel\VelPressureDep.h>

#include <frictionBearing\FlatSliderSimple2d.h>
#include <frictionBearing\FlatSliderSimple3d.h>
#include <frictionBearing\FPBearingPTV.h>
#include <frictionBearing\MultiFP2d.h>
#include <frictionBearing\RJWatsonEQS2d.h>
#include <frictionBearing\RJWatsonEQS3d.h>
#include <frictionBearing\SingleFPSimple2d.h>
#include <frictionBearing\SingleFPSimple3d.h>
#include <frictionBearing\TFP_Bearing.h>
#include <frictionBearing\TFP_Bearing2d.h>
#include <frictionBearing\TPB1D.h>
#include <frictionBearing\TripleFrictionPendulum.h>

#include <HUelements\KikuchiBearing.h>
#include <HUelements\MultipleNormalSpring.h>
#include <HUelements\MultipleShearSpring.h>
#include <HUelements\YamamotoBiaxialHDR.h>



#include "../../taggeds/TaggedObjectWrapper.h"
#include "../../actors/IMovableObjectWrapper.h"
#include "../ElementWrapper.h"
#include "../../materials/section/SectionForceDeformationWrapper.h"
#include "../../materials/uniaxial/UniaxialMaterialWrapper.h"
#include "../../materials/ndmaterial/NDMaterialWrapper.h"
#include "../beamColumn/beamIntegration/BeamIntegrationWrapper.h"
#include "../../matrix/IDWrapper.h"
#include "../../matrix/MatrixWrapper.h"
#include "../../matrix/VectorWrapper.h"


using namespace System;
using namespace OpenSees;
using namespace OpenSees::Elements::BeamIntegrations;
using namespace OpenSees::Materials::Sections;
using namespace OpenSees::Materials::Uniaxials;
using namespace OpenSees::Materials::NDMaterials;

namespace OpenSees {
	namespace Elements {

#pragma region ElastomericBearing 

		public ref class ElastomericBearingBoucWen2dWrapper : ElementWrapper
		{

		public:
			ElastomericBearingBoucWen2dWrapper(int tag, int Nd1, int Nd2,
				double kInit, double qd, double alpha1,
				array<UniaxialMaterialWrapper^>^ theMaterials);
			ElastomericBearingBoucWen2dWrapper(int tag, int Nd1, int Nd2,
				double kInit, double qd, double alpha1,
				array<UniaxialMaterialWrapper^>^ theMaterials,
				VectorWrapper^ y , VectorWrapper^ x ,
				double alpha2 , double mu ,
				double eta , double beta ,
				double gamma , double shearDistI ,
				int addRayleigh , double mass ,
				int maxIter , double tol );
			ElastomericBearingBoucWen2dWrapper() {};
			~ElastomericBearingBoucWen2dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class ElastomericBearingBoucWen3dWrapper : ElementWrapper
		{
		public:
			ElastomericBearingBoucWen3dWrapper(int tag, int Nd1, int Nd2,
				double kInit, double qd, double alpha1,
				array<UniaxialMaterialWrapper^>^ theMaterials, VectorWrapper^ y);
			ElastomericBearingBoucWen3dWrapper(int tag, int Nd1, int Nd2,
				double kInit, double qd, double alpha1,
				array<UniaxialMaterialWrapper^>^ theMaterials,
				VectorWrapper^ y, VectorWrapper^ x,
				double alpha2, double mu,
				double eta, double beta,
				double gamma, double shearDistI,
				int addRayleigh, double mass,
				int maxIter, double tol);
			ElastomericBearingBoucWen3dWrapper() {};
			~ElastomericBearingBoucWen3dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class ElastomericBearingBoucWenMod3dWrapper : ElementWrapper
		{
		public:
			ElastomericBearingBoucWenMod3dWrapper(int tag, int Nd1, int Nd2,
				double kInit, double fy,
				double Gr, double Kbulk, double D1, double D2,
				double ts, double tr, int n, double alpha1);
			ElastomericBearingBoucWenMod3dWrapper(int tag, int Nd1, int Nd2,
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
				int maxIter, double tol);
			ElastomericBearingBoucWenMod3dWrapper() {};
			~ElastomericBearingBoucWenMod3dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class ElastomericBearingPlasticity2dWrapper : ElementWrapper
		{
		public:
			ElastomericBearingPlasticity2dWrapper(int tag, int Nd1, int Nd2,
				double kInit, double qd, double alpha1,
				array<UniaxialMaterialWrapper^>^ theMaterials);
			ElastomericBearingPlasticity2dWrapper(int tag, int Nd1, int Nd2,
				double kInit, double qd, double alpha1,
				array<UniaxialMaterialWrapper^>^ theMaterials,
				VectorWrapper^ y, VectorWrapper^ x,
				double alpha2, double mu,
				double shearDistI,
				int addRayleigh, double mass);
			ElastomericBearingPlasticity2dWrapper() {};
			~ElastomericBearingPlasticity2dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class ElastomericBearingPlasticity3dWrapper : ElementWrapper
		{
		public:
			ElastomericBearingPlasticity3dWrapper(int tag, int Nd1, int Nd2,
				double kInit, double qd, double alpha1,
				array<UniaxialMaterialWrapper^>^ theMaterials, VectorWrapper^ y);
			ElastomericBearingPlasticity3dWrapper(int tag, int Nd1, int Nd2,
				double kInit, double qd, double alpha1,
				array<UniaxialMaterialWrapper^>^ theMaterials,
				VectorWrapper^ y, VectorWrapper^ x,
				double alpha2, double mu,
				double shearDistI,
				int addRayleigh, double mass);
			ElastomericBearingPlasticity3dWrapper() {};
			~ElastomericBearingPlasticity3dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class ElastomericBearingUFRP2dWrapper : ElementWrapper
		{
		public:
			ElastomericBearingUFRP2dWrapper(int tag, int Nd1, int Nd2, double uy,
				double a1, double a2, double a3, double a4, double a5,
				double b, double c, array<UniaxialMaterialWrapper^>^ theMaterials);
			ElastomericBearingUFRP2dWrapper(int tag, int Nd1, int Nd2, double uy,
				double a1, double a2, double a3, double a4, double a5,
				double b, double c, array<UniaxialMaterialWrapper^>^ theMaterials,
				VectorWrapper^ y, VectorWrapper^ x,
				double eta, double beta, double gamma,
				double shearDistI, int addRayleigh, double mass,
				int maxIter, double tol);
			ElastomericBearingUFRP2dWrapper() {};
			~ElastomericBearingUFRP2dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class ElastomericXWrapper : ElementWrapper
		{
		public:
			ElastomericXWrapper(int eleTag, int Nd1, int Nd2, double qd, double alpha, double Gr, double Kbulk,
				double D1, double D2, double ts, double tr, double n, VectorWrapper^ y);
			ElastomericXWrapper(int eleTag, int Nd1, int Nd2, double qd, double alpha, double Gr, double Kbulk,
				double D1, double D2, double ts, double tr, double n, VectorWrapper^ y, VectorWrapper^ x,
				double kc, double PhiM, double ac, double sDratio, double m,
				double cd, double tc, int tag1, int tag2, int tag3, int tag4);
			ElastomericXWrapper() {};
			~ElastomericXWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class HDRWrapper : ElementWrapper
		{
		public:
			HDRWrapper(int tag, int Nd1, int Nd2, double Gr, double Kbulk, double D1, double D2, double ts,
				double tr, int n, double a1, double a2, double a3, double b1, double b2, double b3,
				double c1, double c2, double c3, double c4, VectorWrapper^ y);
			HDRWrapper(int tag, int Nd1, int Nd2, double Gr, double Kbulk, double D1, double D2, double ts,
				double tr, int n, double a1, double a2, double a3, double b1, double b2, double b3,
				double c1, double c2, double c3, double c4, VectorWrapper^ y, VectorWrapper^ x ,
				double kc, double PhiM, double ac, double sDratio, double m,
				double tc);
			HDRWrapper() {};
			~HDRWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class LeadRubberXWrapper : ElementWrapper
		{
		public:
			LeadRubberXWrapper(int eleTag, int Nd1, int Nd2, double qd, double alpha, double Gr, double Kbulk,
				double D1, double D2, double ts, double tr, double n, VectorWrapper^ y);
			LeadRubberXWrapper(int eleTag, int Nd1, int Nd2, double qd, double alpha, double Gr, double Kbulk,
				double D1, double D2, double ts, double tr, double n, VectorWrapper^ y, VectorWrapper^ x,
				double kc, double PhiM, double ac, double sDratio, double m,
				double cd, double tc, double qL, double cL, double kS,
				double aS, int tag1, int tag2, int tag3, int tag4, int tag5);
			LeadRubberXWrapper() {};
			~LeadRubberXWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

#pragma endregion

#pragma region FrictionModel
		public ref class FrictionModelWrapper abstract: TaggedObjectWrapper, IMovableObjectWrapper
		{
		public:
			FrictionModelWrapper() { };
			~FrictionModelWrapper() {
				if (_FrictionModel != 0)
					delete _FrictionModel;
			};

		internal:
			FrictionModel * _FrictionModel;
		private:

		};

		public ref class CoulombWrapper abstract : FrictionModelWrapper
		{
		public:
			CoulombWrapper(int tag, double mu);
			CoulombWrapper() {};
			~CoulombWrapper() {
				if (_FrictionModel != 0)
					delete _FrictionModel;
			};

		internal:
			
		private:

		};

		public ref class VelDependentWrapper abstract : FrictionModelWrapper
		{
		public:
			VelDependentWrapper(int tag, double muSlow, double muFast, double transRate);
			VelDependentWrapper() {};
			~VelDependentWrapper() {
				if (_FrictionModel != 0)
					delete _FrictionModel;
			};

		internal:
			
		private:

		};

		public ref class VelDepMultiLinearWrapper abstract : FrictionModelWrapper
		{
		public:
			VelDepMultiLinearWrapper(int tag,
				 VectorWrapper^ velocityPoints,
				 VectorWrapper^ frictionPoints);
			VelDepMultiLinearWrapper() {};
			~VelDepMultiLinearWrapper() {
				if (_FrictionModel != 0)
					delete _FrictionModel;
			};

		internal:
			
		private:

		};

		public ref class VelNormalFrcDepWrapper abstract : FrictionModelWrapper
		{
		public:
			VelNormalFrcDepWrapper(int tag,
				double aSlow, double nSlow, double aFast, double nFast,
				double alpha0, double alpha1, double alpha2, double maxMuFact);
			VelNormalFrcDepWrapper() {};
			~VelNormalFrcDepWrapper() {
				if (_FrictionModel != 0)
					delete _FrictionModel;
			};

		internal:

		private:

		};

		public ref class VelPressureDepWrapper abstract : FrictionModelWrapper
		{
		public:
			VelPressureDepWrapper(int tag, double muSlow, double muFast0, double A,
				double deltaMu, double alpha, double transRate);
			VelPressureDepWrapper() {};
			~VelPressureDepWrapper() {
				if (_FrictionModel != 0)
					delete _FrictionModel;
			};

		internal:

		private:

		};

#pragma endregion

#pragma region FrictionBearing
		public ref class FlatSliderSimple2dWrapper : ElementWrapper
		{

		public:
			FlatSliderSimple2dWrapper(int tag, int Nd1, int Nd2,
				FrictionModelWrapper^ theFrnMdl, double kInit,
				array<UniaxialMaterialWrapper^>^ theMaterials);
			FlatSliderSimple2dWrapper(int tag, int Nd1, int Nd2,
				FrictionModelWrapper^ theFrnMdl, double kInit,
				array<UniaxialMaterialWrapper^>^ theMaterials,
				VectorWrapper^ y, VectorWrapper^ x,
				double shearDistI,
				int addRayleigh, double mass,
				int maxIter, double tol);
			FlatSliderSimple2dWrapper() {};
			~FlatSliderSimple2dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class FlatSliderSimple3dWrapper : ElementWrapper
		{

		public:
			FlatSliderSimple3dWrapper(int tag, int Nd1, int Nd2,
				FrictionModelWrapper^ theFrnMdl, double kInit,
				array<UniaxialMaterialWrapper^>^ theMaterials);
			FlatSliderSimple3dWrapper(int tag, int Nd1, int Nd2,
				FrictionModelWrapper^ theFrnMdl, double kInit,
				array<UniaxialMaterialWrapper^>^ theMaterials,
				VectorWrapper^ y, VectorWrapper^ x,
				double shearDistI,
				int addRayleigh, double mass,
				int maxIter, double tol);
			FlatSliderSimple3dWrapper() {};
			~FlatSliderSimple3dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class FPBearingPTVWrapper : ElementWrapper
		{

		public:
			FPBearingPTVWrapper(int tag, int Nd1, int Nd2, double MuReference,
				int IsPDependent, double refP, int IsTDependent, double Diffusivity_Steel,
				double Conductivity_Steel, int IsVDependent, double rate_v_mu, double Reff,
				double r_Contact, double kInit,
				UniaxialMaterialWrapper^ theMatA, UniaxialMaterialWrapper^ theMatB,
				UniaxialMaterialWrapper^ theMatC, UniaxialMaterialWrapper^ theMatD);
			FPBearingPTVWrapper(int tag, int Nd1, int Nd2, double MuReference,
				int IsPDependent, double refP, int IsTDependent, double Diffusivity_Steel,
				double Conductivity_Steel, int IsVDependent, double rate_v_mu, double Reff,
				double r_Contact, double kInit,
				UniaxialMaterialWrapper^ theMatA, UniaxialMaterialWrapper^ theMatB,
				UniaxialMaterialWrapper^ theMatC, UniaxialMaterialWrapper^ theMatD,
				VectorWrapper^ x, VectorWrapper^  y,
				double shearDistI,
				int addRayleigh, double mass,
				int maxIter, double tol, int unit);
			FPBearingPTVWrapper() {};
			~FPBearingPTVWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class MultiFP2dWrapper : ElementWrapper
		{

		public:
			MultiFP2dWrapper(int tag,
				int Nd1, int Nd2,
				UniaxialMaterialWrapper^ theFrictionModel,
				UniaxialMaterialWrapper^ theVerticalModel,
				double w0, int axialCase);
			MultiFP2dWrapper(int tag,
				int Nd1, int Nd2,
				int type,
				VectorWrapper^ R,
				VectorWrapper^ h,
				VectorWrapper^ D,
				VectorWrapper^ d,
				VectorWrapper^ mu,
				double Kvert,
				double w0, int axialCase);
			MultiFP2dWrapper() {};
			~MultiFP2dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class RJWatsonEQS2dWrapper : ElementWrapper
		{

		public:
			RJWatsonEQS2dWrapper(int tag, int Nd1, int Nd2,
				FrictionModelWrapper^ theFrnMdl, double kInit,
				array<UniaxialMaterialWrapper^>^ theMaterials);
			RJWatsonEQS2dWrapper(int tag, int Nd1, int Nd2,
				FrictionModelWrapper^ theFrnMdl, double kInit,
				array<UniaxialMaterialWrapper^>^ theMaterials,
				VectorWrapper^  y, VectorWrapper^  x,
				double shearDistI,
				int addRayleigh, double mass,
				int maxIter, double tol,
				double kFactUplift);
			RJWatsonEQS2dWrapper() {};
			~RJWatsonEQS2dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class RJWatsonEQS3dWrapper : ElementWrapper
		{

		public:
			RJWatsonEQS3dWrapper(int tag, int Nd1, int Nd2,
				FrictionModelWrapper^ theFrnMdl, double kInit,
				array<UniaxialMaterialWrapper^>^ theMaterials);
			RJWatsonEQS3dWrapper(int tag, int Nd1, int Nd2,
				FrictionModelWrapper^ theFrnMdl, double kInit,
				array<UniaxialMaterialWrapper^>^ theMaterials,
				VectorWrapper^  y, VectorWrapper^  x,
				double shearDistI,
				int addRayleigh, double mass,
				int maxIter, double tol,
				double kFactUplift);
			RJWatsonEQS3dWrapper() {};
			~RJWatsonEQS3dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class SingleFPSimple2dWrapper : ElementWrapper
		{

		public:
			SingleFPSimple2dWrapper(int tag, int Nd1, int Nd2,
				FrictionModelWrapper^ theFrnMdl, double Reff, double kInit,
				array<UniaxialMaterialWrapper^>^ theMaterials);
			SingleFPSimple2dWrapper(int tag, int Nd1, int Nd2,
				FrictionModelWrapper^ theFrnMdl, double Reff, double kInit,
				array<UniaxialMaterialWrapper^>^ theMaterials,
				VectorWrapper^  y, VectorWrapper^  x,
				double shearDistI , int addRayleigh ,
				int inclVertDisp , double mass ,
				int maxIter , double tol ,
				double kFactUplift );
			SingleFPSimple2dWrapper() {};
			~SingleFPSimple2dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class SingleFPSimple3dWrapper : ElementWrapper
		{

		public:
			SingleFPSimple3dWrapper(int tag, int Nd1, int Nd2,
				FrictionModelWrapper^ theFrnMdl, double Reff, double kInit,
				array<UniaxialMaterialWrapper^>^ theMaterials);
			SingleFPSimple3dWrapper(int tag, int Nd1, int Nd2,
				FrictionModelWrapper^ theFrnMdl, double Reff, double kInit,
				array<UniaxialMaterialWrapper^>^ theMaterials,
				VectorWrapper^  y, VectorWrapper^  x,
				double shearDistI,
				int addRayleigh, int inclVertDisp, double mass,
				int maxIter, double tol,
				double kFactUplift);
			SingleFPSimple3dWrapper() {};
			~SingleFPSimple3dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class TFP_BearingWrapper : ElementWrapper
		{

		public:
			TFP_BearingWrapper(int tag,
				int Nd1, int Nd2,
				array<double>^ r,
				array<double>^ dio,
				array<double>^ di,
				array<double>^ mu,
				array<double>^ h,
				double H0,
				double a,
				double K,
				double vYield);
			TFP_BearingWrapper() {};
			~TFP_BearingWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class TFP_Bearing2dWrapper : ElementWrapper
		{

		public:
			TFP_Bearing2dWrapper(int tag,
				int Nd1, int Nd2,
				array<double>^ r,
				array<double>^ dio,
				array<double>^ di,
				array<double>^ mu,
				array<double>^ h,
				double H0,
				double a,
				double K,
				double vYield);
			TFP_Bearing2dWrapper() {};
			~TFP_Bearing2dWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class TPB1DWrapper : ElementWrapper
		{

		public:
			TPB1DWrapper(int tag,
				int Nd1,
				int Nd2,
				int dir,
				array<double>^ mu,
				array<double>^ R,
				array<double>^ h,
				array<double>^ D,
				array<double>^ d,
				double W);
			TPB1DWrapper() {};
			~TPB1DWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class TripleFrictionPendulumWrapper : ElementWrapper
		{

		public:
			TripleFrictionPendulumWrapper(int tag,
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
				double tol);
			TripleFrictionPendulumWrapper() {};
			~TripleFrictionPendulumWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class KikuchiBearingWrapper : ElementWrapper
		{

		public:
			KikuchiBearingWrapper(int Tag, int Nd1, int Nd2,
				int Shape, double Size, double TotalRubber, double TotalHeight,
				int NMSS, UniaxialMaterialWrapper^ MatMSS, double LimDisp,
				int NMNS, UniaxialMaterialWrapper^ MatMNS, double Lambda,
				VectorWrapper^ OriYp, VectorWrapper^ OriX, double Mass,
				bool IfPDInput, bool IfTilt,
				double AdjCi, double AdjCj,
				bool IfBalance, double LimFo, double LimFi, int NIter);
			KikuchiBearingWrapper() {};
			~KikuchiBearingWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class MultipleNormalSpringWrapper : ElementWrapper
		{

		public:
			MultipleNormalSpringWrapper(int Tag, int Nd1, int Nd2,
				int NDivide,
				UniaxialMaterialWrapper^ Material,
				int Shape,
				double Size,
				double Lambda,
				VectorWrapper^  OriYp, VectorWrapper^  OriX,
				double Mass);
			MultipleNormalSpringWrapper() {};
			~MultipleNormalSpringWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class MultipleShearSpringWrapper : ElementWrapper
		{

		public:
			MultipleShearSpringWrapper(int Tag, int Nd1, int Nd2,
				int NSpring,
				UniaxialMaterialWrapper^ Material,
				double LimDisp,
				VectorWrapper^  OriYp, VectorWrapper^  OriX,
				double Mass);
			MultipleShearSpringWrapper(int Tag, int Nd1, int Nd2,
				array<UniaxialMaterialWrapper^>^ Materials,
				int NSpring,
				double LimDisp,
				VectorWrapper^  OriYp, VectorWrapper^  OriX,
				double Mass);

			MultipleShearSpringWrapper() {};
			~MultipleShearSpringWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};

		public ref class YamamotoBiaxialHDRWrapper : ElementWrapper
		{

		public:
			YamamotoBiaxialHDRWrapper(int Tag, int Nd1, int Nd2, int Tp, double DDo, double DDi, double Hr,
				double Cr, double Cs, VectorWrapper^  OriYp, VectorWrapper^  OriX,
				double Mass );

			YamamotoBiaxialHDRWrapper() {};
			~YamamotoBiaxialHDRWrapper() {
				if (_Element != 0)
					delete _Element;
			};

		internal:

		private:

		};



#pragma endregion

	}

}
