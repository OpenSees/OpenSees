#pragma once
#include <TransientIntegrator.h>
#include <Newmark.h>
#include <AlphaOS.h>
#include <AlphaOS_TP.h>
#include <AlphaOSGeneralized.h>
#include <AlphaOSGeneralized_TP.h>
#include <BackwardEuler.h>
#include <CentralDifference.h>
#include <CentralDifferenceAlternative.h>
#include <CentralDifferenceNoDamping.h>
#include <Collocation.h>
#include <CollocationHSFixedNumIter.h>
#include <CollocationHSIncrLimit.h>
#include <CollocationHSIncrReduct.h>
#include <GeneralizedAlpha.h>
#include <HHT.h>
#include <HHT_TP.h>
#include <HHTExplicit.h>
#include <HHTExplicit_TP.h>
#include <HHTGeneralized.h>
#include <HHTGeneralized_TP.h>
#include <HHTGeneralizedExplicit.h>
#include <HHTGeneralizedExplicit_TP.h>
#include <HHTHSFixedNumIter.h>
#include <HHTHSFixedNumIter_TP.h>
#include <HHTHSIncrLimit.h>
#include <HHTHSIncrLimit_TP.h>
#include <HHTHSIncrReduct.h>
#include <HHTHSIncrReduct_TP.h>
#include <Houbolt.h>
#include <KRAlphaExplicit.h>
#include <KRAlphaExplicit_TP.h>
#include <Newmark1.h>
#include <NewmarkExplicit.h>
#include <NewmarkHSFixedNumIter.h>
#include <NewmarkHSIncrLimit.h>
#include <NewmarkHSIncrReduct.h>
#include <ParkLMS3.h>
#include <PFEMIntegrator.h>
#include <TRBDF2.h>
#include <TRBDF3.h>
#include <WilsonTheta.h>


#include "IncrementalIntegratorWrapper.h"


namespace OpenSees {
	namespace Integrators {
		namespace Transient {

			public ref class TransientIntegratorWrapper abstract : IncrementalIntegratorWrapper
			{
			public:
				TransientIntegratorWrapper();
				~TransientIntegratorWrapper();
			internal:
				TransientIntegrator * _TransientIntegrator;
			private:

			};

			public ref class NewmarkWrapper : TransientIntegratorWrapper
			{
			public:
				NewmarkWrapper(double gamma, double beta, bool disp, bool aflag);
				NewmarkWrapper(double gamma, double beta);
				~NewmarkWrapper();
			internal:

			private:

			};

			public ref class AlphaOSWrapper : TransientIntegratorWrapper
			{
			public:
				AlphaOSWrapper();
				AlphaOSWrapper(double alpha, double beta, double gamma, bool updElemDisp);
				AlphaOSWrapper(double alpha, bool updElemDisp);
				~AlphaOSWrapper();
			internal:

			private:

			};

			public ref class AlphaOS_TPWrapper : TransientIntegratorWrapper
			{
			public:
				AlphaOS_TPWrapper();
				AlphaOS_TPWrapper(double alpha, double beta, double gamma, bool updElemDisp);
				AlphaOS_TPWrapper(double alpha, bool updElemDisp);
				~AlphaOS_TPWrapper();
			internal:

			private:

			};

			public ref class AlphaOSGeneralizedWrapper : TransientIntegratorWrapper
			{
			public:
				AlphaOSGeneralizedWrapper();
				AlphaOSGeneralizedWrapper(double rhoInf, bool updElemDisp);
				AlphaOSGeneralizedWrapper(double alphaI, double alphaF, double beta, double gamma, bool updElemDisp);
				~AlphaOSGeneralizedWrapper();
			internal:

			private:

			};

			public ref class AlphaOSGeneralized_TPWrapper : TransientIntegratorWrapper
			{
			public:
				AlphaOSGeneralized_TPWrapper();
				AlphaOSGeneralized_TPWrapper(double rhoInf, bool updElemDisp);
				AlphaOSGeneralized_TPWrapper(double alphaI, double alphaF, double beta, double gamma, bool updElemDisp);
				~AlphaOSGeneralized_TPWrapper();
			internal:

			private:

			};

			public ref class BackwardEulerWrapper : TransientIntegratorWrapper
			{
			public:
				BackwardEulerWrapper(int eulerOption);
				~BackwardEulerWrapper();
			internal:

			private:

			};

			public ref class CentralDifferenceWrapper : TransientIntegratorWrapper
			{
			public:
				CentralDifferenceWrapper();
				CentralDifferenceWrapper(double alphaM, double betaK, double betaKi, double betaKc);
				~CentralDifferenceWrapper();
			internal:

			private:

			};

			public ref class CentralDifferenceAlternativeWrapper : TransientIntegratorWrapper
			{
			public:
				CentralDifferenceAlternativeWrapper();
				~CentralDifferenceAlternativeWrapper();
			internal:

			private:

			};

			public ref class CentralDifferenceNoDampingWrapper : TransientIntegratorWrapper
			{
			public:
				CentralDifferenceNoDampingWrapper();
				~CentralDifferenceNoDampingWrapper();
			internal:

			private:

			};

			public ref class CollocationWrapper : TransientIntegratorWrapper
			{
			public:
				CollocationWrapper();
				CollocationWrapper(double theta);
				CollocationWrapper(double theta, double beta, double gamma);
				~CollocationWrapper();
			internal:

			private:

			};

			public ref class CollocationHSFixedNumIterWrapper : TransientIntegratorWrapper
			{
			public:
				CollocationHSFixedNumIterWrapper();
				CollocationHSFixedNumIterWrapper(double theta, int polyOrder);
				CollocationHSFixedNumIterWrapper(double theta, double beta, double gamma, int polyOrder);
				~CollocationHSFixedNumIterWrapper();
			internal:

			private:

			};

			public ref class CollocationHSIncrLimitWrapper : TransientIntegratorWrapper
			{
			public:
				CollocationHSIncrLimitWrapper();
				CollocationHSIncrLimitWrapper(double theta, double limit, int normType);
				CollocationHSIncrLimitWrapper(double theta, double beta, double gamma, double limit, int normType);
				~CollocationHSIncrLimitWrapper();
			internal:

			private:

			};

			public ref class CollocationHSIncrReductWrapper : TransientIntegratorWrapper
			{
			public:
				CollocationHSIncrReductWrapper();
				CollocationHSIncrReductWrapper(double theta, double reduct);
				CollocationHSIncrReductWrapper(double theta, double beta, double gamma, double reduct);
				~CollocationHSIncrReductWrapper();
			internal:

			private:

			};

			public ref class GeneralizedAlphaWrapper : TransientIntegratorWrapper
			{
			public:
				GeneralizedAlphaWrapper();
				GeneralizedAlphaWrapper(double alphaM, double alphaF);
				GeneralizedAlphaWrapper(double alphaM, double alphaF, double beta, double gamma);
				~GeneralizedAlphaWrapper();
			internal:

			private:

			};

			public ref class HHTWrapper : TransientIntegratorWrapper
			{
			public:
				HHTWrapper();
				HHTWrapper(double alpha);
				HHTWrapper(double alpha, double beta, double gamma);
				~HHTWrapper();
			internal:

			private:

			};

			public ref class HHT_TPWrapper : TransientIntegratorWrapper
			{
			public:
				HHT_TPWrapper();
				HHT_TPWrapper(double alpha);
				HHT_TPWrapper(double alpha, double beta, double gamma);
				~HHT_TPWrapper();
			internal:

			private:

			};

			public ref class HHTExplicitWrapper : TransientIntegratorWrapper
			{
			public:
				HHTExplicitWrapper();
				HHTExplicitWrapper(double alpha, bool updElemDisp);
				HHTExplicitWrapper(double alpha, double gamma, bool updElemDisp);
				~HHTExplicitWrapper();
			internal:

			private:

			};

			public ref class HHTExplicit_TPWrapper : TransientIntegratorWrapper
			{
			public:
				HHTExplicit_TPWrapper();
				HHTExplicit_TPWrapper(double alpha);
				HHTExplicit_TPWrapper(double alpha, double gamma);
				~HHTExplicit_TPWrapper();
			internal:

			private:

			};

			public ref class HHTGeneralizedWrapper : TransientIntegratorWrapper
			{
			public:
				HHTGeneralizedWrapper();
				HHTGeneralizedWrapper(double rhoInf);
				HHTGeneralizedWrapper(double alphaI, double alphaF, double beta, double gamma);
				~HHTGeneralizedWrapper();
			internal:

			private:

			};

			public ref class HHTGeneralized_TPWrapper : TransientIntegratorWrapper
			{
			public:
				HHTGeneralized_TPWrapper();
				HHTGeneralized_TPWrapper(double rhoInf);
				HHTGeneralized_TPWrapper(double alphaI, double alphaF, double beta, double gamma);
				~HHTGeneralized_TPWrapper();
			internal:

			private:

			};

			public ref class HHTGeneralizedExplicitWrapper : TransientIntegratorWrapper
			{
			public:
				HHTGeneralizedExplicitWrapper();
				HHTGeneralizedExplicitWrapper(double rhoB, double alphaF, bool updElemDisp);
				HHTGeneralizedExplicitWrapper(double alphaI, double alphaF, double beta, double gamma, bool updElemDisp);
				~HHTGeneralizedExplicitWrapper();
			internal:

			private:

			};

			public ref class HHTGeneralizedExplicit_TPWrapper : TransientIntegratorWrapper
			{
			public:
				HHTGeneralizedExplicit_TPWrapper();
				HHTGeneralizedExplicit_TPWrapper(double rhoB, double alphaF);
				HHTGeneralizedExplicit_TPWrapper(double alphaI, double alphaF, double beta, double gamma);
				~HHTGeneralizedExplicit_TPWrapper();
			internal:

			private:

			};

			public ref class HHTHSFixedNumIterWrapper : TransientIntegratorWrapper
			{
			public:
				HHTHSFixedNumIterWrapper();
				HHTHSFixedNumIterWrapper(double rhoInf, int polyOrder, bool updDomFlag);
				HHTHSFixedNumIterWrapper(double alphaI, double alphaF, double beta, double gamma, int polyOrder, bool updDomFlag);
				~HHTHSFixedNumIterWrapper();
			internal:

			private:

			};

			public ref class HHTHSFixedNumIter_TPWrapper : TransientIntegratorWrapper
			{
			public:
				HHTHSFixedNumIter_TPWrapper();
				HHTHSFixedNumIter_TPWrapper(double rhoInf, int polyOrder, bool updDomFlag);
				HHTHSFixedNumIter_TPWrapper(double alphaI, double alphaF, double beta, double gamma, int polyOrder, bool updDomFlag);
				~HHTHSFixedNumIter_TPWrapper();
			internal:

			private:

			};

			public ref class HHTHSIncrLimitWrapper : TransientIntegratorWrapper
			{
			public:
				HHTHSIncrLimitWrapper();
				HHTHSIncrLimitWrapper(double rhoInf, double limit, int normType);
				HHTHSIncrLimitWrapper(double alphaI, double alphaF, double beta, double gamma, double limit, int normType);
				~HHTHSIncrLimitWrapper();
			internal:

			private:

			};

			public ref class HHTHSIncrLimit_TPWrapper : TransientIntegratorWrapper
			{
			public:
				HHTHSIncrLimit_TPWrapper();
				HHTHSIncrLimit_TPWrapper(double rhoInf, double limit, int normType);
				HHTHSIncrLimit_TPWrapper(double alphaI, double alphaF, double beta, double gamma, double limit, int normType);
				~HHTHSIncrLimit_TPWrapper();
			internal:

			private:

			};

			public ref class HHTHSIncrReductWrapper : TransientIntegratorWrapper
			{
			public:
				HHTHSIncrReductWrapper();
				HHTHSIncrReductWrapper(double rhoInf, double reduct);
				HHTHSIncrReductWrapper(double alphaI, double alphaF, double beta, double gamma, double reduct);
				~HHTHSIncrReductWrapper();
			internal:

			private:

			};

			public ref class HouboltWrapper : TransientIntegratorWrapper
			{
			public:
				HouboltWrapper();
				~HouboltWrapper();
			internal:

			private:

			};

			public ref class KRAlphaExplicitWrapper : TransientIntegratorWrapper
			{
			public:
				KRAlphaExplicitWrapper();
				KRAlphaExplicitWrapper(double rhoInf, bool updElemDisp);
				~KRAlphaExplicitWrapper();
			internal:

			private:

			};

			public ref class KRAlphaExplicit_TPWrapper : TransientIntegratorWrapper
			{
			public:
				KRAlphaExplicit_TPWrapper();
				KRAlphaExplicit_TPWrapper(double rhoInf);
				~KRAlphaExplicit_TPWrapper();
			internal:

			private:

			};

			public ref class Newmark1Wrapper : TransientIntegratorWrapper
			{
			public:
				Newmark1Wrapper();
				Newmark1Wrapper(double gamma, double beta, bool disp);
				Newmark1Wrapper(double gamma, double beta, double alphaM, double betaK, double betaKi, double betaKc);
				~Newmark1Wrapper();
			internal:

			private:

			};

			public ref class NewmarkExplicitWrapper : TransientIntegratorWrapper
			{
			public:
				NewmarkExplicitWrapper();
				NewmarkExplicitWrapper(double gamma);
				~NewmarkExplicitWrapper();
			internal:

			private:

			};

			public ref class NewmarkHSFixedNumIterWrapper : TransientIntegratorWrapper
			{
			public:
				NewmarkHSFixedNumIterWrapper();
				NewmarkHSFixedNumIterWrapper(double gamma, double beta, int polyOrder, bool updDomFlag);
				~NewmarkHSFixedNumIterWrapper();
			internal:

			private:

			};

			public ref class NewmarkHSIncrLimitWrapper : TransientIntegratorWrapper
			{
			public:
				NewmarkHSIncrLimitWrapper();
				NewmarkHSIncrLimitWrapper(double gamma, double beta, double limit, int normType);
				~NewmarkHSIncrLimitWrapper();
			internal:

			private:

			};

			public ref class NewmarkHSIncrReductWrapper : TransientIntegratorWrapper
			{
			public:
				NewmarkHSIncrReductWrapper();
				NewmarkHSIncrReductWrapper(double gamma, double beta, double reduct);
				~NewmarkHSIncrReductWrapper();
			internal:

			private:

			};

			public ref class ParkLMS3Wrapper : TransientIntegratorWrapper
			{
			public:
				ParkLMS3Wrapper();
				~ParkLMS3Wrapper();
			internal:

			private:

			};

			public ref class PFEMIntegratorWrapper : TransientIntegratorWrapper
			{
			public:
				PFEMIntegratorWrapper();
				~PFEMIntegratorWrapper();
			internal:

			private:

			};

			public ref class TRBDF2Wrapper : TransientIntegratorWrapper
			{
			public:
				TRBDF2Wrapper();
				~TRBDF2Wrapper();
			internal:

			private:

			};

			public ref class TRBDF3Wrapper : TransientIntegratorWrapper
			{
			public:
				TRBDF3Wrapper();
				~TRBDF3Wrapper();
			internal:

			private:

			};

			public ref class WilsonThetaWrapper : TransientIntegratorWrapper
			{
			public:
				WilsonThetaWrapper();
				WilsonThetaWrapper(double theta);
				~WilsonThetaWrapper();
			internal:

			private:

			};
		}
	}
}
