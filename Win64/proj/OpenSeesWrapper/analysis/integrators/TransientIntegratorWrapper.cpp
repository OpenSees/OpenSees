#include "stdafx.h"
#include "TransientIntegratorWrapper.h"

using namespace OpenSees;
using namespace OpenSees::Integrators;
using namespace OpenSees::Integrators::Transient;

TransientIntegratorWrapper::TransientIntegratorWrapper()
{
	
}

TransientIntegratorWrapper::~TransientIntegratorWrapper()
{
	
}

NewmarkWrapper::NewmarkWrapper(double gamma, double beta, bool disp, bool aflag)
{
	_TransientIntegrator = new Newmark(gamma, beta, disp, aflag);
}

NewmarkWrapper::NewmarkWrapper(double gamma, double beta)
{
	_TransientIntegrator = new Newmark(gamma, beta);
}

NewmarkWrapper::~NewmarkWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}

AlphaOSWrapper::AlphaOSWrapper(double alpha, double beta, double gamma, bool updElemDisp)
{
	_TransientIntegrator = new AlphaOS(alpha, beta, gamma, updElemDisp);
}

AlphaOSWrapper::AlphaOSWrapper(double alpha, bool updElemDisp)
{
	_TransientIntegrator = new AlphaOS(alpha, updElemDisp);
}

AlphaOSWrapper::AlphaOSWrapper()
{
	_TransientIntegrator = new AlphaOS();
}

AlphaOSWrapper::~AlphaOSWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}

AlphaOS_TPWrapper::AlphaOS_TPWrapper(double alpha, double beta, double gamma, bool updElemDisp)
{
	_TransientIntegrator = new AlphaOS_TP(alpha, beta, gamma, updElemDisp);
}

AlphaOS_TPWrapper::AlphaOS_TPWrapper(double alpha, bool updElemDisp)
{
	_TransientIntegrator = new AlphaOS_TP(alpha, updElemDisp);
}

AlphaOS_TPWrapper::AlphaOS_TPWrapper()
{
	_TransientIntegrator = new AlphaOS_TP();
}

AlphaOS_TPWrapper::~AlphaOS_TPWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}

AlphaOSGeneralizedWrapper::AlphaOSGeneralizedWrapper(double alphaI, double alphaF, double beta, double gamma, bool updElemDisp)
{
	_TransientIntegrator = new AlphaOSGeneralized(alphaI, alphaF, beta, gamma, updElemDisp);
}

AlphaOSGeneralizedWrapper::AlphaOSGeneralizedWrapper(double rhoInf, bool updElemDisp)
{
	_TransientIntegrator = new AlphaOSGeneralized(rhoInf, updElemDisp);
}

AlphaOSGeneralizedWrapper::AlphaOSGeneralizedWrapper()
{
	_TransientIntegrator = new AlphaOSGeneralized();
}

AlphaOSGeneralizedWrapper::~AlphaOSGeneralizedWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}

AlphaOSGeneralized_TPWrapper::AlphaOSGeneralized_TPWrapper(double alphaI, double alphaF, double beta, double gamma, bool updElemDisp)
{
	_TransientIntegrator = new AlphaOSGeneralized_TP(alphaI, alphaF, beta, gamma, updElemDisp);
}

AlphaOSGeneralized_TPWrapper::AlphaOSGeneralized_TPWrapper(double rhoInf, bool updElemDisp)
{
	_TransientIntegrator = new AlphaOSGeneralized_TP(rhoInf, updElemDisp);
}

AlphaOSGeneralized_TPWrapper::AlphaOSGeneralized_TPWrapper()
{
	_TransientIntegrator = new AlphaOSGeneralized_TP();
}

AlphaOSGeneralized_TPWrapper::~AlphaOSGeneralized_TPWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}

BackwardEulerWrapper::BackwardEulerWrapper(int eulerOption)
{
	_TransientIntegrator = new BackwardEuler(eulerOption);
}

BackwardEulerWrapper::~BackwardEulerWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}

CentralDifferenceWrapper::CentralDifferenceWrapper()
{
	_TransientIntegrator = new CentralDifference();
}

CentralDifferenceWrapper::CentralDifferenceWrapper(double alphaM, double betaK, double betaKi, double betaKc)
{
	_TransientIntegrator = new CentralDifference(alphaM, betaK, betaKi, betaKc);
}

CentralDifferenceWrapper::~CentralDifferenceWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}

CentralDifferenceAlternativeWrapper::CentralDifferenceAlternativeWrapper()
{
	_TransientIntegrator = new CentralDifferenceAlternative();
}

CentralDifferenceAlternativeWrapper::~CentralDifferenceAlternativeWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}

CentralDifferenceNoDampingWrapper::CentralDifferenceNoDampingWrapper()
{
	_TransientIntegrator = new CentralDifferenceAlternative();
}

CentralDifferenceNoDampingWrapper::~CentralDifferenceNoDampingWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}


CollocationWrapper::CollocationWrapper()
{
	_TransientIntegrator = new Collocation();
}

CollocationWrapper::CollocationWrapper(double theta)
{
	_TransientIntegrator = new Collocation(theta);
}

CollocationWrapper::CollocationWrapper(double theta, double beta, double gamma)
{
	_TransientIntegrator = new Collocation(theta, beta, gamma);
}

CollocationWrapper::~CollocationWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}

CollocationHSFixedNumIterWrapper::CollocationHSFixedNumIterWrapper()
{
	_TransientIntegrator = new CollocationHSFixedNumIter();
}

CollocationHSFixedNumIterWrapper::CollocationHSFixedNumIterWrapper(double theta, int polyOrder)
{
	_TransientIntegrator = new CollocationHSFixedNumIter(theta, polyOrder);
}

CollocationHSFixedNumIterWrapper::CollocationHSFixedNumIterWrapper(double theta, double beta, double gamma, int polyOrder)
{
	_TransientIntegrator = new CollocationHSFixedNumIter(theta, beta, gamma, polyOrder);
}

CollocationHSFixedNumIterWrapper::~CollocationHSFixedNumIterWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}

CollocationHSIncrLimitWrapper::CollocationHSIncrLimitWrapper()
{
	_TransientIntegrator = new CollocationHSIncrLimit();
}

CollocationHSIncrLimitWrapper::CollocationHSIncrLimitWrapper(double theta, double limit, int normType)
{
	_TransientIntegrator = new CollocationHSIncrLimit(theta, limit, normType);
}

CollocationHSIncrLimitWrapper::CollocationHSIncrLimitWrapper(double theta, double beta, double gamma, double limit, int normType)
{
	_TransientIntegrator = new CollocationHSIncrLimit(theta, beta, gamma, limit, normType);
}

CollocationHSIncrLimitWrapper::~CollocationHSIncrLimitWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}

CollocationHSIncrReductWrapper::CollocationHSIncrReductWrapper()
{
	_TransientIntegrator = new CollocationHSIncrReduct();
}

CollocationHSIncrReductWrapper::CollocationHSIncrReductWrapper(double theta, double reduct)
{
	_TransientIntegrator = new CollocationHSIncrReduct(theta, reduct);
}

CollocationHSIncrReductWrapper::CollocationHSIncrReductWrapper(double theta, double beta, double gamma, double reduct)
{
	_TransientIntegrator = new CollocationHSIncrReduct(theta, beta, gamma, reduct);
}

CollocationHSIncrReductWrapper::~CollocationHSIncrReductWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}

GeneralizedAlphaWrapper::GeneralizedAlphaWrapper()
{
	_TransientIntegrator = new GeneralizedAlpha();
}

GeneralizedAlphaWrapper::GeneralizedAlphaWrapper(double alphaM, double alphaF)
{
	_TransientIntegrator = new GeneralizedAlpha(alphaM, alphaF);
}

GeneralizedAlphaWrapper::GeneralizedAlphaWrapper(double alphaM, double alphaF, double beta, double gamma)
{
	_TransientIntegrator = new GeneralizedAlpha(alphaM, alphaF, beta, gamma);
}

GeneralizedAlphaWrapper::~GeneralizedAlphaWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}

HHTWrapper::HHTWrapper()
{
	_TransientIntegrator = new HHT();
}

HHTWrapper::HHTWrapper(double alpha)
{
	_TransientIntegrator = new HHT(alpha);
}

HHTWrapper::HHTWrapper(double alpha, double beta, double gamma)
{
	_TransientIntegrator = new HHT(alpha, beta, gamma);
}

HHTWrapper::~HHTWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}


HHT_TPWrapper::HHT_TPWrapper()
{
	_TransientIntegrator = new HHT_TP();
}

HHT_TPWrapper::HHT_TPWrapper(double alpha)
{
	_TransientIntegrator = new HHT_TP(alpha);
}

HHT_TPWrapper::HHT_TPWrapper(double alpha, double beta, double gamma)
{
	_TransientIntegrator = new HHT_TP(alpha, beta, gamma);
}

HHT_TPWrapper::~HHT_TPWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}


HHTExplicitWrapper::HHTExplicitWrapper()
{
	_TransientIntegrator = new HHTExplicit();
}

HHTExplicitWrapper::HHTExplicitWrapper(double alpha, bool updElemDisp)
{
	_TransientIntegrator = new HHTExplicit(alpha, updElemDisp);
}

HHTExplicitWrapper::HHTExplicitWrapper(double alpha, double gamma, bool updElemDisp)
{
	_TransientIntegrator = new HHTExplicit(alpha, gamma, updElemDisp);
}

HHTExplicitWrapper::~HHTExplicitWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}


HHTExplicit_TPWrapper::HHTExplicit_TPWrapper()
{
	_TransientIntegrator = new HHTExplicit_TP();
}

HHTExplicit_TPWrapper::HHTExplicit_TPWrapper(double alpha)
{
	_TransientIntegrator = new HHTExplicit_TP(alpha);
}

HHTExplicit_TPWrapper::HHTExplicit_TPWrapper(double alpha, double gamma)
{
	_TransientIntegrator = new HHTExplicit_TP(alpha, gamma);
}

HHTExplicit_TPWrapper::~HHTExplicit_TPWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}


HHTGeneralizedWrapper::HHTGeneralizedWrapper()
{
	_TransientIntegrator = new HHTGeneralized();
}

HHTGeneralizedWrapper::HHTGeneralizedWrapper(double rhoInf)
{
	_TransientIntegrator = new HHTGeneralized(rhoInf);
}

HHTGeneralizedWrapper::HHTGeneralizedWrapper(double alphaI, double alphaF, double beta, double gamma)
{
	_TransientIntegrator = new HHTGeneralized(alphaI, alphaF, beta, gamma);
}

HHTGeneralizedWrapper::~HHTGeneralizedWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}


HHTGeneralized_TPWrapper::HHTGeneralized_TPWrapper()
{
	_TransientIntegrator = new HHTGeneralized_TP();
}

HHTGeneralized_TPWrapper::HHTGeneralized_TPWrapper(double rhoInf)
{
	_TransientIntegrator = new HHTGeneralized_TP(rhoInf);
}

HHTGeneralized_TPWrapper::HHTGeneralized_TPWrapper(double alphaI, double alphaF, double beta, double gamma)
{
	_TransientIntegrator = new HHTGeneralized_TP(alphaI, alphaF, beta, gamma);
}

HHTGeneralized_TPWrapper::~HHTGeneralized_TPWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}


HHTGeneralizedExplicitWrapper::HHTGeneralizedExplicitWrapper()
{
	_TransientIntegrator = new HHTGeneralizedExplicit();
}

HHTGeneralizedExplicitWrapper::HHTGeneralizedExplicitWrapper(double rhoB, double alphaF, bool updElemDisp)
{
	_TransientIntegrator = new HHTGeneralizedExplicit(rhoB, alphaF, updElemDisp);
}

HHTGeneralizedExplicitWrapper::HHTGeneralizedExplicitWrapper(double alphaI, double alphaF, double beta, double gamma, bool updElemDisp)
{
	_TransientIntegrator = new HHTGeneralizedExplicit(alphaI, alphaF, beta, gamma, updElemDisp);
}

HHTGeneralizedExplicitWrapper::~HHTGeneralizedExplicitWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}


HHTGeneralizedExplicit_TPWrapper::HHTGeneralizedExplicit_TPWrapper()
{
	_TransientIntegrator = new HHTGeneralizedExplicit_TP();
}

HHTGeneralizedExplicit_TPWrapper::HHTGeneralizedExplicit_TPWrapper(double rhoB, double alphaF)
{
	_TransientIntegrator = new HHTGeneralizedExplicit_TP(rhoB, alphaF);
}

HHTGeneralizedExplicit_TPWrapper::HHTGeneralizedExplicit_TPWrapper(double alphaI, double alphaF, double beta, double gamma)
{
	_TransientIntegrator = new HHTGeneralizedExplicit_TP(alphaI, alphaF, beta, gamma);
}

HHTGeneralizedExplicit_TPWrapper::~HHTGeneralizedExplicit_TPWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}


HHTHSFixedNumIterWrapper::HHTHSFixedNumIterWrapper()
{
	_TransientIntegrator = new HHTHSFixedNumIter();
}

HHTHSFixedNumIterWrapper::HHTHSFixedNumIterWrapper(double rhoInf, int polyOrder, bool updDomFlag)
{
	_TransientIntegrator = new HHTHSFixedNumIter(rhoInf, polyOrder, updDomFlag);
}

HHTHSFixedNumIterWrapper::HHTHSFixedNumIterWrapper(double alphaI, double alphaF, double beta, double gamma, int polyOrder, bool updDomFlag)
{
	_TransientIntegrator = new HHTHSFixedNumIter(alphaI, alphaF, beta, gamma, polyOrder, updDomFlag);
}

HHTHSFixedNumIterWrapper::~HHTHSFixedNumIterWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}


HHTHSFixedNumIter_TPWrapper::HHTHSFixedNumIter_TPWrapper()
{
	_TransientIntegrator = new HHTHSFixedNumIter_TP();
}

HHTHSFixedNumIter_TPWrapper::HHTHSFixedNumIter_TPWrapper(double rhoInf, int polyOrder, bool updDomFlag)
{
	_TransientIntegrator = new HHTHSFixedNumIter_TP(rhoInf, polyOrder, updDomFlag);
}

HHTHSFixedNumIter_TPWrapper::HHTHSFixedNumIter_TPWrapper(double alphaI, double alphaF, double beta, double gamma, int polyOrder, bool updDomFlag)
{
	_TransientIntegrator = new HHTHSFixedNumIter_TP(alphaI, alphaF, beta, gamma, polyOrder, updDomFlag);
}

HHTHSFixedNumIter_TPWrapper::~HHTHSFixedNumIter_TPWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}


HHTHSIncrLimitWrapper::HHTHSIncrLimitWrapper()
{
	_TransientIntegrator = new HHTHSIncrLimit();
}

HHTHSIncrLimitWrapper::HHTHSIncrLimitWrapper(double rhoInf, double limit, int normType)
{
	_TransientIntegrator = new HHTHSIncrLimit(rhoInf, limit, normType);
}

HHTHSIncrLimitWrapper::HHTHSIncrLimitWrapper(double alphaI, double alphaF, double beta, double gamma, double limit, int normType)
{
	_TransientIntegrator = new HHTHSIncrLimit(alphaI, alphaF, beta, gamma, limit, normType);
}

HHTHSIncrLimitWrapper::~HHTHSIncrLimitWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}


HHTHSIncrLimit_TPWrapper::HHTHSIncrLimit_TPWrapper()
{
	_TransientIntegrator = new HHTHSIncrLimit_TP();
}

HHTHSIncrLimit_TPWrapper::HHTHSIncrLimit_TPWrapper(double rhoInf, double limit, int normType)
{
	_TransientIntegrator = new HHTHSIncrLimit_TP(rhoInf, limit, normType);
}

HHTHSIncrLimit_TPWrapper::HHTHSIncrLimit_TPWrapper(double alphaI, double alphaF, double beta, double gamma, double limit, int normType)
{
	_TransientIntegrator = new HHTHSIncrLimit_TP(alphaI, alphaF, beta, gamma, limit, normType);
}

HHTHSIncrLimit_TPWrapper::~HHTHSIncrLimit_TPWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}


HHTHSIncrReductWrapper::HHTHSIncrReductWrapper()
{
	_TransientIntegrator = new HHTHSIncrReduct();
}

HHTHSIncrReductWrapper::HHTHSIncrReductWrapper(double rhoInf, double reduct)
{
	_TransientIntegrator = new HHTHSIncrReduct(rhoInf, reduct);
}

HHTHSIncrReductWrapper::HHTHSIncrReductWrapper(double alphaI, double alphaF, double beta, double gamma, double reduct)
{
	_TransientIntegrator = new HHTHSIncrReduct(alphaI, alphaF, beta, gamma, reduct);
}

HHTHSIncrReductWrapper::~HHTHSIncrReductWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}


HouboltWrapper::HouboltWrapper()
{
	_TransientIntegrator = new Houbolt();
}

HouboltWrapper::~HouboltWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}


KRAlphaExplicitWrapper::KRAlphaExplicitWrapper()
{
	_TransientIntegrator = new KRAlphaExplicit();
}

KRAlphaExplicitWrapper::KRAlphaExplicitWrapper(double rhoInf, bool updElemDisp)
{
	_TransientIntegrator = new KRAlphaExplicit(rhoInf, updElemDisp);
}

KRAlphaExplicitWrapper::~KRAlphaExplicitWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}


KRAlphaExplicit_TPWrapper::KRAlphaExplicit_TPWrapper()
{
	_TransientIntegrator = new KRAlphaExplicit_TP();
}

KRAlphaExplicit_TPWrapper::KRAlphaExplicit_TPWrapper(double rhoInf)
{
	_TransientIntegrator = new KRAlphaExplicit_TP(rhoInf);
}

KRAlphaExplicit_TPWrapper::~KRAlphaExplicit_TPWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}


Newmark1Wrapper::Newmark1Wrapper()
{
	_TransientIntegrator = new Newmark1();
}

Newmark1Wrapper::Newmark1Wrapper(double gamma, double beta, bool disp)
{
	_TransientIntegrator = new Newmark1(gamma, beta, disp);
}

Newmark1Wrapper::Newmark1Wrapper(double gamma, double beta, double alphaM, double betaK, double betaKi, double betaKc)
{
	_TransientIntegrator = new Newmark1(gamma, beta, alphaM, betaK, betaKi, betaKc);
}

Newmark1Wrapper::~Newmark1Wrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}


NewmarkExplicitWrapper::NewmarkExplicitWrapper()
{
	_TransientIntegrator = new NewmarkExplicit();
}

NewmarkExplicitWrapper::NewmarkExplicitWrapper(double gamma)
{
	_TransientIntegrator = new NewmarkExplicit(gamma);
}

NewmarkExplicitWrapper::~NewmarkExplicitWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}


NewmarkHSFixedNumIterWrapper::NewmarkHSFixedNumIterWrapper()
{
	_TransientIntegrator = new NewmarkHSFixedNumIter();
}

NewmarkHSFixedNumIterWrapper::NewmarkHSFixedNumIterWrapper(double gamma, double beta, int polyOrder, bool updDomFlag)
{
	_TransientIntegrator = new NewmarkHSFixedNumIter(gamma,beta, polyOrder, updDomFlag);
}

NewmarkHSFixedNumIterWrapper::~NewmarkHSFixedNumIterWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}


NewmarkHSIncrLimitWrapper::NewmarkHSIncrLimitWrapper()
{
	_TransientIntegrator = new NewmarkHSIncrLimit();
}

NewmarkHSIncrLimitWrapper::NewmarkHSIncrLimitWrapper(double gamma, double beta, double limit, int normType)
{
	_TransientIntegrator = new NewmarkHSIncrLimit(gamma, beta, limit, normType);
}

NewmarkHSIncrLimitWrapper::~NewmarkHSIncrLimitWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}


NewmarkHSIncrReductWrapper::NewmarkHSIncrReductWrapper()
{
	_TransientIntegrator = new NewmarkHSIncrReduct();
}

NewmarkHSIncrReductWrapper::NewmarkHSIncrReductWrapper(double gamma, double beta, double reduct)
{
	_TransientIntegrator = new NewmarkHSIncrReduct(gamma, beta, reduct);
}

NewmarkHSIncrReductWrapper::~NewmarkHSIncrReductWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}


ParkLMS3Wrapper::ParkLMS3Wrapper()
{
	_TransientIntegrator = new ParkLMS3();
}

ParkLMS3Wrapper::~ParkLMS3Wrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}


PFEMIntegratorWrapper::PFEMIntegratorWrapper()
{
	_TransientIntegrator = new PFEMIntegrator();
}

PFEMIntegratorWrapper::~PFEMIntegratorWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}


TRBDF2Wrapper::TRBDF2Wrapper()
{
	_TransientIntegrator = new TRBDF2();
}

TRBDF2Wrapper::~TRBDF2Wrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}


TRBDF3Wrapper::TRBDF3Wrapper()
{
	_TransientIntegrator = new TRBDF3();
}

TRBDF3Wrapper::~TRBDF3Wrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}


WilsonThetaWrapper::WilsonThetaWrapper()
{
	_TransientIntegrator = new TRBDF3();
}

WilsonThetaWrapper::WilsonThetaWrapper(double theta)
{
	_TransientIntegrator = new WilsonTheta(theta);
}

WilsonThetaWrapper::~WilsonThetaWrapper()
{
	if (_TransientIntegrator != 0)
		delete _TransientIntegrator;
}