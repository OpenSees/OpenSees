#include "stdafx.h"
#include "ULBeamColumnWrapper.h"

using namespace System::Runtime::InteropServices;
using namespace OpenSees;
using namespace OpenSees::Elements;
using namespace OpenSees::Materials::YieldSurfaces;


Elastic2dGNLWrapper::Elastic2dGNLWrapper(int tag, double A, double E, double I, int Nd1, int Nd2,
	bool islinear, double rho)
{
	_Element = new Elastic2dGNL(tag, A, E, I, Nd1, Nd2,
		islinear, rho);
}

Inelastic2DYS01Wrapper::Inelastic2DYS01Wrapper(int tag, double A, double E, double Iz,
	int Nd1, int Nd2,
	YieldSurface_BCWrapper^ ysEnd1, YieldSurface_BCWrapper^ ysEnd2,
	int rf_algo, bool islinear, double rho)
{
	_Element = new Inelastic2DYS01(tag, A, E, Iz,
		Nd1, Nd2,
		ysEnd1->_YieldSurface_BC, ysEnd2->_YieldSurface_BC,
		rf_algo, islinear, rho);
}

Inelastic2DYS02Wrapper::Inelastic2DYS02Wrapper(int tag, double A, double E, double Iz,
	int Nd1, int Nd2,
	YieldSurface_BCWrapper^ ysEnd1, YieldSurface_BCWrapper^ ysEnd2,
	CyclicModelWrapper^ cycModel,
	double del_p_max,
	double Alpha, double Beta, int rf_algo,
	bool islinear, double rho)
{
	_Element = new Inelastic2DYS02(tag, A, E, Iz,
		Nd1, Nd2,
		ysEnd1->_YieldSurface_BC, ysEnd2->_YieldSurface_BC,
		cycModel->_CyclicModel,
		del_p_max,
		Alpha, Beta, rf_algo,
		islinear, rho);
}

Inelastic2DYS03Wrapper::Inelastic2DYS03Wrapper(int tag, double a_ten, double a_com, double e,
	double iz_pos, double iz_neg, int Nd1, int Nd2,
	YieldSurface_BCWrapper^ ysEnd1, YieldSurface_BCWrapper^ ysEnd2,
	int rf_algo, bool islinear, double rho)
{
	_Element = new Inelastic2DYS03(tag, a_ten, a_com, e,
		iz_pos, iz_neg, Nd1, Nd2,
		ysEnd1->_YieldSurface_BC, ysEnd2->_YieldSurface_BC,
		rf_algo, islinear, rho);
}




