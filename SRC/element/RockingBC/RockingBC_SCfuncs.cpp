/* ****************************************************************** **
**    RockingBC Shear calc functions  			    				  **
** ****************************************************************** */

// Written by: Evangelos Avgenakis

#include "RockingBC_SCfuncs.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
// #include <Eigen/Dense>

void Dt_calc(const Vector& P, double& d, Vector& dddP)
{
	double p{ P(0) };
	double q{ P(1) };

	double a = SC_A(p);
	double b = SC_B(p);
	double c = SC_C(p);

	double da_dp = SC_DA(p);
	double db_dp = SC_DB(p);
	double dc_dp = SC_DC(p);

	d = a*pow(1.0 - pow(q, b), c);

	double dddp{ 0.0 };
	double dddq{ 0.0 };
	if (q > 0 && q < 1)
	{
		dddp = a*(std::log(1 - pow(q,b))*pow(1 - pow(q, b),c)*dc_dp - log(q)*c*pow(1 - pow(q, b),c-1)*db_dp*pow(q, b)) + pow(1 - pow(q, b),c)*da_dp;
		dddq = -a*b*c*pow(q,b-1)*pow(1 - pow(q, b),c-1);
	}
	else if (q == 1)
	{
		dddp = a*(-std::log(q)*c*pow(1 - pow(q, b),c-1)*db_dp*pow(q, b)) + pow(1 - pow(q, b),c)*da_dp;
		dddq = -a*b*c*pow(q, b-1)*pow(1 - pow(q, b),c-1);
	}
	else
	{
		throw;
	}
	dddP(0) = dddp;
	dddP(1) = dddq;

	return;

}

void Rt_calc(const Vector& P, double& th, Vector& dthdP)
{
	double p{ P(0) };
	double q{ P(1) };

	double d = SC_D(p);
	double e = SC_E(p);
	double f = SC_F(p);

	double dd_dp = SC_DD(p);
	double de_dp = SC_DE(p);
	double df_dp = SC_DF(p);

	th = d*pow(1 - q, e) + f;

	double dthdp{ 0.0 };
	double dthdq{ 0.0 };
	if (q > 0 && q < 1)
	{
		dthdp = pow(1 - q,e)*dd_dp + df_dp + std::log(1 - q)*d*pow(1 - q,e)*de_dp;
		dthdq = -d*e*pow(1 - q, e - 1);
	}
	else if (q == 1)
	{
		dthdp = df_dp;
		dthdq = -d*e*pow(1 - q, e - 1);
	}
	else
	{
		throw;
	}
	dthdP(0) = dthdp;
	dthdP(1) = dthdq;

	return;

}

void se_shear_1der(const Vector& Youter, Vector& Ut, Matrix& dUt_dYouter)
{
	double d{ 0.0 };
	double th{ 0.0 };
	static Vector dddP = Vector(2);
	static Vector dthdP = Vector(2);
	static Vector dp_dYouter = Vector(2);
	static Vector dq_dYouter = Vector(2);
	static Vector P(2);
	static Matrix dP_dYouter(2,2);

	if ((Youter(0) + Youter(1)) / 2.0 <= 0) {
		double q = 1. + (Youter(0) + Youter(1)) / 2.;
		double p = (1. + Youter(0)) / (1. + Youter(1));
		P(0) = p;
		P(1) = q;

		dp_dYouter(0) = 1. / (1. + Youter(1));
		dp_dYouter(1) = -(1. + Youter(0)) / (1. + Youter(1)) / (1. + Youter(1));
		dq_dYouter(0) = 1. / 2.;
		dq_dYouter(1) = 1. / 2.;
		dP_dYouter(0,0) = dp_dYouter(0);
		dP_dYouter(0,1) = dp_dYouter(1);
		dP_dYouter(1,0) = dq_dYouter(0);
		dP_dYouter(1,1) = dq_dYouter(1);

		Dt_calc(P,d,dddP);
		Rt_calc(P,th,dthdP);
	}
	else {

		double q = 1. - (Youter(0) + Youter(1)) / 2.;
		double p = (1. - Youter(1)) / (1. - Youter(0));
		P(0) = p;
		P(1) = q;

		dq_dYouter(0) = -1. / 2.;
		dq_dYouter(1) = -1. / 2.;
		dp_dYouter(0) = (1. - Youter(1)) / (1. - Youter(0)) / (1. - Youter(0));
		dp_dYouter(1) = -1. / (1. - Youter(0));
		dP_dYouter(0, 0) = dp_dYouter(0);
		dP_dYouter(0, 1) = dp_dYouter(1);
		dP_dYouter(1, 0) = dq_dYouter(0);
		dP_dYouter(1, 1) = dq_dYouter(1);

		Dt_calc(P, d, dddP);
		Rt_calc(P, th, dthdP);
		d = -d;
		dddP(0) = -dddP(0);
		dddP(1) = -dddP(1);
	}

	Ut(0) = d;
	Ut(1) = th;
	static Matrix dUt_dP(2, 2);
	dUt_dP(0,0) = dddP(0);
	dUt_dP(0,1) = dddP(1);
	dUt_dP(1,0) = dthdP(0);
	dUt_dP(1,1) = dthdP(1);
	dUt_dYouter = dUt_dP*dP_dYouter;

	return;

}

