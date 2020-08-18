/* ****************************************************************** **
**    RockingBC Semi-infinite strip functions						  **
** ****************************************************************** */

// Written by: Evangelos Avgenakis

#include "RockingBC_SISfuncs.h"

const double SISfunclim{ 1.0e-15 };

double YMXLOGYMX(double y, double p)
{
	double ymx{ p - y };
	if (std::fabs(ymx)<SISfunclim)
	{
		return 0.0;
	}
	else
	{
		return ymx*std::log(std::fabs(ymx));
	}
}

double OMXYLOGOMXYOXY(double yp)
{
	if (std::fabs(yp)<SISfunclim)
	{
		return -1.0;
	}
	else if (std::fabs(yp - 1) < SISfunclim)
	{
		return 0.0;
	}
	else
	{
		return (1.0 - yp)*std::log1p(-yp) / (yp);
	}
}

double J2(double yp)
{
	if (std::fabs(yp)<SISfunclim)
	{
		return 0.5;
	}
	else if (std::fabs(yp - 1) < SISfunclim)
	{
		return 1.0;
	}
	else
	{
		return (OMXYLOGOMXYOXY(yp) + 1.0) / (yp);
	}
}

double OMXATANYMOOXMO(double y, double p)
{
	if (std::fabs(y-1)<SISfunclim)
	{
		return 0.0;
	}
	else
	{ 
		return (1.0 - y)*std::atan((p - 1.0) / (y - 1.0));
	}
}

double OMYLOGSQ(double y, double p)
{
	if (std::fabs(p - 1)<SISfunclim)
	{
		return 0.0;
	}
	else
	{
		return (1.0 - p)*std::log((y - 1.0)*(y - 1.0) + (p - 1.0)*(p - 1.0));
	}
}

double I_FA(double y, double p)
{
	double FA1 = 2 * YMXLOGYMX(y, p);
	double FA2 = -OMXYLOGOMXYOXY(y*p)*p / 3. * (2 * y*y * p*p + 5 * y*p - 1);
	double FA3 = OMXYLOGOMXYOXY(-y*p)*p / 3. * (1 + y*p)*(2 * y*p - 1);
	double FA4 = 4. / 3. * y*p*p;

	return FA1 + FA2 + FA3 + FA4;
}

double J_FA(double y, double p)
{
	double FA1 = +YMXLOGYMX(y, p)*(p + y);
	double FA2 = -p*p / 6. * (OMXYLOGOMXYOXY(y*p) + YMXLOGYMX(y*p, 1.)*(3. * y*p + 7.) + J2(y*p));
	double FA3 = -p*p / 6. * (OMXYLOGOMXYOXY(-y*p) + YMXLOGYMX(-y*p, 1.)*(3. * y*p + 1.) + J2(-y*p));
	double FA4 = y*p*p*p + p*p / 3. - p*y;

	return FA1 + FA2 + FA3 + FA4;
}

double pImJ_FA(double y, double p)
{
	double FA1 = (p - y)*YMXLOGYMX(y, p);
	double FA2 = -OMXYLOGOMXYOXY(y*p)*p*p / 3. * (2. * y*y * p*p + 5. * y*p - 3. / 2.) + p*p / 6. * YMXLOGYMX(y*p, 1.)*(3 * y*p + 7) + p*p / 6. * J2(y*p);
	double FA3 = OMXYLOGOMXYOXY(-y*p)*p*p / 3. * ((1 + y*p)*(2 * y*p - 1) + 1. / 2.) + p*p / 6. * YMXLOGYMX(-y*p, 1.)*(3 * y*p + 1) + p*p / 6. * J2(-y*p);
	double FA4 = 1. / 3. * y*p*p*p - p*p / 3. + p*y;

	return FA1 + FA2 + FA3 + FA4;
}

double pImJ_FA_nochecks(double y, double p)
{
	double FA1 = (p - y)*(p - y) * log(abs(p - y));
	double FA2 = -(1. - y*p)*log1p(-y*p) / (y*p)*p*p / 3. * (2. * y*y * p*p + 5 * y*p - 3. / 2.) + p*p / 6. * (1 - y*p)*log(1 - y*p)*(3. * y*p + 7.) + p*p / 6. * ((1. - y*p)*log1p(-y*p) / (y*p) + 1.) / (y*p);
	double FA3 = -(1. + y*p)*log1p(y*p) / (y*p)*p*p / 3. * ((1 + y*p)*(2 * y*p - 1) + 1. / 2.) + p*p / 6. * (1 + y*p)*log(1 + y*p)*(3. * y*p + 1.) + p*p / 6. * ((1. + y*p)*log1p(y*p) / (y*p) - 1.) / (y*p);
	double FA4 = 1. / 3. * y*p*p*p - p*p / 3. + p*y;

	return FA1 + FA2 + FA3 + FA4;
}

double I_FB(double y, double p)
{
	static const double log2{ std::log(2) };
	static const double pi{ std::atan(1) * 4 };

	double FB1 = -4. * OMXATANYMOOXMO(y, p) - 2. * OMYLOGSQ(y, p);
	double FB2 = 4. * OMXATANYMOOXMO(-y, -p) + 2. * OMYLOGSQ(-y, -p);

	double FB3 = 3. / 2. * p*p * ((1 + y)*YMXLOGYMX(-y, 1.) - (1 - y)*YMXLOGYMX(y, 1.));
	double FB4 = y*((1 + p)*(1 + p) * YMXLOGYMX(-p, 1.) + (1 - p)*(1 - p) * YMXLOGYMX(p, 1.));

	double FB5 = std::log((y - 1)*(y - 1) + 4)*p / 4. * (3 * p*y*y - 6 * p*y + 3 * p - 8);
	double FB6 = -std::log((y + 1)*(y + 1) + 4)*p / 4. * (3 * p*y*y + 6 * p*y + 3 * p + 8);
	double FB7 = +std::log((p - 1)*(p - 1) + 4)*(2 - y / 2. - 2 * p + 3. / 2. * y*p - 3. / 2. * y*p*p + y*p*p*p / 2.);
	double FB8 = -std::log((p + 1)*(p + 1) + 4)*(2 + y / 2. + 2 * p + 3. / 2. * y*p + 3. / 2. * y*p*p + y*p*p*p / 2.);

	double FB9 = std::atan(y / 2. - 1. / 2.)*p*(3 * p + 2)*(y - 1);
	double FB10 = -std::atan(y / 2. + 1. / 2.)*p*(3 * p - 2)*(y + 1);
	double FB11 = -std::atan(p / 2. - 1. / 2.)*((3 * y + 1)*(-5 + 2 * p - p*p) + 8 * (1 + y));
	double FB12 = -std::atan(p / 2. + 1. / 2.)*((3 * y - 1)*(5 + 2 * p + p*p) + 8 * (1 - y));

	double FB13 = +2 * p*(6. * log2 - pi) + 3. * y*p*p * (2 * log2 + pi + 1);

	return FB1 + FB2 + FB3 + FB4 + FB5 + FB6 + FB7 + FB8 + FB9 + FB10 + FB11 + FB12 + FB13;
}

double J_FB(double y, double p)
{
	static const double log2{ std::log(2) };
	static const double pi{ std::atan(1) * 4 };

	double FB1 = -4. * OMXATANYMOOXMO(y, p) + (1 - y)*OMYLOGSQ(p, y) - (p + 1)*OMYLOGSQ(y, p);
	double FB2 = -4. * OMXATANYMOOXMO(-y, -p) + (1 + y)*OMYLOGSQ(-p, -y) + (p - 1)*OMYLOGSQ(-y, -p);

	double FB3 = p*p*p * ((1 + y)*YMXLOGYMX(-y, 1.) - (1 - y)*YMXLOGYMX(y, 1.));
	double FB4 = y / 4. * ((3. * p - 1)*(1 + p)*(1 + p) * YMXLOGYMX(-p, 1.) + (3 * p + 1)*(1 - p)*(1 - p) * YMXLOGYMX(p, 1.));

	double FB5 = std::log((y - 1)*(y - 1) + 4.)*p*p * (p*y*y - 2 * p*y + p - 2.) / 2.;
	double FB6 = -std::log((y + 1)*(y + 1) + 4.)*p*p * (p*y*y + 2 * p*y + p + 2.) / 2.;
	double FB7 = std::log((p - 1)*(p - 1) + 4.)*(-1. / 3. - p*p + 15. / 8. * y + 3. / 4. * y*p*p - y*p*p*p + 3. / 8. * y*pow(p,4));
	double FB8 = std::log((p + 1)*(p + 1) + 4.)*(-1. / 3. - p*p - 15. / 8. * y - 3. / 4. * y*p*p - y*p*p*p - 3. / 8. * y*pow(p,4));

	double FB9 = std::atan(y / 2. - 1. / 2.)*p*p * (2 * p + 1)*(y - 1);
	double FB10 = -std::atan(y / 2. + 1. / 2.)*p*p * (2 * p - 1)*(y + 1);
	double FB11 = -std::atan(p / 2. - 1. / 2.)*((y + 1. / 3.)*(-13. + 3 * p*p - 2 * p*p*p) + 8. * (1 + y));
	double FB12 = -std::atan(p / 2. + 1. / 2.)*((y - 1. / 3.)*(-13. + 3 * p*p + 2 * p*p*p) - 8. * (1 - y));

	double FB13 = 2. * y*(2. * log2 + pi + 1.)*p*p*p + (6. * log2 - pi + 2. / 3.)*p*p - (2. * y*p);

	return FB1 + FB2 + FB3 + FB4 + FB5 + FB6 + FB7 + FB8 + FB9 + FB10 + FB11 + FB12 + FB13;

}

double pImJ_FB(double y, double p)
{
	static const double log2{ std::log(2) };
	static const double pi{ std::atan(1) * 4 };

	double FB1 = 4. * (1 - p)*OMXATANYMOOXMO(y, p) - (1. - y)*OMYLOGSQ(p, y) + (1 - p)*OMYLOGSQ(y, p);
	double FB2 = 4. * (1 + p)*OMXATANYMOOXMO(-y, -p) - (1. + y)*OMYLOGSQ(-p, -y) + (1 + p)*OMYLOGSQ(-y, -p);

	double FB3 = p*p*p * ((1 + y)*YMXLOGYMX(-y, 1) - (1 - y)*YMXLOGYMX(y, 1)) / 2.;
	double FB4 = y*(pow((1 + p),3) * YMXLOGYMX(-p, 1) - pow((1 - p),3) * YMXLOGYMX(p, 1)) / 4.;

	double FB5 = log((y - 1)*(y - 1) + 4)*p*p / 4. * (p*y*y - 2 * p*y + p - 4);
	double FB6 = -log((y + 1)*(y + 1) + 4)*p*p / 4. * (p*y*y + 2 * p*y + p + 4);
	double FB7 = log((p - 1)*(p - 1) + 4)*(+1. / 3. - 15. / 8. * y + 2. * p - y*p / 2. - p*p + 3. / 4 * y*p*p - 1. / 2. * y*p*p*p + 1. / 8. * y*p*p*p*p);
	double FB8 = -log((p + 1)*(p + 1) + 4)*(-1. / 3. - 15. / 8. * y + 2. * p + y*p / 2. + p*p + 3. / 4 * y*p*p + 1. / 2. * y*p*p*p + 1. / 8. * y*p*p*p*p);

	double FB9 = atan(y / 2. - 1. / 2.)*p*p * (1 + p)*(y - 1);
	double FB10 = atan(y / 2. + 1. / 2.)*p*p * (1 - p)*(y + 1);
	double FB11 = atan(p / 2. - 1. / 2.)*(1 - p)*(2. * p - 15. * y + 6. * p*y - 3. * p*p * y - p*p + 11.) / 3.;
	double FB12 = -atan(p / 2. + 1. / 2.)*(1 + p)*(-2. * p + 15. * y + 6. * p*y + 3. * p*p * y - p*p + 11.) / 3.;

	double FB13 = p*p * (6. * log2 - pi - 2. / 3.) + y*p*p*p * (2. * log2 + pi + 1.) + (2 * y*p);

	return FB1 + FB2 + FB3 + FB4 + FB5 + FB6 + FB7 + FB8 + FB9 + FB10 + FB11 + FB12 + FB13;
}

double pImJ_FB_nochecks(double y, double p)
{
	static const double log2{ std::log(2) };
	static const double pi{ std::atan(1) * 4 };

	double FB1 = 4. * (1 - p)*(1 - y)*atan((p - 1) / (y - 1)) + ((1 - p)*(1 - p) - (1 - y)*(1 - y))*log((y - 1)*(y - 1) + (p - 1)*(p - 1));
	double FB2 = 4. * (1 + p)*(1 + y)*atan((p + 1) / (y + 1)) + ((1 + p)*(1 + p) - (1 + y)*(1 + y))*log((y + 1)*(y + 1) + (p + 1)*(p + 1));

	double FB3 = p*p*p * ((1 + y)*(1 + y) * log(1 + y) - (1 - y)*(1 - y) * log(1 - y)) / 2.;
	double FB4 = y*(pow((1 + p),4) * log(1 + p) - pow((1 - p),4) * log(1 - p)) / 4.;

	double FB5 = log((y - 1)*(y - 1) + 4)*p*p / 4. * (p*y*y - 2 * p*y + p - 4);
	double FB6 = -log((y + 1)*(y + 1) + 4)*p*p / 4. * (p*y*y + 2 * p*y + p + 4);
	double FB7 = log((p - 1)*(p - 1) + 4)*(+1. / 3. - 15. / 8. * y + 2. * p - y*p / 2. - p*p + 3. / 4 * y*p*p - 1. / 2. * y*p*p*p + 1. / 8. * y*p*p*p*p);
	double FB8 = -log((p + 1)*(p + 1) + 4)*(-1. / 3. - 15. / 8. * y + 2. * p + y*p / 2. + p*p + 3. / 4 * y*p*p + 1. / 2. * y*p*p*p + 1. / 8. * y*p*p*p*p);

	double FB9 = atan(y / 2. - 1. / 2.)*p*p * (1 + p)*(y - 1);
	double FB10 = atan(y / 2. + 1. / 2.)*p*p * (1 - p)*(y + 1);
	double FB11 = atan(p / 2. - 1. / 2.)*(1 - p)*(2. * p - 15. * y + 6. * p*y - 3. * p*p * y - p*p + 11.) / 3.;
	double FB12 = -atan(p / 2. + 1. / 2.)*(1 + p)*(-2. * p + 15. * y + 6. * p*y + 3. * p*p * y - p*p + 11.) / 3.;

	double FB13 = p*p * (6. * log2 - pi - 2. / 3.) + y*p*p*p * (2. * log2 + pi + 1.) + (2 * y*p);

	return FB1 + FB2 + FB3 + FB4 + FB5 + FB6 + FB7 + FB8 + FB9 + FB10 + FB11 + FB12 + FB13;
}

double I_FP(double y, double p)
{
	static const double q22{ -1.23991592 }, q42{ 1.08897876 }, q44{ -9.726553 }, q62{ -1.50465044 }, q64{ 18.273236 }, q66{ -38.99971412 }, q82{ 0.74180336 }, q84{ -9.64366612 }, q86{ 22.03387365 }, q88{ -13.05630027 };
	static const double q31{ 0.74952005 }, q51{ -0.08175407 }, q53{ 5.10578057 }, q71{ 0.04700608 }, q73{ -2.9709584 }, q75{ 9.15391675 };

	double FP1 = ((q88*pow(y,8)) / 9 + (q86*pow(y,6)) / 9 + (q84*pow(y,4)) / 9 + (q82*y*y) / 9 - q82 / 27 - q84 / 45 - q86 / 63 - q88 / 81)*pow(p,9) + ((q86*pow(y,8)) / 7 + (q66*pow(y,6)) / 7 + (q64*pow(y,4)) / 7 + (q62*y*y) / 7 - q62 / 21 - q64 / 35 - q66 / 49 - q86 / 63)*pow(p,7) + ((q84*pow(y,8)) / 5 + (q64*pow(y,6)) / 5 + (q44*pow(y,4)) / 5 + (q42*y*y) / 5 - q42 / 15 - q44 / 25 - q64 / 35 - q84 / 45)*pow(p,5) + ((q82*pow(y,8)) / 3 + (q62*pow(y,6)) / 3 + (q42*pow(y,4)) / 3 + (q22*y*y) / 3 - q22 / 9 - q42 / 15 - q62 / 21 - q82 / 27)*p*p*p + ((-q82 / 3 - q84 / 5 - q86 / 7 - q88 / 9)*pow(y,8) + (-q62 / 3 - q64 / 5 - q66 / 7 - q86 / 9)*pow(y,6) + (-q42 / 3 - q44 / 5 - q64 / 7 - q84 / 9)*pow(y,4) + (-q22 / 3 - q42 / 5 - q62 / 7 - q82 / 9)*y*y + q22 / 9 + (2 * q42) / 15 + q44 / 25 + (2 * q62) / 21 + (2 * q64) / 35 + q66 / 49 + (2 * q82) / 27 + (2 * q84) / 45 + (2 * q86) / 63 + q88 / 81)*p;
	double FP2 = ((q71*y) / 8 - (pow(y,7) * (3 * q71 + (9 * q73) / 5 + (9 * q75) / 7)) / 8 + (q73*y*y*y) / 8 + (q75*pow(y,5)) / 8)*pow(p,8) + ((q51*y) / 6 - (pow(y,5) * ((7 * q51) / 3 + (7 * q53) / 5 + (7 * q75) / 9)) / 6 + (q53*y*y*y) / 6 + (q75*pow(y,7)) / 6)*pow(p,6) + ((q31*y) / 4 - (y*y*y * ((5 * q31) / 3 + (5 * q53) / 7 + (5 * q73) / 9)) / 4 + (q53*pow(y,5)) / 4 + (q73*pow(y,7)) / 4)*pow(p,4) + ((q31*y*y*y) / 2 - (y*((3 * q31) / 5 + (3 * q51) / 7 + q71 / 3)) / 2 + (q51*pow(y,5)) / 2 + (q71*pow(y,7)) / 2)*p*p;

	return FP1 + FP2;
}

double J_FP(double y, double p)
{
	static const double q22{ -1.23991592 }, q42{ 1.08897876 }, q44{ -9.726553 }, q62{ -1.50465044 }, q64{ 18.273236 }, q66{ -38.99971412 }, q82{ 0.74180336 }, q84{ -9.64366612 }, q86{ 22.03387365 }, q88{ -13.05630027 };
	static const double q31{ 0.74952005 }, q51{ -0.08175407 }, q53{ 5.10578057 }, q71{ 0.04700608 }, q73{ -2.9709584 }, q75{ 9.15391675 };

	double FP1 = ((q88*pow(y,8)) / 10 + (q86*pow(y,6)) / 10 + (q84*pow(y,4)) / 10 + (q82*y*y) / 10 - q82 / 30 - q84 / 50 - q86 / 70 - q88 / 90)*pow(p,10) + ((q86*pow(y,8)) / 8 + (q66*pow(y,6)) / 8 + (q64*pow(y,4)) / 8 + (q62*y*y) / 8 - q62 / 24 - q64 / 40 - q66 / 56 - q86 / 72)*pow(p,8) + ((q84*pow(y,8)) / 6 + (q64*pow(y,6)) / 6 + (q44*pow(y,4)) / 6 + (q42*y*y) / 6 - q42 / 18 - q44 / 30 - q64 / 42 - q84 / 54)*pow(p,6) + ((q82*pow(y,8)) / 4 + (q62*pow(y,6)) / 4 + (q42*pow(y,4)) / 4 + (q22*y*y) / 4 - q22 / 12 - q42 / 20 - q62 / 28 - q82 / 36)*pow(p,4) + (q22 / 18 + q42 / 15 + q44 / 50 + q62 / 21 + q64 / 35 + q66 / 98 + q82 / 27 + q84 / 45 + q86 / 63 + q88 / 162 - (q22*y*y) / 6 - (q42*y*y) / 10 - (q42*pow(y,4)) / 6 - (q44*pow(y,4)) / 10 - (q62*y*y) / 14 - (q62*pow(y,6)) / 6 - (q64*pow(y,4)) / 14 - (q64*pow(y,6)) / 10 - (q66*pow(y,6)) / 14 - (q82*y*y) / 18 - (q84*pow(y,4)) / 18 - (q82*pow(y,8)) / 6 - (q84*pow(y,8)) / 10 - (q86*pow(y,6)) / 18 - (q86*pow(y,8)) / 14 - (q88*pow(y,8)) / 18)*p*p;
	double FP2 = (y*p*p*p * (p*p - 1)*(63 * q31 + 45 * q51 + 35 * q71 - 105 * q31*y*y - 105 * q51*pow(y,4) - 105 * q71*pow(y,6) + 45 * q51*p*p + 35 * q71*p*p + 35 * q71*pow(p,4) - 105 * q51*pow(y,4) * p*p + 45 * q53*y*y * p*p - 63 * q53*pow(y,4) * p*p + 35 * q73*y*y * p*p - 105 * q71*pow(y,6) * p*p + 35 * q73*y*y * pow(p,4) - 105 * q71*pow(y,6) * pow(p,4) - 63 * q73*pow(y,6) * p*p - 63 * q73*pow(y,6) * pow(p,4) + 35 * q75*pow(y,4) * pow(p,4) - 45 * q75*pow(y,6) * pow(p,4))) / 315;

	return FP1 + FP2;
}

double I_FP_alt(double y, double p)
{
	static const double q22{ -1.23991592 }, q42{ 1.08897876 }, q44{ -9.726553 }, q62{ -1.50465044 }, q64{ 18.273236 }, q66{ -38.99971412 }, q82{ 0.74180336 }, q84{ -9.64366612 }, q86{ 22.03387365 }, q88{ -13.05630027 };
	static const double q31{ 0.74952005 }, q51{ -0.08175407 }, q53{ 5.10578057 }, q71{ 0.04700608 }, q73{ -2.9709584 }, q75{ 9.15391675 };
	static const double a1 = -q82 / 27 - q84 / 45 - q86 / 63 - q88 / 81;
	static const double a2 = -q62 / 21 - q64 / 35 - q66 / 49 - q86 / 63;
	static const double a3 = -q42 / 15 - q44 / 25 - q64 / 35 - q84 / 45;
	static const double a4 = -q22 / 9 - q42 / 15 - q62 / 21 - q82 / 27;
	static const double a5 = -q82 / 3 - q84 / 5 - q86 / 7 - q88 / 9;
	static const double a6 = -q62 / 3 - q64 / 5 - q66 / 7 - q86 / 9;
	static const double a7 = -q42 / 3 - q44 / 5 - q64 / 7 - q84 / 9;
	static const double a8 = -q22 / 3 - q42 / 5 - q62 / 7 - q82 / 9;
	static const double a9 = +q22 / 9 + (2 * q42) / 15 + q44 / 25 + (2 * q62) / 21 + (2 * q64) / 35 + q66 / 49 + (2 * q82) / 27 + (2 * q84) / 45 + (2 * q86) / 63 + q88 / 81;
	static const double a10 = 3 * q71 + (9 * q73) / 5 + (9 * q75) / 7;
	static const double a11 = (7 * q51) / 3 + (7 * q53) / 5 + (7 * q75) / 9;
	static const double a12 = (5 * q31) / 3 + (5 * q53) / 7 + (5 * q73) / 9;
	static const double a13 = (3 * q31) / 5 + (3 * q51) / 7 + q71 / 3;
	
	double FP1 = ((q88*pow(y,8)) / 9 + (q86*pow(y,6)) / 9 + (q84*pow(y,4)) / 9 + (q82*y*y) / 9 + a1)*pow(p,9) + ((q86*pow(y,8)) / 7 + (q66*pow(y,6)) / 7 + (q64*pow(y,4)) / 7 + (q62*y*y) / 7 +a2)*pow(p,7) + ((q84*pow(y,8)) / 5 + (q64*pow(y,6)) / 5 + (q44*pow(y,4)) / 5 + (q42*y*y) / 5 +a3)*pow(p,5) + ((q82*pow(y,8)) / 3 + (q62*pow(y,6)) / 3 + (q42*pow(y,4)) / 3 + (q22*y*y) / 3 +a4)*p*p*p + ((a5)*pow(y,8) + (a6)*pow(y,6) + (a7)*pow(y,4) + (a8)*y*y +a9)*p;
	double FP2 = ((q71*y) / 8 - (pow(y,7) * (a10)) / 8 + (q73*y*y*y) / 8 + (q75*pow(y,5)) / 8)*pow(p,8) + ((q51*y) / 6 - (pow(y,5) * (a11)) / 6 + (q53*y*y*y) / 6 + (q75*pow(y,7)) / 6)*pow(p,6) + ((q31*y) / 4 - (y*y*y * (a12)) / 4 + (q53*pow(y,5)) / 4 + (q73*pow(y,7)) / 4)*pow(p,4) + ((q31*y*y*y) / 2 - (y*(a13)) / 2 + (q51*pow(y,5)) / 2 + (q71*pow(y,7)) / 2)*p*p;

	return FP1 + FP2;
}

double pImJ_FP(double y, double p)
{
	static const double q22{ -1.23991592 }, q42{ 1.08897876 }, q44{ -9.726553 }, q62{ -1.50465044 }, q64{ 18.273236 }, q66{ -38.99971412 }, q82{ 0.74180336 }, q84{ -9.64366612 }, q86{ 22.03387365 }, q88{ -13.05630027 };
	static const double q31{ 0.74952005 }, q51{ -0.08175407 }, q53{ 5.10578057 }, q71{ 0.04700608 }, q73{ -2.9709584 }, q75{ 9.15391675 };

	double FP1 = -p*p * ((p*p * q22) / 36 - q42 / 15 - q44 / 50 - q62 / 21 - q64 / 35 - q66 / 98 - q82 / 27 - q84 / 45 - q86 / 63 - q88 / 162 - q22 / 18 + (p*p * q42) / 60 + (pow(p,4) * q42) / 90 + (pow(p,4) * q44) / 150 + (p*p * q62) / 84 + (pow(p,4) * q64) / 210 + (pow(p,6) * q62) / 168 + (pow(p,6) * q64) / 280 + (pow(p,6) * q66) / 392 + (p*p * q82) / 108 + (pow(p,4) * q84) / 270 + (pow(p,8) * q82) / 270 + (pow(p,6) * q86) / 504 + (pow(p,8) * q84) / 450 + (pow(p,8) * q86) / 630 + (pow(p,8) * q88) / 810 + (q22*y*y) / 6 + (q42*y*y) / 10 + (q42*pow(y,4)) / 6 + (q44*pow(y,4)) / 10 + (q62*y*y) / 14 + (q62*pow(y,6)) / 6 + (q64*pow(y,4)) / 14 + (q64*pow(y,6)) / 10 + (q66*pow(y,6)) / 14 + (q82*y*y) / 18 + (q84*pow(y,4)) / 18 + (q82*pow(y,8)) / 6 + (q84*pow(y,8)) / 10 + (q86*pow(y,6)) / 18 + (q86*pow(y,8)) / 14 + (q88*pow(y,8)) / 18 - (p*p * q22*y*y) / 12 - (p*p * q42*pow(y,4)) / 12 - (pow(p,4) * q42*y*y) / 30 - (pow(p,4) * q44*pow(y,4)) / 30 - (p*p * q62*pow(y,6)) / 12 - (pow(p,6) * q62*y*y) / 56 - (pow(p,4) * q64*pow(y,6)) / 30 - (pow(p,6) * q64*pow(y,4)) / 56 - (pow(p,6) * q66*pow(y,6)) / 56 - (p*p * q82*pow(y,8)) / 12 - (pow(p,8) * q82*y*y) / 90 - (pow(p,4) * q84*pow(y,8)) / 30 - (pow(p,8) * q84*pow(y,4)) / 90 - (pow(p,6) * q86*pow(y,8)) / 56 - (pow(p,8) * q86*pow(y,6)) / 90 - (pow(p,8) * q88*pow(y,8)) / 90);
	double FP2 = -(p*p*p * y*(756 * q31 + 540 * q51 + 420 * q71 - 378 * p*p * q31 - 180 * pow(p,4) * q51 - 105 * pow(p,6) * q71 - 1260 * q31*y*y - 1260 * q51*pow(y,4) - 1260 * q71*pow(y,6) + 630 * p*p * q31*y*y + 270 * p*p * q53*y*y - 378 * p*p * q53*pow(y,4) + 420 * pow(p,4) * q51*pow(y,4) - 180 * pow(p,4) * q53*y*y + 252 * pow(p,4) * q53*pow(y,4) + 210 * p*p * q73*y*y - 378 * p*p * q73*pow(y,6) - 105 * pow(p,6) * q73*y*y + 140 * pow(p,4) * q75*pow(y,4) + 315 * pow(p,6) * q71*pow(y,6) - 180 * pow(p,4) * q75*pow(y,6) + 189 * pow(p,6) * q73*pow(y,6) - 105 * pow(p,6) * q75*pow(y,4) + 135 * pow(p,6) * q75*pow(y,6))) / 7560;

	return FP1 + FP2;
}

double I_calc(double y, double r)
{
	static const double pi{ std::atan(1) * 4 };
	static const double A = -1 / pi;
	static const double B = -0.19532775;

	return A*I_FA(y, r) + B*I_FB(y, r) + I_FP(y, r);
}

double J_calc(double y, double r)
{
	static const double pi{ std::atan(1.) * 4 };
	static const double A = -1 / pi;
	static const double B = -0.19532775;

	return A*J_FA(y, r) + B*J_FB(y, r) + J_FP(y, r);
}

double FAa(double y, double p)
{
	static const double pi{ std::atan(1) * 4 };
	static const double A = -1 / pi;
	return A*2.0*std::log(std::fabs(p - y));
}

double dFAa_dp(double y, double p)
{
	static const double pi{ std::atan(1) * 4 };
	static const double A = -1 / pi;
	return -A*2.0/(y-p);
}

double FA(double y, double p)
{
	double FA1 = 2. * std::log(fabs(p - y));
	double FA2 = 2. * std::log1p(-y*p)*(-1. + y*p + y*y * p*p);
	double FA3 = 2. * std::log1p(y*p)*(-y*p - y*y * p*p);
	double FA4 = 2. * (2 * y*p + 1);

	return FA1 + FA2 + FA3 + FA4;
}

double D_FA(double y, double p)
{
	double FA1 = 2. / (y - p);
	double FA2 = 2. * std::log1p(-y*p)*(2 * y*p*p + p) + (2 * p*(y*y * p*p + y*p - 1)) / (y*p - 1);
	double FA3 = -2. * std::log1p(y*p)*(2 * y*p*p + p) - (2 * p*(y*y * p*p + y*p)) / (y*p + 1);
	double FA4 = 4. * p;

	return FA1 + FA2 + FA3 + FA4;
}

double FB(double y, double p)
{
	static const double log2{ std::log(2) };
	static const double pi{ std::atan(1) * 4 };

	double FB1 = 2. * log((1. - y)*(1. - y) + (1. - p)*(1. - p));
	double FB2 = 2. * log((1. + y)*(1. + y) + (1. + p)*(1. + p));

	double FB3 = 3. * p*((1. + y)*YMXLOGYMX(-y, 1.) - (1. - y)*YMXLOGYMX(y, 1.));
	double FB4 = 3. * y*((1. + p)*YMXLOGYMX(-p, 1.) - (1. - p)*YMXLOGYMX(p, 1.));

	double FB5 = (3. / 2. * p - 3 * y*p + 3. / 2. * p*y*y - 2.)*log((1. - y)*(1. - y) + 4.);
	double FB6 = (-3. / 2. * p - 3 * y*p - 3. / 2. * p*y*y - 2.)*log((1. + y)*(1. + y) + 4.);
	double FB7 = (3. / 2. * y - 3 * p*y + 3. / 2. * y*p*p - 2.)*log((1. - p)*(1. - p) + 4.);
	double FB8 = (-3. / 2. * y - 3 * p*y - 3. / 2. * y*p*p - 2.)*log((1. + p)*(1. + p) + 4.);

	double FB9 = 2. * (1. - y + 3. * p - 3. * y*p)*atan((1. - y) / 2.);
	double FB10 = 2. * (1. + y - 3. * p - 3. * y*p)*atan((1. + y) / 2.);
	double FB11 = 2. * (1. - p + 3. * y - 3. * y*p)*atan((1. - p) / 2.);
	double FB12 = 2. * (1. + p - 3. * y - 3. * y*p)*atan((1. + p) / 2.);

	double FB13 = (6. * y*p*(pi + 2. * log2 + 1.) - 2. * (pi - 6. * log2 - 2.));

	return FB1 + FB2 + FB3 + FB4 + FB5 + FB6 + FB7 + FB8 + FB9 + FB10 + FB11 + FB12 + FB13;
}

double D_FB(double y, double p)
{
	static const double log2{ std::log(2) };
	static const double pi{ std::atan(1) * 4 };

	double FB1 = (2. * (2. * y - 2.)) / ((y - 1.)*(y - 1.) + (p - 1.)*(p - 1.));
	double FB2 = (2. * (2. * y + 2.)) / ((y + 1.)*(y + 1.) + (p + 1.)*(p + 1.));

	double FB3 = 6. * p*(YMXLOGYMX(-y, 1.) + YMXLOGYMX(y, 1.) + 1.);
	double FB4 = 3. * YMXLOGYMX(-p, 1.)*(p + 1.) - 3. * YMXLOGYMX(p, 1.)*(1. - p);

	double FB5 = ((2. * y - 2.)*((3. * p*y*y) / 2. - 3. * p*y + (3. * p) / 2. - 2.)) / ((y - 1)*(y - 1) + 4.) - log((y - 1)*(y - 1) + 4.)*(3. * p - 3. * y*p);
	double FB6 = -log((y + 1)*(y + 1) + 4.)*(3. * p + 3. * y*p) - ((2. * y + 2.)*((3. * p*y*y) / 2. + 3. * p*y + (3. * p) / 2. + 2.)) / ((y + 1)*(y + 1) + 4.);
	double FB7 = log((p - 1)*(p - 1) + 4.)*((3. * p*p) / 2. - 3. * p + 3. / 2.);
	double FB8 = -log((p + 1)*(p + 1) + 4.)*((3. * p*p) / 2. + 3. * p + 3. / 2.);

	double FB9 = atan(y / 2. - 1. / 2.)*(6. * p + 2.) + (4. * (3. * p + 1.)*(y - 1.)) / (y*y - 2. * y + 5.);
	double FB10 = -atan(y / 2. + 1. / 2.)*(6. * p - 2.) - (4. * (3. * p - 1.)*(y + 1.)) / (y*y + 2. * y + 5.);
	double FB11 = atan(p / 2. - 1. / 2.)*(6. * p - 6.);
	double FB12 = -atan(p / 2. + 1. / 2.)*(6. * p + 6.);

	double FB13 = 6. * p*(2. * log2 + pi + 1.);

	return FB1 + FB2 + FB3 + FB4 + FB5 + FB6 + FB7 + FB8 + FB9 + FB10 + FB11 + FB12 + FB13;
}

double FP(double y, double p)
{
	static const double q22{ -1.23991592 }, q42{ 1.08897876 }, q44{ -9.726553 }, q62{ -1.50465044 }, q64{ 18.273236 }, q66{ -38.99971412 }, q82{ 0.74180336 }, q84{ -9.64366612 }, q86{ 22.03387365 }, q88{ -13.05630027 };
	static const double q31{ 0.74952005 }, q51{ -0.08175407 }, q53{ 5.10578057 }, q71{ 0.04700608 }, q73{ -2.9709584 }, q75{ 9.15391675 };

	double FP1 = q22 / 9 + (2 * q42) / 15 + q44 / 25 + (2 * q62) / 21 + (2 * q64) / 35 + q66 / 49 + (2 * q82) / 27 + (2 * q84) / 45 + (2 * q86) / 63 + q88 / 81 - (y*y + p*p)*(q22 / 3 + q42 / 5 + q62 / 7 + q82 / 9) - (pow(y,4) + pow(p,4))*(q42 / 3 + q44 / 5 + q64 / 7 + q84 / 9) - (pow(y,6) + pow(p,6))*(q62 / 3 + q64 / 5 + q66 / 7 + q86 / 9) - (pow(y,8) + pow(p,8))*(q82 / 3 + q84 / 5 + q86 / 7 + q88 / 9) + q22*y*y * p*p + q44*pow(y,4) * pow(p,4) + q66*pow(y,6) * pow(p,6) + q88*pow(y,8) * pow(p,8) + q42*y*y * p*p * (y*y + p*p) + q62*y*y * p*p * (pow(y,4) + pow(p,4)) + q64*pow(y,4) * pow(p,4) * (y*y + p*p) + q82*y*y * p*p * (pow(y,6) + pow(p,6)) + q84*pow(y,4) * pow(p,4) * (pow(y,4) + pow(p,4)) + q86*pow(y,6) * pow(p,6) * (y*y + p*p);
	double FP2 = q53*y*y*y * p*p*p * (y*y + p*p) - y*y*y * p*p*p * ((5 * q31) / 3 + (5 * q53) / 7 + (5 * q73) / 9) - pow(y,5) * pow(p,5) * ((7 * q51) / 3 + (7 * q53) / 5 + (7 * q75) / 9) - pow(y,7) * pow(p,7) * (3 * q71 + (9 * q73) / 5 + (9 * q75) / 7) - y*p*((3 * q31) / 5 + (3 * q51) / 7 + q71 / 3) + q73*y*y*y * p*p*p * (pow(y,4) + pow(p,4)) + q75*pow(y,5) * pow(p,5) * (y*y + p*p) + q31*y*p*(y*y + p*p) + q51*y*p*(pow(y,4) + pow(p,4)) + q71*y*p*(pow(y,6) + pow(p,6));
	
	return FP1 + FP2;
}

double D_FP(double y, double p)
{
	static const double q22{ -1.23991592 }, q42{ 1.08897876 }, q44{ -9.726553 }, q62{ -1.50465044 }, q64{ 18.273236 }, q66{ -38.99971412 }, q82{ 0.74180336 }, q84{ -9.64366612 }, q86{ 22.03387365 }, q88{ -13.05630027 };
	static const double q31{ 0.74952005 }, q51{ -0.08175407 }, q53{ 5.10578057 }, q71{ 0.04700608 }, q73{ -2.9709584 }, q75{ 9.15391675 };

	double FP1 = 2 * q22*y*p*p - (2 * q42*y) / 5 - (2 * q62*y) / 7 - (2 * q82*y) / 9 - (4 * q42*y*y*y) / 3 - (4 * q44*y*y*y) / 5 - 2 * q62*pow(y,5) - (4 * q64*y*y*y) / 7 - (6 * q64*pow(y,5)) / 5 - (6 * q66*pow(y,5)) / 7 - (4 * q84*y*y*y) / 9 - (8 * q82*pow(y,7)) / 3 - (8 * q84*pow(y,7)) / 5 - (2 * q86*pow(y,5)) / 3 - (8 * q86*pow(y,7)) / 7 - (8 * q88*pow(y,7)) / 9 - (2 * q22*y) / 3 + 2 * q42*y*pow(p,4) + 2 * q62*y*pow(p,6) + 2 * q82*y*pow(p,8) + 4 * q42*y*y*y * p*p + 4 * q44*y*y*y * pow(p,4) + 6 * q62*pow(y,5) * p*p + 4 * q64*y*y*y * pow(p,6) + 6 * q64*pow(y,5) * pow(p,4) + 6 * q66*pow(y,5) * pow(p,6) + 8 * q82*pow(y,7) * p*p + 4 * q84*y*y*y * pow(p,8) + 8 * q84*pow(y,7) * pow(p,4) + 6 * q86*pow(y,5) * pow(p,8) + 8 * q86*pow(y,7) * pow(p,6) + 8 * q88*pow(y,7) * pow(p,8);
	double FP2 = 2 * q31*y*y * p - p*((3 * q31) / 5 + (3 * q51) / 7 + q71 / 3) + 4 * q51*pow(y,4) * p + 6 * q71*pow(y,6) * p + q31*p*(y*y + p*p) + q51*p*(pow(y,4) + pow(p,4)) + q71*p*(pow(y,6) + pow(p,6)) + 2 * q53*pow(y,4) * p*p*p + 4 * q73*pow(y,6) * p*p*p + 2 * q75*pow(y,6) * pow(p,5) - 3 * y*y * p*p*p * ((5 * q31) / 3 + (5 * q53) / 7 + (5 * q73) / 9) - 5 * pow(y,4) * pow(p,5) * ((7 * q51) / 3 + (7 * q53) / 5 + (7 * q75) / 9) - 7 * pow(y,6) * pow(p,7) * (3 * q71 + (9 * q73) / 5 + (9 * q75) / 7) + 3 * q53*y*y * p*p*p * (y*y + p*p) + 3 * q73*y*y * p*p*p * (pow(y,4) + pow(p,4)) + 5 * q75*pow(y,4) * pow(p,5) * (y*y + p*p);

	return FP1 + FP2;
}

double I_FAb(double y, double p)
{
	double FA2 = -OMXYLOGOMXYOXY(y*p)*p / 3. * (2 * y*y * p*p + 5 * y*p - 1);
	double FA3 = OMXYLOGOMXYOXY(-y*p)*p / 3. * (1 + y*p)*(2 * y*p - 1);
	double FA4 = 4. / 3. * y*p*p+2.*(p-y);

	return FA2 + FA3 + FA4;
}

double J_FAb(double y, double p)
{
	double FA2 = -p*p / 6. * (OMXYLOGOMXYOXY(y*p) + YMXLOGYMX(y*p, 1.)*(3. * y*p + 7.) + J2(y*p));
	double FA3 = -p*p / 6. * (OMXYLOGOMXYOXY(-y*p) + YMXLOGYMX(-y*p, 1.)*(3. * y*p + 1.) + J2(-y*p));
	double FA4 = y*p*p*p + 5.*p*p/6.;

	return FA2 + FA3 + FA4;
}

double Ib_calc(double y, double r)
{
	static const double pi{ std::atan(1) * 4 };
	static const double A = -1 / pi;
	static const double B = -0.19532775;

	return A*I_FAb(y, r) + B*I_FB(y, r) + I_FP(y, r);
}

double Jb_calc(double y, double r)
{
	static const double pi{ std::atan(1.) * 4 };
	static const double A = -1 / pi;
	static const double B = -0.19532775;

	return A*J_FAb(y, r) + B*J_FB(y, r) + J_FP(y, r);
}

double pImJ_calc(double y, double r)
{
	static const double pi{ std::atan(1.) * 4 };
	static const double A = -1 / pi;
	static const double B = -0.19532775;

	return A*pImJ_FA(y, r) + B*pImJ_FB(y, r) + pImJ_FP(y, r);
}

void Imat_calc(const Vector& Y, const Vector& R, Matrix& Imat)
{
	for (size_t i = 0; i != Y.Size(); i++)
	{
		for (size_t j = 0; j != R.Size(); j++)
		{
			Imat(i,j) = I_calc(Y[i], R[j]);
		}
	}
	return;
}

void Jmat_calc(const Vector& Y, const Vector& R, Matrix& Jmat)
{
	for (size_t i = 0; i != Y.Size(); i++)
	{
		for (size_t j = 0; j != R.Size(); j++)
		{
			Jmat(i,j) = J_calc(Y[i], R[j]);
		}
	}
	return;
}

void Imatb_calc(const Vector& Y, const Vector& R, Matrix& Imat)
{
	for (size_t i = 0; i != Y.Size(); i++)
	{
		for (size_t j = 0; j != R.Size(); j++)
		{
			Imat(i, j) = Ib_calc(Y[i], R[j]);
		}
	}
	return;
}

void Jmatb_calc(const Vector& Y, const Vector& R, Matrix& Jmat)
{
	for (size_t i = 0; i != Y.Size(); i++)
	{
		for (size_t j = 0; j != R.Size(); j++)
		{
			Jmat(i, j) = Jb_calc(Y[i], R[j]);
		}
	}
	return;
}

void pImJmat_calc(const Vector& Y, const Vector& R, Matrix& pImJmat)
{
	for (size_t i = 0; i != Y.Size(); i++)
	{
		for (size_t j = 0; j != R.Size(); j++)
		{
			pImJmat(i, j) = pImJ_calc(Y[i], R[j]);
		}
	}
	return;
}

void Im1_calc(const Vector& Y, Vector& Im1)
{
	for (size_t i = 0; i != Y.Size(); i++)
	{
		Im1(i) =I_calc(Y[i], -1.0);
	}
	return;
}

void Jm1_calc(const Vector& Y, Vector& Jm1)
{
	for (size_t i = 0; i != Y.Size(); i++)
	{
		Jm1(i) =J_calc(Y[i], -1.0);
	}
	return;
}

void Im1b_calc(const Vector& Y, Vector& Im1)
{
	for (size_t i = 0; i != Y.Size(); i++)
	{
		Im1(i) = Ib_calc(Y[i], -1.0);
	}
	return;
}

void Jm1b_calc(const Vector& Y, Vector& Jm1)
{
	for (size_t i = 0; i != Y.Size(); i++)
	{
		Jm1(i) = Jb_calc(Y[i], -1.0);
	}
	return;
}

void Usgm_trapz(const Vector& Yw, Matrix& Usgm)
{
	Matrix CC (Yw.Size(),Yw.Size());
	for (int i = 0; i != Yw.Size(); i++)
	{
		if (i > 0) {
			CC(i - 1,i) += -1.0 / (Yw(i - 1) - Yw(i));
			CC(i,i) += 1.0 / (Yw(i - 1) - Yw(i));
		}
		if (i < Yw.Size() - 1) {
			CC(i,i) += 1.0 / (Yw(i) - Yw(i + 1));
			CC(i + 1,i) += -1.0 / (Yw(i) - Yw(i + 1));
		}
	}

	Matrix Imat(Yw.Size(), Yw.Size());
	Matrix Jmat(Yw.Size(), Yw.Size());
	Vector Im1(Yw.Size());
	Vector Jm1(Yw.Size());
	Imat_calc(Yw, Yw, Imat);
	Jmat_calc(Yw, Yw, Jmat);
	Im1_calc(Yw, Im1);
	Jm1_calc(Yw, Jm1);

	Matrix Us(Yw.Size(), Yw.Size());
	for (size_t i = 0; i != Yw.Size(); i++)
	{
		for (size_t k = 0; k != Yw.Size(); k++)
		{
			Us(k,i) = (Yw(i) * Imat(k,i) - Jmat(k,i)) - Im1(k) * Yw(i) + Jm1(k);
		}
	}

	Usgm = Us* CC;

	return;
}

void triangle_dispslope_disps(const Vector& R, const Vector& Y, Matrix& U, Matrix& dU_dR)
{
	Matrix Imat(Y.Size(), R.Size());
	Matrix Jmat(Y.Size(), R.Size());
	Vector Im1(Y.Size());
	Vector Jm1(Y.Size());
	Imat_calc(Y, R, Imat);
	Jmat_calc(Y, R, Jmat);
	Im1_calc(Y, Im1);
	Jm1_calc(Y, Jm1);

	for (size_t i = 0; i != R.Size(); i++)
	{
		for (size_t k = 0; k != Y.Size(); k++)
		{
			U(k,i) = (R(i) * Imat(k,i) - Jmat(k,i)) - Im1(k) * R(i) + Jm1(k);
			dU_dR(k,i) = Imat(k,i) - Im1(k);
		}
	}

	return;
}

void triangle_dispslope_disps_givenMat1(const Vector& R, const Vector& Y, const Vector& Im1, const Vector& Jm1, Matrix& U, Matrix& dU_dR)
{
	Matrix Imat(Y.Size(), R.Size());
	Matrix Jmat(Y.Size(), R.Size());
	Imat_calc(Y, R, Imat);
	Jmat_calc(Y, R, Jmat);

	for (size_t i = 0; i != R.Size(); i++)
	{
		for (size_t k = 0; k != Y.Size(); k++)
		{
			U(k, i) = (R(i) * Imat(k, i) - Jmat(k, i)) - Im1(k) * R(i) + Jm1(k);
			dU_dR(k, i) = Imat(k, i) - Im1(k);
		}
	}

	return;
}

void triangle_dispslope_disps_2(const Vector& R, const Vector& Y, const Vector& Im1, const Vector& Jm1, Matrix& U, Matrix& dU_dR)
{
	Matrix pImJmat(Y.Size(), R.Size());
	Matrix Imat(Y.Size(), R.Size());
	pImJmat_calc(Y, R, pImJmat);
	Imat_calc(Y, R, Imat);

	for (size_t i = 0; i != R.Size(); i++)
	{
		for (size_t k = 0; k != Y.Size(); k++)
		{
			U(k, i) = pImJmat(k, i) - Im1(k) * R(i) + Jm1(k);
			dU_dR(k, i) = Imat(k, i) - Im1(k);
		}
	}

	return;
}

void UNM_trapz(const Vector& R2, const Vector& R1, const Vector& Y, Matrix& U)
{
	Matrix Imata(Y.Size(), R1.Size());
	Matrix Jmata(Y.Size(), R1.Size());
	Matrix Imatb(Y.Size(), R2.Size());
	Matrix Jmatb(Y.Size(), R2.Size());
	Vector Im1(Y.Size());
	Imat_calc(Y, R1, Imata);
	Jmat_calc(Y, R1, Jmata);
	Imat_calc(Y, R2, Imatb);
	Jmat_calc(Y, R2, Jmatb);
	Im1_calc(Y, Im1);

	U = Matrix(Y.Size(), R2.Size());
	for (size_t i = 0; i != R2.Size(); i++)
	{
		for (size_t k = 0; k != Y.Size(); k++)
		{
			U(k, i) = (R2(i) * Imatb(k, i) - Jmatb(k, i)) - (R1(i) * Imata(k, i) - Jmata(k, i)) - Im1(k) * (R2(i)-R1(i));
		}
	}

	return;
}

void UNM_rect(const Vector& R, const Vector& Y, Matrix& U)
{
	Matrix Imat(Y.Size(), R.Size());
	Vector Im1(Y.Size());
	Imat_calc(Y, R, Imat);
	Im1_calc(Y, Im1);

	U = Matrix(Y.Size(), R.Size());
	for (size_t i = 0; i != R.Size(); i++)
	{
		for (size_t k = 0; k != Y.Size(); k++)
		{
			U(k, i) = Imat(k, i) - Im1(k);
		}
	}

	return;
}

void UNM_calc(const Vector& Yw, Matrix& UN, Matrix& UM)
{
	Vector R1(Yw.Size() - 1);
	for (int i = 0; i < Yw.Size() - 1; i++) {
		R1(i) = Yw(i);
	}
	Vector R2(Yw.Size() - 1);
	for (int i = 0; i < Yw.Size() - 1; i++) {
		R2(i) = Yw(i+1);
	}

	Matrix Utr{};
	Matrix Ur{};
	UNM_trapz(R2, R1, Yw, Utr);
	UNM_rect(Yw, Yw, Ur);

	Matrix Ur1(Ur.noRows(), Ur.noCols() -1);
	for (int i = 0; i < Ur.noRows(); i++) {
		for (int j = 0; j < Ur.noCols() - 1; j++) {
			Ur1(i, j) = Ur(i, j);
		}
	}
	Matrix Ur2(Ur.noRows(), Ur.noCols() - 1);
	for (int i = 0; i < Ur.noRows(); i++) {
		for (int j = 0; j < Ur.noCols() - 1; j++) {
			Ur2(i, j) = Ur(i, j+1);
		}
	}

	UN = Matrix(Yw.Size(), Yw.Size() - 1);
	UM = Matrix(Yw.Size(), Yw.Size() - 1);
	for (size_t i = 0; i != Yw.Size()-1; i++)
	{
		for (size_t k = 0; k != Yw.Size(); k++)
		{
			UN(k, i) = 6. * (Yw[i + 1] + Yw[i]) / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) * Utr(k, i) - 2. * (2. * Yw[i + 1] + Yw[i]) / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) * Ur1(k, i) - 2. * (Yw[i + 1] + 2. * Yw[i]) / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) * Ur2(k, i);
			UM(k, i) = -12. / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) * Utr(k, i) + 6. / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) * (Ur1(k, i) + Ur2(k, i));
		}
	}
	return;
}

void UNMb_trapz(const Vector& R2, const Vector& R1, const Vector& Y, Matrix& U)
{
	Matrix Imata(Y.Size(), R1.Size());
	Matrix Jmata(Y.Size(), R1.Size());
	Matrix Imatb(Y.Size(), R2.Size());
	Matrix Jmatb(Y.Size(), R2.Size());
	Vector Im1(Y.Size());
	Imatb_calc(Y, R1, Imata);
	Jmatb_calc(Y, R1, Jmata);
	Imatb_calc(Y, R2, Imatb);
	Jmatb_calc(Y, R2, Jmatb);
	Im1b_calc(Y, Im1);

	U = Matrix(Y.Size(), R2.Size());
	for (size_t i = 0; i != R2.Size(); i++)
	{
		for (size_t k = 0; k != Y.Size(); k++)
		{
			U(k, i) = (R2(i) * Imatb(k, i) - Jmatb(k, i)) - (R1(i) * Imata(k, i) - Jmata(k, i)) - Im1(k) * (R2(i) - R1(i));
		}
	}

	return;
}

void UNMb_rect(const Vector& R, const Vector& Y, Matrix& U)
{
	Matrix Imat(Y.Size(), R.Size());
	Vector Im1(Y.Size());
	Imatb_calc(Y, R, Imat);
	Im1b_calc(Y, Im1);

	U = Matrix(Y.Size(), R.Size());
	for (size_t i = 0; i != R.Size(); i++)
	{
		for (size_t k = 0; k != Y.Size(); k++)
		{
			U(k, i) = Imat(k, i) - Im1(k);
		}
	}

	return;
}

void UNMb_calc(const Vector& Yw, Matrix& UN, Matrix& UM)
{
	Vector R1(Yw.Size() - 1);
	for (int i = 0; i < Yw.Size() - 1; i++) {
		R1(i) = Yw(i);
	}
	Vector R2(Yw.Size() - 1);
	for (int i = 0; i < Yw.Size() - 1; i++) {
		R2(i) = Yw(i + 1);
	}

	Matrix Utr{};
	Matrix Ur{};
	UNMb_trapz(R2, R1, Yw, Utr);
	UNMb_rect(Yw, Yw, Ur);
	//Matrix Ur1 = Ur.leftCols(Yw.Size() - 1);
	//Matrix Ur2 = Ur.rightCols(Yw.Size() - 1);

	Matrix Ur1(Ur.noRows(), Ur.noCols() - 1);
	for (int i = 0; i < Ur.noRows(); i++) {
		for (int j = 0; j < Ur.noCols() - 1; j++) {
			Ur1(i, j) = Ur(i, j);
		}
	}
	Matrix Ur2(Ur.noRows(), Ur.noCols() - 1);
	for (int i = 0; i < Ur.noRows(); i++) {
		for (int j = 0; j < Ur.noCols() - 1; j++) {
			Ur2(i, j) = Ur(i, j + 1);
		}
	}

	UN = Matrix(Yw.Size(), Yw.Size() - 1);
	UM = Matrix(Yw.Size(), Yw.Size() - 1);
	for (size_t i = 0; i != Yw.Size() - 1; i++)
	{
		for (size_t k = 0; k != Yw.Size(); k++)
		{
			UN(k, i) = 6. * (Yw[i + 1] + Yw[i]) / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) * Utr(k, i) - 2. * (2. * Yw[i + 1] + Yw[i]) / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) * Ur1(k, i) - 2. * (Yw[i + 1] + 2. * Yw[i]) / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) * Ur2(k, i);
			UM(k, i) = -12. / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) * Utr(k, i) + 6. / (Yw[i + 1] - Yw[i]) / (Yw[i + 1] - Yw[i]) * (Ur1(k, i) + Ur2(k, i));
		}
	}
	return;
}