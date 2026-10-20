// Code Developed by: Shiye Wang, Tongji University & EPFL (shiyewang@tongji.edu.cn, shiye.wang@epfl.ch, 465317046@qq.com)

#include <math.h>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <CFSTsteel.h>
#include <float.h>
#include <Vector.h>
#include <classTags.h>
#include <string>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <elementAPI.h>
#include <OPS_Globals.h>

static int numCFSTsteel = 0;

void* OPS_CFSTsteel(void)
{
	if (numCFSTsteel == 0) {
		numCFSTsteel++;
		OPS_Error("CFSTsteel uniaxial material model - Code by Shiye Wang", 1);
	}
	UniaxialMaterial* theMaterial = 0;

	const int N_TAGS = 1;
	const int N_YIELD_PROPERTIES = 2;
	const int N_LOCALBUCKLING_PROPERTIES = 5;
	const int N_ENERGY_PROPERTIES = 5;
	const int N_DETERIORATION_PROPERTIES = 5;
	const int N_ISO_PROPERTIES = 2;
	const int N_PARAM_PER_BACK = 2;
	const int MAX_BACKSTRESSES = 8;
	const int BACKSTRESS_SPACE = MAX_BACKSTRESSES * N_PARAM_PER_BACK;

	const char* inputInstructions =
        "Want: uniaxialMaterial CFSTsteel tag E0 fy b2 sigma_cm a sigma_rm b3 "
        "Lambda_t Lambda_cp Lambda_a Lambda_yp Lambda_r "
        "c_t c_cp c_a c_yp c_r Qinf b N <C1 gamma1 ... CN gammaN>";
    if (OPS_GetNumRemainingInputArgs() < 21) {
        opserr << "WARNING CFSTsteel: insufficient arguments\n" << inputInstructions << endln;
        return 0;
    }

	int nInputsToRead;
	int nBackstresses[1];
	int materialTag[1];
	double yiProps[2];
	double lbProps[5];
	double enProps[5];
	double deProps[5];
	double isProps[2];
	double backstressProps[BACKSTRESS_SPACE];
	std::vector<double> cK;
	std::vector<double> gammaK;

	nInputsToRead = N_TAGS;
	if (OPS_GetIntInput(&nInputsToRead, materialTag) != 0) {
		opserr << "WARNING invalid uniaxialMaterial CFSTsteel tag" << endln;
		return 0;
	}

	nInputsToRead = N_YIELD_PROPERTIES;
	if (OPS_GetDoubleInput(&nInputsToRead, yiProps) != 0) {
		opserr << inputInstructions << endln;
		return 0;
	}

	nInputsToRead = N_LOCALBUCKLING_PROPERTIES;
	if (OPS_GetDoubleInput(&nInputsToRead, lbProps) != 0) {
		opserr << inputInstructions << endln;
		return 0;
	}

	nInputsToRead = N_ENERGY_PROPERTIES;
	if (OPS_GetDoubleInput(&nInputsToRead, enProps) != 0) {
		opserr << inputInstructions << endln;
		return 0;
	}

	nInputsToRead = N_DETERIORATION_PROPERTIES;
	if (OPS_GetDoubleInput(&nInputsToRead, deProps) != 0) {
		opserr << inputInstructions << endln;
		return 0;
	}

	nInputsToRead = N_ISO_PROPERTIES;
	if (OPS_GetDoubleInput(&nInputsToRead, isProps) != 0) {
		opserr << inputInstructions << endln;
		return 0;
	}

	nInputsToRead = 1;
	if (OPS_GetIntInput(&nInputsToRead, nBackstresses) != 0) {
		opserr << "WARNING N must be an integer" <<
			inputInstructions << endln;
		return 0;
	}
	if (nBackstresses[0] < 0 || nBackstresses[0] > MAX_BACKSTRESSES) {
		opserr << "WARNING CFSTsteel: N must be between zero and " <<
			MAX_BACKSTRESSES << endln <<
			inputInstructions << endln;
		return 0;
	}

	nInputsToRead = 2 * nBackstresses[0];
	if (OPS_GetNumRemainingInputArgs() != nInputsToRead ||
        (nInputsToRead > 0 && OPS_GetDoubleInput(&nInputsToRead, backstressProps) != 0)) {
		opserr << inputInstructions << endln;
		return 0;
	}

    for (double v : yiProps) if (!std::isfinite(v)) return 0;
    for (double v : lbProps) if (!std::isfinite(v)) return 0;
    for (double v : enProps) if (!std::isfinite(v)) return 0;
    for (double v : deProps) if (!std::isfinite(v)) return 0;
    for (double v : isProps) if (!std::isfinite(v)) return 0;
    bool valid = yiProps[0] > 0.0 && yiProps[1] > 0.0 &&
        lbProps[0] > 0.0 && lbProps[1] < 0.0 && lbProps[2] > 0.0 &&
        lbProps[3] > lbProps[1] && lbProps[3] <= 0.0 && lbProps[4] > 0.0 &&
        isProps[0] >= 0.0 && isProps[1] >= 0.0;
    for (double v : enProps) valid = valid && v > 0.0;
    for (double v : deProps) valid = valid && v > 0.0;
    for (int i = 0; i < nBackstresses[0]; ++i)
        valid = valid && std::isfinite(backstressProps[2*i]) &&
            std::isfinite(backstressProps[2*i+1]) &&
            backstressProps[2*i] >= 0.0 && backstressProps[2*i+1] > 0.0;
    if (!valid) {
        opserr << "WARNING CFSTsteel: invalid parameter range\n" << inputInstructions << endln;
        return 0;
    }

	for (int i = 0; i < nBackstresses[0]; ++i) {
		cK.push_back(backstressProps[2 * i]);
		gammaK.push_back(backstressProps[1 + 2 * i]);
	}

	theMaterial = new CFSTsteel(materialTag[0],
		yiProps[0], yiProps[1],
		lbProps[0], lbProps[1], lbProps[2], lbProps[3], lbProps[4],
		enProps[0], enProps[1], enProps[2], enProps[3], enProps[4],
		deProps[0], deProps[1], deProps[2], deProps[3], deProps[4],
		isProps[0], isProps[1],
		cK, gammaK);

	return theMaterial;
}

CFSTsteel::CFSTsteel(int tag,
	double E0, double Fy,
	double b2, double siglb, double a, double sigr, double b3,
	double Lambda_y, double Lambda_lb, double Lambda_a, double Lambda_yneg, double Lambda_re,
	double C_y, double C_lb, double C_a, double C_yneg, double C_re,
	double qInf, double b,
	std::vector<double> cK, std::vector<double> gammaK)
	: UniaxialMaterial(tag, MAT_TAG_CFSTsteel),
	E0(E0), Fy(Fy),
	b2(b2), siglb(siglb), a(a), sigr(sigr), b3(b3),
	Lambda_y(Lambda_y), Lambda_lb(Lambda_lb), Lambda_a(Lambda_a), Lambda_yneg(Lambda_yneg), Lambda_re(Lambda_re),
	C_y(C_y), C_lb(C_lb), C_a(C_a), C_yneg(C_yneg), C_re(C_re),
	qInf(qInf), b(b),
	cK(cK), gammaK(gammaK),
	strainConverged(0.), strainTrial(0.), strainPEqConverged(0.), strainPEqTrial(0.), stressConverged(0.), stressTrial(0.),
	stiffnessInitial(E0), stiffnessConverged(E0), stiffnessTrial(E0), flowDirection(0.), plasticLoading(false)
{
	nBackstresses = cK.size();
	for (int i = 0; i < nBackstresses; ++i) {
		alphaKTrial.push_back(0.);
		alphaKConverged.push_back(0.);
	}
	this->revertToStart();
}

CFSTsteel::CFSTsteel() :
    CFSTsteel(0, 1.0, 1.0, 0.01, -1.1, 1.0, -0.1, 1.0,
              1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
              0.0, 1.0, std::vector<double>(), std::vector<double>())
{
}

CFSTsteel::~CFSTsteel(void)
{

}

double
CFSTsteel::getInitialTangent(void)
{
	return E0;
}

int CFSTsteel::setTrialStrain(double trialStrain, double strainRate)
{
    if (!std::isfinite(trialStrain)) return -1;
	this->revertToLastCommit();
	eps = trialStrain;
	double deps = eps - eps_1;
	double strainIncrement = eps - strainConverged;

	double Eb_pos_j = Eb_pos_j_1;
	double sigy_project_pos_j = sigy_project_pos_j_1;
	double siglb_project_j = siglb_project_j_1;
	double Eb_neg_j = Eb_neg_j_1;
	double sigy_project_neg_j = sigy_project_neg_j_1;
	double epslb_j = epslb_j_1;
	double siglb_j = siglb_j_1;
	double a_j = a_j_1;
	double sigr_j = sigr_j_1;
	double E_r_j = E_r_j_1;

	if (eps >= eps_1) {
		Di = 1;
	}
	else {
		Di = -1;
	}

	if ((Di_1 / Di < 0)) {
		Reversal_Flag = 1.0;
		epsreversal = eps_1;
		sigreversal = sig_1;
	}
	else {
		Reversal_Flag = 0.0;
	}

	if (Energy_Flag == 0.0) {
	if (Excursion_Flag == 1.0) {
		if (eps - eps_1 >= 0.0) {

			beta_y_j = pow((Energy_Excrsn / (Ref_Energy_y - Energy_total)), C_y);
			Eb_pos_j = Eb_pos_j_1 * (1.0 - beta_y_j);
		}
		else {

			beta_lb_j = pow((Energy_Excrsn / (Ref_Energy_lb - Energy_total)), C_lb);
			beta_a_j = pow((Energy_Excrsn / (Ref_Energy_a - Energy_total)), C_a);
			beta_yneg_j = pow((Energy_Excrsn / (Ref_Energy_yneg - Energy_total)), C_yneg);
			beta_re_j = pow((Energy_Excrsn / (Ref_Energy_re - Energy_total)), C_re);
			sigy_project_neg_j = sigy_project_neg_j_1 * (1.0 - beta_yneg_j);
			Eb_neg_j = Eb_neg_j_1 * (1.0 - beta_yneg_j);
			siglb_project_j = siglb_project_j_1 * (1.0 - beta_lb_j);
			sigr_j = sigr_j_1 * (1.0 - beta_lb_j);
			a_j = a_j_1 * (1.0 - beta_a_j);
			E_r_j = E_r_j_1 * (1.0 - beta_re_j);
		}
	}
	else {
		if (eps - eps_1 >= 0.0) {
			beta_y_j = beta_y_j_1;
			Eb_pos_j = Eb_pos_j_1;
		}
		else {
			beta_lb_j = beta_lb_j_1;
			beta_a_j = beta_a_j_1;
			beta_yneg_j = beta_yneg_j_1;
			a_j = a_j_1;
			sigy_project_neg_j = sigy_project_neg_j_1;
			Eb_neg_j = Eb_neg_j_1;
			siglb_project_j = siglb_project_j_1;
			sigr_j = sigr_j_1;
			E_r_j = E_r_j_1;
		}
	}
	}

	if (Energy_Flag == 1) {

		beta_y_j = beta_y_j_1;
		beta_lb_j = beta_lb_j_1;
		beta_a_j = beta_a_j_1;
		beta_yneg_j = beta_yneg_j_1;
		beta_re_j = beta_re_j_1;
	}

	if ((fabs(sig) >= Fy || fabs(eps) > fabs(epsy_neg)) && Yield_Flag == 0.0) {
		Yield_Flag = 1.0;
	}

	if ((eps_max_Flag == 1.0 && deps < 0.0 && Yield_Flag == 1.0) || (Buckling_flag == 0.0 && Reversal_Flag == 1.0 && deps < 0.0 && lbstage == 2.0 && reloading_Flag != 1.0)) {
		epsdiatance = epsreversal - sigreversal / E0;
	}

	if (lbstage == 1.0 && Yield_Flag == 1.0 && Reversal_Flag == 1.0 && deps < 0.0 && reloading_Flag != 1.0) {
		epsdiatance = epsreversal - sigreversal / E0;
	}

	if (Yield_Flag != 0.0) {
		if ((reloading_Flag == 1.0) && (deps / deps_re > 0.0) && (kon == 2)) {
			kon = 6.0;
		}
	}

	if ((eps_max_Flag == 1.0) && (Reversal_Flag == 1.0) && (eps_max > 0.02)) {
		if (fabs(sig) < 0.1 * Fy) {
			Damage_Flag = 1.0;
		}
	}

	if (kon == 0 || kon == 10) {
		if (fabs(deps) < 10.0 * DBL_EPSILON)
		{
			e = E0;
			kon = 10;
			return 0;
		}
		else
		{
			if (deps < 0.0)
			{
				kon = 2;
			}
			else {
				kon = 1;
			}
		}
	}

	if ((kon == 1) && (deps > 0.0)) {
		kon = 1;
		returnMapping(strainIncrement);
		calculateStiffness();
		sig = stressTrial;
		e = stiffnessTrial;
	}

	if ((kon == 2) && (deps > 0.0) && (Yield_Flag == 0.0)) {
		returnMapping(strainIncrement);
		calculateStiffness();
		if (eps < 0.0) {
			kon = 2;
			sig = sig_1 + deps * E0;
			e = E0;
		}
		else {
			kon = 1;
			sig = stressTrial;
			e = stiffnessTrial;
		}
	}

	if ((kon == 1 || kon == 2) && (deps < 0.0) && (Yield_Flag == 0.0)) {
		returnMapping(strainIncrement);
		calculateStiffness();
		if (eps > 0.0) {
			kon = 1;
			sig = stressTrial;
			e = stiffnessTrial;
		}
		else {
			kon = 2;
			sig = sig_1 + deps * E0;
			e = E0;
		}
	}

	if ((kon == 2 || kon == 1 || kon == 5 || kon == 6) && (deps < 0.0) && (Yield_Flag == 1.0)) {
		kon = 2;
		returnMapping(strainIncrement);
		calculateStiffness();
		sig_Trial = sig_1 + deps * E0;
		epsy_neg = sigy_project_neg_j / (E0 - Eb_neg_j);
		double sig_eps_y = sigr_j + siglb_project_j * exp(a_j * (epsy_neg));
		if (sig_eps_y <= (E0 * epsy_neg)) {
			double JD_1 = sigy_project_neg_j / (E0 - Eb_neg_j);
			double JD_2 = Newtoniteration(&CFSTsteel::fxfunction1, &CFSTsteel::fdxfunction1, 0.0);
			lbstrength = sigr_j + siglb_project_j * exp(a_j * JD_2);
			if (((JD_1 + epsdiatance) < eps) && (eps <= epsreversal)) {
				lbstage = 1.0;

				double sig_Elastic = E0 * (eps - epsreversal) + sigreversal;
				double sig_Hardening = Eb_neg_j * (eps - epsdiatance) + sigy_project_neg_j;
				if (sig_Elastic <= sig_Hardening) {
					sig = sig_Hardening;
					e = Eb_neg_j;
				}
				else {
					sig = sig_Elastic;
					e = E0;
				}

			}
			if (((JD_2 + epsdiatance) < eps) && (eps <= (JD_1 + epsdiatance))) {
				lbstage = 1.0;
				sig = Eb_neg_j * (eps - epsdiatance) + sigy_project_neg_j;
				e = Eb_neg_j;
			}
			if (eps <= (JD_2 + epsdiatance)) {
				lbstage = 2.0;
				sig = sigr_j + siglb_project_j * exp(a_j * (eps - epsdiatance));
				e = (a_j)*siglb_project_j * exp(a_j * (eps - epsdiatance));
			}
			curve_epslb = JD_2 + epsdiatance;
			siglb_j = Eb_neg_j * JD_2 + sigy_project_neg_j;
			epslb_j = JD_2 + epsdiatance;
		}
		else {
			double JD_1 = Newtoniteration(&CFSTsteel::fxfunction2, &CFSTsteel::fdxfunction2, 0.0);
			lbstrength = sigr_j + siglb_project_j * exp(a_j * JD_1);
			if (((JD_1 + epsdiatance) < eps) && (eps <= epsreversal)) {
				lbstage = 1.0;

				double sig_Elastic = E0 * (eps - epsreversal) + sigreversal;
				double sig_Softening = sigr_j + siglb_project_j * exp(a_j * (eps - epsdiatance));
				if (sig_Elastic <= sig_Softening) {
					sig = sig_Softening;
					e = (a_j)*siglb_project_j * exp(a_j * (eps - epsdiatance));
				}
				else {
					sig = sig_Elastic;
					e = E0;
				}

			}
			if (eps <= (JD_1 + epsdiatance)) {
				lbstage = 2.0;
				sig = sigr_j + siglb_project_j * exp(a_j * (eps - epsdiatance));
				e = (a_j)*siglb_project_j * exp(a_j * (eps - epsdiatance));
			}
			curve_epslb = JD_1 + epsdiatance;
			siglb_j = sigr_j + siglb_project_j * exp(a_j * JD_1);
			epslb_j = JD_1 + epsdiatance;
		}
		if (sig >= sig_Trial) {
			sig = sig;
			e = e;
		}
		else {
			sig = sig_Trial;
			e = E0;
		}
	}

	if ((kon == 5 || kon == 2) && (deps > 0.0) && (Yield_Flag == 1.0)) {
		kon = 5;
		sig_Trial = E_r_j * (eps - epsreversal) + sigreversal;
		returnMapping(strainIncrement);
		calculateStiffness();
		double JD_3 = Newtoniteration(&CFSTsteel::fxfunction3, &CFSTsteel::fdxfunction3, 0.0);
		if (lbstage == 1.0) {
			if (epsreversal <= eps) {
				if (sig_Trial < sig_iso) {
					sig = sig_Trial;
					e = E_r_j;
				}
				else {
					sig = sig_iso;
					e = stiffnessTrial;
				}
			}
		}
		else if (lbstage == 2.0) {
			if ((epsreversal <= eps) && (eps < JD_3)) {
				sig = E_r_j * (eps - epsreversal) + sigreversal;
				e = E_r_j;
			}
			if (JD_3 <= eps) {
				sig = sig_1 + deps * (-b3 * a_j * siglb_project_j * exp(a_j * (eps - epsdiatance)));
				e = -b3 * a_j * siglb_project_j * exp(a_j * (eps - epsdiatance));
			}
			if (sig < sig_iso) {
				sig = sig;
				e = e;
			}
			else {
				sig = sig_iso;
				e = stiffnessTrial;
			}
		}
	}

	if ((kon == 6) && (deps > 0.0) && (Yield_Flag == 1.0)) {
		kon = 6;
		returnMapping(strainIncrement);
		calculateStiffness();
		sig_Trial = E0 * (eps - epsreversal) + sigreversal;
		if (sig_Trial < sig_re) {
			sig = sig_Trial;
			e = E0;
		}
		else {
			sig = sig_re;
			e = e_re;
			epsreversal = epsreversal_re;
			sigreversal = sigreversal_re;
			Minus_Flag = Minus_Flag_re;
			kon = kon_re;
			lbstage = lbstage_re;
			kon = 5;
		}
	}

	Energy_total = Energy_total + (sig + sig_1) * 0.5 * (eps - eps_1);

	if ((sig / sig_1) < 0.0) {
		Energy_Excrsn = Energy_total - Energy_Excrsni_1;
		Energy_Excrsni_1 = Energy_total;
		Excursion_Flag = 1;
	}
	else {
		Excursion_Flag = 0.0;
	}

	if (Excursion_Flag == 1.0) {
		if (eps_max >= 0.08) {
			Energy_Flag = 1.0;
		}
		if (Damage_Flag == 1.0) {
			Energy_Flag = 1.0;
		}
	}

	if ((Energy_total >= Ref_Energy_y) || (Energy_total >= Ref_Energy_lb) || (Energy_total >= Ref_Energy_a) || (Energy_total >= Ref_Energy_yneg) || (Energy_total >= Ref_Energy_re)) {
		Energy_Flag = 1.0;
	}
	if ((beta_y_j > 1) || (beta_lb_j > 1) || (beta_a_j > 1) || (beta_yneg_j > 1) || (beta_re_j > 1)) {
		Energy_Flag = 1.0;
	}

	if ((deps < 0) && lbstage == 2.0) {
		Bucklingdegree = eps - curve_epslb;
	}
	if ((deps > 0) && lbstage == 2.0) {
		Bucklingdegree = Bucklingdegree + deps;
	}
	if (Bucklingdegree > 0) {
		Bucklingdegree = 0;
	}
	if (deps > 0) {
		if (Bucklingdegree == 0) {
			Buckling_flag = 0;
		}
		else {
			Buckling_flag = 1;
		}
	}

	if (Energy_Flag == 1) {

	}

	if (eps_max > eps) {
		eps_max = eps_max;
		eps_max_Flag = 0.0;
	}
	else {
		eps_max = eps;
		eps_max_Flag = 1.0;
	}

	if (Yield_Flag != 0.0) {
		if ((e != E0) && (e != E_r_j)) {
			kon_re = kon;
			eps_re = eps;
			sig_re = sig;
			e_re = e;
			deps_re = eps - eps_1;
			epsreversal_re = epsreversal;
			sigreversal_re = sigreversal;
			lbstage_re = lbstage;
			Minus_Flag_re = Minus_Flag;
		}
		if ((e == E0) || (e == E_r_j)) {
			reloading_Flag = 1.0;
		}
		if (fabs(e / e_1) <= 0.9) {
			reloading_Flag = 0.0;
		}
	}

	Eb_pos_j_1 = Eb_pos_j;
	sigy_project_pos_j_1 = sigy_project_pos_j;
	Eb_neg_j_1 = Eb_neg_j;
	E_r_j_1 = E_r_j;
	sigy_project_neg_j_1 = sigy_project_neg_j;
	epslb_j_1 = epslb_j;
	siglb_project_j_1 = siglb_project_j;
	siglb_j_1 = siglb_j;
	a_j_1 = a_j;
	eps_2 = eps_1;
	eps_1 = eps;
	sig_2 = sig_1;
	sig_1 = sig;
	e_2 = e_1;
	e_1 = e;
	Di_2 = Di_1;
	Di_1 = Di;
	kon_2 = kon_1;
	kon_1 = kon;
	Minus_Flag_2 = Minus_Flag_1;
	Minus_Flag_1 = Minus_Flag;

	beta_y_j_1 = beta_y_j;
	beta_lb_j_1 = beta_lb_j;
	beta_a_j_1 = beta_a_j;
	beta_yneg_j_1 = beta_yneg_j;
	beta_re_j_1 = beta_re_j;
	Minus_1 = Minus;

	return 0;
}

double CFSTsteel::fxfunction1(double x)
{
	return ((sigr_j_1 + siglb_project_j_1 * exp(a_j_1 * x)) - (Eb_neg_j_1 * x + sigy_project_neg_j_1));
}
double CFSTsteel::fdxfunction1(double x)
{
	return (a_j_1 * siglb_project_j_1 * exp(a_j_1 * x) - Eb_neg_j_1);
}
double CFSTsteel::fxfunction2(double x)
{
	return ((sigr_j_1 + siglb_project_j_1 * exp(a_j_1 * x)) - E0 * x);
}
double CFSTsteel::fdxfunction2(double x)
{
	return (a_j_1 * siglb_project_j_1 * exp(a_j_1 * x) - E0);
}
double CFSTsteel::fxfunction3(double x)
{
	return (-0.7 * (sigr_j_1 + siglb_project_j_1 * exp(a_j_1 * (x - epsdiatance))) - (E_r_j_1 * (x - epsreversal) + sigreversal));
}
double CFSTsteel::fdxfunction3(double x)
{
	return (-0.7 * (a_j_1)*siglb_project_j_1 * exp(a_j_1 * (x - epsdiatance)) - E_r_j_1);
}

void CFSTsteel::returnMapping(double strainIncrement) {

	bool converged = true;
	int iterationNumber = 0;
	double sigmaY1 = 0.;
	double dit = 0.;
	double plasticStrainIncrement = 0.;
	double aux = 0.;
	double alpha = 0.;
	double sy = 0.;
	double stressRadius = 0.;
	double phi = 0.;
	double ePEq = strainPEqConverged;

	for (int i = 0; i < nBackstresses; ++i) {
		alpha += alphaKConverged[i];
	}
	sigmaY1 = qInf * (1. - exp(-b * ePEq));
	sy = Fy + sigmaY1;
	stressTrial = stressConverged + E0 * strainIncrement;
	stressRadius = stressTrial - alpha;
	phi = pow(stressRadius, 2) - pow(sy, 2);

	if (phi > RETURN_MAP_TOL) {
		converged = false;
	}
	while ((!converged) && (iterationNumber < MAXIMUM_ITERATIONS)) {
		iterationNumber++;

		aux = E0;
		for (int i = 0; i < nBackstresses; ++i) {
			aux = aux + sgn<double>(stressRadius) * cK[i] - gammaK[i] * alphaKTrial[i];
		}

		dit = 2. * stressRadius * aux +
			2. * sy * qInf * b * exp(-b * ePEq);
		plasticStrainIncrement = phi / dit;

		if (abs(plasticStrainIncrement) > abs(stressTrial / E0)) {
			plasticStrainIncrement = sgn<double>(plasticStrainIncrement) * 0.95 * abs(stressTrial / E0);
		}

		ePEq = ePEq + abs(plasticStrainIncrement);
		stressTrial = stressTrial - E0 * plasticStrainIncrement;
		sigmaY1 = qInf * (1. - exp(-b * ePEq));
		sy = Fy + sigmaY1;

		alpha = 0.;
		for (int i = 0; i < nBackstresses; ++i) {
			alphaKTrial[i] = sgn<double>(stressRadius) * cK[i] / gammaK[i] -
				(sgn<double>(stressRadius) * cK[i] / gammaK[i] - alphaKConverged[i]) * exp(-gammaK[i] * (ePEq - strainPEqConverged));
			alpha += alphaKTrial[i];
		}

		stressRadius = stressTrial - alpha;
		phi = pow(stressRadius, 2) - pow(sy, 2);
		if (abs(phi) < RETURN_MAP_TOL) {
			converged = true;
		}
	}
	sig_iso = sy + alpha;

	if (!converged) {
		opserr << "WARNING: return mapping in CFSTsteel does not converge!" << endln;
		opserr << "\tStrain increment = " << strainIncrement << endln;
		opserr << "\tExiting with phi = " << phi << " > " << RETURN_MAP_TOL << endln;
	}

	if (iterationNumber == 0) {
		plasticLoading = false;
	}
	else {
		plasticLoading = true;
	}

	flowDirection = sgn<double>(stressRadius);
	strainPEqTrial = ePEq;
	return;
}

void CFSTsteel::calculateStiffness() {
	if (!plasticLoading) {
		stiffnessTrial = E0;
	}
	else {
		double sigmaY1 = qInf * (1. - exp(-b * strainPEqTrial));
		double plasticModulus = 0.;

		plasticModulus = b * (qInf - sigmaY1);
		for (int i = 0; i < nBackstresses; ++i) {
			plasticModulus += gammaK[i] *
				(cK[i] / gammaK[i] - flowDirection * alphaKTrial[i]);
		}
		stiffnessTrial = (E0 * plasticModulus) /
			(E0 + plasticModulus);
	}
	return;
}

double CFSTsteel::getStrain(void)
{
	return eps;
}

double CFSTsteel::getStress(void)
{
	return sig;
}

double CFSTsteel::getTangent(void)
{
	return e;
}

int CFSTsteel::commitState(void)
{
	cEb_pos_j_1 = Eb_pos_j_1;
	csigy_project_pos_j_1 = sigy_project_pos_j_1;
	cEb_neg_j_1 = Eb_neg_j_1;
	cE_r_j_1 = E_r_j_1;
	csigy_project_neg_j_1 = sigy_project_neg_j_1;
	cepslb_j_1 = epslb_j_1;
	csiglb_j_1 = siglb_j_1;
	ca_j_1 = a_j_1;
	csigr_j_1 = sigr_j_1;
	cslope_lb_j_1 = slope_lb_j_1;
	cslope_r_j_1 = slope_r_j_1;
	csiglb_project_j_1 = siglb_project_j_1;
	csigr_project_j_1 = sigr_project_j_1;
	cDi = Di;
	cDi_1 = Di_1;
	cDi_2 = Di_2;
	cepsreversal = epsreversal;
	csigreversal = sigreversal;
	cepsreversal_1 = epsreversal_1;
	csigreversal_1 = sigreversal_1;
	cepsreversal_2 = epsreversal_2;
	csigreversal_2 = sigreversal_2;
	csig = sig;
	csig_1 = sig_1;
	csig_2 = sig_2;
	ceps = eps;
	ceps_1 = eps_1;
	ceps_2 = eps_2;
	ckon = kon;
	ckon_1 = kon_1;
	ckon_2 = kon_2;
	clbstage = lbstage;
	cbeta_y_j_1 = beta_y_j_1;
	cbeta_lb_j_1 = beta_lb_j_1;
	cbeta_a_j_1 = beta_a_j_1;
	cbeta_yneg_j_1 = beta_yneg_j_1;
	cbeta_re_j_1 = beta_re_j_1;
	ccurve_epslb = curve_epslb;
	ccurve_epsr = curve_epsr;
	cbeta_y_j = beta_y_j;
	cbeta_lb_j = beta_lb_j;
	cbeta_a_j = beta_a_j;
	cbeta_yneg_j = beta_yneg_j;
	cbeta_re_j = beta_re_j;

	cExcursion_Flag = Excursion_Flag;
	cReversal_Flag = Reversal_Flag;
	cEnergy_Excrsni_1 = Energy_Excrsni_1;
	cEnergy_Flag = Energy_Flag;
	cYield_Flag = Yield_Flag;
	cBucklingdegree = Bucklingdegree;
	cBuckling_flag = Buckling_flag;
	cStiffness_neg_Flag = Stiffness_neg_Flag;
	cEnergy_Excrsn = Energy_Excrsn;
	cEnergy_total = Energy_total;
	ce = e;
	ce_1 = e_1;
	ce_2 = e_2;

	ceps_inflection = eps_inflection;
	csig_inflection = sig_inflection;
	ce_inflection = e_inflection;
	cepsdiatance = epsdiatance;

	cMinus = Minus;
	cMinus_1 = Minus_1;
	cMinus_Flag = Minus_Flag;
	cMinus_Flag_1 = Minus_Flag_1;
	cMinus_Flag_2 = Minus_Flag_2;
	cStiffness_neg = Stiffness_neg;
	csig_iso = sig_iso;
	csig_lb_re = sig_lb_re;

	ckon_re = kon_re;
	ceps_re = eps_re;
	csig_re = sig_re;
	ce_re = e_re;
	cdeps_re = deps_re;
	cepsreversal_re = epsreversal_re;
	csigreversal_re = sigreversal_re;
	creloading_Flag = reloading_Flag;

	csig_Trial = sig_Trial;
	ceps_max = eps_max;
	ceps_max_Flag = eps_max_Flag;
	clbstage_re = lbstage_re;
	cMinus_Flag_re = Minus_Flag_re;
	csig_boundary = sig_boundary;
	calpha_boundary = alpha_boundary;
	cDamage_Flag = Damage_Flag;
	clbstrength = lbstrength;

	strainConverged = eps;
	strainPEqConverged = strainPEqTrial;
	stressConverged = stressTrial;
	alphaKConverged = alphaKTrial;
	stiffnessConverged = stiffnessTrial;

	return 0;
}

int CFSTsteel::revertToLastCommit(void)
{
	Eb_pos_j_1 = cEb_pos_j_1;
	sigy_project_pos_j_1 = csigy_project_pos_j_1;
	Eb_neg_j_1 = cEb_neg_j_1;
	E_r_j_1 = cE_r_j_1;
	sigy_project_neg_j_1 = csigy_project_neg_j_1;
	epslb_j_1 = cepslb_j_1;
	siglb_j_1 = csiglb_j_1;
	a_j_1 = ca_j_1;
	sigr_j_1 = csigr_j_1;
	slope_lb_j_1 = cslope_lb_j_1;
	slope_r_j_1 = cslope_r_j_1;
	siglb_project_j_1 = csiglb_project_j_1;
	sigr_project_j_1 = csigr_project_j_1;
	Di = cDi;
	Di_1 = cDi_1;
	Di_2 = cDi_2;
	epsreversal = cepsreversal;
	sigreversal = csigreversal;
	epsreversal_1 = cepsreversal_1;
	sigreversal_1 = csigreversal_1;
	epsreversal_2 = cepsreversal_2;
	sigreversal_2 = csigreversal_2;
	sig = csig;
	sig_1 = csig_1;
	sig_2 = csig_2;
	eps = ceps;
	eps_1 = ceps_1;
	eps_2 = ceps_2;
	kon = ckon;
	kon_1 = ckon_1;
	kon_2 = ckon_2;
	lbstage = clbstage;
	beta_y_j_1 = cbeta_y_j_1;
	beta_lb_j_1 = cbeta_lb_j_1;
	beta_a_j_1 = cbeta_a_j_1;
	beta_yneg_j_1 = cbeta_yneg_j_1;
	beta_re_j_1 = cbeta_re_j_1;
	beta_y_j = cbeta_y_j;
	beta_lb_j = cbeta_lb_j;
	beta_a_j = cbeta_a_j;
	beta_yneg_j = cbeta_yneg_j;
	beta_re_j = cbeta_re_j;
	curve_epslb = ccurve_epslb;
	sig_iso = csig_iso;
	sig_lb_re = csig_lb_re;

	Excursion_Flag = cExcursion_Flag;
	Reversal_Flag = cReversal_Flag;
	Energy_Excrsni_1 = cEnergy_Excrsni_1;
	Energy_Flag = cEnergy_Flag;
	Yield_Flag = cYield_Flag;
	Buckling_flag = cBuckling_flag;
	Stiffness_neg_Flag = cStiffness_neg_Flag;
	Bucklingdegree = cBucklingdegree;

	eps_inflection = ceps_inflection;
	sig_inflection = csig_inflection;
	e_inflection = ce_inflection;
	epsdiatance = cepsdiatance;

	Minus = cMinus;
	Minus_1 = cMinus_1;
	Minus_Flag = cMinus_Flag;
	Minus_Flag_1 = cMinus_Flag_1;
	Minus_Flag_2 = cMinus_Flag_2;
	Stiffness_neg = cStiffness_neg;

	kon_re = ckon_re;
	eps_re = ceps_re;
	sig_re = csig_re;
	e_re = ce_re;
	deps_re = cdeps_re;
	epsreversal_re = cepsreversal_re;
	sigreversal_re = csigreversal_re;
	reloading_Flag = creloading_Flag;

	sig_Trial = csig_Trial;
	eps_max = ceps_max;
	eps_max_Flag = ceps_max_Flag;
	Minus_Flag_re = cMinus_Flag_re;
	lbstage_re = clbstage_re;
	lbstage_re = clbstage_re;
	alpha_boundary = calpha_boundary;
	Damage_Flag = cDamage_Flag;
	lbstrength = clbstrength;

	Energy_Excrsn = cEnergy_Excrsn;
	Energy_total = cEnergy_total;

	e = ce;
	e_1 = ce_1;
	e_2 = ce_2;

	eps = strainConverged;
    curve_epsr = ccurve_epsr;
    sig_boundary = csig_boundary;
	strainPEqTrial = strainPEqConverged;
	stressTrial = stressConverged;
	alphaKTrial = alphaKConverged;
	stiffnessTrial = stiffnessConverged;

	return 0;
}

int CFSTsteel::revertToStart(void)
{
    slope_lb = slope_r = sigr_project = 0.0;
    siglb_project_j = csiglb_project_j = exp_point = cexp_point = 0.0;
    curve_epsr = ccurve_epsr = 0.0;
    Minus_Flag_re = cMinus_Flag_re = 0.0;
    strainTrial = strainPEqTrial = stressTrial = flowDirection = 0.0;
    stiffnessTrial = E0;
    plasticLoading = false;

	Eb_pos = b * E0;
	epsy_pos = Fy / E0;
	sigy_project_pos = Fy - Eb_pos * epsy_pos;

	Eb_neg = b2 * E0;
	epsy_neg = -Fy / E0;
	sigy_project_neg = -Fy - Eb_neg * epsy_neg;

	epslb = (siglb - sigy_project_neg) / Eb_neg;
	siglb_project = (siglb - sigr) / exp(a * (epslb));

	if (fabs(siglb) <= fabs(Fy)) {
		epslb = siglb / E0;
	}
	else {
		epslb = epsy_neg + (siglb + Fy) / Eb_neg;
	}

	Ref_Energy_y = Lambda_y * Fy;
	Ref_Energy_lb = Lambda_lb * Fy;
	Ref_Energy_a = Lambda_a * Fy;
	Ref_Energy_yneg = Lambda_yneg * Fy;
	Ref_Energy_re = Lambda_re * Fy;

	Eb_pos_j_1 = Eb_pos;
	sigy_project_pos_j_1 = sigy_project_pos;
	Eb_neg_j_1 = Eb_neg;
	E_r_j_1 = E0;
	sigy_project_neg_j_1 = sigy_project_neg;
	epslb_j_1 = epslb;
	siglb_j_1 = siglb;
	a_j_1 = a;
	sigr_j_1 = sigr;
	slope_lb_j_1 = slope_lb;
	slope_r_j_1 = slope_r;
	siglb_project_j_1 = siglb_project;
	sigr_project_j_1 = sigr_project;

	cEb_pos_j_1 = Eb_pos;
	csigy_project_pos_j_1 = sigy_project_pos;
	cEb_neg_j_1 = Eb_neg;
	cE_r_j_1 = E0;
	csigy_project_neg_j_1 = sigy_project_neg;
	cepslb_j_1 = epslb;
	csiglb_j_1 = siglb;
	ca_j_1 = a;
	csigr_j_1 = sigr;
	cslope_lb_j_1 = slope_lb;
	cslope_r_j_1 = slope_r;
	csiglb_project_j_1 = siglb_project;
	csigr_project_j_1 = sigr_project;

	Di = 0;
	cDi = 0;
	Di_1 = 0;
	cDi_1 = 0;
	Di_2 = 0;
	cDi_2 = 0;
	epsreversal = cepsreversal = 0;
	sigreversal = csigreversal = 0;
	epsreversal_1 = cepsreversal_1 = 0;
	sigreversal_1 = csigreversal_1 = 0;
	epsreversal_2 = cepsreversal_2 = 0;
	sigreversal_2 = csigreversal_2 = 0;
	sig = 0;
	csig = 0;
	sig_1 = 0;
	csig_1 = 0;
	sig_2 = 0;
	csig_2 = 0;
	eps = 0;
	ceps = 0;
	eps_1 = 0;
	ceps_1 = 0;
	eps_2 = 0;
	ceps_2 = 0;
	kon = ckon = 0;
	kon_1 = ckon_1 = 0;
	kon_2 = ckon_2 = 0;
	lbstage = clbstage = 1;
	beta_y_j_1 = cbeta_y_j_1 = 0;
	beta_lb_j_1 = cbeta_lb_j_1 = 0;
	beta_a_j_1 = cbeta_a_j_1 = 0;
	beta_yneg_j_1 = cbeta_yneg_j_1 = 0;
	beta_re_j_1 = cbeta_re_j_1 = 0;
	beta_y_j = cbeta_y_j = 0;
	beta_lb_j = cbeta_lb_j = 0;
	beta_a_j = cbeta_a_j = 0;
	beta_yneg_j = cbeta_yneg_j = 0;
	beta_re_j = cbeta_re_j = 0;
	curve_epslb = ccurve_epslb = 0;

	Excursion_Flag = cExcursion_Flag = 0;
	Reversal_Flag = cReversal_Flag = 0;
	Energy_Excrsni_1 = cEnergy_Excrsni_1 = 0;
	Energy_Flag = cEnergy_Flag = 0;
	Yield_Flag = cYield_Flag = 0;
	Bucklingdegree = cBucklingdegree = 0;
	Buckling_flag = cBuckling_flag = 0;
	Stiffness_neg_Flag = cStiffness_neg_Flag = 0;

	eps_inflection = ceps_inflection = 0;
	sig_inflection = csig_inflection = 0;
	e_inflection = ce_inflection = 0;
	epsdiatance = cepsdiatance = 0;

	Minus = cMinus = 0;
	Minus_1 = cMinus_1 = 0;
	Minus_Flag = cMinus_Flag = 0;
	Minus_Flag_1 = cMinus_Flag_1 = 0;
	Minus_Flag_2 = cMinus_Flag_2 = 0;
	Stiffness_neg = cStiffness_neg = 0;
	sig_iso = csig_iso = 0;
	sig_lb_re = csig_lb_re = 0;

	kon_re = ckon_re = 0;
	eps_re = ceps_re = 0;
	sig_re = csig_re = 0;
	e_re = ce_re = 0;
	deps_re = cdeps_re = 0;
	epsreversal_re = cepsreversal_re = 0;
	sigreversal_re = csigreversal_re = 0;
	reloading_Flag = creloading_Flag = 0;
	sig_Trial = csig_Trial = 0;
	eps_max = ceps_max = 0;
	eps_max_Flag = ceps_max_Flag = 0;

	sig_Trial = csig_Trial = 0;
	eps_max = ceps_max = 0;
	eps_max_Flag = ceps_max_Flag = 0;
	lbstage_re = clbstage_re = 0;
	Minus_Flag_re = cMinus_Flag_re;
	sig_boundary = csig_boundary = 0;
	alpha_boundary = calpha_boundary = 0;
	Damage_Flag = cDamage_Flag = 0;
	lbstrength = clbstrength = Fy;

	Energy_Excrsn = cEnergy_Excrsn = 0;
	Energy_total = cEnergy_total = 0;

	e = ce = E0;
	e_1 = ce_1 = E0;
	e_2 = ce_2 = E0;

	strainConverged = 0.;
	strainPEqConverged = 0.;
	stressConverged = 0.;
	stiffnessConverged = E0;
	for (int i = 0; i < nBackstresses; ++i) {
		alphaKConverged[i] = 0.;
		alphaKTrial[i] = 0.;
	}

	for (int i = 0; i < nBackstresses; ++i) {
		alpha_boundary += (cK[i] / gammaK[i]);
	}
	sig_boundary = Fy + qInf + alpha_boundary;

	return 0;
}

UniaxialMaterial*
CFSTsteel::getCopy(void)
{
    CFSTsteel* copy = new CFSTsteel(*this);
    copy->setDbTag(0);
    return copy;
}

std::vector<double*> CFSTsteel::stateVariables()
{
    return {
        &Fy, &E0, &b2, &siglb,
        &a, &sigr, &b3, &Lambda_y,
        &Lambda_lb, &Lambda_a, &Lambda_yneg, &Lambda_re,
        &C_y, &C_lb, &C_a, &C_yneg,
        &C_re, &Eb_pos, &Eb_neg, &epsy_pos,
        &epsy_neg, &sigy_project_pos, &sigy_project_neg, &epslb,
        &slope_lb, &slope_r, &siglb_project, &sigr_project,
        &Eb_pos_j_1, &cEb_pos_j_1, &sigy_project_pos_j_1, &csigy_project_pos_j_1,
        &Eb_neg_j_1, &cEb_neg_j_1, &sigy_project_neg_j_1, &csigy_project_neg_j_1,
        &epslb_j_1, &cepslb_j_1, &siglb_j_1, &csiglb_j_1,
        &a_j_1, &ca_j_1, &sigr_j_1, &csigr_j_1,
        &slope_lb_j_1, &cslope_lb_j_1, &slope_r_j_1, &cslope_r_j_1,
        &siglb_project_j_1, &csiglb_project_j_1, &siglb_project_j, &csiglb_project_j,
        &sigr_project_j_1, &csigr_project_j_1, &Di, &cDi,
        &Di_1, &cDi_1, &Di_2, &cDi_2,
        &epsreversal, &cepsreversal, &sigreversal, &csigreversal,
        &epsreversal_1, &cepsreversal_1, &sigreversal_1, &csigreversal_1,
        &epsreversal_2, &cepsreversal_2, &sigreversal_2, &csigreversal_2,
        &sig, &csig, &sig_1, &csig_1,
        &sig_2, &csig_2, &eps, &ceps,
        &eps_1, &ceps_1, &eps_2, &ceps_2,
        &kon, &ckon, &kon_1, &ckon_1,
        &kon_2, &ckon_2, &lbstage, &clbstage,
        &beta_y_j_1, &cbeta_y_j_1, &beta_lb_j_1, &cbeta_lb_j_1,
        &beta_a_j_1, &cbeta_a_j_1, &beta_yneg_j_1, &cbeta_yneg_j_1,
        &beta_re_j_1, &cbeta_re_j_1, &curve_epslb, &ccurve_epslb,
        &curve_epsr, &ccurve_epsr, &beta_y_j, &beta_lb_j,
        &beta_a_j, &beta_yneg_j, &beta_re_j, &cbeta_y_j,
        &cbeta_lb_j, &cbeta_a_j, &cbeta_yneg_j, &cbeta_re_j,
        &Excursion_Flag, &cExcursion_Flag, &Reversal_Flag, &cReversal_Flag,
        &Energy_Excrsni_1, &cEnergy_Excrsni_1, &Energy_Flag, &cEnergy_Flag,
        &Energy_Excrsn, &cEnergy_Excrsn, &Energy_total, &cEnergy_total,
        &Ref_Energy_y, &Ref_Energy_lb, &Ref_Energy_a, &Ref_Energy_yneg,
        &Ref_Energy_re, &e, &ce, &e_1,
        &ce_1, &e_2, &ce_2, &Yield_Flag,
        &cYield_Flag, &Bucklingdegree, &cBucklingdegree, &Buckling_flag,
        &cBuckling_flag, &Stiffness_neg_Flag, &cStiffness_neg_Flag, &eps_inflection,
        &ceps_inflection, &sig_inflection, &csig_inflection, &e_inflection,
        &ce_inflection, &Minus, &cMinus, &Minus_1,
        &cMinus_1, &Minus_Flag, &cMinus_Flag, &Minus_Flag_1,
        &cMinus_Flag_1, &Minus_Flag_2, &cMinus_Flag_2, &epsdiatance,
        &cepsdiatance, &E_r_j_1, &cE_r_j_1, &Stiffness_neg,
        &cStiffness_neg, &sig_iso, &csig_iso, &exp_point,
        &cexp_point, &sig_lb_re, &csig_lb_re, &kon_re,
        &ckon_re, &eps_re, &ceps_re, &sig_re,
        &csig_re, &e_re, &ce_re, &deps_re,
        &cdeps_re, &epsreversal_re, &cepsreversal_re, &sigreversal_re,
        &csigreversal_re, &reloading_Flag, &creloading_Flag, &sig_Trial,
        &csig_Trial, &eps_max, &ceps_max, &eps_max_Flag,
        &ceps_max_Flag, &lbstage_re, &clbstage_re, &Minus_Flag_re,
        &cMinus_Flag_re, &sig_boundary, &csig_boundary, &alpha_boundary,
        &calpha_boundary, &Damage_Flag, &cDamage_Flag, &lbstrength,
        &clbstrength, &qInf, &b, &stiffnessInitial,
        &strainConverged, &strainTrial, &strainPEqConverged, &strainPEqTrial,
        &stressConverged, &stressTrial, &stiffnessConverged, &stiffnessTrial,
        &flowDirection,
    };
}

int CFSTsteel::sendSelf(int commitTag, Channel& theChannel)
{
    const std::vector<double*> state = stateVariables();

    Vector data(4 + static_cast<int>(state.size()) + 4 * 8);
    data.Zero();
    data(0) = 1;
    data(1) = getTag();
    data(2) = nBackstresses;
    data(3) = plasticLoading ? 1.0 : 0.0;
    for (unsigned i = 0; i < state.size(); ++i) data(4+i) = *state[i];
    const int start = 4 + static_cast<int>(state.size());
    for (int i = 0; i < nBackstresses; ++i) {
        data(start+i) = cK[i];
        data(start+8+i) = gammaK[i];
        data(start+16+i) = alphaKConverged[i];
        data(start+24+i) = alphaKTrial[i];
    }
    if (theChannel.sendVector(getDbTag(), commitTag, data) < 0) {
        opserr << "CFSTsteel::sendSelf - failed to send state" << endln;
        return -1;
    }
    return 0;
}

int CFSTsteel::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker&)
{
    const std::vector<double*> state = stateVariables();
    const int size = 4 + static_cast<int>(state.size()) + 4 * 8;
    Vector data(size);
    if (theChannel.recvVector(getDbTag(), commitTag, data) < 0 || data.Size() != size) {
        opserr << "CFSTsteel::recvSelf - failed to receive state" << endln;
        return -1;
    }
    for (int i = 0; i < data.Size(); ++i) if (!std::isfinite(data(i))) return -1;
    if (data(0) != 1.0 || data(2) < 0 || data(2) > 8 ||
        std::floor(data(2)) != data(2) || (data(3) != 0.0 && data(3) != 1.0)) {
        opserr << "CFSTsteel::recvSelf - invalid state header" << endln;
        return -1;
    }
    setTag(static_cast<int>(data(1)));
    nBackstresses = static_cast<int>(data(2));
    plasticLoading = data(3) != 0.0;
    for (unsigned i = 0; i < state.size(); ++i) *state[i] = data(4+i);
    cK.resize(nBackstresses);
    gammaK.resize(nBackstresses);
    alphaKConverged.resize(nBackstresses);
    alphaKTrial.resize(nBackstresses);
    const int start = 4 + static_cast<int>(state.size());
    for (int i = 0; i < nBackstresses; ++i) {
        cK[i] = data(start+i);
        gammaK[i] = data(start+8+i);
        alphaKConverged[i] = data(start+16+i);
        alphaKTrial[i] = data(start+24+i);
    }
    return 0;
}

void CFSTsteel::Print(OPS_Stream& s, int flag)
{
	if (flag != OPS_PRINT_PRINTMODEL_JSON) {
		s << "CFSTsteel tag: " << this->getTag() << endln;
		s << "  E0: " << E0 << " ";
		s << "  Fy: " << Fy << " ";
		s << "   b2: " << b2 << " ";
		s << "  siglb: " << siglb << " ";
		s << "  a: " << a << " ";
		s << " sigr: " << sigr << " ";
		s << " b3: " << b3 << " ";
		s << "  Lambda_y: " << Lambda_y << " ";
		s << "  Lambda_lb: " << Lambda_lb << " ";
		s << "  Lambda_a: " << Lambda_a << " ";
		s << "  Lambda_yneg: " << Lambda_yneg << " ";
		s << "  Lambda_re: " << Lambda_re << " ";
		s << "  C_y: " << C_y << " ";
		s << "  C_lb: " << C_lb << " ";
		s << "  C_a: " << C_a << " ";
		s << "  C_yneg: " << C_yneg << " ";
		s << "  C_re: " << C_re << " ";
		s << "  qInf: " << qInf << " ";
		s << "   b: " << b << " ";
		for (int i = 0; i < nBackstresses; ++i) {
			s << "  C" << (i + 1) << ": " << cK[i] << " ";
			s << "gam" << (i + 1) << ": " << gammaK[i] << " ";
		}
	}

	if (flag == 25000) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"CFSTsteel\", ";
		s << "\"E0\": " << E0 << ", ";
		s << "\"Fy\": " << Fy << ", ";
		s << "\"b2\": " << b2 << ", ";
		s << "\"siglb\": " << siglb << ", ";
		s << "\"a\": " << a << ", ";
		s << "\"sigr\": " << sigr << ", ";
		s << "\"b3\": " << b3 << ", ";
		s << "\"Lambda_y\": " << Lambda_y << ", ";
		s << "\"Lambda_lb\": " << Lambda_lb << ", ";
		s << "\"Lambda_a\": " << Lambda_a << ", ";
		s << "\"Lambda_yneg\": " << Lambda_yneg << ", ";
		s << "\"Lambda_re\": " << Lambda_re << ", ";
		s << "\"C_y\": " << C_y << ", ";
		s << "\"C_lb\": " << C_lb << ", ";
		s << "\"C_a\": " << C_a << ", ";
		s << "\"C_yneg\": " << C_yneg << ", ";
		s << "\"C_re\": " << C_re << ", ";
		s << "\"qInf\": " << qInf << ", ";
		s << "\"b\": " << b << ", ";
		for (int i = 0; i < nBackstresses; ++i) {
			s << "\"C\": " << cK[i] << ", ";
			s << "\"gam\": " << gammaK[i] << ", ";
		}
	}
}

template <typename T> int CFSTsteel::sgn(T val) {
	return (T(0) < val) - (val < T(0));
}
