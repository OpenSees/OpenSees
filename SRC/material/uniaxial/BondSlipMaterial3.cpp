/* Written by: Mohammad Salehi (mohammad.salehi@rice.edu)
** Created: 2020
** Description: The source code for BondSlip material model used in model BarSlip2
**
**
** Reference:
**
** Mohammad Salehi, Petros Sideris, and Reginald DesRoches (2022)
** “Numerical modeling of repaired reinforced concrete bridge columns”
** Engineering Structures, 253: 113801
*/

#include <stdlib.h>
#include <BondSlipMaterial3.h>
#include <Channel.h>
#include <math.h>
#include <float.h>
#include <iostream>
#include <elementAPI.h>

void*
OPS_BondSlipMaterial3(void)
{
	// pointer to a uniaxial material that will be returned
	UniaxialMaterial* theMaterial = 0;

	// get Input Values
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs != 6) {
		opserr << "WARNING! Invalid number of args in uniaxialMaterial BondSlip" << endln;
		opserr << "want: uniaxialMaterial BondSlip tag? s_o? s_1? tau_bo? tau_fo? k_o?" << endln;
		return 0;
	}

	int tag;
	double dData[5];

	int numData = 1;
	if (OPS_GetIntInput(&numData, &tag) != 0) {
		opserr << "WARNING! Invalid tags for uniaxialMaterial BondSlip" << endln;
		return 0;
	}

	numData = 5;
	if (OPS_GetDoubleInput(&numData, dData) != 0) {
		opserr << "WARNING! Invalid double parameters for uniaxialMaterial BondSlip " << tag << endln;
		return 0;
	}

	// allocate material model
	theMaterial = new BondSlipMaterial3(tag, dData[0], dData[1], dData[2], dData[3], dData[4]);

	if (theMaterial == 0) {
		opserr << "WARNING! Could not create uniaxialMaterial of type BondSlip\n";
		return 0;
	}

	return theMaterial;
}


// Constructors
BondSlipMaterial3::BondSlipMaterial3(int tag, double so, double s1, double to, double t1, double ko)
	:UniaxialMaterial(tag, MAT_TAG_BondSlip),
	s_o(so), s_1(s1), tau_o(to), tau_1(t1),
	k_o(ko), E_o(0.0), E_of(0.0),
	k_init(0.0), PI(3.14159265358979323846),
	s_commit(0.0), s_trial(0.0), tau_commit(0.0), tau_trial(0.0),
	k_tan(0.0), 
	tau_or_pos(to), tau_1r_pos(t1), tau_f_pos(0.1 * t1), tau_fr_pos(0.1 * t1),
	tau_or_neg(to), tau_1r_neg(t1), tau_f_neg(0.1 * t1), tau_fr_neg(0.1 * t1),
	s_max(0.0), s_min(0.0), s_max_commit(0.0), s_min_commit(0.0), s_max_env(0.0), s_min_env(0.0),
	E_d(0.0), E_f_pos(0.0), E_f_neg(0.0)
{
	if (s_o > s_1) {
		opserr << "ERROR! BondSlipMaterial3::BondSlipMaterial3 - material: " << this->getTag()
			<< "\ns_o > s_1, so switched them!" << endln;
		s_o = s1;
		s_1 = so;
	}

	if (tau_o < tau_1) {
		opserr << "ERROR! BondSlipMaterial3::BondSlipMaterial3 - material: " << this->getTag()
			<< "\ntau_o < tau_1, so switched them!" << endln;
		tau_o = t1;
		tau_1 = to;
	}

	// Calculate constant values
	k_init = 6.0 * tau_o / (1.2 * s_o);			// initial stiffness

	E_o = tau_o * 0.8 * s_o + (s_1 - 1.2 * s_o) * (tau_o + tau_1) / 2.0;		// approximated energy dissipation under monotonic response
	E_of = tau_1 * s_1;

	if (k_o < k_init) {
		opserr << "ERROR! BondSlipMaterial3::BondSlipMaterial3 - material: " << this->getTag()
			<< "\nk_o < k_init, so used k_o = k_init" << endln;

		k_o = k_init;
	}
}

BondSlipMaterial3::BondSlipMaterial3()
	:UniaxialMaterial(0, MAT_TAG_BondSlip),
	s_o(0.0), s_1(0.0), tau_o(0.0), tau_1(0.0),
	k_o(0.0), E_o(0.0), E_of(0.0),
	k_init(0.0), PI(3.14159265358979323846),
	s_commit(0.0), s_trial(0.0), tau_commit(0.0), tau_trial(0.0),
	k_tan(0.0),
	tau_or_pos(0.0), tau_1r_pos(0.0), tau_f_pos(0.0), tau_fr_pos(0.0),
	tau_or_neg(0.0), tau_1r_neg(0.0), tau_f_neg(0.0), tau_fr_neg(0.0),
	s_max(0.0), s_min(0.0), s_max_commit(0.0), s_min_commit(0.0), s_max_env(0.0), s_min_env(0.0),
	E_d(0.0), E_f_pos(0.0), E_f_neg(0.0)
{
	// does nothing
}

// Destructor
BondSlipMaterial3::~BondSlipMaterial3()
{
	// does nothing
}

// Method to Copy Material
UniaxialMaterial*
BondSlipMaterial3::getCopy(void)
{
	BondSlipMaterial3* theCopy = new BondSlipMaterial3(this->getTag(), s_o, s_1, tau_o, tau_1, k_o);

	theCopy->s_commit = s_commit;
	theCopy->s_trial = s_trial;
	theCopy->tau_commit = tau_commit;
	theCopy->tau_trial = tau_trial;
	theCopy->k_tan = k_tan;
	theCopy->s_max = s_max;
	theCopy->s_min = s_min;
	theCopy->s_max_env = s_max_env;
	theCopy->s_min_env = s_min_env;
	theCopy->s_max_commit = s_max_commit;
	theCopy->s_min_commit = s_min_commit;
	theCopy->E_d = E_d;
	theCopy->E_f_pos = E_f_pos;
	theCopy->E_f_neg = E_f_neg;

	return theCopy;
}

// Methods to Set the State of Material
int
BondSlipMaterial3::setTrialStrain(double strain, double strainRate)
{

	// Set Trial Strain
	s_trial = strain;

	double ds = s_trial - s_commit;

	// Compute Trial Bond Stress and Corresponding Tangent
	tau_trial = tau_commit + k_o * ds;
	k_tan = k_o;

	double tau_env, k_env, s_o_env;		// envelope bond stress and tangent stiffness
	double xi, xi_1, dxi2ds;			// for envelope calculations
	double k_r = 1.0E-6 * k_o;
	double expFact = 1.0 / (1.0 - exp(-6.0));

	if (ds >= 0.0) {	// positive slip increase
		
		// calculate slip corresponding to peak bond stress
		s_o_env = 1.2 * (s_o + s_max_env);

		if (s_o_env > 0.9 * s_1) {
			s_o_env = 0.9 * s_1;
			s_max_env = s_o_env / 1.2 - s_o;
		}

		xi = s_trial / s_o_env;
		xi_1 = s_1 / s_o_env;
		dxi2ds = 1.0 / s_o_env;	// derivative of xi with respect to s

		// determine upper envelope bond stress for trial slip
		if (s_trial < 0.0) {
			k_env = k_r;	// to avoid zero stiffness
			tau_env = tau_fr_pos + k_env * s_trial;
		}
		else if (s_trial <= s_o_env) {
			if (s_max_env > 0.0) {		// hyperbolic form
				double th = tanh(3.0 * (s_trial - s_max_env) / fmin(s_max_env, s_o));
				double tho = tanh(3.0 * (s_o_env - s_max_env) / fmin(s_max_env, s_o));
				double w_r = (1.0 + th) / (1.0 + tho);
				double dw2ds = (3.0 / fmin(s_max_env, s_o)) * (1.0 - th * th) / (1.0 + tho);

				//tau_env = tau_fr_pos + w_r * (tau_or_pos - tau_fr_pos) * sinf(PI * xi * (1.0 - xi / 2.0));
				tau_env = tau_fr_pos + w_r * (tau_or_pos - tau_fr_pos) * (1.0 - exp(-6.0 * xi)) * expFact;
				//k_env = dxi2ds * (w_r * (tau_or_pos - tau_fr_pos) * cosf(PI * xi * (1.0 - xi / 2.0)) * PI * (1.0 - xi)) + dw2ds * (tau_or_pos - tau_fr_pos) * sinf(PI * xi * (1.0 - xi / 2.0));
				k_env = dxi2ds * (w_r * (tau_or_pos - tau_fr_pos) * 6.0 * exp(-6.0 * xi) * expFact) + dw2ds * (tau_or_pos - tau_fr_pos) * (1.0 - exp(-6 * xi)) * expFact;
			}
			else if (s_min_env < 0.0) {
				double w_r = (s_trial < 0.2 * s_o) ? (s_trial / (0.2 * s_o)) : 1.0;
				double dw2ds = (s_trial < 0.2 * s_o) ? (1.0 / (0.2 * s_o)) : 0.0;

				//tau_env = tau_fr_pos + w_r * (tau_or_pos - tau_fr_pos) * sinf(PI * xi * (1.0 - xi / 2.0));
				tau_env = tau_fr_pos + w_r * (tau_or_pos - tau_fr_pos) * (1.0 - exp(-6.0 * xi)) * expFact;
				//k_env = dxi2ds * (w_r * (tau_or_pos - tau_fr_pos) * cosf(PI * xi * (1.0 - xi / 2.0)) * PI * (1.0 - xi)) + dw2ds * (tau_or_pos - tau_fr_pos) * sinf(PI * xi * (1.0 - xi / 2.0));
				k_env = dxi2ds * (w_r * (tau_or_pos - tau_fr_pos) * 6.0 * exp(-6.0 * xi) * expFact) + dw2ds * (tau_or_pos - tau_fr_pos) * (1.0 - exp(-6 * xi)) * expFact;
			}
			else {
				//tau_env = tau_or_pos * sinf(PI * xi * (1.0 - xi / 2.0));
				tau_env = tau_or_pos * (1.0 - exp(-6 * xi)) * expFact;
				//k_env = dxi2ds * (tau_or_pos * cosf(PI * xi * (1.0 - xi / 2.0)) * PI * (1.0 - xi));
				k_env = dxi2ds * (tau_or_pos * 6.0 * exp(-6.0 * xi) * expFact);
			}
		}
		else if (s_trial <= s_1) {
			tau_env = (xi - 1.0) * (tau_1r_pos - tau_or_pos) / (xi_1 - 1.0) + tau_or_pos;
			k_env = dxi2ds * (tau_1r_pos - tau_or_pos) / (xi_1 - 1.0);
		}
		else {
			k_env = k_r;	// to avoid zero stiffness
			tau_env = tau_1r_pos + k_env * (s_trial - s_1);
		}

		// check if envelope stress has been exceeded
		if (tau_trial >= tau_env) {
			tau_trial = tau_env;
			k_tan = k_env;
		}
	}
	else {	// negative slip increase

		// calculate slip corresponding to peak bond stress
		s_o_env = 1.2 * (s_o - s_min_env);

		if (s_o_env > 0.9 * s_1) {
			s_o_env = 0.9 * s_1;
			s_min_env = s_o - s_o_env / 1.2;
		}

		xi = -s_trial / s_o_env;
		xi_1 = s_1 / s_o_env;
		dxi2ds = 1.0 / s_o_env;	// derivative of xi with respect to s

		// determine lower envelope bond stress for trial slip
		if (s_trial > 0.0) {
			k_env = k_r;	// to avoid zero stiffness
			tau_env = -tau_fr_neg + k_env * s_trial;
		}
		else if (s_trial >= -s_o_env) {
			if (s_min_env < 0.0) {		// hyperbolic form
				double th = tanh(3.0 * (-s_trial + s_min_env) / fmin(-s_min_env, s_o));
				double tho = tanh(3.0 * (s_o_env + s_min_env) / fmin(-s_min_env, s_o));
				double w_r = (1.0 + th) / (1.0 + tho);
				double dw2ds = (3.0 / fmin(-s_min_env, s_o)) * (1.0 - th * th) / (1.0 + tho);

				//tau_env = -1.0 * (tau_fr_neg + w_r * (tau_or_neg - tau_fr_neg) * sinf(PI * xi * (1.0 - xi / 2.0)));
				tau_env = -1.0 * (tau_fr_neg + w_r * (tau_or_neg - tau_fr_neg) * (1.0 - exp(-6.0 * xi)) * expFact);
				//k_env = dxi2ds * (w_r * (tau_or_neg - tau_fr_neg) * cosf(PI * xi * (1.0 - xi / 2.0)) * PI * (1.0 - xi)) + dw2ds * (tau_or_neg - tau_fr_neg) * sinf(PI * xi * (1.0 - xi / 2.0));
				k_env = dxi2ds * (w_r * (tau_or_neg - tau_fr_neg) * 6.0 * exp(-6.0 * xi) * expFact) + dw2ds * (tau_or_neg - tau_fr_neg) * (1.0 - exp(-6.0 * xi)) * expFact;
			}
			else if (s_max_env > 0.0) {
				double w_r = (s_trial > -0.2 * s_o) ? (-s_trial / (0.2 * s_o)) : 1.0;
				double dw2ds = (s_trial > -0.2 * s_o) ? (1.0 / (0.2 * s_o)) : 0.0;

				//tau_env = -1.0 * (tau_fr_neg + w_r * (tau_or_neg - tau_fr_neg) * sinf(PI * xi * (1.0 - xi / 2.0)));
				tau_env = -1.0 * (tau_fr_neg + w_r * (tau_or_neg - tau_fr_neg) * (1.0 - exp(-6.0 * xi)) * expFact);
				//k_env = dxi2ds * (w_r * (tau_or_neg - tau_fr_neg) * cosf(PI * xi * (1.0 - xi / 2.0)) * PI * (1.0 - xi)) + dw2ds * (tau_or_neg - tau_fr_neg) * sinf(PI * xi * (1.0 - xi / 2.0));
				k_env = dxi2ds * (w_r * (tau_or_neg - tau_fr_neg) * 6.0 * exp(-6.0 * xi) * expFact) + dw2ds * (tau_or_neg - tau_fr_neg) * (1.0 - exp(-6.0 * xi)) * expFact;
			}
			else {
				//tau_env = -tau_or_neg * sinf(PI * xi * (1.0 - xi / 2.0));
				tau_env = -tau_or_neg * (1.0 - exp(-6.0 * xi)) * expFact;
				//k_env = dxi2ds * (tau_or_neg * cosf(PI * xi * (1.0 - xi / 2.0)) * PI * (1.0 - xi));
				k_env = dxi2ds * (tau_or_neg * 6.0 * exp(-6.0 * xi) * expFact);
			}
		}
		else if (s_trial >= -s_1) {
			tau_env = -1.0 * ((xi - 1.0) * (tau_1r_neg - tau_or_neg) / (xi_1 - 1.0) + tau_or_neg);
			k_env = dxi2ds * (tau_1r_neg - tau_or_neg) / (xi_1 - 1.0);
		}
		else {
			k_env = k_r;	// to avoid zero stiffness
			tau_env = -1.0 * (tau_1r_neg + k_env * (-s_trial - s_1));
		}

		// check if envelope stress has been exceeded
		if (tau_trial <= tau_env) {
			tau_trial = tau_env;
			k_tan = k_env;
		}
	}

	return 0;
}

int
BondSlipMaterial3::commitState(void)
{
	// update maximum, minimum, and accumulated slip
	s_max = fmax(s_trial, s_max);
	s_min = fmin(s_trial, s_min);

	// update dissipated energy
	double ds = s_trial - s_commit;

	if (s_trial * tau_trial >= 0.0)		// not only frictional
		if (0.5 * (tau_trial + tau_commit) > tau_fr_pos) {
			E_d += 0.5 * (tau_trial + tau_commit - tau_fr_pos) * ds;
			E_f_neg += tau_fr_pos * ds;
			E_f_pos += tau_fr_pos * ds;
		}
		else if (0.5 * (tau_trial + tau_commit) < -tau_fr_neg) {
			E_d += 0.5 * (tau_trial + tau_commit + tau_fr_neg) * ds;
			E_f_neg -= tau_fr_neg * ds;
			E_f_pos -= tau_fr_neg * ds;
		}
		else {
			E_d += 0.5 * 0.5 * (tau_trial + tau_commit) * ds;
			E_f_neg += 0.5 * (tau_trial + tau_commit) * ds;
			E_f_pos += 0.5 * (tau_trial + tau_commit) * ds;
		}
	else {
		E_d += 0.5 * 0.5 * (tau_trial + tau_commit) * ds;
		E_f_neg += 0.5 * (tau_trial + tau_commit) * ds;
		E_f_pos += 0.5 * (tau_trial + tau_commit) * ds;
	}

	E_f_pos = fmax(0.0, E_f_pos);
	E_f_neg = fmax(0.0, E_f_neg);
	E_d = fmax(0.0, E_d);

	// determine maximum stress damage factor
	double d_o = fmin(1.0 - exp(-1.2 * pow(E_d / E_o, 1.1)), 0.9);		// limit maximum damage factor to 0.9

	// determine residual stress damage factor
	double d_1 = d_o / (2.0 - d_o);

	// determine reduced maximum and resdiual stress limits and maximum/minimum slip values used to obtain envelopes
	if (s_trial > s_commit && tau_trial >= 0.0) {	// positive slip direction
		tau_or_neg = (1.0 - d_o) * tau_o;
		tau_1r_neg = (1.0 - d_1) * tau_1;

		s_min_env = s_min;
	}
	else if (s_trial < s_commit && tau_trial <= 0.0) {	// negative slip direction
		tau_or_pos = (1.0 - d_o) * tau_o;
		tau_1r_pos = (1.0 - d_1) * tau_1;

		s_max_env = s_max;
	}

	// determine reduced friction stress
	if (s_max > s_max_commit) {
		double tf2t1 = tau_fr_neg / tau_1r_neg + 1.8 * (s_max - s_max_commit) / s_1;

		if (tf2t1 > 1.0)
			tf2t1 = 1.0;

		tau_f_neg = tf2t1 * tau_1r_neg;

		E_f_neg = 0.0;		// reset to zero after first reversal in current maximum slip range
	}

	if (s_min < s_min_commit) {
		double tf2t1 = tau_fr_pos / tau_1r_pos + 1.8 * (-s_min + s_min_commit) / s_1;

		if (tf2t1 > 1)
			tf2t1 = 1.0;

		tau_f_pos = tf2t1 * tau_1r_pos;

		E_f_pos = 0.0;		// reset to zero after first reversal in current maximum slip range
	}

	if (s_trial > s_commit && tau_trial >= 0.0) {	// positive slip direction
		double d_f_neg = fmin(1.0 - exp(-1.2 * pow(E_f_neg / E_of, 0.67)), 0.9);	// limit maximum damage factor to 0.9
		tau_fr_neg = (1.0 - d_f_neg) * tau_f_neg;
	}
	else if (s_trial < s_commit && tau_trial <= 0.0) {	// negative slip direction
		double d_f_pos = fmin(1.0 - exp(-1.2 * pow(E_f_pos / E_of, 0.67)), 0.9);	// limit maximum damage factor to 0.9
		tau_fr_pos = (1.0 - d_f_pos) * tau_f_pos;
	}

	// commit state variables
	s_commit = s_trial;
	s_max_commit = s_max;
	s_min_commit = s_min;
	tau_commit = tau_trial;

	return 0;
}


int
BondSlipMaterial3::revertToLastCommit(void)
{
	s_trial = s_commit;
	tau_trial = tau_commit;

	return 0;
}


int
BondSlipMaterial3::revertToStart(void)
{
	s_trial = s_commit = 0.0;
	tau_trial = tau_commit = 0.0;
	tau_1r_pos = tau_1r_neg = tau_o;
	tau_1r_pos = tau_1r_neg = tau_1;
	tau_f_pos = tau_f_neg = tau_fr_pos = tau_fr_neg = 0.1 * tau_1;
	k_tan = k_init;
	s_max = s_min = s_max_commit = s_min_commit = s_max_env = s_min_env = 0.0;
	E_d = E_f_neg = E_f_pos = 0.0;

	return 0;
}

// Methods to Get Response Variables
double
BondSlipMaterial3::getStrain(void)
{
	return s_trial;
}

double
BondSlipMaterial3::getStress(void)
{
	return tau_trial;
}

double
BondSlipMaterial3::getTangent(void)
{
	return k_tan; // (k_tan > 0.0 ? k_tan : 1.0e-6 * k_o);
}

double
BondSlipMaterial3::getInitialTangent(void)
{
	return k_init;
}

// Methods to Do Parallel Processing
int
BondSlipMaterial3::sendSelf(int cTag, Channel& theChannel)
{
	opserr << "ERROR! BondSlipMaterial3::sendSelf - material: " << this->getTag() << " - incapable of parallel processing";
	return -1;
}

int
BondSlipMaterial3::recvSelf(int cTag, Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	opserr << "ERROR! BondSlipMaterial3::recvSelf - material: " << this->getTag() << " - incapable of parallel processing";
	return -1;
}

// Methods to Obtain Information
void
BondSlipMaterial3::Print(OPS_Stream& s, int flag)
{
	s << "BondSlip tag: " << this->getTag() << endln;
	s << "  so:    " << s_o << endln;
	s << "  s1:    " << s_1 << endln;
	s << "  tau_o: " << tau_o << endln;
	s << "  tau_1: " << tau_1 << endln;
	s << "  k_o:   " << k_o << endln;
}