/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.0 $
// $Date: 2042-06-14 11:29:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ASDSteel1DMaterial.cpp,v $

// Massimo Petracca - ASDEA Software, Italy
//
// A Simple and robust plastic-damage model for concrete and masonry
//

#include <ASDSteel1DMaterial.h>
#include <Channel.h>
#include <OPS_Globals.h>
#include <Information.h>
#include <Parameter.h>
#include <elementAPI.h>
#include <Element.h>
#include <MaterialResponse.h>
#include <cmath>
#include <algorithm>
#include <limits>
#include <string>
#include <sstream>
#include <iomanip>

// anonymous namespace for utilities
namespace {

	/**
	Converts a string into a vector of doubles using whitespace as delimiter
	*/
	bool string_to_double(const std::string& text, double& num) {
		num = 0.0;
		try {
			num = std::stod(text);
			return true;
		}
		catch (...) {
			return false;
		}
	}

	inline double sign(double x) { return x == 0.0 ? 0.0 : (x > 0.0 ? 1.0 : -1.0); }

}

void* OPS_ASDSteel1DMaterial()
{
	// some kudos
	static bool first_done = false;
	if (!first_done) {
		opserr << "Using ASDSteel1D - Developed by: Alessia Casalucci, Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
		first_done = true;
	}
	static const char* msg = "uniaxialMaterial ASDSteel1D $tag $E $sy $su $eu";

	// check arguments
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 5) {
		opserr <<
			"nDMaterial ASDSteel1D Error: Few arguments (< 5).\n" << msg << "\n";
		return nullptr;
	}

	// numData
	int numData = 1;

	// data
	int tag;
	double E;
	double sy;
	double su;
	double eu;

	// get tag
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "nDMaterial ASDSteel1D Error: invalid 'tag'.\n";
		return nullptr;
	}

	// get steel base arguments
	auto lam_get_dparam = [&numData](double* val, const char* valname) -> bool {
		if (OPS_GetDouble(&numData, val) != 0) {
			opserr << "nDMaterial ASDSteel1D Error: invalid '" << valname << "'.\n" << msg << "\n";
			return false;
		}
		if (*val <= 0.0) {
			opserr << "nDMaterial ASDSteel1D Error: invalid value for '" << valname << "' (" << *val << "). It should be strictly positive.\n" << msg << "\n";
			return false;
		}
		return true;
	};
	if (!lam_get_dparam(&E, "E")) return nullptr;
	if (!lam_get_dparam(&sy, "sy")) return nullptr;
	if (!lam_get_dparam(&su, "su")) return nullptr;
	if (!lam_get_dparam(&eu, "eu")) return nullptr;
	if (sy >= su) {
		opserr << "nDMaterial ASDSteel1D Error: invalid value for 'su' (" << su << "). It should be larger than sy.\n" << msg << "\n";
		return nullptr;
	}

	// obtain chaboche params from E, sy, su, eu
	// we want to use 2 hardening functions as per chaboche model.
	// so that the initial slope is close to E and the the stress apporaches su at eu
	double dy = su - sy;
	double H1 = E / 1000.0 / eu * dy / 40.0;
	double gamma1 = H1 / dy;
	double H2 = H1 * 50;
	double gamma2 = gamma1 * 50;
	double alpha = 0.9;
	ASDSteel1DMaterial::InputParameters params;
	params.E = E;
	params.sy = sy;
	params.H1 = H1 * alpha;
	params.gamma1 = gamma1;
	params.H2 = H2 * (1.0 - alpha);
	params.gamma2 = gamma2;

	// create the material
	UniaxialMaterial* instance = new ASDSteel1DMaterial(
		// tag
		tag, params);
		// base steel args
		//E, sy, H1*alpha, gamma1, H2*(1.0-alpha), gamma2
		// others..
	//);
	if (instance == nullptr) {
		opserr << "UniaxialMaterial ASDSteel1D Error: failed to allocate a new material.\n";
		return nullptr;
	}
	return instance;
}

void ASDSteel1DMaterial::StateVariablesSteel::commit(const  ASDSteel1DMaterial::InputParameters& params)
{
	// state variables
	alpha1_commit = alpha1;
	alpha2_commit = alpha2;
	lambda_commit = lambda;
	strain_commit = strain;
	stress_commit = stress;
}

void ASDSteel1DMaterial::StateVariablesSteel::revertToLastCommit(const  ASDSteel1DMaterial::InputParameters& params)
{
	// state variables
	alpha1 = alpha1_commit;
	alpha2 = alpha2_commit;
	lambda = lambda_commit;
	strain = strain_commit;
	stress = stress_commit;
}

void ASDSteel1DMaterial::StateVariablesSteel::revertToStart(const  ASDSteel1DMaterial::InputParameters& params)
{
	// state variables
	alpha1 = 0.0;
	alpha1_commit = 0.0;
	alpha2 = 0.0;
	alpha2_commit = 0.0;
	lambda = 0.0;
	lambda_commit = 0.0;

	// strain, stress and tangent
	strain = 0.0;
	strain_commit = 0.0;
	stress = 0.0;
	stress_commit = 0.0;
	C = params.E;
}
void ASDSteel1DMaterial::StateVariablesSteel::sendSelf(int& counter, Vector& ddata) {
	ddata(counter++) = alpha1;
	ddata(counter++) = alpha1_commit;
	ddata(counter++) = alpha2;
	ddata(counter++) = alpha2_commit;
	ddata(counter++) = lambda;
	ddata(counter++) = lambda_commit;
}
void ASDSteel1DMaterial::StateVariablesSteel::recvSelf(int& counter, Vector& ddata) {
	alpha1 = ddata(counter++);
	alpha1_commit = ddata(counter++);
	alpha2 = ddata(counter++);
	alpha2_commit = ddata(counter++);
	lambda = ddata(counter++);
	lambda_commit = ddata(counter++);
}

ASDSteel1DMaterial::ASDSteel1DMaterial(
	int _tag,
	const InputParameters& _params)
	: UniaxialMaterial(_tag, MAT_TAG_ASDSteel1DMaterial)
	, params(_params)
{
	// intialize C as C0
	C = getInitialTangent();
}

ASDSteel1DMaterial::ASDSteel1DMaterial()
	: UniaxialMaterial(0, MAT_TAG_ASDSteel1DMaterial)
{
}

ASDSteel1DMaterial::~ASDSteel1DMaterial()
{
}

int ASDSteel1DMaterial::setTrialStrain(double v, double r)
{
	// retval
	int retval = 0;

	// compute real response
	steel.strain = v;
	retval = computeBaseSteel(steel);
	if (retval < 0) return retval;

	// todo: homogenize
	strain = steel.strain;
	stress = steel.stress;
	C = steel.C;

	// done
	return retval;

}

double ASDSteel1DMaterial::getStress(void)
{
	return stress;
}

double ASDSteel1DMaterial::getTangent(void)
{
	return C;
}

double ASDSteel1DMaterial::getInitialTangent(void)
{
	return params.E;
}

double ASDSteel1DMaterial::getStrain(void)
{
	return strain;
}

int ASDSteel1DMaterial::commitState(void)
{
	// compute energy
	energy += 0.5 * (stress_commit + stress) * (strain - strain_commit);

	// store committed variables
	steel.commit(params);

	// state variables
	strain_commit = strain;
	stress_commit = stress;

	// done
	return 0;
}

int ASDSteel1DMaterial::revertToLastCommit(void)
{
	// restore converged values
	steel.revertToLastCommit(params);

	// state variables
	strain = strain_commit;
	stress = stress_commit;

	// done
	return 0;
}

int ASDSteel1DMaterial::revertToStart(void)
{
	// state variables
	steel.revertToStart(params);

	// strain, stress and tangent
	strain = 0.0;
	strain_commit = 0.0;
	stress = 0.0;
	stress_commit = 0.0;
	C = getInitialTangent();

	// output variables
	energy = 0.0;

	// done
	return 0;
}

UniaxialMaterial* ASDSteel1DMaterial::getCopy(void)
{
	// we can safely use the default copy-constructor according to the member variables we're using
	return new ASDSteel1DMaterial(*this);
}

void ASDSteel1DMaterial::Print(OPS_Stream& s, int flag)
{
	s << "ASDSteel1D Material, tag: " << this->getTag() << "\n";
}

int ASDSteel1DMaterial::sendSelf(int commitTag, Channel &theChannel)
{
	// aux
	int counter;

	// send DBL data
	Vector ddata(InputParameters::NDATA + 1*StateVariablesSteel::NDATA + 7);
	counter = 0;
	ddata(counter++) = static_cast<double>(getTag());
	ddata(counter++) = params.E;
	ddata(counter++) = params.sy;
	ddata(counter++) = params.H1;
	ddata(counter++) = params.H2;
	ddata(counter++) = params.gamma1;
	ddata(counter++) = params.gamma2;
	steel.sendSelf(counter, ddata);
	ddata(counter++) = strain;
	ddata(counter++) = strain_commit;
	ddata(counter++) = stress;
	ddata(counter++) = stress_commit;
	ddata(counter++) = C;
	ddata(counter++) = energy;
	if (theChannel.sendVector(getDbTag(), commitTag, ddata) < 0) {
		opserr << "ASDSteel1DMaterial::sendSelf() - failed to send DBL data\n";
		return -1;
	}

	// done
	return 0;
}

int ASDSteel1DMaterial::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	// aux
	int counter;

	// recv DBL data
	Vector ddata(InputParameters::NDATA + 1 * StateVariablesSteel::NDATA + 7);
	if (theChannel.recvVector(getDbTag(), commitTag, ddata) < 0) {
		opserr << "ASDSteel1DMaterial::recvSelf() - failed to receive DBL data\n";
		return -1;
	}
	counter = 0;
	setTag(ddata(counter++));
	params.E = ddata(counter++);
	params.sy = ddata(counter++);
	params.H1 = ddata(counter++);
	params.H2 = ddata(counter++);
	params.gamma1 = ddata(counter++);
	params.gamma2 = ddata(counter++);
	steel.recvSelf(counter, ddata);
	strain = ddata(counter++);
	strain_commit = ddata(counter++);
	stress = ddata(counter++);
	stress_commit = ddata(counter++);
	C = ddata(counter++);
	energy = ddata(counter++);

	// done
	return 0;
}

int ASDSteel1DMaterial::setParameter(const char** argv, int argc, Parameter& param)
{
	// default
	return -1;
}

int ASDSteel1DMaterial::updateParameter(int parameterID, Information& info)
{
	switch (parameterID) {
		// default
	default:
		return -1;
	}
}

Response* ASDSteel1DMaterial::setResponse(const char** argv, int argc, OPS_Stream& output)
{
	// utils
	auto make_resp = [&output, this](int rid, const Vector& v, const std::vector<std::string>* labels = nullptr) -> MaterialResponse* {
		output.tag("UniaxialMaterialOutput");
		output.attr("matType", getClassType());
		output.attr("matTag", getTag());
		if (labels) {
			for (const auto& item : (*labels))
				output.tag("ResponseType", item.c_str());
		}
		MaterialResponse* resp = new MaterialResponse(this, rid, v);
		output.endTag();
		return resp;
	};

	// labels
	static std::vector<std::string> lb_eqpl_strain = { "PLE" };

	// all outputs are 1D
	static Vector out1(1);

	// check specific responses
	if (argc > 0) {
		// 1000 - base steel output
		if (strcmp(argv[0], "equivalentPlasticStrain") == 0 || strcmp(argv[0], "EquivalentPlasticStrain") == 0) {
			out1(0) = steel.lambda;
			return make_resp(1001, out1, &lb_eqpl_strain);
		}
	}

	// otherwise return base-class response
	return UniaxialMaterial::setResponse(argv, argc, output);
}

int ASDSteel1DMaterial::getResponse(int responseID, Information& matInformation)
{
	// all outputs are 1D
	static Vector out1(1);

	switch (responseID) {
		// 1000 - base steel output
	case 1001:
		out1(0) = steel.lambda;
		return matInformation.setVector(out1);
	default:
		break;
	}
	return UniaxialMaterial::getResponse(responseID, matInformation);
}

double ASDSteel1DMaterial::getEnergy(void)
{
	return energy;
}

int ASDSteel1DMaterial::computeBaseSteel(ASDSteel1DMaterial::StateVariablesSteel& sv)
{
	// return value
	int retval = 0;

	// settings
	constexpr int MAX_ITER = 100;
	constexpr double F_REL_TOL = 1.0e-6;
	constexpr double L_ABS_TOL = 1.0e-8;
	// base steel response
	sv.alpha1 = sv.alpha1_commit;
	sv.alpha2 = sv.alpha2_commit;
	sv.lambda = sv.lambda_commit;
	// elastic predictor
	double dstrain = sv.strain - sv.strain_commit;
	double sigma = sv.stress_commit + params.E * dstrain;
	double tangent = params.E;
	// plastic utilities
	auto lam_rel_stress = [&sv, &sigma]() -> double {
		return sigma - sv.alpha1 - sv.alpha2;
	};
	auto lam_yield_function = [this, &lam_rel_stress]() -> double {
		return std::abs(lam_rel_stress()) - params.sy;
	};
	auto lam_yield_derivative = [this, &lam_rel_stress, &sv](double dlambda) -> double {
		// plastic flow direction
		double sg = sign(lam_rel_stress());
		// d stress / d lambda
		double dsigma = -params.E * sg;
		// d backstress / d lambda
		double dalpha1 = params.H1 * sg - params.gamma1 * sv.alpha1;
		double dalpha2 = params.H2 * sg - params.gamma2 * sv.alpha2;
		return sg * (dsigma - dalpha1 - dalpha2);
	};
	auto lam_yield_update = [this, &sigma, &lam_rel_stress, &sv](double dlambda, double delta_lambda) {
		// plastic flow direction
		double sg = sign(lam_rel_stress());
		// update stress
		sigma -= sg * dlambda * params.E;
		// update backstress
		sv.alpha1 = sg * params.H1 / params.gamma1 - (sg * params.H1 / params.gamma1 - sv.alpha1_commit) * std::exp(-params.gamma1 * delta_lambda);
		sv.alpha2 = sg * params.H2 / params.gamma2 - (sg * params.H2 / params.gamma2 - sv.alpha2_commit) * std::exp(-params.gamma2 * delta_lambda);
	};
	// plastic corrector
	double F = lam_yield_function();
	if (F > 0.0) {
		double delta_lambda = 0.0;
		double dlambda = 0.0;
		bool converged = false;
		for (int niter = 0; niter < MAX_ITER; ++niter) {
			// form tangent
			double dF = lam_yield_derivative(dlambda);
			if (dF == 0.0) break;
			// solve for dlambda
			dlambda = -F / dF;
			delta_lambda += dlambda;
			// update plastic multiplier increment and sigma
			lam_yield_update(dlambda, delta_lambda);
			// update residual
			F = lam_yield_function();
			// check convergence
			if (std::abs(F) < F_REL_TOL * params.sy && std::abs(dlambda) < L_ABS_TOL) {
				converged = true;
				// update plastic multiplier
				sv.lambda += delta_lambda;
				// compute tangent
				double sg = sign(lam_rel_stress());
				double PE =
					params.gamma1 * (params.H1 / params.gamma1 - sg * sv.alpha1) +
					params.gamma2 * (params.H2 / params.gamma2 - sg * sv.alpha2);
				tangent = (params.E * PE) / (params.E + PE);
				break;
			}
		}
		if (!converged)
			retval = -1;
	}
	// accept solution
	sv.stress = sigma;
	sv.C = tangent;

	// done
	return retval;
}

