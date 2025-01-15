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
	enum ErrorCodes {
		EC_IMPLEX_Error_Control = -10
	};
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
	static const char* msg = "uniaxialMaterial ASDSteel1D $tag $E $sy $su $eu <-implex>  <-implexControl $implexErrorTolerance $implexTimeReductionLimit>  ";

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
	bool implex = false;
	bool implex_control = false;
	double implex_error_tolerance = 0.05;  //to set
	double implex_time_redution_limit = 0.01; //to set

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
	auto lam_optional_double = [&numData](const char* variable, double& value) -> bool {
		if (OPS_GetNumRemainingInputArgs() > 0) {
			if (OPS_GetDouble(&numData, &value) < 0) {
				opserr << "nDMaterial ASDSteel1D Error: failed to get '" << variable << "'.\n";
				return false;
			}
		}
		else {
			opserr << "nDMaterial ASDSteel1D Error: '" << variable << "' requested but not provided.\n";
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

	// parse optional arguments
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* value = OPS_GetString();
		if (strcmp(value, "-implex") == 0) {
			implex = true;
		}
		else if (strcmp(value, "-implexControl") == 0) {
			implex_control = true;
			if (OPS_GetNumRemainingInputArgs() < 2) {
				opserr << "nDMaterial ASDSteel1D Error: '-implexControl' given without the next 2 arguments $implexErrorTolerance $implexTimeReductionLimit.\n";
				return nullptr;
			}
			if (!lam_optional_double("implexErrorTolerance", implex_error_tolerance))
				return nullptr;
			if (!lam_optional_double("implexTimeReductionLimit", implex_time_redution_limit))
				return nullptr;
		}
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
	params.implex = implex;
	params.implex_control = implex_control;
	params.implex_error_tolerance = implex_error_tolerance;
	params.implex_time_redution_limit = implex_time_redution_limit;

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
	// store the previously committed variables for next move from n to n - 1
	lambda_commit_old = lambda_commit;
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
	lambda_commit_old = 0.0;
	sg_commit = 0.0;

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
	ddata(counter++) = lambda_commit_old;
	ddata(counter++) = sg_commit;
}
void ASDSteel1DMaterial::StateVariablesSteel::recvSelf(int& counter, Vector& ddata) {
	alpha1 = ddata(counter++);
	alpha1_commit = ddata(counter++);
	alpha2 = ddata(counter++);
	alpha2_commit = ddata(counter++);
	lambda = ddata(counter++);
	lambda_commit = ddata(counter++);
	lambda_commit_old = ddata(counter++);
	sg_commit = ddata(counter++);
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

	// save dT
	if (!params.dtime_is_user_defined) {
		dtime_n = ops_Dt;
		if (!commit_done) {
			dtime_0 = dtime_n;
			dtime_n_commit = dtime_n;
		}
	}

	// compute real response
	steel.strain = v;
	//retval = computeBaseSteel(steel, params.implex);
	//if (retval < 0) return retval;
	if (params.implex) {
		if (params.implex_control) {
			// initial state
			double stress_implicit = 0.0;
			retval = computeBaseSteel(steel, false); // implicit solution
			if (retval < 0) return retval;
			stress_implicit = steel.stress;

			// explicit solution
			retval = computeBaseSteel(steel, params.implex); // Implex solution
			if (retval < 0) return retval;

			// Implex error
			params.implex_error = std::abs(stress_implicit - steel.stress)/params.sy;
			if (params.implex_error > params.implex_error_tolerance) {
				if (dtime_n >= params.implex_time_redution_limit * dtime_0) {
					return EC_IMPLEX_Error_Control;
				}
			}
		} else {
			retval = computeBaseSteel(steel, true); // Implex solution
			if (retval < 0) return retval;
		}
	}
	else {
		retval = computeBaseSteel(steel, false);
	}

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
	// implicit stage
	if (params.implex) {
		int retval;
		retval = computeBaseSteel(steel, false);
		if (retval < 0) return retval;
		// todo: homogenize
		strain = steel.strain;
		stress = steel.stress;
		C = steel.C;
		/*
		// IMPL-EX error
		params.implex_error = std::abs(steel.stress - stress);
		if (params.implex_control && params.implex_error > params.implex_error_tolerance) {
			return -1;   // avevamo detto quello implicito?
		}*/
	}

	// compute energy
	energy += 0.5 * (stress_commit + stress) * (strain - strain_commit);

	// store committed variables
	steel.commit(params);

	// state variables
	strain_commit = strain;
	stress_commit = stress;

	// implex
	dtime_n_commit = dtime_n;
	commit_done = true;

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

	// implex
	dtime_n = dtime_n_commit;

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

	// implex
	dtime_n = 0.0;
	dtime_n_commit = 0.0;
	dtime_0 = 0.0;
	

	commit_done = false;

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
	Vector ddata(InputParameters::NDATA + 1*StateVariablesSteel::NDATA + 10);
	counter = 0;
	ddata(counter++) = static_cast<double>(getTag());
	ddata(counter++) = params.E;
	ddata(counter++) = params.sy;
	ddata(counter++) = params.H1;
	ddata(counter++) = params.H2;
	ddata(counter++) = params.gamma1;
	ddata(counter++) = params.gamma2;
	ddata(counter++) = static_cast<double>(params.implex);
	ddata(counter++) = params.implex_error_tolerance;
	ddata(counter++) = params.implex_time_redution_limit;
	ddata(counter++) = static_cast<int>(params.implex_control);
	ddata(counter++) = static_cast<double>(params.dtime_is_user_defined);
	ddata(counter++) = params.implex_error;
	steel.sendSelf(counter, ddata);
	ddata(counter++) = dtime_n;
	ddata(counter++) = dtime_n_commit;
	ddata(counter++) = dtime_0;
	ddata(counter++) = static_cast<double>(commit_done);
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
	Vector ddata(InputParameters::NDATA + 1 * StateVariablesSteel::NDATA + 10);
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
	params.implex = static_cast<bool>(ddata(counter++));
	params.implex_error_tolerance = ddata(counter++);
	params.implex_time_redution_limit = ddata(counter++);
	params.implex_control = static_cast<bool>(ddata(counter++));
	params.dtime_is_user_defined = static_cast<bool>(ddata(counter++));
	params.implex_error = ddata(counter++);
	steel.recvSelf(counter, ddata);
	dtime_n = ddata(counter++);
	dtime_n_commit = ddata(counter++);
	dtime_0 = ddata(counter++);
	commit_done = static_cast<bool>(ddata(counter++));
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
	case 2000:
		dtime_n = info.theDouble;
		params.dtime_is_user_defined = true;
		return 0;
	case 2001:
		dtime_n_commit = info.theDouble;
		params.dtime_is_user_defined = true;
		return 0;
	case 2002:
		dtime_0 = info.theDouble;
		params.dtime_is_user_defined = true;
		return 0;

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
	static std::vector<std::string> lb_implex_error = { "Error" };

	// all outputs are 1D
	static Vector out1(1);

	// check specific responses
	if (argc > 0) {
		// 1000 - base steel output
		if (strcmp(argv[0], "equivalentPlasticStrain") == 0 || strcmp(argv[0], "EquivalentPlasticStrain") == 0) {
			out1(0) = steel.lambda;
			return make_resp(1001, out1, &lb_eqpl_strain);
		}
		// 3000 - implex error
		//if (strcmp(argv[0], "implexError") == 0 || strcmp(argv[0], "ImplexError") == 0) {
			//return make_resp(3000, getImplexError(), &lb_implex_error);
		//}
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
	//case 3000: return matInformation.setVector(getImplexError());
	default:
		break;
	}
	return UniaxialMaterial::getResponse(responseID, matInformation);
}

double ASDSteel1DMaterial::getEnergy(void)
{
	return energy;
}

int ASDSteel1DMaterial::computeBaseSteel(ASDSteel1DMaterial::StateVariablesSteel& sv, bool do_implex)
{
	// return value
	int retval = 0;

	// time factor for explicit extrapolation
	double time_factor = 1.0;
	if (params.implex && do_implex && (dtime_n_commit > 0.0))
		time_factor = dtime_n / dtime_n_commit;

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
	double sg = 0.0; // plastic flow direction
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
	if (params.implex && do_implex) {
		// extrapolate lambda
		//  xn + time_factor * (xn - xnn);
		double delta_lambda = time_factor * (sv.lambda_commit - sv.lambda_commit_old);
		sv.lambda = sv.lambda_commit + delta_lambda;
		// extrapolate plastic flow direction
		sg = sv.sg_commit;
		// update stress
		sigma -= sg * delta_lambda * params.E;
		// update backstress
		sv.alpha1 = sg * params.H1 / params.gamma1 - (sg * params.H1 / params.gamma1 - sv.alpha1_commit) * std::exp(-params.gamma1 * delta_lambda);
		sv.alpha2 = sg * params.H2 / params.gamma2 - (sg * params.H2 / params.gamma2 - sv.alpha2_commit) * std::exp(-params.gamma2 * delta_lambda);
	}
	else {
		// standard implicit evaluation of lambda
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
					sg = sign(lam_rel_stress());
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
	}
	// accept solution
	sv.stress = sigma;
	sv.C = tangent;

	// save real plastic flow direction, if mp.implex and !do_implex -> called from commit
	if (params.implex && !do_implex) {
		// save it in implex mode during implicit phase
		sv.sg_commit = sg;
	}

	// done
	return retval;
}
/*
const Vector& ASDSteel1DMaterial::getImplexError() const
{
	static Vector d(1); 
	d(0) = params.implex_error; 
	return d; 
}*/