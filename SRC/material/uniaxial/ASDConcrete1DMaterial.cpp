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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ASDConcrete1DMaterial.cpp,v $

// Massimo Petracca - ASDEA Software, Italy
//
// A Simple and robust plastic-damage model for concrete and masonry
//

#include <ASDConcrete1DMaterial.h>
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
		EC_Generic = -1,
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
	bool string_to_list_of_doubles(const std::string& text, char sep, std::vector<double>& out) {
		if (out.size() > 0) out.clear();
		std::size_t start = 0, end = 0;
		double value;
		while (true) {
			end = text.find(sep, start);
			if (end == std::string::npos) {
				if (start < text.size()) {
					if (!string_to_double(text.substr(start), value))
						return false;
					out.push_back(value);
				}
				break;
			}
			std::string subs = text.substr(start, end - start);
			if (subs.size() > 0) {
				if (!string_to_double(subs, value))
					return false;
				out.push_back(value);
			}
			start = end + 1;
		}
		return true;
	}

	// Heavyside function
	inline double Heavyside(double X) { return X > 0.0 ? 1.0 : (X < 0.0 ? 0.0 : 0.5); }

	// Macauley function
	inline double Macauley(double X) { return X > 0.0 ? X : 0.0; }

	/**
	global parameters storage
	*/
	class GlobalParameters {
	private:
		double max_error = 0.0;
		double avg_error = 0.0;
		int avg_counter = 0;
	private:
		GlobalParameters() = default;
		GlobalParameters(const GlobalParameters&) = delete;
		GlobalParameters& operator = (const GlobalParameters&) = delete;
	public:
		static GlobalParameters& instance() {
			static GlobalParameters _instance;
			return _instance;
		}
		inline double getMaxError() const {
			return max_error;
		}
		inline void setMaxError(double x) {
			max_error = x;
		}
		inline double getAverageError() {
			if (avg_counter > 0) {
				avg_error /= static_cast<double>(avg_counter);
				avg_counter = 0;
			}
			return avg_error;
		}
		inline void accumulateAverageError(double x) {
			avg_error += x;
			++avg_counter;
		}
		inline void setAverageError(double x) {
			avg_error = x;
			avg_counter = 0;
		}
	};

	double bezier3(double xi,
		double x0, double x1, double x2,
		double y0, double y1, double y2)
	{
		double A = x0 - 2.0 * x1 + x2;
		double B = 2.0 * (x1 - x0);
		double C = x0 - xi;
		if (fabs(A) < 1.0e-12) {
			x1 = x1 + 1.0E-6 * (x2 - x0);
			A = x0 - 2.0 * x1 + x2;
			B = 2.0 * (x1 - x0);
			C = x0 - xi;
		}
		if (A == 0.0)
			return 0.0;

		double D = B * B - 4.0 * A * C;
		double t = (sqrt(D) - B) / (2.0 * A);

		return (y0 - 2.0 * y1 + y2) * t * t + 2.0 * (y1 - y0) * t + y0;
	}

}

void* OPS_ASDConcrete1DMaterial()
{
	// some kudos
	static bool first_done = false;
	if (!first_done) {
		opserr << "Using ASDConcrete1D - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
		first_done = true;
	}

	// check arguments
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 2) {
		opserr <<
			"nDMaterial ASDConcrete1D Error: Few arguments (< 2).\n"
			"nDMaterial ASDConcrete1D $tag $E "
			"<-fc $fc> <-ft $ft> "
			"<-Te $Te -Ts $Ts <-Td $Td>> <-Ce $Ce -Cs $Cs <-Cd $Cd>> "
			"<-implex> <-implexControl $implexErrorTolerance $implexTimeReductionLimit> <-implexAlpha $alpha>"
			"<-eta $eta> <-tangent> <-autoRegularization $lch_ref>\n";
		return nullptr;
	}

	// numData
	int numData = 1;

	// data
	int tag;
	double E;
	bool implex = false;
	bool implex_control = false;
	double implex_error_tolerance = 0.05;
	double implex_time_redution_limit = 0.01;
	double implex_alpha = 1.0;
	double eta = 0.0;
	bool tangent = false;
	bool auto_regularization = false;
	double lch_ref = 1.0;
	std::vector<double> Te, Ts, Td, Ce, Cs, Cd;

	// get tag
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "nDMaterial ASDConcrete1D Error: invalid 'tag'.\n";
		return nullptr;
	}

	// get Elasticity arguments
	if (OPS_GetDouble(&numData, &E) != 0) {
		opserr << "nDMaterial ASDConcrete1D Error: invalid 'E'.\n";
		return nullptr;
	}
	if (E <= 0.0) {
		opserr << "nDMaterial ASDConcrete1D Error: invalid value for 'E' (" << E << "). It should be strictly positive.\n";
		return nullptr;
	}

	// utilities (code re-use)
	auto lam_optional_int = [&numData](const char* variable, int& value) -> bool {
		if (OPS_GetNumRemainingInputArgs() > 0) {
			if (OPS_GetInt(&numData, &value) < 0) {
				opserr << "nDMaterial ASDConcrete1D Error: failed to get '" << variable << "'.\n";
				return false;
			}
		}
		else {
			opserr << "nDMaterial ASDConcrete1D Error: '" << variable << "' requested but not provided.\n";
			return false;
		}
		return true;
	};
	auto lam_optional_double = [&numData](const char* variable, double& value) -> bool {
		if (OPS_GetNumRemainingInputArgs() > 0) {
			if (OPS_GetDouble(&numData, &value) < 0) {
				opserr << "nDMaterial ASDConcrete1D Error: failed to get '" << variable << "'.\n";
				return false;
			}
		}
		else {
			opserr << "nDMaterial ASDConcrete1D Error: '" << variable << "' requested but not provided.\n";
			return false;
		}
		return true;
	};
	auto lam_optional_list = [&numData](const char* variable, std::vector<double>& value) -> bool {
		// first try expanded list like {*}$the_list,
		// also used in python like *the_list
		value.clear();
		while (OPS_GetNumRemainingInputArgs() > 0) {
			double item;
			auto old_num_rem = OPS_GetNumRemainingInputArgs();
			if (OPS_GetDoubleInput(&numData, &item) < 0) {
				auto new_num_rem = OPS_GetNumRemainingInputArgs();
				if (new_num_rem < old_num_rem)
					OPS_ResetCurrentInputArg(-1);
				break;
			}
			value.push_back(item);
		}
		// try Tcl list (it's a string after all...)
		if (value.size() == 0 && OPS_GetNumRemainingInputArgs() > 0) {
			std::string list_string = OPS_GetString();
			if (!string_to_list_of_doubles(list_string, ' ', value)) {
				opserr << "nDMaterial ASDConcrete1D Error: cannot parse the '" << variable << "' list.\n";
				return false;
			}
		}
		return true;
	};

	double fc;
	double ft;
	bool have_fc = false;
	bool have_ft = false;
	bool have_lch_ref = false;

	// optional parameters
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* value = OPS_GetString();
		if (strcmp(value, "-fc") == 0) {
			if (!lam_optional_double("fc", fc))
				return nullptr;
			have_fc = true;
		}
		else if (strcmp(value, "-ft") == 0) {
			if (!lam_optional_double("ft", ft))
				return nullptr;
			have_ft = true;
		}
		else if (strcmp(value, "-implex") == 0) {
			implex = true;
		}
		else if (strcmp(value, "-implexControl") == 0) {
			implex_control = true;
			if (OPS_GetNumRemainingInputArgs() < 2) {
				opserr << "nDMaterial ASDConcrete1D Error: '-implexControl' given without the next 2 arguments $implexErrorTolerance $implexTimeReductionLimit.\n";
				return nullptr;
			}
			if (!lam_optional_double("implexErrorTolerance", implex_error_tolerance))
				return nullptr;
			if (!lam_optional_double("implexTimeReductionLimit", implex_time_redution_limit))
				return nullptr;
		}
		else if (strcmp(value, "-implexAlpha") == 0) {
			if (!lam_optional_double("alpha", implex_alpha))
				return nullptr;
		}
		else if (strcmp(value, "-eta") == 0) {
			if (!lam_optional_double("eta", eta))
				return nullptr;
		}
		else if (strcmp(value, "-tangent") == 0) {
			tangent = true;
		}
		else if (strcmp(value, "-autoRegularization") == 0) {
			auto_regularization = true;
			if (OPS_GetNumRemainingInputArgs() < 1) {
				opserr << "nDMaterial ASDConcrete1D Error: '-autoRegularization' given without the next 1 argument $lch_ref.\n";
				return nullptr;
			}
			if (!lam_optional_double("lch_ref", lch_ref))
				return nullptr;
			have_lch_ref = true;
		}
		else if (strcmp(value, "-Te") == 0) {
			if (!lam_optional_list("Te", Te))
				return nullptr;
		}
		else if (strcmp(value, "-Ts") == 0) {
			if (!lam_optional_list("Ts", Ts))
				return nullptr;
		}
		else if (strcmp(value, "-Td") == 0) {
			if (!lam_optional_list("Td", Td))
				return nullptr;
		}
		else if (strcmp(value, "-Ce") == 0) {
			if (!lam_optional_list("Ce", Ce))
				return nullptr;
		}
		else if (strcmp(value, "-Cs") == 0) {
			if (!lam_optional_list("Cs", Cs))
				return nullptr;
		}
		else if (strcmp(value, "-Cd") == 0) {
			if (!lam_optional_list("Cd", Cd))
				return nullptr;
		}
	}

	// Set a default value of tension strength if none specified
	if (have_fc && !have_ft)
		ft = 0.1 * fc;

	if (have_fc) {
		double ec = 2 * fc / E;
		double Gt = 0.073 * pow(fc, 0.18);
		double Gc = 2 * Gt * (fc * fc) / (ft * ft);


		if (!have_lch_ref) {
			//
			// _get_lch_ref from ASDConcrete1D_MakeLaws.py
			//

			// min lch for tension
			double et_el = ft / E;
			double Gt_min = 0.5 * ft * et_el;
			double hmin_t = 0.01 * Gt / Gt_min;

			// min lch for compression
			double ec1 = fc / E;
			double ec_pl = (ec - ec1) * 0.4 + ec1;
			double Gc_min = 0.5 * fc * (ec - ec_pl);
			double hmin_c = 0.01 * Gc / Gc_min;

			lch_ref = std::min(hmin_c, hmin_t);
		}

		//
		// _make_tension from ASDConcrete1D_MakeLaws.py
		//

		Gt = Gt / lch_ref;

		double f0 = 0.9 * ft;
		double f1 = ft;
		double e0 = f0 / E;
		double e1 = 1.5 * f1 / E;
		double ep = e1 - f1 / E;
		double f2 = 0.2 * ft;
		double f3 = 1.0e-3 * ft;
		double w2 = Gt / ft;
		double w3 = 5.0 * w2;
		double e2 = w2 + f2 / E + ep;
		if (e2 <= e1)
			e2 = 1.001 * e1;
		double e3 = w3 + f3 / E + ep;
		if (e3 <= e2)
			e3 = 1.001 * e2;
		double e4 = 10.0 * e3;
		Te.resize(6); Te = { 0.0, e0, e1, e2, e3, e4 };
		Ts.resize(6); Ts = { 0.0, f0, f1, f2, f3, f3 };
		Td.resize(6); Td = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
		double Tpl[6] = { 0.0, 0.0, ep, 0.9 * e2, 0.8 * e3, 0.8 * e3 };

		for (int i = 2; i < 6; i++) {
			double xi = Te[i];
			double si = Ts[i];
			double xipl = Tpl[i];
			double xipl_max = xi - si / E;
			xipl = std::min(xipl, xipl_max);
			double qi = (xi - xipl) * E;
			Td[i] = 1.0 - si / qi;
		}

		//
		// _make_compression from ASDConcrete1D_MakeLaws.py
		//

		Gc = Gc / lch_ref;

		double fc0 = 0.5 * fc;
		double ec0 = fc0 / E;
		double ec1 = fc / E;
		double fcr = 0.1 * fc;
		double ec_pl = (ec - ec1) * 0.4 + ec1;
		double Gc1 = 0.5 * fc * (ec - ec_pl);
		double Gc2 = std::max(0.01 * Gc1, Gc - Gc1);
		double ecr = ec + 2.0 * Gc2 / (fc + fcr);
		const int nc = 10;
		Ce.resize(nc + 3); Ce[0] = 0.0; Ce[1] = ec0;
		Cs.resize(nc + 3); Cs[0] = 0.0; Cs[1] = fc0;
		double Cpl[nc + 3]; Cpl[0] = 0.0; Cpl[1] = 0.0;
		double dec = (ec - ec0) / (nc - 1);
		for (int i = 0; i < nc - 1; i++) {
			double iec = ec0 + (i + 1) * dec;
			Ce[i + 2] = iec;
			Cs[i + 2] = bezier3(iec, ec0, ec1, ec, fc0, fc, fc);
			Cpl[i + 2] = Cpl[i + 1] + 0.7 * (iec - Cpl[i + 1]);
		}
		Ce[nc + 1] = ecr;
		Cs[nc + 1] = fcr;
		Cpl[nc + 1] = Cpl[nc] + 0.7 * (ecr - Cpl[nc]);
		Ce[nc + 2] = ecr + ec0;
		Cs[nc + 2] = fcr;
		Cpl[nc + 2] = Cpl[nc + 1];
		Cd.resize(nc + 3); Cd[0] = 0.0; Cd[1] = 0.0;
		for (int i = 2; i < nc + 3; i++) {
			double xi = Ce[i];
			double si = Cs[i];
			double xipl = Cpl[i];
			double xipl_max = xi - si / E;
			xipl = std::min(xipl, xipl_max);
			double qi = (xi - xipl) * E;
			Cd[i] = 1.0 - si / qi;
		}
	}

	// check lists
	if (Te.size() < 1) {
		opserr << "nDMaterial ASDConcrete1D Error: 'Te' list is empty. At least 1 non-zero value should be provided.\n";
		return nullptr;
	}
	if (Ts.size() != Te.size()) {
		opserr << "nDMaterial ASDConcrete1D Error: 'Te' (size = " <<
			static_cast<int>(Te.size()) << ") and 'Ts' (size = " <<
			static_cast<int>(Ts.size()) << ") lists should have the same size.\n";
		return nullptr;
	}
	if (Td.size() == 0) {
		Td.resize(Te.size(), 0.0);
	}
	else if (Td.size() != Te.size()) {
		opserr << "nDMaterial ASDConcrete1D Error: 'Te' (size = " <<
			static_cast<int>(Te.size()) << ") and 'Td' (size = " <<
			static_cast<int>(Td.size()) << ") lists should have the same size.\n";
		return nullptr;
	}
	if (Ce.size() < 1) {
		opserr << "nDMaterial ASDConcrete1D Error: 'Tc' list is empty. At least 1 non-zero value should be provided.\n";
		return nullptr;
	}
	if (Cs.size() != Ce.size()) {
		opserr << "nDMaterial ASDConcrete1D Error: 'Ce' (size = " <<
			static_cast<int>(Ce.size()) << ") and 'Cs' (size = " <<
			static_cast<int>(Cs.size()) << ") lists should have the same size.\n";
		return nullptr;
	}
	if (Cd.size() == 0) {
		Cd.resize(Ce.size(), 0.0);
	}
	else if (Cd.size() != Ce.size()) {
		opserr << "nDMaterial ASDConcrete1D Error: 'Ce' (size = " <<
			static_cast<int>(Ce.size()) << ") and 'Cd' (size = " <<
			static_cast<int>(Cd.size()) << ") lists should have the same size.\n";
		return nullptr;
	}

	// build the hardening laws
	ASDConcrete1DMaterial::HardeningLaw HT(tag, ASDConcrete1DMaterial::HardeningLawType::Tension, E, Te, Ts, Td);
	if (!HT.isValid()) {
		opserr << "nDMaterial ASDConcrete1D Error: Tensile hardening law is not valid.\n";
		return nullptr;
	}
	ASDConcrete1DMaterial::HardeningLaw HC(tag, ASDConcrete1DMaterial::HardeningLawType::Compression, E, Ce, Cs, Cd);
	if (!HC.isValid()) {
		opserr << "nDMaterial ASDConcrete1D Error: Compressive hardening law is not valid.\n";
		return nullptr;
	}

	// create the material
	UniaxialMaterial* instance = new ASDConcrete1DMaterial(
		tag,
		E, eta,
		implex, implex_control, implex_error_tolerance, implex_time_redution_limit, implex_alpha,
		tangent, auto_regularization, lch_ref,
		HT, HC);
	if (instance == nullptr) {
		opserr << "UniaxialMaterial ASDConcrete1D Error: failed to allocate a new material.\n";
		return nullptr;
	}
	return instance;
}

ASDConcrete1DMaterial::HardeningLaw::HardeningLaw(
	int tag, HardeningLawType type,
	double E,
	const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& d)
	: m_tag(tag)
	, m_type(type)
{
	// initial checks
	if (!(x.size() > 0 && x.size() == y.size() && x.size() == d.size())) {
		opserr << "ASDConcrete1D Fatal Error: HardeningLaw::c-tor - found incompatible sizes.\n";
		return;
	}
	// fill it
	m_points.resize(x.size());
	// make sure they are positive values
	double xmax = 0.0;
	double ymax = 0.0;
	for (std::size_t i = 0; i < x.size(); ++i) {
		auto& pi = m_points[i];
		pi.x = std::abs(x[i]);
		pi.y = std::abs(y[i]);
		pi.d = std::min(1.0, std::max(0.0, d[i]));
		xmax = std::max(xmax, pi.x);
		ymax = std::max(ymax, pi.y);
	}
	if (xmax == 0.0) {
		opserr << "ASDConcrete1D Fatal Error: HardeningLaw::c-tor - max(X) == 0 " << xmax << "\n";
		return;
	}
	if (ymax == 0.0) {
		opserr << "ASDConcrete1D Fatal Error: HardeningLaw::c-tor - max(Y) == 0 " << ymax << "\n";
		return;
	}
	// check first values and, if needed, append (0,0,0)
	if (m_points[0].x > 0.0) {
		if (m_points[0].y == 0.0) {
			// first x is > 0, y cannot be 0, so we set it to E*x
			m_points[0].y = E * m_points[0].x;
		}
		else {
			// ok, the user did not define the first point at (0,0), we make it
			m_points.insert(m_points.begin(), HardeningLawPoint());
		}
	}
	else {
		// first x is 0, set y to zero as well
		m_points[0].y = 0.0;
	}
	// make sure the first damage is 0!
	m_points[0].d = 0.0;
	// make sure the slope of the first 2 points is E (note p1 is now 0,0)
	m_points[1].y = E * m_points[1].x;
	// define tolerances: find the min(dx or dy) > 0
	double dxmin = xmax;
	double dymin = ymax;
	for (std::size_t i = 1; i < m_points.size(); ++i) {
		const auto& pi = m_points[i];
		const auto& pold = m_points[i - 1];
		double dx = std::abs(pi.x - pold.x);
		double dy = std::abs(pi.y - pold.y);
		if (dx > 0.0)
			dxmin = std::min(dx, dxmin);
		if (dy > 0.0)
			dymin = std::min(dy, dymin);
	}
	m_xtolerance = 1.0e-6 * dxmin;
	m_ytolerance = 1.0e-6 * dymin;
	// make valid
	m_valid = true;
	// make a first adjustment
	adjust();
	// compute fracture energy. If the user selects the auto-regularization,
	// this will store the real (non-regularized) fracture energy
	computeFractureEnergy();
	// store it
	HardeningLawStorage::instance().store(*this);
}

void ASDConcrete1DMaterial::HardeningLaw::regularize(double lch, double lch_ref)
{
	// quick return if not valid
	if (!m_valid)
		return;
	// quick return if un-bounded (inf fracture energy), or invalid lch
	double lch_scale = lch > 0.0 ? lch_ref / lch : 0.0;
	if (!m_fracture_energy_is_bounded || lch_scale <= 0.0 || lch_scale == 1.0)
		return;
	// back to original
	deRegularize();
	// the initial fracture energy has been computed in the full constructor, and we
	// are back to it after deRegularize...
	// compute the required specific fracture energy
	double gnew = m_fracture_energy * lch_scale;
	// compute the minimum fracture energy (in case lch is too large)
	const auto& peak = m_points[m_softening_begin];
	double gmin = (peak.y * peak.x / 2.0) * 1.01; // make it 1% larger to ensure a monotonically increasing abscissa
	gnew = std::max(gnew, gmin);
	// iteratively scale the curve until the gnew
	double tol = 1.0e-3 * gnew;
	double dscale = gnew / m_fracture_energy;
	double E = m_points[1].y / m_points[1].x;
	double x0 = peak.x;
	double scale = dscale;
	static constexpr int max_iter = 10;
	for (int iter = 0; iter < max_iter; ++iter) {
		// scale points after the peak
		for (std::size_t i = m_softening_begin + 1; i < m_points.size(); ++i) {
			auto& pi = m_points[i];
			// save the plastic-to-inelastic ratio before
			double xi_inel = std::max(pi.x - pi.y / E, 0.0);
			double xi_pl = pi.x - pi.q / E;
			double xi_ratio = xi_inel > 0.0 ? xi_pl / xi_inel : 0.0;
			// update the total strain
			pi.x = x0 + (pi.x - x0) * dscale;
			// update damage to keep the same ratio
			xi_inel = std::max(pi.x - pi.y / E, 0.0);
			xi_pl = xi_inel * xi_ratio;
			pi.q = E * (pi.x - xi_pl);
			if (pi.q > 0)
				pi.d = 1.0 - pi.y / pi.q;
		}
		// update g and check
		adjust();
		computeFractureEnergy();
		if (std::abs(m_fracture_energy - gnew) < tol)
			break;
		// update scale
		dscale = gnew / m_fracture_energy;
		scale *= dscale;
	}
	// re-adjust.
	adjust();
}

void ASDConcrete1DMaterial::HardeningLaw::deRegularize()
{
	auto source = HardeningLawStorage::instance().recover(m_tag, m_type);
	if (source)
		*this = *source;
}

ASDConcrete1DMaterial::HardeningLawPoint ASDConcrete1DMaterial::HardeningLaw::evaluateAt(double x) const
{
	// quick return
	if (!m_valid)
		return HardeningLawPoint();
	// search for x
	double x1, x2, y1, y2, q1, q2;
	bool found = false;
	for (std::size_t i = 1; i < m_points.size(); ++i) {
		const auto& p1 = m_points[i - 1];
		const auto& p2 = m_points[i];
		if (x <= p2.x + m_xtolerance) {
			x1 = p1.x;
			x2 = p2.x;
			y1 = p1.y;
			y2 = p2.y;
			q1 = p1.q;
			q2 = p2.q;
			found = true;
			break;
		}
	}
	if (!found) {
		// we're beyond last point
		x1 = m_points.back().x;
		x2 = x;
		double span = x - x1;
		// interpolate last tangent if positive, otherwise keep constant
		y1 = m_points.back().y;
		double tangent = (y1 - m_points[m_points.size() - 2].y) / (x1 - m_points[m_points.size() - 2].x);
		y2 = tangent > 0.0 ? y1 + span * tangent : y1;
		// interpolate last tangent if positive, otherwise keep constant
		q1 = m_points.back().q;
		tangent = (q1 - m_points[m_points.size() - 2].q) / (x1 - m_points[m_points.size() - 2].x);
		q2 = tangent > 0.0 ? q1 + span * tangent : q1;
	}
	// interpolate
	double xspan = x2 - x1;
	double xratio = xspan > 0.0 ? (x - x1) / xspan : 0.0;
	double y = std::max(m_ytolerance, y1 + (y2 - y1) * xratio);
	double q = std::max(m_ytolerance, q1 + (q2 - q1) * xratio);
	double d = 1.0 - y / q;
	// done
	return HardeningLawPoint(x, y, d, q);
}

double ASDConcrete1DMaterial::HardeningLaw::computeMaxStress() const
{
	double smax = 0.0;
	for (const auto& p : m_points) {
		smax = std::max(smax, p.y);
	}
	return smax;
}

int ASDConcrete1DMaterial::HardeningLaw::serializationDataSize() const
{
	// number of points (variable, 4 components each)
	int np = static_cast<int>(m_points.size());
	// number of fixed data
	int nn = 10;
	// we need to save 2 copies
	return (nn + np * 4) * 2;
}

void ASDConcrete1DMaterial::HardeningLaw::serialize(Vector& data, int& pos)
{
	// internal
	auto lam = [&data, &pos](HardeningLaw& x) {
		data(pos++) = static_cast<double>(x.m_tag);
		data(pos++) = static_cast<double>(static_cast<int>(x.m_type));
		data(pos++) = static_cast<double>(x.m_points.size());
		data(pos++) = x.m_fracture_energy;
		data(pos++) = static_cast<double>(x.m_fracture_energy_is_bounded);
		data(pos++) = static_cast<double>(x.m_softening_begin);
		data(pos++) = static_cast<double>(x.m_softening_end);
		data(pos++) = static_cast<double>(x.m_valid);
		data(pos++) = x.m_xtolerance;
		data(pos++) = x.m_ytolerance;
		for (const auto& p : x.m_points) {
			data(pos++) = p.x;
			data(pos++) = p.y;
			data(pos++) = p.d;
			data(pos++) = p.q;
		}
	};
	// save the current
	lam(*this);
	// save the original
	HardeningLaw original(*this);
	original.deRegularize();
	lam(original);
}

void ASDConcrete1DMaterial::HardeningLaw::deserialize(Vector& data, int& pos)
{
	// internal
	auto lam = [&data, &pos](HardeningLaw& x) {
		x.m_tag = static_cast<int>(data(pos++));
		x.m_type = static_cast<HardeningLawType>(static_cast<int>(data(pos++)));
		x.m_points.resize(static_cast<std::size_t>(data(pos++)));
		x.m_fracture_energy = data(pos++);
		x.m_fracture_energy_is_bounded = static_cast<bool>(data(pos++));
		x.m_softening_begin = static_cast<std::size_t>(data(pos++));
		x.m_softening_end = static_cast<std::size_t>(data(pos++));
		x.m_valid = static_cast<bool>(data(pos++));
		x.m_xtolerance = data(pos++);
		x.m_ytolerance = data(pos++);
		for (auto& p : x.m_points) {
			p.x = data(pos++);
			p.y = data(pos++);
			p.d = data(pos++);
			p.q = data(pos++);
		}
	};
	// recover the current
	lam(*this);
	// recover the original
	HardeningLaw original;
	lam(original);
	// store if necessary
	HardeningLawStorage::instance().store(original);
}

void ASDConcrete1DMaterial::HardeningLaw::adjust()
{
	// quick return
	if (!m_valid)
		return;
	// get initial tangent
	double E = m_points[1].y / m_points[1].x;
	// check all points
	for (std::size_t i = 1; i < m_points.size(); ++i) {
		const auto& pold = m_points[i - 1];
		auto& pi = m_points[i];
		double xi = pi.x;
		double xold = pold.x;
		double yi = pi.y;
		double yold = pold.y;
		double di = pi.d;
		double dold = pold.d;
		// check: strictly monotonic x
		if (xi <= xold)
			xi += m_xtolerance;
		// check yi is not exactly 0.0
		if (yi < m_ytolerance)
			yi = m_ytolerance;
		// check tangent stifness is not higher than the initial
		double Ei = (yi - yold) / (xi - xold); // denom > 0 < --strictly monotonic x
		if (Ei > E)
			yi = yold + (xi - xold) * E;
		// check damage
		// current plastic strain cannot be lower then the old one
		double eepd_old = dold < 1.0 ? xold - yold / ((1.0 - dold) * E) : xold;
		double eepd = di < 1.0 ? xi - yi / ((1.0 - di) * E) : -1.0;
		if (eepd < eepd_old)
			di = std::min(1.0, std::max(0.0, 1.0 - yi / E / (xi - eepd_old)));
		// damage cannot decrease
		if (di < dold)
			di = dold;
		// last check: make sure (1-d)*E is not < tangent (Ei)
		Ei = (yi - yold) / (xi - xold);
		double Ed = (1.0 - di) * E;
		if (Ei > Ed)
			yi = Ed * (xi - eepd_old);
		// update
		pi.x = xi;
		pi.y = yi;
		pi.d = di;
		pi.q = yi / (1.0 - di);
	}
}

void ASDConcrete1DMaterial::HardeningLaw::computeFractureEnergy()
{
	// initialize as un-bounded
	m_fracture_energy = 0.0;
	m_fracture_energy_is_bounded = false;
	m_softening_begin = 0;
	m_softening_end = 0;
	// quick return
	if (!m_valid)
		return;
	// find the first point where the slope is negative
	std::size_t pos1 = 0;
	bool pos1_found = false;
	for (std::size_t i = 1; i < m_points.size(); ++i) {
		const auto& p1 = m_points[i - 1];
		const auto& p2 = m_points[i];
		double k = (p2.y - p1.y) / (p2.x - p1.x);
		if (k < 0.0) {
			pos1 = i - 1;
			pos1_found = true;
			break;
		}
	}
	// exit if not found -> infinite energy
	if (!pos1_found)
		return;
	// find the last point where the slope is negative
	std::size_t pos2 = 0;
	bool pos2_found = false;
	for (std::size_t i = pos1 + 1; i < m_points.size(); ++i) {
		const auto& p1 = m_points[i - 1];
		const auto& p2 = m_points[i];
		double k = (p2.y - p1.y) / (p2.x - p1.x);
		if (k >= 0.0) {
			pos2 = i - 1;
			pos2_found = true;
			break;
		}
	}
	// if the last point is not found, it means
	// that we need to extend the last portion
	double g_add = 0.0;
	if (!pos2_found) {
		const auto& p2 = m_points.back();
		if (p2.y > 0.0) {
			const auto& p1 = m_points[m_points.size() - 2];
			double k = (p2.y - p1.y) / (p2.x - p1.x);
			double x3 = p2.x - p2.y / k;
			g_add = p2.y * (x3 - p2.x) / 2.0;
		}
		pos2 = m_points.size() - 1;
	}
	// now compute g
	double g = 0.0;
	// first add the initial triangle based on the unloading
	// stiffness at the pos1 point
	double E = m_points[1].y / m_points[1].x;
	double Ed = (1.0 - m_points[pos1].d) * E;
	g += (std::pow(m_points[pos1].y, 2) / Ed / 2.0);
	// then add other components
	for (std::size_t i = pos1 + 1; i <= pos2; ++i) {
		const auto& p1 = m_points[i - 1];
		const auto& p2 = m_points[i];
		g += (p2.x - p1.x) * (p1.y + p2.y) / 2.0;
	}
	// finally add the final one
	g += g_add;
	// done. store values now
	m_fracture_energy = g;
	m_fracture_energy_is_bounded = true;
	m_softening_begin = pos1;
	m_softening_end = pos2;
}

ASDConcrete1DMaterial::HardeningLawStorage& ASDConcrete1DMaterial::HardeningLawStorage::instance()
{
	static HardeningLawStorage _instance;
	return _instance;
}

void ASDConcrete1DMaterial::HardeningLawStorage::store(const HardeningLaw& hl)
{
	if (hl.type() == HardeningLawType::Tension) {
		auto& item = m_tension[hl.tag()];
		if (item == nullptr)
			item = std::make_shared<HardeningLaw>(hl);
	}
	else {
		auto& item = m_compression[hl.tag()];
		if (item == nullptr)
			item = std::make_shared<HardeningLaw>(hl);
	}
}

ASDConcrete1DMaterial::HardeningLawStorage::PointerType ASDConcrete1DMaterial::HardeningLawStorage::recover(int tag, HardeningLawType type)
{
	if (type == HardeningLawType::Tension) {
		auto it = m_tension.find(tag);
		if (it != m_tension.end())
			return it->second;
	}
	else {
		auto it = m_compression.find(tag);
		if (it != m_compression.end())
			return it->second;
	}
	return nullptr;
}

ASDConcrete1DMaterial::ASDConcrete1DMaterial(
	int _tag,
	double _E,
	double _eta,
	bool _implex,
	bool _implex_control,
	double _implex_error_tolerance,
	double _implex_time_reduction_limit,
	double _implex_alpha,
	bool _tangent,
	bool _auto_regularize,
	double _lch_ref,
	const HardeningLaw& _ht,
	const HardeningLaw& _hc)
	: UniaxialMaterial(_tag, MAT_TAG_ASDConcrete1DMaterial)
	, E(_E)
	, eta(_eta)
	, implex(_implex)
	, implex_control(_implex_control)
	, implex_error_tolerance(_implex_error_tolerance)
	, implex_time_redution_limit(_implex_time_reduction_limit)
	, implex_alpha(_implex_alpha)
	, tangent(_tangent)
	, auto_regularize(_auto_regularize)
	, lch_ref(_lch_ref)
	, ht(_ht)
	, hc(_hc)
{
	// intialize C as C0
	C = getInitialTangent();

	// initialize PT_commit as eye(6)*0.5
	PT_commit = 0.5;
}

ASDConcrete1DMaterial::ASDConcrete1DMaterial()
	: UniaxialMaterial(0, MAT_TAG_ASDConcrete1DMaterial)
{
}

ASDConcrete1DMaterial::~ASDConcrete1DMaterial()
{
}

int ASDConcrete1DMaterial::setTrialStrain(double v, double r)
{
	// return value
	int retval = 0;

	// get characteristic length and perform regularization
	if (!regularization_done) {
		if (ops_TheActiveElement)
			lch = ops_TheActiveElement->getCharacteristicLength();
		regularization_done = true;
		if (auto_regularize) {
			ht.regularize(lch, lch_ref);
			hc.regularize(lch, lch_ref);
		}
	}

	// save dT
	if (!dtime_is_user_defined) {
		dtime_n = ops_Dt;
		if (!commit_done) {
			dtime_0 = dtime_n;
			dtime_n_commit = dtime_n;
		}
	}

	// if the user requested the tangent matrix
	// and not the IMPL-EX (in IMPL-EX the tangent coincides with the secant) ...
	if (tangent && !implex) {
		// numerical tangent tensor
		double Cnum = 0.0;
		// strain perturbation parameter
		double PERT = (ht.strainTolerance() + hc.strainTolerance()) / 2.0;
		// compute the forward perturbed solution and store in Cnum
		strain = v + PERT;
		retval = compute(true, false);
		if (retval < 0) return retval;
		Cnum = stress;
		// compute unperturbed solution
		strain = v;
		retval = compute(true, false);
		if (retval < 0) return retval;
		Cnum = (Cnum - stress) / PERT;
		// save tangent
		C = Cnum;
	}
	else {
		strain = v;
		if (implex) {
			if (implex_control) {
				// implicit solution
				double aux = PT_commit;
				retval = compute(false, false);
				if (retval < 0) return retval;
				double dt_implicit = dt_bar;
				double dc_implicit = dc_bar;
				PT_commit = aux;
				// standard call 
				retval = compute(true, true);
				if (retval < 0) return retval;
				double dt = dt_bar;
				double dc = dc_bar;
				implex_error = std::max(std::abs(dt - dt_implicit), std::abs(dc - dc_implicit));
				if (implex_error > implex_error_tolerance) {
					if (dtime_n >= implex_time_redution_limit * dtime_0) {
						retval = EC_IMPLEX_Error_Control;
					}
				}
			}
			else {
				// standard call 
				retval = compute(true, true);
			}
		}
		else {
			// standard call 
			retval = compute(true, true);
		}
	}

	// done
	return retval;
}

double ASDConcrete1DMaterial::getStress(void)
{
	return stress;
}

double ASDConcrete1DMaterial::getTangent(void)
{
	return C;
}

double ASDConcrete1DMaterial::getInitialTangent(void)
{
	return E;
}

double ASDConcrete1DMaterial::getStrain(void)
{
	return strain;
}

int ASDConcrete1DMaterial::commitState(void)
{
	// implicit stage
	if (implex) {
		// compute implex error here always
		double dt = dt_bar;
		double dc = dc_bar;
		// implicit solution
		compute(false, false);
		// compute implex error here always
		double dt_implicit = dt_bar;
		double dc_implicit = dc_bar;
		implex_error = std::max(std::abs(dt - dt_implicit), std::abs(dc - dc_implicit));
		GlobalParameters::instance().setMaxError(std::max(implex_error, GlobalParameters::instance().getMaxError()));
		GlobalParameters::instance().accumulateAverageError(implex_error);
	}
	// compute energy
	energy += 0.5 * (stress_commit + stress) * (strain - strain_commit);
	// store the previously committed variables for next move from n to n - 1
	xt_commit_old = xt_commit;
	xc_commit_old = xc_commit;
	// store committed variables
	xt_commit = xt;
	xc_commit = xc;
	strain_commit = strain;
	stress_commit = stress;
	stress_eff_commit = stress_eff;
	dtime_n_commit = dtime_n;
	// done
	commit_done = true;
	return 0;
}

int ASDConcrete1DMaterial::revertToLastCommit(void)
{
	// restore converged values
	xt = xt_commit;
	xc = xc_commit;
	strain = strain_commit;
	stress = stress_commit;
	stress_eff = stress_eff_commit;
	dtime_n = dtime_n_commit;
	// done
	return 0;
}

int ASDConcrete1DMaterial::revertToStart(void)
{
	// State variables
	xt = 0.0;
	xt_commit = 0.0;
	xt_commit_old = 0.0;
	xc = 0.0;
	xc_commit = 0.0;
	xc_commit_old = 0.0;

	// Time step
	dtime_n = 0.0;
	dtime_n_commit = 0.0;
	dtime_0 = 0.0;
	dtime_is_user_defined = false;

	// Commit flag
	commit_done = false;

	// IMPL-EX error
	implex_error = 0.0;

	// Strain, Stress and Tangent
	strain = 0.0;
	strain_commit = 0.0;
	stress = 0.0;
	stress_commit = 0.0;
	stress_eff = 0.0;
	stress_eff_commit = 0.0;
	C = getInitialTangent();
	PT_commit = 0.5;

	// Output variables
	dt_bar = 0.0;
	dc_bar = 0.0;
	energy = 0.0;

	// Done
	return 0;
}

UniaxialMaterial* ASDConcrete1DMaterial::getCopy(void)
{
	// we can safely use the default copy-constructor according to the member variables we're using
	return new ASDConcrete1DMaterial(*this);
}

void ASDConcrete1DMaterial::Print(OPS_Stream& s, int flag)
{
	s << "ASDConcrete1D Material, tag: " << this->getTag() << "\n";
}

int ASDConcrete1DMaterial::sendSelf(int commitTag, Channel &theChannel)
{
	// aux
	int counter;

	// variable DBL data size
	int nv_dbl = 28 +
		ht.serializationDataSize() +
		hc.serializationDataSize();

	// send INT data
	static ID idata(9);
	counter = 0;
	idata(counter++) = getTag();
	idata(counter++) = static_cast<int>(implex);
	idata(counter++) = static_cast<int>(implex_control);
	idata(counter++) = static_cast<int>(tangent);
	idata(counter++) = static_cast<int>(auto_regularize);
	idata(counter++) = static_cast<int>(regularization_done);
	idata(counter++) = static_cast<int>(dtime_is_user_defined);
	idata(counter++) = static_cast<int>(commit_done);
	idata(counter++) = nv_dbl;
	if (theChannel.sendID(getDbTag(), commitTag, idata) < 0) {
		opserr << "ASDConcrete1DMaterial::sendSelf() - failed to send INT data\n";
		return -1;
	}

	// send DBL data
	Vector ddata(nv_dbl);
	counter = 0;
	ddata(counter++) = E;
	ddata(counter++) = eta;
	ddata(counter++) = implex_error_tolerance;
	ddata(counter++) = implex_time_redution_limit;
	ddata(counter++) = implex_alpha;
	ddata(counter++) = lch;
	ddata(counter++) = lch_ref;
	ddata(counter++) = xt;
	ddata(counter++) = xt_commit;
	ddata(counter++) = xt_commit_old;
	ddata(counter++) = xc;
	ddata(counter++) = xc_commit;
	ddata(counter++) = xc_commit_old;
	ddata(counter++) = dtime_n;
	ddata(counter++) = dtime_n_commit;
	ddata(counter++) = dtime_0;
	ddata(counter++) = implex_error;
	ddata(counter++) = PT_commit;
	ddata(counter++) = strain;
	ddata(counter++) = strain_commit;
	ddata(counter++) = stress;
	ddata(counter++) = stress_commit;
	ddata(counter++) = stress_eff;
	ddata(counter++) = stress_eff_commit;
	ddata(counter++) = C;
	ddata(counter++) = dt_bar;
	ddata(counter++) = dc_bar;
	ddata(counter++) = energy;
	ht.serialize(ddata, counter);
	hc.serialize(ddata, counter);
	if (theChannel.sendVector(getDbTag(), commitTag, ddata) < 0) {
		opserr << "ASDConcrete1DMaterial::sendSelf() - failed to send DBL data\n";
		return -1;
	}

	// done
	return 0;
}

int ASDConcrete1DMaterial::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	// aux
	int counter;

	// recv INT data
	static ID idata(9);
	if (theChannel.recvID(getDbTag(), commitTag, idata) < 0) {
		opserr << "ASDConcrete1DMaterial::recvSelf() - failed to receive INT data\n";
		return -1;
	}
	counter = 0;
	setTag(idata(counter++));
	implex = static_cast<bool>(idata(counter++));
	implex_control = static_cast<bool>(idata(counter++));
	tangent = static_cast<bool>(idata(counter++));
	auto_regularize = static_cast<bool>(idata(counter++));
	regularization_done = static_cast<bool>(idata(counter++));
	dtime_is_user_defined = static_cast<bool>(idata(counter++));
	commit_done = static_cast<bool>(idata(counter++));
	int nv_dbl = idata(counter++);

	// recv DBL data
	Vector ddata(nv_dbl);
	if (theChannel.recvVector(getDbTag(), commitTag, ddata) < 0) {
		opserr << "ASDConcrete1DMaterial::recvSelf() - failed to receive DBL data\n";
		return -1;
	}
	counter = 0;
	E = ddata(counter++);
	eta = ddata(counter++);
	implex_error_tolerance = ddata(counter++);
	implex_time_redution_limit = ddata(counter++);
	implex_alpha = ddata(counter++);
	lch = ddata(counter++);
	lch_ref = ddata(counter++);
	xt = ddata(counter++);
	xt_commit = ddata(counter++);
	xt_commit_old = ddata(counter++);
	xc = ddata(counter++);
	xc_commit = ddata(counter++);
	xc_commit_old = ddata(counter++);
	dtime_n = ddata(counter++);
	dtime_n_commit = ddata(counter++);
	dtime_0 = ddata(counter++);
	implex_error = ddata(counter++);
	PT_commit = ddata(counter++);
	strain = ddata(counter++);
	strain_commit = ddata(counter++);
	stress = ddata(counter++);
	stress_commit = ddata(counter++);
	stress_eff = ddata(counter++);
	stress_eff_commit = ddata(counter++);
	C = ddata(counter++);
	dt_bar = ddata(counter++);
	dc_bar = ddata(counter++);
	energy = ddata(counter++);
	ht.deserialize(ddata, counter);
	hc.deserialize(ddata, counter);

	// done
	return 0;
}

int ASDConcrete1DMaterial::setParameter(const char** argv, int argc, Parameter& param)
{
	// 1000 - elasticity & mass
	if (strcmp(argv[0], "E") == 0) {
		param.setValue(E);
		return param.addObject(1000, this);
	}

	// 2000 - time
	if (strcmp(argv[0], "dTime") == 0) {
		param.setValue(dtime_n);
		return param.addObject(2000, this);
	}
	if (strcmp(argv[0], "dTimeCommit") == 0) {
		param.setValue(dtime_n_commit);
		return param.addObject(2001, this);
	}
	if (strcmp(argv[0], "dTimeInitial") == 0) {
		param.setValue(dtime_0);
		return param.addObject(2002, this);
	}

	// 3000 - globals
	if (strcmp(argv[0], "implexError") == 0 || strcmp(argv[0], "ImplexError") == 0) {
		param.setValue(GlobalParameters::instance().getMaxError());
		return param.addObject(3000, this);
	}
	if (strcmp(argv[0], "avgImplexError") == 0 || strcmp(argv[0], "AvgImplexError") == 0) {
		param.setValue(GlobalParameters::instance().getAverageError());
		return param.addObject(3001, this);
	}

	// default
	return -1;
}

int ASDConcrete1DMaterial::updateParameter(int parameterID, Information& info)
{
	switch (parameterID) {
		// 1000 - elasticity & mass
	case 1000:
		E = info.theDouble;
		return 0;

		// 2000 - time
	case 2000:
		dtime_n = info.theDouble;
		dtime_is_user_defined = true;
		return 0;
	case 2001:
		dtime_n_commit = info.theDouble;
		dtime_is_user_defined = true;
		return 0;
	case 2002:
		dtime_0 = info.theDouble;
		dtime_is_user_defined = true;
		return 0;

	case 3000:
		GlobalParameters::instance().setMaxError(info.theDouble);
		return 0;
	case 3001:
		GlobalParameters::instance().setAverageError(info.theDouble);
		return 0;

		// default
	default:
		return -1;
	}
}

Response* ASDConcrete1DMaterial::setResponse(const char** argv, int argc, OPS_Stream& output)
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
	static std::vector<std::string> lb_damage = { "d+", "d-" };
	static std::vector<std::string> lb_eqpl_strain = { "PLE+", "PLE-" };
	static std::vector<std::string> lb_tot_strain = { "TE+", "TE-" };
	static std::vector<std::string> lb_cw = { "cw" };
	static std::vector<std::string> lb_crackpattern = { "C1x", "C1y", "C1z",   "C2x", "C2y", "C2z",   "C3x", "C3y", "C3z" };
	static std::vector<std::string> lb_implex_error = { "Error" };
	static std::vector<std::string> lb_time = { "dTime", "dTimeCommit", "dTimeInitial" };
	static std::vector<std::string> lb_crack_strain = { "CS+", "LchRef" };
	static std::vector<std::string> lb_crush_strain = { "CS-", "LchRef" };
	static Vector Cinfo(2);

	// check specific responses
	if (argc > 0) {
		// 1000 - compressive hardening variables
		if (strcmp(argv[0], "Ce") == 0)
			return make_resp(1000, getHardeningLawVector(HardeningLawType::Compression, HardeningLawPointComponent::TotalStrain));
		if (strcmp(argv[0], "Cs") == 0)
			return make_resp(1001, getHardeningLawVector(HardeningLawType::Compression, HardeningLawPointComponent::NominalStress));
		if (strcmp(argv[0], "Cq") == 0)
			return make_resp(1002, getHardeningLawVector(HardeningLawType::Compression, HardeningLawPointComponent::EffectiveStress));
		// 1100 - tensile hardening variables
		if (strcmp(argv[0], "Te") == 0)
			return make_resp(1100, getHardeningLawVector(HardeningLawType::Tension, HardeningLawPointComponent::TotalStrain));
		if (strcmp(argv[0], "Ts") == 0)
			return make_resp(1101, getHardeningLawVector(HardeningLawType::Tension, HardeningLawPointComponent::NominalStress));
		if (strcmp(argv[0], "Tq") == 0)
			return make_resp(1102, getHardeningLawVector(HardeningLawType::Tension, HardeningLawPointComponent::EffectiveStress));
		// 2000 - damage variables
		if (strcmp(argv[0], "damage") == 0 || strcmp(argv[0], "Damage") == 0)
			return make_resp(2000, getDamage(), &lb_damage);
		if (strcmp(argv[0], "equivalentPlasticStrain") == 0 || strcmp(argv[0], "EquivalentPlasticStrain") == 0)
			return make_resp(2001, getEquivalentPlasticStrain(), &lb_eqpl_strain);
		if (strcmp(argv[0], "equivalentTotalStrain") == 0 || strcmp(argv[0], "EquivalentTotalStrain") == 0)
			return make_resp(2002, getStrainMeasure(), &lb_tot_strain);
		if (strcmp(argv[0], "cw") == 0 || strcmp(argv[0], "crackWidth") == 0 || strcmp(argv[0], "CrackWidth") == 0)
			return make_resp(2003, getCrackWidth(), &lb_cw);
		if (strcmp(argv[0], "crackStrain") == 0 || strcmp(argv[0], "CrackStrain") == 0) {
			if (argc > 2 && strcmp(argv[1], "-lchRef") == 0) {
				double lch_ref = 0.0;
				if (string_to_double(argv[2], lch_ref)) {
					Cinfo(1) = lch_ref;
					Cinfo(0) = getCrackWidth()(0) / lch_ref;
					return make_resp(2004, Cinfo, &lb_crack_strain);
				}
			}
		}
		if (strcmp(argv[0], "crushStrain") == 0 || strcmp(argv[0], "CrushStrain") == 0) {
			if (argc > 2 && strcmp(argv[1], "-lchRef") == 0) {
				double lch_ref = 0.0;
				if (string_to_double(argv[2], lch_ref)) {
					Cinfo(1) = lch_ref;
					Cinfo(0) = getCrushWidth()(0) / lch_ref;
					return make_resp(2005, Cinfo, &lb_crush_strain);
				}
			}
		}
		// 3000 - implex error
		if (strcmp(argv[0], "implexError") == 0 || strcmp(argv[0], "ImplexError") == 0) {
			return make_resp(3000, getImplexError(), &lb_implex_error);
		}
		// 4000 - internal time
		if (strcmp(argv[0], "time") == 0 || strcmp(argv[0], "Time") == 0) {
			return make_resp(4000, getTimeIncrements(), &lb_time);
		}
	}

	// otherwise return base-class response
	return UniaxialMaterial::setResponse(argv, argc, output);
}

int ASDConcrete1DMaterial::getResponse(int responseID, Information& matInformation)
{
	switch (responseID) {
		// 1000 - compressive hardening variables
	case 1000: return matInformation.setVector(getHardeningLawVector(HardeningLawType::Compression, HardeningLawPointComponent::TotalStrain));
	case 1001: return matInformation.setVector(getHardeningLawVector(HardeningLawType::Compression, HardeningLawPointComponent::NominalStress));
	case 1002: return matInformation.setVector(getHardeningLawVector(HardeningLawType::Compression, HardeningLawPointComponent::EffectiveStress));
		// 1100 - tensile hardening variables
	case 1100: return matInformation.setVector(getHardeningLawVector(HardeningLawType::Tension, HardeningLawPointComponent::TotalStrain));
	case 1101: return matInformation.setVector(getHardeningLawVector(HardeningLawType::Tension, HardeningLawPointComponent::NominalStress));
	case 1102: return matInformation.setVector(getHardeningLawVector(HardeningLawType::Tension, HardeningLawPointComponent::EffectiveStress));
		// 2000 - damage variables
	case 2000: return matInformation.setVector(getDamage());
	case 2001: return matInformation.setVector(getEquivalentPlasticStrain());
	case 2002: return matInformation.setVector(getStrainMeasure());
	case 2003: return matInformation.setVector(getCrackWidth());
	case 2004:
		if (matInformation.theVector && matInformation.theVector->Size() == 2) {
			double lch_ref = matInformation.theVector->operator()(1);
			matInformation.theVector->operator()(0) = getCrackWidth()(0) / lch_ref;
			return 0;
		}
		break;
	case 2005:
		if (matInformation.theVector && matInformation.theVector->Size() == 2) {
			double lch_ref = matInformation.theVector->operator()(1);
			matInformation.theVector->operator()(0) = getCrushWidth()(0) / lch_ref;
			return 0;
		}
		break;
		// 3000 - implex error
	case 3000: return matInformation.setVector(getImplexError());
		// 4000 - internal time
	case 4000: return matInformation.setVector(getTimeIncrements());
	default:
		break;
	}
	return UniaxialMaterial::getResponse(responseID, matInformation);
}

double ASDConcrete1DMaterial::getEnergy(void)
{
	return energy;
}

int ASDConcrete1DMaterial::compute(bool do_implex, bool do_tangent)
{
	// get committed variables
	xt = xt_commit;
	xc = xc_commit;
	stress = stress_commit;
	stress_eff = stress_eff_commit;

	// time factor for explicit extrapolation
	double time_factor = 1.0;
	if (implex && do_implex && (dtime_n_commit > 0.0))
		time_factor = dtime_n / dtime_n_commit * implex_alpha;

	// compute rate coefficients
	double rate_coeff_1 = 0.0;
	double rate_coeff_2 = 1.0;
	if ((dtime_n > 0.0) && (eta > 0.0)) {
		rate_coeff_1 = eta / (eta + dtime_n);
		rate_coeff_2 = dtime_n / (eta + dtime_n);
	}

	// compute elastic effective stress: SEFFn = C0 : (En - En-1)
	double dStrain = strain - strain_commit;
	stress_eff += E * dStrain;

	// compute stress split
	double PT = (implex && do_implex) ? PT_commit : Heavyside(stress_eff);
	double PC = 1.0 - PT;
	double ST = PT * stress_eff;
	double SC = PC * stress_eff;

	// compute committed hardening variables
	HardeningLawPoint pt = ht.evaluateAt(xt);
	HardeningLawPoint pc = hc.evaluateAt(xc);

	// temporary clone of old equivalent plastic strains
	double xt_pl = pt.plasticStrain(E);
	double xc_pl = pc.plasticStrain(E);

	// compute new trial equivalent strain measures
	if (implex && do_implex) {
		// extrapolated equivalent strain measures (explicit)
		xt = xt_commit + time_factor * (xt_commit - xt_commit_old);
		xc = xc_commit + time_factor * (xc_commit - xc_commit_old);
	}
	else {
		// compute trial strain measures (implicit)
		double xt_trial = ST / E + xt_pl;
		double xc_trial = -SC / E + xc_pl;
		// update hardening variables
		if (xt_trial > pt.x) 
			xt = rate_coeff_1 * pt.x + rate_coeff_2 * xt_trial;
		if (xc_trial > pc.x) 
			xc = rate_coeff_1 * pc.x + rate_coeff_2 * xc_trial;
	}
	pt = ht.evaluateAt(xt);
	pc = hc.evaluateAt(xc);

	// compute plastic damage
	double seff_eq_t = (pt.x - xt_pl) * E;
	double dt_plastic = seff_eq_t > 0.0 ? 1.0 - pt.q / seff_eq_t : 0.0;
	double seff_eq_c = (pc.x - xc_pl) * E;
	double dc_plastic = seff_eq_c > 0.0 ? 1.0 - pc.q / seff_eq_c : 0.0;

	// update effective stress
	stress_eff = (1.0 - dt_plastic) * ST + (1 - dc_plastic) * SC;

	// update nominal stress
	dt_bar = pt.d + dt_plastic - pt.d * dt_plastic;
	dc_bar = pc.d + dc_plastic - pc.d * dc_plastic;
	stress = (1.0 - dt_bar) * ST + (1.0 - dc_bar) * SC;

	// secant matrix
	if (do_tangent) {
		double W = 1.0 - dt_bar * PT - dc_bar * PC;
		C = W * E;
	}

	// save real PT and R, if mp.implex and !do_implex -> called from commit
	if (implex && !do_implex) {
		// save it in implex mode during implicit phase
		PT_commit = PT;
	}

	// done
	return 0;
}

Vector ASDConcrete1DMaterial::getHardeningLawVector(HardeningLawType ltype, HardeningLawPointComponent c) const
{
	Vector r;
	const HardeningLaw& h = ltype == HardeningLawType::Tension ? ht : hc;
	r.resize(static_cast<int>(h.points().size()));
	for (std::size_t i = 0; i < h.points().size(); ++i) {
		const HardeningLawPoint& p = h.points()[i];
		switch (c)
		{
		case ASDConcrete1DMaterial::HardeningLawPointComponent::TotalStrain:
			r(static_cast<int>(i)) = p.totalStrain();
			break;
		case ASDConcrete1DMaterial::HardeningLawPointComponent::EffectiveStress:
			r(static_cast<int>(i)) = p.effectiveStress();
			break;
		case ASDConcrete1DMaterial::HardeningLawPointComponent::NominalStress:
			r(static_cast<int>(i)) = p.stress();
			break;
		default:
			break;
		}
	}
	return r;
}

const Vector& ASDConcrete1DMaterial::getStrainMeasure() const
{
	static Vector d(2);
	d(0) = xt;
	d(1) = xc;
	return d;
}

const Vector& ASDConcrete1DMaterial::getDamage() const
{
	static Vector d(2);
	const Vector& x = getStrainMeasure();
	d(0) = ht.evaluateAt(x(0)).crackingDamage();
	d(1) = hc.evaluateAt(x(1)).crackingDamage();
	return d;
}

const Vector& ASDConcrete1DMaterial::getEquivalentPlasticStrain() const
{
	static Vector d(2);
	const Vector& x = getStrainMeasure();
	d(0) = ht.evaluateAt(x(0)).plasticStrain(E);
	d(1) = hc.evaluateAt(x(1)).plasticStrain(E);
	return d;
}

const Vector& ASDConcrete1DMaterial::getCrackWidth() const
{
	static Vector d(1);
	d.Zero();
	if (ht.hasStrainSoftening()) {
		double e0 = ht.strainAtOnsetOfCrack();
		const Vector& x = getStrainMeasure();
		d(0) = std::max(x(0) - e0, 0.0) * lch;
	}
	return d;
}

const Vector& ASDConcrete1DMaterial::getCrushWidth() const
{
	static Vector d(1);
	d.Zero();
	if (hc.hasStrainSoftening()) {
		double e0 = hc.strainAtOnsetOfCrack();
		const Vector& x = getStrainMeasure();
		d(0) = std::max(x(1) - e0, 0.0) * lch;
	}
	return d;
}

const Vector& ASDConcrete1DMaterial::getImplexError() const
{
	static Vector d(1);
	d(0) = implex_error;
	return d;
}

const Vector& ASDConcrete1DMaterial::getTimeIncrements() const
{
	static Vector d(3);
	d(0) = dtime_n;
	d(1) = dtime_n_commit;
	d(2) = dtime_0;
	return d;
}

