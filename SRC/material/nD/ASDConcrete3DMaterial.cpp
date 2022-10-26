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
                                                                        
// $Revision: 1.4 $
// $Date: 2020-04-19 23:01:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ASDConcrete3DMaterial.h,v $

// Massimo Petracca, Guido Camata - ASDEA Software, Italy
//
// A Simple and robust plastic-damage model for concrete and masonry
//

#include <ASDConcrete3DMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <Element.h>
#include <algorithm>
#include <limits>
#include <string>
#include <sstream>
#include <iomanip>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#endif

#define _FFMT(X) std::setw(10) << std::setprecision(3) << X

// anonymous namespace for utilities
namespace {

	enum ErrorCodes {
		EC_Generic = -1,
		EC_IMPLEX_Error_Control = -10,
		EC_Eigen_Error = -1000
	};

	/**
	Converts a string into a vector of doubles using whitespace as delimiter
	*/
	bool string_to_list_of_doubles(const std::string& text, char sep, std::vector<double>& out) {
		auto to_double = [](const std::string& text, double& num) -> bool {
			num = 0.0;
			try {
				num = std::stod(text);
				return true;
			}
			catch (...) {
				return false;
			}
		};
		if (out.size() > 0) out.clear();
		std::size_t start = 0, end = 0;
		double value;
		while (true) {
			end = text.find(sep, start);
			if (end == std::string::npos) {
				if (start < text.size()) {
					if (!to_double(text.substr(start), value))
						return false;
					out.push_back(value);
				}
				break;
			}
			std::string subs = text.substr(start, end - start);
			if (subs.size() > 0) {
				if (!to_double(subs, value))
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

	/*
	Taken from Eigen3 in Matrix class. We need it because we have to access eigenvectors.
	#TODO: Should we extend it in Matrix class?
	*/
	int Eigen3(Vector& d, Matrix& v) {
		if (v.noRows() != 3 || v.noCols() != 3 || d.Size() != 3)
			return -1;

		int     rot, its, i, j, k;
		double  g, h, aij, sm, thresh, t, c, s, tau;

		static Vector  a(3);
		static Vector  b(3);
		static Vector  z(3);

		static const double tol = 1.0e-08;

		//.... move array into one-d arrays
		a(0) = v(0, 1);
		a(1) = v(1, 2);
		a(2) = v(2, 0);

		for (i = 0; i < 3; i++) {
			d(i) = v(i, i);
			b(i) = v(i, i);
			z(i) = 0.0;

			for (j = 0; j < 3; j++)
				v(i, j) = 0.0;

			v(i, i) = 1.0;

		} //end for i

		rot = 0;
		its = 0;

		sm = fabs(a(0)) + fabs(a(1)) + fabs(a(2));

		while (sm > tol) {
			//.... set convergence test and threshold
			if (its < 3)
				thresh = 0.011 * sm;
			else
				thresh = 0.0;

			//.... perform sweeps for rotations
			for (i = 0; i < 3; i++) {

				j = (i + 1) % 3;
				k = (j + 1) % 3;

				aij = a(i);

				g = 100.0 * fabs(aij);

				if (fabs(d(i)) + g != fabs(d(i)) ||
					fabs(d(j)) + g != fabs(d(j))) {

					if (fabs(aij) > thresh) {

						a(i) = 0.0;
						h = d(j) - d(i);

						if (fabs(h) + g == fabs(h))
							t = aij / h;
						else {
							//t = 2.0 * sign(h/aij) / ( fabs(h/aij) + sqrt(4.0+(h*h/aij/aij)));
							double hDIVaij = h / aij;
							if (hDIVaij > 0.0)
								t = 2.0 / (hDIVaij + sqrt(4.0 + (hDIVaij * hDIVaij)));
							else
								t = -2.0 / (-hDIVaij + sqrt(4.0 + (hDIVaij * hDIVaij)));
						}

						//.... set rotation parameters

						c = 1.0 / sqrt(1.0 + t * t);
						s = t * c;
						tau = s / (1.0 + c);

						//.... rotate diagonal terms

						h = t * aij;
						z(i) = z(i) - h;
						z(j) = z(j) + h;
						d(i) = d(i) - h;
						d(j) = d(j) + h;

						//.... rotate off-diagonal terms

						h = a(j);
						g = a[k];
						a(j) = h + s * (g - h * tau);
						a(k) = g - s * (h + g * tau);

						//.... rotate eigenvectors

						for (k = 0; k < 3; k++) {
							g = v(k, i);
							h = v(k, j);
							v(k, i) = g - s * (h + g * tau);
							v(k, j) = h + s * (g - h * tau);
						} // end for k

						rot = rot + 1;

					} // end if fabs > thresh 
				} //else
				else
					a(i) = 0.0;

			}  // end for i

			//.... update the diagonal terms
			for (i = 0; i < 3; i++) {
				b(i) = b(i) + z(i);
				d(i) = b(i);
				z(i) = 0.0;
			} // end for i

			its += 1;

			sm = fabs(a(0)) + fabs(a(1)) + fabs(a(2));

		} //end while sm
		
		// sort in descending order (unrolled bubble sort)
		auto sortij = [&d, &v](int i, int j) {
			if (d(i) < d(j)) {
				std::swap(d(i), d(j));
				for (int k = 0; k < 3; ++k)
					std::swap(v(k, i), v(k, j));
			}
		};
		sortij(0, 1);
		sortij(1, 2);
		sortij(0, 1);

		// done
		return 0;
	}

	/*
	computes the vj x vj x vj x vj tensor in voigt notation
	*/
	void computePjj(const Matrix& V, int j, Matrix& pjj) {
		double A1 = V(2, j) * V(2, j);
		double A2 = V(1, j) * V(1, j);
		double A3 = V(0, j) * V(0, j);
		double A4 = V(2, j);
		double A5 = V(0, j);
		double A6 = V(1, j);
		double A7 = 2.0 * A1 * A5 * A6;
		double A8 = 2.0 * A3 * A4 * A6;
		double A9 = 2.0 * A2 * A4 * A5;
		double A10 = A1 * A2;
		pjj(0, 0) = A3 * A3;
		pjj(0, 1) = A2 * A3;
		pjj(0, 2) = A1 * A3;
		pjj(0, 3) = 2.0 * A3 * A5 * A6;
		pjj(0, 4) = A8;
		pjj(0, 5) = 2.0 * A3 * A4 * A5;
		pjj(1, 0) = A2 * A3;
		pjj(1, 1) = A2 * A2;
		pjj(1, 2) = A10;
		pjj(1, 3) = 2.0 * A2 * A5 * A6;
		pjj(1, 4) = 2.0 * A2 * A4 * A6;
		pjj(1, 5) = A9;
		pjj(2, 0) = A1 * A3;
		pjj(2, 1) = A10;
		pjj(2, 2) = A1 * A1;
		pjj(2, 3) = A7;
		pjj(2, 4) = 2.0 * A1 * A4 * A6;
		pjj(2, 5) = 2.0 * A1 * A4 * A5;
		pjj(3, 0) = A3 * A5 * A6;
		pjj(3, 1) = A2 * A5 * A6;
		pjj(3, 2) = A1 * A5 * A6;
		pjj(3, 3) = 2.0 * A2 * A3;
		pjj(3, 4) = A9;
		pjj(3, 5) = A8;
		pjj(4, 0) = A3 * A4 * A6;
		pjj(4, 1) = A2 * A4 * A6;
		pjj(4, 2) = A1 * A4 * A6;
		pjj(4, 3) = A9;
		pjj(4, 4) = 2.0 * A1 * A2;
		pjj(4, 5) = A7;
		pjj(5, 0) = A3 * A4 * A5;
		pjj(5, 1) = A2 * A4 * A5;
		pjj(5, 2) = A1 * A4 * A5;
		pjj(5, 3) = A8;
		pjj(5, 4) = A7;
		pjj(5, 5) = 2.0 * A1 * A3;
	}

}

void *OPS_ASDConcrete3DMaterial(void)
{
	// some kudos
	static bool first_done = false;
	if (!first_done) {
		opserr << "Using ASDConcrete3DM - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
		first_done = true;
	}

	// check arguments
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 3) {
		opserr << 
			"nDMaterial ASDConcrete3D Error: Few arguments (< 3).\n"
			"nDMaterial ASDConcrete3D $tag $E $v "
			"-Te $Te -Ts $Ts <-Td $Td> -Cs $Cs <-Cd $Cd> "
			"<-rho $rho> "
			"<-implex> <-implexControl $implexErrorTolerance $implexTimeReductionLimit>"
			"<-crackPlanes $nct $ncc $smoothingAngle>"
			"<-eta $eta> <-tangent> <-autoRegularization>\n";
		return nullptr;
	}

	// numData
	int numData = 1;

	// data
	int tag;
	double E;
	double v;
	double rho = 0.0;
	bool implex = false;
	bool implex_control = false;
	double implex_error_tolerance = 0.05;
	double implex_time_redution_limit = 0.01;
	double eta = 0.0;
	bool tangent = false;
	bool auto_regularization = false;
	std::vector<double> Te, Ts, Td, Ce, Cs, Cd;
	int nct = 0;
	int ncc = 0;
	double smoothing_angle = 45.0;
	
	// get tag
	if (OPS_GetInt(&numData, &tag) != 0)  {
		opserr << "nDMaterial ASDConcrete3D Error: invalid 'tag'.\n";
		return nullptr;
	}

	// get Elasticity arguments
	if (OPS_GetDouble(&numData, &E) != 0) {
		opserr << "nDMaterial ASDConcrete3D Error: invalid 'E'.\n";
		return nullptr;
	}
	if (E <= 0.0) {
		opserr << "nDMaterial ASDConcrete3D Error: invalid value for 'E' (" << E << "). It should be strictly positive.\n";
		return nullptr;
	}
	if (OPS_GetDouble(&numData, &v) != 0) {
		opserr << "nDMaterial ASDConcrete3D Error: invalid 'v'.\n";
		return nullptr;
	}

	// utilities (code re-use)
	auto lam_optional_int = [&numData](const char* variable, int& value) -> bool {
		if (OPS_GetNumRemainingInputArgs() > 0) {
			if (OPS_GetInt(&numData, &value) < 0) {
				opserr << "nDMaterial ASDConcrete3D Error: failed to get '" << variable << "'.\n";
				return false;
			}
		}
		else {
			opserr << "nDMaterial ASDConcrete3D Error: '" << variable << "' requested but not provided.\n";
			return false;
		}
		return true;
	};
	auto lam_optional_double = [&numData](const char* variable, double& value) -> bool {
		if (OPS_GetNumRemainingInputArgs() > 0) {
			if (OPS_GetDouble(&numData, &value) < 0) {
				opserr << "nDMaterial ASDConcrete3D Error: failed to get '" << variable << "'.\n";
				return false;
			}
		}
		else {
			opserr << "nDMaterial ASDConcrete3D Error: '" << variable << "' requested but not provided.\n";
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
				opserr << "nDMaterial ASDConcrete3D Error: cannot parse the '" << variable << "' list.\n";
				return false;
			}
		}
		return true;
	};

	// optional parameters
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* value = OPS_GetString();
		if (strcmp(value, "-rho") == 0) {
			if (!lam_optional_double("rho", rho))
				return nullptr;
		}
		else if (strcmp(value, "-implex") == 0) {
			implex = true;
		}
		else if (strcmp(value, "-implexControl") == 0) {
			implex_control = true;
			if (OPS_GetNumRemainingInputArgs() < 2) {
				opserr << "nDMaterial ASDConcrete3D Error: '-implexControl' given without the next 2 arguments $implexErrorTolerance $implexTimeReductionLimit.\n";
				return nullptr;
			}
			if (!lam_optional_double("implexErrorTolerance", implex_error_tolerance))
				return nullptr;
			if (!lam_optional_double("implexTimeReductionLimit", implex_time_redution_limit))
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
		else if (strcmp(value, "-crackPlanes") == 0) {
			implex_control = true;
			if (OPS_GetNumRemainingInputArgs() < 3) {
				opserr << "nDMaterial ASDConcrete3D Error: '-crackPlanes' given without the next 3 arguments $nct $ncc and $smoothingAngle.\n";
				return nullptr;
			}
			if (!lam_optional_int("nct", nct))
				return nullptr;
			if (!lam_optional_int("ncc", ncc))
				return nullptr;
			if (!lam_optional_double("$smoothingAngle", smoothing_angle))
				return nullptr;
		}
	}

	// check lists
	if (Te.size() < 1) {
		opserr << "nDMaterial ASDConcrete3D Error: 'Te' list is empty. At least 1 non-zero value should be provided.\n";
		return nullptr;
	}
	if (Ts.size() != Te.size()) {
		opserr << "nDMaterial ASDConcrete3D Error: 'Te' (size = " <<
			static_cast<int>(Te.size()) << ") and 'Ts' (size = " <<
			static_cast<int>(Ts.size()) << ") lists should have the same size.\n";
		return nullptr;
	}
	if (Td.size() == 0) {
		Td.resize(Te.size(), 0.0);
	}
	else if (Td.size() != Te.size()) {
		opserr << "nDMaterial ASDConcrete3D Error: 'Te' (size = " <<
			static_cast<int>(Te.size()) << ") and 'Td' (size = " <<
			static_cast<int>(Td.size()) << ") lists should have the same size.\n";
		return nullptr;
	}
	if (Ce.size() < 1) {
		opserr << "nDMaterial ASDConcrete3D Error: 'Tc' list is empty. At least 1 non-zero value should be provided.\n";
		return nullptr;
	}
	if (Cs.size() != Ce.size()) {
		opserr << "nDMaterial ASDConcrete3D Error: 'Ce' (size = " <<
			static_cast<int>(Ce.size()) << ") and 'Cs' (size = " <<
			static_cast<int>(Cs.size()) << ") lists should have the same size.\n";
		return nullptr;
	}
	if (Cd.size() == 0) {
		Cd.resize(Ce.size(), 0.0);
	}
	else if (Cd.size() != Ce.size()) {
		opserr << "nDMaterial ASDConcrete3D Error: 'Ce' (size = " <<
			static_cast<int>(Ce.size()) << ") and 'Cd' (size = " <<
			static_cast<int>(Cd.size()) << ") lists should have the same size.\n";
		return nullptr;
	}

	// build the hardening laws
	ASDConcrete3DMaterial::HardeningLaw HT(tag, ASDConcrete3DMaterial::HardeningLawType::Tension, E, Te, Ts, Td);
	if (!HT.isValid()) {
		opserr << "nDMaterial ASDConcrete3D Error: Tensile hardening law is not valid.\n";
		return nullptr;
	}
	ASDConcrete3DMaterial::HardeningLaw HC(tag, ASDConcrete3DMaterial::HardeningLawType::Compression, E, Ce, Cs, Cd);
	if (!HC.isValid()) {
		opserr << "nDMaterial ASDConcrete3D Error: Compressive hardening law is not valid.\n";
		return nullptr;
	}

	// create the material
	NDMaterial* instance = new ASDConcrete3DMaterial(
		tag, 
		E, v, rho, eta, 
		implex, implex_control, implex_error_tolerance, implex_time_redution_limit,
		tangent, auto_regularization,
		HT, HC,
		nct, ncc, smoothing_angle);
	if (instance == nullptr) {
		opserr << "nDMaterial ASDConcrete3D Error: failed to allocate a new material.\n";
		return nullptr;
	}
	return instance;
}

int ASDConcrete3DMaterial::StressDecomposition::compute(const Vector& S)
{
	// copy S in V
	V(0, 0) = S(0);
	V(1, 1) = S(1);
	V(2, 2) = S(2);
	V(0, 1) = V(1, 0) = S(3);
	V(1, 2) = V(2, 1) = S(4);
	V(0, 2) = V(2, 0) = S(5);

	// solve
	if (Eigen3(Si, V) < 0)
		return EC_Eigen_Error;

	// construct matrices PT and PC
	static Matrix pjj(6, 6);
	PT.Zero();
	PC.Zero();

	// method 1: compute PT = sum(H(sj)*pjj) --- PC = I - PT
	//for (int j = 0; j < 3; j++) {
	//	double hsj = Heavyside(Si(j));
	//	if (hsj > 0.0) {
	//		computePjj(V, j, pjj);
	//		PT.addMatrix(1.0, pjj, hsj);
	//	}
	//}
	//for (int i = 0; i < 6; ++i)
	//	PC(i, i) = 1.0;
	//PC.addMatrix(1.0, PT, -1.0);

	// method 2: compute PT = sum(H(sj)*pjj) --- PC = sum(H(-sj)*pjj)
	// then PO = I - PT - PC
	// finally PT += PO/2 --- PC += PO/2
	for (int j = 0; j < 3; j++) {
		computePjj(V, j, pjj);
		double hsj = Heavyside(Si(j));
		if(hsj > 0.0)
			PT.addMatrix(1.0, pjj, hsj);
		hsj = Heavyside(-Si(j));
		if (hsj > 0.0)
			PC.addMatrix(1.0, pjj, hsj);
	}
	static Matrix PO(6, 6); // PO = I - PT - PC
	PO.addMatrix(0.0, PT, -1.0);
	PO.addMatrix(1.0, PC, -1.0);
	for (int i = 0; i < 6; ++i)
		PO(i, i) += 1.0;
	PT.addMatrix(1.0, PO, 0.5);
	PC.addMatrix(1.0, PO, 0.5);

	// compute positive and negative parts of the stress
	ST.addMatrixVector(0.0, PT, S, 1.0);
	SC.addMatrixVector(0.0, PC, S, 1.0);

	// done
	return 0;
}

ASDConcrete3DMaterial::HardeningLaw::HardeningLaw(
	int tag, HardeningLawType type,
	double E,
	const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& d)
	: m_tag(tag)
	, m_type(type)
{
	// initial checks
	if (!(x.size() > 0 && x.size() == y.size() && x.size() == d.size())) {
		opserr << "ASDConcrete3D Fatal Error: HardeningLaw::c-tor - found incompatible sizes.\n";
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
		opserr << "ASDConcrete3D Fatal Error: HardeningLaw::c-tor - max(X) == 0 " << xmax << "\n";
		return;
	}
	if (ymax == 0.0) {
		opserr << "ASDConcrete3D Fatal Error: HardeningLaw::c-tor - max(Y) == 0 " << ymax << "\n";
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
	// adjust
	adjust();
	// compute fracture energy
	computeFractureEnergy();
	// store
	HardeningLawStorage::instance().store(*this);
}

void ASDConcrete3DMaterial::HardeningLaw::regularize(double lch)
{
	// back to original
	deRegularize();
	// quick return
	if (!m_valid)
		return;
	// check lch (< 0 -> invalid, 1 -> no-op)
	if (lch <= 0.0 || lch == 1.0)
		return;
	// return if un-bounded (inf fracture energy)
	if (!m_fracture_energy_is_bounded)
		return;
	// the initial fracture energy has been computed in the full constructor.
	// so this method should be called only once!
	// compute the required specific fracture energy
	double gnew = m_fracture_energy / lch;
	// compute the minimum fracture energy (in case lch is too large)
	const auto& peak = m_points[m_softening_begin];
	double gmin = (peak.y * peak.x / 2.0) * 1.01; // make it 1% larger to ensure a monotonically increasing abscissa
	gnew = std::max(gnew, gmin);
	// iteratively scale the curve until the gnew
	double tol = 1.0e-3 * gnew;
	double dscale = gnew / m_fracture_energy;
	double x0 = peak.x;
	double scale = dscale;
	static constexpr int max_iter = 10;
	for (int iter = 0; iter < max_iter; ++iter) {
		// scale points after the peak
		for (std::size_t i = m_softening_begin + 1; i < m_points.size(); ++i) {
			auto& pi = m_points[i];
			pi.x = x0 + (pi.x - x0) * dscale;
		}
		// update g and check
		computeFractureEnergy();
		if (std::abs(m_fracture_energy - gnew) < tol)
			break;
		// update scale
		dscale = gnew / m_fracture_energy;
		scale *= dscale;
	}
	// re-adjust
	adjust();
}

void ASDConcrete3DMaterial::HardeningLaw::deRegularize()
{
	auto source = HardeningLawStorage::instance().recover(m_tag, m_type);
	if (source)
		*this = *source;
}

ASDConcrete3DMaterial::HardeningLawPoint ASDConcrete3DMaterial::HardeningLaw::evaluateAt(double x) const
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
		// interpolate last tangent for y (it can go negative..., it will be set to a minimum later on)
		y1 = m_points.back().y;
		double tangent = (y1 - m_points[m_points.size() - 2].y) / (x1 - m_points[m_points.size() - 2].x);
		y2 = y1 + span * tangent;
		// interpolate last tangent if positive, otherwise keep constant
		q1 = m_points.back().q;
		tangent = (q1 - m_points[m_points.size() - 2].q) / (x1 - m_points[m_points.size() - 2].x);
		q2 = tangent > 0.0 ? q1 + span * tangent : q1;
	}
	// interpolate
	double  xspan = x2 - x1;
	double xratio = xspan > 0.0 ? (x - x1) / xspan : 0.0;
	double y = std::max(m_ytolerance, y1 + (y2 - y1) * xratio);
	double q = std::max(m_ytolerance, q1 + (q2 - q1) * xratio);
	double d = 1.0 - y / q;
	// done
	return HardeningLawPoint(x, y, d, q);
}

void ASDConcrete3DMaterial::HardeningLaw::adjust()
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

void ASDConcrete3DMaterial::HardeningLaw::computeFractureEnergy()
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

ASDConcrete3DMaterial::HardeningLawStorage& ASDConcrete3DMaterial::HardeningLawStorage::instance()
{
	static HardeningLawStorage _instance;
	return _instance;
}

void ASDConcrete3DMaterial::HardeningLawStorage::store(const HardeningLaw& hl)
{
	if (hl.type() == HardeningLawType::Tension) {
		m_tension[hl.tag()] = std::make_shared<HardeningLaw>(hl);
	}
	else {
		m_compression[hl.tag()] = std::make_shared<HardeningLaw>(hl);
	}
}

ASDConcrete3DMaterial::HardeningLawStorage::PointerType ASDConcrete3DMaterial::HardeningLawStorage::recover(int tag, HardeningLawType type)
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

ASDConcrete3DMaterial::CrackPlanesStorage& ASDConcrete3DMaterial::CrackPlanesStorage::instance()
{
	static CrackPlanesStorage _instance;
	return _instance;
}

ASDConcrete3DMaterial::CrackPlanesStorage::Vector3ListPointer ASDConcrete3DMaterial::CrackPlanesStorage::get(int n90)
{
	// return nullptr for 0
	if (n90 < 1)
		return nullptr;
	// get it or make it
	auto& item = m_map[n90];
	if (item == nullptr) {
		// make it
		item = std::make_shared<Vector3List>();
		Vector3List& normals = *item;
		double dangle = M_PI / 2.0 / static_cast<double>(n90); // angle increment
		std::size_t N = static_cast<std::size_t>(n90) * 2; // number of subdivision per linear range [0 - 180[
		std::size_t M = N * (N - 1) + 1; // number of total unique samples
		normals.resize(M);
		std::size_t counter = 0;
		for (std::size_t j = 0; j < N; ++j) {
			double ay = static_cast<double>(j)* dangle;
			std::size_t Ni = j == 0 ? 1 : N;
			for (std::size_t i = 0; i < Ni; ++i) {
				double ax = static_cast<double>(i)* dangle;
				Vector3& nij = normals[counter++];
				nij.x = std::cos(ax) * std::sin(ay);
				nij.y = std::sin(ax) * std::sin(ay);
				nij.z = std::cos(ay);
			}
		}
	}
	// done
	return item;
}

ASDConcrete3DMaterial::CrackPlanes::CrackPlanes(int n90)
	: m_n90(std::max(0, n90))
	, m_normals(CrackPlanesStorage::instance().get(m_n90))
{
	std::size_t num = m_normals ? m_normals->size() : 1;
	m_equivalent_strain.resize(num, 0.0);
}

void ASDConcrete3DMaterial::CrackPlanes::setCurrentNormal(const Vector3& ni)
{
	// always reset to 0
	m_closest_normal_loc = 0;
	// search only if necessary
	if (m_normals) {
		m_current_normal = ni;
		double max_dot = -1.0;
		const auto& normals = *m_normals;
		for (std::size_t i = 0; i < normals.size(); ++i) {
			double adot = std::abs(m_current_normal.dot(normals[i]));
			if (adot > max_dot) {
				max_dot = adot;
				m_closest_normal_loc = i;
			}
		}
	}
}

double ASDConcrete3DMaterial::CrackPlanes::getCurrentEquivalentStrain() const
{
	if (m_closest_normal_loc < m_equivalent_strain.size())
		return m_equivalent_strain[m_closest_normal_loc];
	return 0.0;
}

void ASDConcrete3DMaterial::CrackPlanes::updateCurrentEquivalentStrain(double x, double smooth_angle)
{
	if (m_closest_normal_loc < m_equivalent_strain.size()) {
		// smooth only if necessary
		if (m_normals) {
			double sig = std::max(1.0e-6, smooth_angle) / 2.3546;
			double sden = 2.0 * sig * sig;
			const auto& normals = *m_normals;
			for (std::size_t i = 0; i < normals.size(); ++i) {
				double adot = std::abs(m_current_normal.dot(normals[i]));
				double angle = std::acos(adot);
				double X1 = std::exp(-std::pow(angle, 2) / sden);
				double X2 = std::exp(-std::pow(angle + M_PI, 2) / sden);
				double X3 = std::exp(-std::pow(angle - M_PI, 2) / sden);
				double XX = std::max(X1, std::max(X2, X3));
				double EEQ = m_equivalent_strain[i];
				EEQ = std::max(EEQ, x * XX);
				m_equivalent_strain[i] = EEQ;
			}
		}
		// impose the closest always
		double EEQ_closest = m_equivalent_strain[m_closest_normal_loc];
		EEQ_closest = std::max(EEQ_closest, x);
		m_equivalent_strain[m_closest_normal_loc] = EEQ_closest;
	}
}

void ASDConcrete3DMaterial::CrackPlanes::reset()
{
	for (std::size_t i = 0; i < m_equivalent_strain.size(); ++i)
		m_equivalent_strain[i] = 0.0;
}

ASDConcrete3DMaterial::ASDConcrete3DMaterial(
	int _tag,
	double _E,
	double _v,
	double _rho,
	double _eta,
	bool _implex,
	bool _implex_control,
	double _implex_error_tolerance,
	double _implex_time_reduction_limit,
	bool _tangent,
	bool _auto_regularize,
	const HardeningLaw& _ht,
	const HardeningLaw& _hc,
	int _nct,
	int _ncc,
	double _smoothing_angle)
	: NDMaterial(_tag, ND_TAG_ASDConcrete3DMaterial)
	, E(_E)
	, v(_v)
	, rho(_rho)
	, eta(_eta)
	, implex(_implex)
	, implex_control(_implex_control)
	, implex_error_tolerance(_implex_error_tolerance)
	, implex_time_redution_limit(_implex_time_reduction_limit)
	, tangent(_tangent)
	, auto_regularize(_auto_regularize)
	, ht(_ht)
	, hc(_hc)
	, nct(std::max(0, _nct))
	, ncc(std::max(0, _ncc))
	, smoothing_angle(std::abs(_smoothing_angle) * M_PI / 180.0)
	, svt(nct)
	, svt_commit(nct)
	, svt_commit_old(nct)
	, svc(ncc)
	, svc_commit(ncc)
	, svc_commit_old(ncc)
{
	// Initialize the committed PT (for implex) to Identity/2
	PT_commit.Zero();
	for (int i = 0; i < 6; ++i)
		PT_commit(i, i) = 0.5;
}

ASDConcrete3DMaterial::ASDConcrete3DMaterial()
	: NDMaterial(0, ND_TAG_ASDConcrete3DMaterial)
{
}

ASDConcrete3DMaterial::~ASDConcrete3DMaterial()
{ 
}

double ASDConcrete3DMaterial::getRho(void)
{
	return rho;
}

int ASDConcrete3DMaterial::setTrialStrain(const Vector& v)
{
	// return value
	int retval = 0;

	// get charcteristic length and perform regularization
	if (auto_regularize && !regularization_done) {
		if (ops_TheActiveElement) {
			double lch = ops_TheActiveElement->getCharacteristicLength();
			ht.regularize(lch);
			hc.regularize(lch);
			regularization_done = true;
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
		static Matrix Cnum(6, 6);
		// strain perturbation parameter
		double PERT = (ht.strainTolerance() + hc.strainTolerance()) / 2.0;
		// compute the forward perturbed solution and store in Cnum columns
		for (int j = 0; j < 6; ++j) {
			strain = v;
			strain(j) += PERT;
			retval = compute(true, false);
			if (retval < 0) return retval;
			for (int i = 0; i < 6; ++i)
				Cnum(i, j) = stress(i);
		}
		// compute unperturbed solution
		strain = v;
		retval = compute(true, false);
		if (retval < 0) return retval;
		for (int j = 0; j < 6; ++j) {
			for (int i = 0; i < 6; ++i)
				Cnum(i, j) = (Cnum(i, j) - stress(i)) / PERT;
		}
		// save tangent
		C = Cnum;
	}
	else {
		strain = v;
		if (implex) {
			if (implex_control) {
				// implicit solution 
				static Matrix aux(6, 6);
				aux = PT_commit;
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
					if (dtime_n >= implex_time_redution_limit * dtime_0)
						retval = EC_IMPLEX_Error_Control;
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

int ASDConcrete3DMaterial::setTrialStrain(const Vector& v, const Vector& /*r*/)
{
	return setTrialStrain(v);
}

int ASDConcrete3DMaterial::setTrialStrainIncr(const Vector& v)
{
	static Vector aux(6);
	aux = strain;
	aux.addVector(1.0, v, 1.0);
	return setTrialStrain(aux);
}

int ASDConcrete3DMaterial::setTrialStrainIncr(const Vector& v, const Vector& /*r*/)
{
	return setTrialStrainIncr(v);
}

const Vector &ASDConcrete3DMaterial::getStrain(void)
{
	return strain;
}

const Vector &ASDConcrete3DMaterial::getStress(void)
{
	return stress;
}

const Matrix &ASDConcrete3DMaterial::getTangent(void)
{
	return C;
}

const Matrix &ASDConcrete3DMaterial::getInitialTangent(void)
{
	static Matrix D(6, 6);
	D.Zero();
	double mu2 = E / (1.0 + v);
	double lam = v * mu2 / (1.0 - 2.0 * v);
	double mu = 0.50 * mu2;
	mu2 += lam;
	D(0, 0) = D(1, 1) = D(2, 2) = mu2;
	D(0, 1) = D(1, 0) = D(0, 2) = D(2, 0) = D(1, 2) = D(2, 1) = lam;
	D(3, 3) = D(4, 4) = D(5, 5) = mu;
	return D;
}

int ASDConcrete3DMaterial::commitState(void)
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
	}
	// store the previously committed variables for next move from n to n - 1
	svt_commit_old = svt_commit;
	svc_commit_old = svc_commit;
	// store committed variables
	svt_commit = svt;
	svc_commit = svc;
	strain_commit = strain;
	stress_eff_commit = stress_eff;
	dtime_n_commit = dtime_n;
	// done
	commit_done = true;
	return 0;
}

int ASDConcrete3DMaterial::revertToLastCommit(void)
{
	// restore converged values
	svt = svt_commit;
	svc = svc_commit;
	strain = strain_commit;
	stress_eff = stress_eff_commit;
	dtime_n = dtime_n_commit;
	// done
	return 0;
}

int ASDConcrete3DMaterial::revertToStart(void)
{
	// State variables
	svt.reset();
	svt_commit.reset();
	svt_commit_old.reset();
	svc.reset();
	svc_commit.reset();
	svc_commit_old.reset();

	// Time step
	dtime_n = 0.0;
	dtime_n_commit = 0.0;
	dtime_0 = 0.0;
	dtime_is_user_defined = false;

	// Commit flag
	commit_done = false;

	// Initialize the committed PT (for implex) to Identity/2
	PT_commit.Zero();
	for (int i = 0; i < 6; ++i)
		PT_commit(i, i) = 0.5;

	// IMPL-EX error
	implex_error = 0.0;

	// Strain, Stress and Tangent
	strain.Zero();
	strain_commit.Zero();
	stress.Zero();
	stress_eff.Zero();
	stress_eff_commit.Zero();
	C.Zero();

	// Output variables
	dt_bar = 0.0;
	dc_bar = 0.0;

	// Done
	return 0;
}

NDMaterial * ASDConcrete3DMaterial::getCopy(void)
{
	// we can safely use the default copy-constructor according to the member variables we're using
	return new ASDConcrete3DMaterial(*this);
}

NDMaterial* ASDConcrete3DMaterial::getCopy(const char* code)
{
	if (strcmp(code, "ThreeDimensional") == 0)
		return getCopy();
	return NDMaterial::getCopy(code);
}

const char* ASDConcrete3DMaterial::getType(void) const
{
	return "ThreeDimensional";
}

int ASDConcrete3DMaterial::getOrder(void) const
{
	return 6;
}

void ASDConcrete3DMaterial::Print(OPS_Stream &s, int flag)
{
	s << "ASDConcrete3D Material, tag: " << this->getTag() << "\n";
}

int ASDConcrete3DMaterial::sendSelf(int commitTag, Channel &theChannel)
{
	return -1; // todo
}

int ASDConcrete3DMaterial::recvSelf(int commitTag, Channel & theChannel, FEM_ObjectBroker & theBroker)
{
	return -1; // todo
}

int ASDConcrete3DMaterial::setParameter(const char** argv, int argc, Parameter& param)
{
	return -1; //todo
}

int ASDConcrete3DMaterial::updateParameter(int parameterID, Information& info)
{
	return -1; //todo
}

Response* ASDConcrete3DMaterial::setResponse(const char** argv, int argc, OPS_Stream& output)
{
	// utils
	auto lam_begin = [&output, this]() {
		output.tag("NdMaterialOutput");
		output.attr("matType", getClassType());
		output.attr("matTag", getTag());
	};
	auto lam_end = [&output]() {
		output.endTag(); // NdMaterialOutput
	};

	// check specific responses
	if (argc > 0) {
		//todo
	}

	// iotherwise return base-class response
	return NDMaterial::setResponse(argv, argc, output);
}

int ASDConcrete3DMaterial::getResponse(int responseID, Information& matInformation)
{
	return -1; // todo
}

int ASDConcrete3DMaterial::compute(bool do_implex, bool do_tangent)
{
	// get committed variables
	svt = svt_commit;
	svc = svc_commit;
	stress_eff = stress_eff_commit;

	// time factor for explicit extrapolation
	double time_factor = 1.0;
	if (implex && do_implex && (dtime_n_commit > 0.0))
		time_factor = dtime_n / dtime_n_commit;

	// compute rate coefficients
	double rate_coeff_1 = 0.0;
	double rate_coeff_2 = 1.0;
	if ((dtime_n > 0.0) && (eta > 0.0)) {
		rate_coeff_1 = eta / (eta + dtime_n);
		rate_coeff_2 = dtime_n / (eta + dtime_n);
	}

	// compute elastic effective stress: SEFFn = C0 : (En - En-1)
	static Vector dStrain(6);
	dStrain = strain;
	dStrain.addVector(1.0, strain_commit, -1.0);
	stress_eff.addMatrixVector(1.0, getInitialTangent(), dStrain, 1.0);

	// compute stress split
	static StressDecomposition D;
	if (implex && do_implex) {
		// take eigenvectors from known (n-1) effective stress.
		// we need them for the crack planes' normals
		D.compute(stress_eff_commit);
		// update the ST and SC with the committed PT and PC, but with current effective stress
		// don't update eigenvalues, we don't need them
		D.ST.addMatrixVector(0.0, D.PT, stress_eff, 1.0);
		D.SC.addMatrixVector(0.0, D.PC, stress_eff, 1.0);
	}
	else {
		D.compute(stress_eff);
	}

	// update normals
	Vector3 Tnormal(D.V(0, 0), D.V(1, 0), D.V(2, 0));
	Vector3 Cnormal(D.V(0, 2), D.V(1, 2), D.V(2, 2));
	svt.setCurrentNormal(Tnormal);
	svc.setCurrentNormal(Cnormal);

	// compute committed hardening variables
	HardeningLawPoint pt = ht.evaluateAt(svt.getCurrentEquivalentStrain());
	HardeningLawPoint pc = hc.evaluateAt(svc.getCurrentEquivalentStrain());

	// tempoary clone of old equivalent plastic strains
	double xt_pl = pt.plasticStrain(E);
	double xc_pl = pc.plasticStrain(E);

	// compute new trial equivalent strain measures
	double xt_trial, xc_trial;
	if (implex && do_implex) {
		// extrapolated equivalent strain measures (explicit)
		// already in equivalent total stress space
		xt_trial = svt_commit.getCurrentEquivalentStrain() + 
			time_factor * (svt_commit.getCurrentEquivalentStrain() - svt_commit_old.getCurrentEquivalentStrain());
		xc_trial = svc_commit.getCurrentEquivalentStrain() +
			time_factor * (svc_commit.getCurrentEquivalentStrain() - svc_commit_old.getCurrentEquivalentStrain());
	}
	else {
		xt_trial = tensileCriterion(D) + xt_pl;
		xc_trial = compressiveCriterion(D) + xc_pl;
	}

	// update hardening variables
	if (xt_trial > pt.x) {
		double xt_new = rate_coeff_1 * pt.x + rate_coeff_2 * xt_trial;
		pt = ht.evaluateAt(xt_new);
		svt.updateCurrentEquivalentStrain(pt.x, smoothing_angle);
	}
	if (xc_trial > pc.x) {
		double xc_new = rate_coeff_1 * pc.x + rate_coeff_2 * xc_trial;
		pc = hc.evaluateAt(xc_new);
		svc.updateCurrentEquivalentStrain(pc.x, smoothing_angle);
	}

	// compute plastic damage
	double seff_eq_t = (pt.x - xt_pl) * E;
	double dt_plastic = seff_eq_t > 0.0 ? 1.0 - pt.q / seff_eq_t : 0.0;
	double seff_eq_c = (pc.x - xc_pl) * E;
	double dc_plastic = seff_eq_c > 0.0 ? 1.0 - pc.q / seff_eq_c : 0.0;

	// update effective stress
	stress_eff.addVector(0.0, D.ST, 1.0 - dt_plastic);
	stress_eff.addVector(1.0, D.SC, 1.0 - dc_plastic);

	// update nominal stress
	dt_bar = pt.d + dt_plastic - pt.d * dt_plastic;
	dc_bar = pc.d + dc_plastic - pc.d * dc_plastic;
	stress.addVector(0.0, D.ST, 1.0 - dt_bar);
	stress.addVector(1.0, D.SC, 1.0 - dc_bar);

	// tangent matrix
	if (do_tangent) {
		static Matrix W(6, 6);
		W.Zero();
		for (int i = 0; i < 6; ++i)
			W(i, i) = 1.0;
		W.addMatrix(1.0, D.PT, -dt_bar);
		W.addMatrix(1.0, D.PC, -dc_bar);
		C.addMatrixProduct(0.0, W, getInitialTangent(), 1.0);
	}

	// done
	return 0;
}

double ASDConcrete3DMaterial::lublinerCriterion(const StressDecomposition& D, double ft, double fc, double k1, double scale)
{
	double fb = 1.16 * fc;
	double Kc = 2.0 / 3.0;
	double gamma = 3.0 * (1.0 - Kc) / (2.0 * Kc - 1.0);
	double alpha = (fb - fc) / (2.0 * fb - fc);
	double I1 = D.Si(0) + D.Si(1) + D.Si(2);
	static Vector Sdev(3);
	for (int i = 0; i < 3; i++)
		Sdev(i) = D.Si(i) - I1 / 3.0;
	double J2 = 0.5 * (Sdev ^ Sdev);
	double beta = fc / ft * (1.0 - alpha) - (1.0 + alpha);
	double smax_alg = D.Si(0);
	double smax = Macauley(smax_alg);
	double smin = Macauley(-smax_alg);
	return (1.0 / (1.0 - alpha) * (alpha * I1 + std::sqrt(3.0 * J2) + k1 * beta * smax - gamma * smin)) * scale / E;
}

double ASDConcrete3DMaterial::tensileCriterion(const StressDecomposition& D)
{
	// skip if maximum principal stress is not strictly positive
	if (D.Si(0) < ht.stressTolerance())
		return 0.0;

	// compute equivalent strain measure
	// - Lubliner criterion
	// - normalized to fc
	// - assume ratio ft/fc = 0.1 (don't take the ratio from the hardening laws!)
	return lublinerCriterion(D, 0.1, 1.0, 1.0, 0.1);
}

double ASDConcrete3DMaterial::compressiveCriterion(const StressDecomposition& D)
{
	// skip if minimum principal stress is not strictly negative
	if (D.Si(2) > -hc.stressTolerance())
		return 0.0;

	// compute equivalent strain measure
	// - Lubliner criterion
	// - normalized to fc
	// - assume ratio ft/fc = 0.1 (don't take the ratio from the hardening laws!)
	return lublinerCriterion(D, 0.1, 1.0, 0.0, 1.0);
}

