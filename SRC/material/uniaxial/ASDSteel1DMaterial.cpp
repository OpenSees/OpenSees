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
// $Date: 2025-01-03 11:29:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ASDSteel1DMaterial.cpp,v $

// Alessia Casalucci - ASDEA Software, Italy
//
// todo...
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
#include <array>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define ASD_STEEL1D__VERBOSE

#ifdef ASD_STEEL1D__VERBOSE
#define asd_print_full(X) opserr  << __func__ << " (Line " << __LINE__ << "): " << X << "\n"
#define asd_print(X) opserr << X << "\n"
#else
#define asd_print_full(X)
#define asd_print(X)
#endif // ASD_STEEL1D__VERBOSE

// anonymous namespace for utilities
namespace {

	inline double sign(double x) { return x == 0.0 ? 0.0 : (x > 0.0 ? 1.0 : -1.0); }

	/**
	A simple 2D vector
	*/
	class V2D
	{
	public:
		double x = 0.0;
		double y = 0.0;
	public:
		V2D() = default;
		V2D(double _x, double _y) : x(_x), y(_y) {}
		inline double normalize() {
			double n = x * x + y * y;
			if (n > 0.0) {
				n = std::sqrt(n);
				x /= n;
				y /= n;
			}
			return n;
		}
		inline double operator()(int i) const {
			if (i == 0) return x;
			return y;
		}
	};
	inline V2D operator + (const V2D& a, const V2D& b) {
		return V2D(a.x + b.x, a.y + b.y);
	}
	inline V2D operator - (const V2D& a, const V2D& b) {
		return V2D(a.x - b.x, a.y - b.y);
	}
	inline V2D operator * (double a, const V2D& b) {
		return V2D(a * b.x, a * b.y);
	}
	inline V2D operator * (const V2D& b, double a) {
		return a * b;		
	}

	/**
	A simple 2D basis matrix
	*/
	class M2D
	{
	public:
		V2D dx;
		V2D dy;
	public:
		M2D() = default;
		M2D(const V2D& x, const V2D& y) : dx(x), dy(y) {}
		inline double operator()(int i, int j) const {
			if (i == 0) return dx(j);
			return dy(j);
		}
	};
	
	/**
	A quaternion optimized for the 2D case
	*/
	class Q2D
	{
	public:
		double w = 0.0;
		double z = 0.0;
	public:
		Q2D() = default;
		Q2D(double _w, double _z) : w(_w), z(_z) {}
	public:
		inline Q2D conjugate()const { return Q2D(w, -z); }
		inline double toRotationVector()const {
			double ww = w;
			double zz = z;
			if (w < 0.0) {
				ww = -ww;
				zz = -zz;
			}
			double vNorm = std::sqrt(zz * zz);
			if (vNorm == 0.0) return 0.0;
			double mult = vNorm < ww ? (2.0 / vNorm * std::asin(vNorm)) : (2.0 / vNorm * std::acos(ww));
			return zz * mult;
		}
		inline V2D rotateVector(const V2D& a) const {
			V2D b(-2.0 * z * a.y, 2.0 * z * a.x);
			V2D c(-z * b.y, z * b.x);
			return a + b * w + c;
		}
		inline void normalize() {
			double n = w * w + z * z;
			if (n > 0.0) {
				n = std::sqrt(n);
				w /= n;
				z /= n;
			}
		}
		static Q2D identity() {
			return Q2D(1.0, 0.0);
		}
		static Q2D fromRotationVector(double axis) {
			double n = std::abs(axis);
			if (n == 0.0)
				return Q2D::identity();
			if (n != 1.0)
				axis /= n;
			double halfAngle = n * 0.5;
			double q0 = std::cos(halfAngle);
			double s = std::sin(halfAngle);
			Q2D result(q0, axis * s);
			result.normalize();
			return result;
		}
		static Q2D fromRotationMatrix(const M2D& m) {
			double tr = m.dx.x + m.dy.y + 1.0;
			if ((tr > m.dx.x) && (tr > m.dy.y) && (tr > 1.0)) {
				double S = std::sqrt(tr + 1.0) * 2.0;
				Q2D Q(0.25 * S, (m.dy.x - m.dx.y) / S);
				Q.normalize();
				return Q;
			}
			else {
				double S = std::sqrt(1.0 + 1.0 - m.dx.x - m.dy.y) * 2.0;
				Q2D Q((m.dy.x - m.dx.y) / S, 0.25 * S);
				Q.normalize();
				return Q;
			}
		}
	};
	inline Q2D operator * (const Q2D& a, const Q2D& b) {
		return Q2D(a.w * b.w - a.z * b.z, a.w * b.z + a.z * b.w);
	}

	/**
	Global storage to avoid many dynamic allocation/deallocations.
	*/
	class Globals
	{
	private:
		Globals() = default;
		Globals(const Globals&) = delete;
		Globals& operator = (const Globals&) = delete;

	public:
		static Globals& instance() {
			static Globals _instance;
			return _instance;
		}

	public:

		// for section
		Vector section_stress = Vector(3);
		Vector section_strain = Vector(3);
		Matrix section_tangent = Matrix(3, 3);

		// for element
		Vector element_RHS = Vector(6);
		Vector element_RHS_local = Vector(6);
		Matrix element_LHS = Matrix(6, 6);
		Matrix element_LHS_local = Matrix(6, 6);
		Matrix element_T = Matrix(6, 6);
		Vector element_U = Vector(6);
		Vector element_U_commit = Vector(6);
		Vector element_U_local = Vector(6);
		Vector element_U_local_commit = Vector(6);
		Matrix element_B = Matrix(3, 6);
		Matrix element_TtK = Matrix(6, 6);
		Matrix element_BtC = Matrix(6, 3);
		Vector element_dU = Vector(6);

		// for RVE
		Vector rve_dU = Vector(7);
		Vector rve_R = Vector(7);
		Matrix rve_K = Matrix(7, 7);
		Vector rve_Kcr = Vector(7);
		Vector rve_Krc = Vector(7);
		Vector rve_Kinv_Kcr = Vector(7);
		Vector rve_dU_el = Vector(10);
		Vector rve_R_el = Vector(10);
		Matrix rve_K_el = Matrix(10, 10);
		Vector rve_Kcr_el = Vector(10);
		Vector rve_Krc_el = Vector(10);
		Vector rve_Kinv_Kcr_el = Vector(10);
		double rve_Krr = 0.0;
		double rve_Kred = 0.0;
		// node positions (to be set at each setTrialStrain)
		std::array<V2D, 4> rve_nodes = { V2D(0.0, 0.0), V2D(0.0, 0.0), V2D(0.0, 0.0), V2D(0.0, 0.0) };
		inline void setRVENodes(double lch) {
			double dy = lch / 3.0;
			for (int i = 0; i < 4; ++i) {
				auto& ip = rve_nodes[i];
				ip.y = dy * static_cast<double>(i); // from bottom to top, bottom at 0,0
				ip.x = ip.y * 1.0e-3; // add imperfection
			}
		}
		std::array<V2D, 5> rve_nodes_el = { V2D(0.0, 0.0), V2D(0.0, 0.0), V2D(0.0, 0.0), V2D(0.0, 0.0), V2D(0.0, 0.0) };
		inline void setRVENodes_el(double lch, double length_el) {
			double dy = lch / 3.0;
			double y_el = (length_el - lch);
			for (int i = 1; i < 5; ++i) {
				auto& ip = rve_nodes_el[i];
				ip.y = dy * static_cast<double>(i-1) + y_el; // from bottom to top, bottom at 0,0
				ip.x = ip.y * 1.0e-3; // add imperfection
			}
			rve_nodes_el[0].y = 0.0;
			rve_nodes_el[0].x = 0.0;
		}
	};

	/**
	Basic steel component.
	*/
	class SteelComponent
	{
	public:
		using param_t = ASDSteel1DMaterial::InputParameters;
		SteelComponent() = default;

		inline int commitState() {
			// store the previously committed variables for next move from n to n - 1
			lambda_commit_old = lambda_commit;
			// state variables
			alpha1_commit = alpha1;
			alpha2_commit = alpha2;
			lambda_commit = lambda;
			strain_commit = strain;
			stress_commit = stress;
			// done 
			return 0;
		}
		inline void revertToLastCommit() {
			// state variables
			alpha1 = alpha1_commit;
			alpha2 = alpha2_commit;
			lambda = lambda_commit;
			strain = strain_commit;
			stress = stress_commit;
		}
		inline void revertToStart() {
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
		}
		inline int compute(const param_t& params, bool do_implex, double time_factor, double _strain, double& sigma, double& tangent) {
			int retval = 0;
			// settings
			constexpr int MAX_ITER = 1000;
			constexpr double F_REL_TOL = 1.0e-6;
			constexpr double L_ABS_TOL = 1.0e-8;
			//constexpr double K_alpha = 0.8;  
			// base steel response
			alpha1 = alpha1_commit;
			alpha2 = alpha2_commit;
			lambda = lambda_commit;
			// elastic predictor
			strain = _strain;
			double dstrain = strain - strain_commit;
			sigma = stress_commit + params.E * dstrain;
			tangent = params.E;
			// plastic utilities
			double sg = 0.0; // plastic flow direction
			auto lam_rel_stress = [this, &sigma]() -> double {
				return sigma - alpha1 - alpha2;
				};
			auto lam_yield_function = [&params, &lam_rel_stress]() -> double {
				return std::abs(lam_rel_stress()) - params.sy;
				};
			auto lam_yield_derivative = [this, &lam_rel_stress, &params](double dlambda) -> double {
				// plastic flow direction
				double sg = sign(lam_rel_stress());
				// d stress / d lambda
				double dsigma = -params.E * sg;
				// d backstress / d lambda
				double dalpha1 = params.H1 * sg - params.gamma1 * alpha1;
				double dalpha2 = params.H2 * sg - params.gamma2 * alpha2;
				return sg * (dsigma - dalpha1 - dalpha2);
				};
			auto lam_yield_update = [this, &sigma, &lam_rel_stress, &params](double dlambda, double delta_lambda) {
				// plastic flow direction
				double sg = sign(lam_rel_stress());
				// update stress
				sigma -= sg * dlambda * params.E;
				// update backstress
				alpha1 = sg * params.H1 / params.gamma1 - (sg * params.H1 / params.gamma1 - alpha1_commit) * std::exp(-params.gamma1 * delta_lambda);
				alpha2 = sg * params.H2 / params.gamma2 - (sg * params.H2 / params.gamma2 - alpha2_commit) * std::exp(-params.gamma2 * delta_lambda);
				};
			// plastic corrector
			if (params.implex && do_implex) {
				// extrapolate lambda
				//  xn + time_factor * (xn - xnn);
				double delta_lambda = time_factor * (lambda_commit - lambda_commit_old);
				lambda = lambda_commit + delta_lambda;
				// extrapolate plastic flow direction
				sg = sg_commit;
				// update stress
				sigma -= sg * delta_lambda * params.E;
				// update backstress
				alpha1 = sg * params.H1 / params.gamma1 - (sg * params.H1 / params.gamma1 - alpha1_commit) * std::exp(-params.gamma1 * delta_lambda);
				alpha2 = sg * params.H2 / params.gamma2 - (sg * params.H2 / params.gamma2 - alpha2_commit) * std::exp(-params.gamma2 * delta_lambda);
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
							lambda += delta_lambda;
							// compute tangent
							sg = sign(lam_rel_stress());
							double PE =
								params.gamma1 * (params.H1 / params.gamma1 - sg * alpha1) +
								params.gamma2 * (params.H2 / params.gamma2 - sg * alpha2);
							tangent = params.K_alpha*((params.E * PE) / (params.E + PE))+(1-params.K_alpha)*params.E;  //tangent correction (K_alpha =1 -> computed tangent, K_alpha=0 -> tangent=E)
							break;
						}
					}
					if (!converged) {
						asd_print("dlambda !converged");
						retval = -1;
					}
				}
			}
			// accept solution
			stress = sigma;
			// save real plastic flow direction, if mp.implex and !do_implex -> called from commit
			if (params.implex && !do_implex) {
				// save it in implex mode during implicit phase
				sg_commit = sg;
			}
			// done
			return retval;
		}

	public:
		// state variables - backstresses
		double alpha1 = 0.0;
		double alpha1_commit = 0.0;
		double alpha2 = 0.0;
		double alpha2_commit = 0.0;
		// state variables - plastic multiplier
		double lambda = 0.0;
		double lambda_commit = 0.0;
		double lambda_commit_old = 0.0; // for implex
		double sg_commit = 0.0; // plastic flow dir for implex
		// strain, stress and tangent
		double strain = 0.0;
		double strain_commit = 0.0;
		double stress = 0.0;
		double stress_commit = 0.0;
		// methods
		static constexpr int NDATA = 12;
	};

	/**
	Section component.
	*/
	template<int NFiber>
	class SectionComponent
	{
	private:
		SectionComponent() = delete;
	};

	/**
	Section component with 1 fiber.
	*/
	template<>
	class SectionComponent<1>
	{
	public:
		SteelComponent fib;

	public:
		SectionComponent() = default;

		inline int commitState() {
			return fib.commitState();
		}
		inline void revertToLastCommit() {
			fib.revertToLastCommit();
		}
		inline void revertToStart() {
			fib.revertToStart();
		}
		inline int compute(
			const  ASDSteel1DMaterial::InputParameters& params, const Vector& strain,
			bool do_implex, double time_factor,
			Vector& stress, Matrix& tangent) {
			/**
			A section with 1 fiber is used for the central element. In this RVE model the
			central element has always zero curvature. However, with only 1 fiber the bending
			stiffness will be zero, so we use closed-form evaluation of the bending stiffness.
			Note also that M = EI*curvature for newton iterations!
			*/
			double constexpr poiss = 0.3; // built-in poisson ratio
			double constexpr kshear = 0.9; // built-in shear correction factor for circular sections
			// area
			double area = M_PI * params.radius * params.radius;
			// compute elastic shear response
			double G = params.E / (2.0 * (1.0 + poiss));
			double GAk = G * area * kshear;
			tangent(2, 2) = GAk;
			stress(2) = GAk * strain(2);
			// compute elastic M response (will be zero upon convergence for this model)
			double I = M_PI * std::pow(params.radius, 4) / 4.0;
			double EI = params.E * I;
			stress(1) = EI * strain(1);
			tangent(1, 1) = EI;
			// nonlinear P reponse
			double fiber_stress;
			double fiber_tangent;
			int retval = fib.compute(params, do_implex, time_factor, strain(0), fiber_stress, fiber_tangent);
			double fiber_force = fiber_stress * area;
			stress(0) = fiber_force;
			tangent(0, 0) = fiber_tangent * area;
			// zero terms
			tangent(0, 1) = tangent(1, 0) = tangent(0, 2) = tangent(2, 0) = tangent(1, 2) = tangent(2, 1) = 0.0;
			// done
			return retval;
		}
	};

	/**
	Section component with 3 fibers.
	*/
	template<>
	class SectionComponent<3>
	{
	public:
		std::array<SteelComponent, 3> fibers;
		static constexpr std::array<double, 3> positions = { -1.0 / 2.0, 0.0, 1.0 / 2.0 };
		static constexpr std::array<double, 3> weights = { 1.0 / 4.0, 1.0 / 2.0, 1.0 / 4.0 };

	public:
		SectionComponent() = default;
		inline int commitState() {
			int retval = 0;
			for (auto& item : fibers) {
				int fiber_retval = item.commitState();
				if (fiber_retval != 0)
					retval = fiber_retval;
			}
			return retval;
		}
		inline void revertToLastCommit() {
			for (auto& item : fibers) item.revertToLastCommit();
		}
		inline void revertToStart() {
			for (auto& item : fibers) item.revertToStart();
		}
		inline int compute(
			const  ASDSteel1DMaterial::InputParameters& params, const Vector& strain,
			bool do_implex, double time_factor,
			Vector& stress, Matrix& tangent) {
			/**
			A section with 3 fibers.
			*/
			double constexpr poiss = 0.3; // built-in poisson ratio
			double constexpr kshear = 0.9; // built-in shear correction factor for circular sections
			// area
			double area = M_PI * params.radius * params.radius;
			// zero P-M components for integration
			stress.Zero();
			tangent.Zero();
			// compute elastic shear response
			double G = params.E / (2.0 * (1.0 + poiss));
			double GAk = G * area * kshear;
			tangent(2, 2) = GAk;
			stress(2) = GAk * strain(2);
			// integrate nonlinear P-M response
			for (int i = 0; i < fibers.size(); ++i) {
				double fiber_pos = positions[i] * params.radius;
				double fiber_area = weights[i] * area;
				double fiber_strain = strain(0) - strain(1) * fiber_pos;
				double fiber_stress;
				double fiber_tangent;
				int retval = fibers[i].compute(params, do_implex, time_factor, fiber_strain, fiber_stress, fiber_tangent);
				if (retval != 0)
					return retval;
				double fiber_force = fiber_stress * fiber_area;
				stress(0) += fiber_force;
				stress(1) -= fiber_force * fiber_pos;
				
				double EtA = fiber_tangent * fiber_area;
				double ya = fiber_pos * EtA;
				tangent(0, 0) += EtA;
				tangent(0, 1) -= ya;
				tangent(1, 0) -= ya;
				tangent(1, 1) += ya * fiber_pos;
			}
			// done
			return 0;
		}
	};
	template<>
	class SectionComponent<0>
	{
	public:
		SteelComponent fib;

	public:
		SectionComponent() = default;
		inline int commitState() {
			return 0;
		}
		inline void revertToLastCommit() {
		}
		inline void revertToStart() {
		}

		inline int compute(
			const  ASDSteel1DMaterial::InputParameters& params, const Vector& strain,
			bool do_implex, double time_factor,
			Vector& stress, Matrix& tangent) {
			/**
			A section with 0 fiber is used for the elastic element.
			*/
			int retval = 0.0;
			double constexpr poiss = 0.3; // built-in poisson ratio
			double constexpr kshear = 0.9; // built-in shear correction factor for circular sections
			// area
			double area = M_PI * params.radius * params.radius;
			// compute elastic shear response
			double G = params.E / (2.0 * (1.0 + poiss));
			double GAk = G * area * kshear;
			tangent(2, 2) = GAk;
			stress(2) = GAk * strain(2);
			// compute elastic M response (will be zero upon convergence for this model)
			double I = M_PI * std::pow(params.radius, 4) / 4.0;
			double EI = params.E * I;
			stress(1) = EI * strain(1);
			tangent(1, 1) = EI;
			// compute elastic shear response
			double EA = params.E * area;
			stress(0) = EA * strain(0);
			tangent(0, 0) = EA;
			// zero terms
			tangent(0, 1) = tangent(1, 0) = tangent(0, 2) = tangent(2, 0) = tangent(1, 2) = tangent(2, 1) = 0.0;
			// done
			return retval;
		}
	};


	/**
	Element Traits:
	- 4 nodes - 3 elements RVE
	- total number of DOFs = 12
	- free DOFs = 7
	- bottom node - Ux = 0, Uy = macro_strain*L, Rz = 0
	- top node - Uy = 0, Rz = 0

	logic:
	for each dof in DOFs:
		if dof < 7
			U_dof = URVE(dof)
		else if dof == 8:
			U_dof = RY
		else
			U_dof = 0.0

	*/

	constexpr int EType_Bot = 1;
	constexpr int EType_Mid = 2;
	constexpr int EType_Top = 3;
	constexpr int EType_El = 0;
	template<int EType>
	class ETypeTraits
	{
	};
	template<>
	class ETypeTraits<EType_Bot>
	{
	public:
		static constexpr std::array<int, 6> DOFs = {8,7,9,    0,1,2 };
		static constexpr int N1 = 1;
		static constexpr int N2 = 2;
	};
	template<>
	class ETypeTraits<EType_Mid>
	{
	public:
		static constexpr std::array<int, 6> DOFs = { 0,1,2,    3,4,5 };
		static constexpr int N1 = 2;
		static constexpr int N2 = 3;
	};
	template<>
	class ETypeTraits<EType_Top>
	{
	public:
		static constexpr std::array<int, 6> DOFs = { 3,4,5,    6,10,11 };
		static constexpr int N1 = 3;
		static constexpr int N2 = 4;
	};
	template<>
	class ETypeTraits<EType_El>
	{
	public:
		static constexpr std::array<int, 6> DOFs = { 12,13,14,    8,7,9 };
		static constexpr int N1 = 0;
		static constexpr int N2 = 1;
	};
	template<int EType>
	inline void getElementDisplacementVector(bool flag, const Vector& G, Vector& L) {
		// G = global vector, size = 8 -> full = 12, semi-full 8 (7 free + 1 Uy imposed)
		// L = local vector , size = 6
		if (flag) {
			for (int i = 0; i < 6; ++i) {
				int gdof = ETypeTraits<EType>::DOFs[i];
				if (gdof < 10 ) {
					L(i) = G(gdof);

				}
				else if (gdof == 13) {
					L(i) = G(10);
				}
				else {
					L(i) = 0.0;
				}
			}
		}
		else {
			for (int i = 0; i < 6; ++i) {
				int gdof = ETypeTraits<EType>::DOFs[i];
				L(i) = gdof < 8 ? G(gdof) : 0.0;
			}
		}
	}
	template<int EType>
	inline void assembleRHS(bool flag, const Vector& L, Vector& G) {
		// G = global vector, size = 7 (-> full = 12, semi-full 8 (7 free + 1 Uy imposed) only free)
		// L = local vector , size = 6
		if (flag) {
			for (int i = 0; i < 6; ++i) {
				int gdof = ETypeTraits<EType>::DOFs[i];
				if (gdof < 10) {
					G(gdof) += L(i);
				}
			}
		}
		else {
			for (int i = 0; i < 6; ++i) {
				int gdof = ETypeTraits<EType>::DOFs[i];
				if (gdof < 7) {
					G(gdof) += L(i);
				}
			}
		}		
	}
	template<int EType>
	inline void assembleLHS(bool flag, const Matrix& M, Matrix& G) {
		// G = global matrix, size = 7x7 -> full = 12x12, semi-full 8 (7 free + 1 Uy imposed)
		// L = local  vector , size = 6
		// M = local  matrix , size = 6x6
		if (flag) {
			for (int i = 0; i < 6; ++i) {
				int idof = ETypeTraits<EType>::DOFs[i];
				if (idof < 10) {
					for (int j = 0; j < 6; ++j) {
						int jdof = ETypeTraits<EType>::DOFs[j];
						if (jdof < 10) {
							G(idof, jdof) += M(i, j);
						}
					}
				}
			}
		}
		else {
			for (int i = 0; i < 6; ++i) {
				int idof = ETypeTraits<EType>::DOFs[i];
				if (idof < 7) {
					for (int j = 0; j < 6; ++j) {
						int jdof = ETypeTraits<EType>::DOFs[j];
						if (jdof < 7) {
							G(idof, jdof) += M(i, j);
						}
					}
				}
			}
		}
		
	}
	template<int EType>
	inline void getRetainedComponents(bool flag, const Matrix& M, double& Krr, Vector& Kcr, Vector& Krc) {
		if (flag) {
			for (int i = 0; i < 6; ++i) {
				int idof = ETypeTraits<EType>::DOFs[i];
				for (int j = 0; j < 6; ++j) {
					int jdof = ETypeTraits<EType>::DOFs[j];
					if (idof == 13) {
						if (jdof == 13) {
							Krr = M(i, j);
						}
						else if (jdof < 10) { //it is used to assemble krc and kcr for the bottom element -> jdof will be always < 3
							Krc(jdof) = M(i, j);
						}
					}
					else if (jdof == 13 && idof < 10) {
						Kcr(idof) = M(j, i);
					}
				}
			}
		}
		else {
			for (int i = 0; i < 6; ++i) {
				int idof = ETypeTraits<EType>::DOFs[i];
				for (int j = 0; j < 6; ++j) {
					int jdof = ETypeTraits<EType>::DOFs[j];
					if (idof == 7) {
						if (jdof == 7) {
							Krr = M(i, j);
						}
						else if (jdof < 7) { //it is used to assemble krc and kcr for the bottom element -> jdof will be always < 3
							Krc(jdof) = M(i, j);
						}
					}
					else if (jdof == 7 && idof < 7) {
						Kcr(idof) = M(j, i);
					}
				}
			}
		}
	}

	inline M2D computeLocalFrameInternal(const V2D& X1, const V2D& X2, Q2D& Qi, V2D& C) {
		C = (X1 + X2) * 0.5;
		V2D dx = X2 - X1;
		double L = dx.normalize();
		V2D dy(-dx.y, dx.x);
		M2D Ri(dx, dy);
		Qi = Q2D::fromRotationMatrix(Ri);
		return Ri;
	}
	
	inline void computeLocalFrame(const V2D& X1, const V2D& X2, Q2D& Qi, V2D& C) {
		computeLocalFrameInternal(X1, X2, Qi, C);
	}
	inline void computeLocalFrame(const V2D& X1, const V2D& X2, Q2D& Qi, V2D& C, Matrix& T) {
		M2D Ri = computeLocalFrameInternal(X1, X2, Qi, C);
		//T matrix : considering P_t @ R
		for (int i = 0; i < 2;++i) {
			for (int j = 0; j < 2;++j) {
				T(i, j) = T(i + 3, j + 3) = Ri(i, j) / 2.0;
				T(i, j + 3) = T(i + 3, j) = -Ri(i, j) / 2.0;
			}
		}
		T(2,2) = T(5,5) = 1.0;
	}

	inline void computeBmatrix(const V2D& X1, const V2D& X2, double& L,Matrix& B) {
		V2D dx = X2 - X1;
		L = dx.normalize();
		B(0,0) = B(1,2) = B(2,1) = -1.0 / L;
		B(0,3) = B(1,5) = B(2,4) = 1.0 / L;
		B(2,2) = B(2,5) = -1.0 / 2.0;
	}

	class RVEStateVariables
	{
	public:
		Vector UG = Vector(8);
		Vector UG_commit = Vector(8);
		Vector UG_el = Vector(11);
		Vector UG_el_commit = Vector(11);
	};

	/**
	Element component.
	*/
	template<int NFiber, int EType>
	class ElementComponent
	{
	public:
		SectionComponent<NFiber> section;
		Vector UL_commit = Vector(6);
		std::array<Q2D, 2> qn = { Q2D::identity(), Q2D::identity()};
		std::array<Q2D, 2> qn_commit = { Q2D::identity(), Q2D::identity() };
		V2D rn = V2D(0.0, 0.0);
		V2D rn_commit = V2D(0.0, 0.0);
		
		ElementComponent() = default;
		
		inline int compute(bool flag, const RVEStateVariables& rve, const ASDSteel1DMaterial::InputParameters& params, bool do_implex, double time_factor ) {

			auto& globals = Globals::instance(); 
			auto& U = globals.element_U;
			auto& U_commit = globals.element_U_commit;
			auto& UL = globals.element_U_local;
			auto& dU = globals.element_dU;
			auto& B = globals.element_B;
			auto& T = globals.element_T;
			auto& RHS = globals.element_RHS;
			auto& RHSL = globals.element_RHS_local;
			auto& LHS = globals.element_LHS;
			auto& LHSL = globals.element_LHS_local;
			auto& BtC = globals.element_BtC;
			auto& TtK = globals.element_TtK;
			auto& stress_s = globals.section_stress;
			auto& strain_s = globals.section_strain;
			auto& tangent_s = globals.section_tangent;
			auto& rve_nodes = globals.rve_nodes;
			auto& rve_nodes_el = globals.rve_nodes_el;
			int inode;
			int jnode;
			if (flag) {
				inode = ETypeTraits<EType>::N1;
				jnode = ETypeTraits<EType>::N2;
			}
			else {
				inode = ETypeTraits<EType>::N1-1;
				jnode = ETypeTraits<EType>::N2-1;
			}
			
			
			if (flag) {
				getElementDisplacementVector<EType>(flag, rve.UG_el, U);
				getElementDisplacementVector<EType>(flag, rve.UG_el_commit, U_commit);
			}
			else {

				getElementDisplacementVector<EType>(flag, rve.UG, U);
				getElementDisplacementVector<EType>(flag, rve.UG_commit, U_commit);
			}
			

			// undeformed nodes
			V2D node_i;
			V2D node_j;
			if (flag) {
				node_i = rve_nodes_el[inode];
				node_j = rve_nodes_el[jnode];
			}
			else {
				node_i = rve_nodes[inode];
				node_j = rve_nodes[jnode];
			}
			
			
			//deformed node position
			V2D deformed_node_i = node_i +V2D(U(0), U(1));
			V2D deformed_node_j = node_j +V2D(U(3), U(4));

			V2D deformed_node_i_old = node_i + V2D(U_commit(0), U_commit(1));
			V2D deformed_node_j_old = node_j + V2D(U_commit(3), U_commit(4));
			
			Q2D Q, Q0;

			V2D C;
			V2D C0;
			double L = 0.0;
			
			if (do_implex) {
				computeLocalFrame(deformed_node_i_old, deformed_node_j_old, Q, C, T);  //deformed old nodes
			}
			else {
				computeLocalFrame(deformed_node_i, deformed_node_j, Q, C, T);  //deformed nodes
			}
			computeLocalFrame(node_i, node_j, Q0, C0);  //rve nodes

			if (!do_implex) {
				V2D rn_current = V2D(U(2), U(5)); // rotation at current iteration
				V2D dtheta = rn_current - rn; // iterative correction
				rn = rn_current; // save current for next iteration
				for (int i = 0; i < 2; ++i) {
					qn[i] = Q2D::fromRotationVector(dtheta(i)) * qn[i]; // increment quaterion from previous iteration
				}

				//position vectors in local reference frame
				V2D rf_node1_position = Q0.rotateVector(node_i - C0);
				V2D rf_node2_position = Q0.rotateVector(node_j - C0);

				//deformed position vectors in local current frame
				V2D cf_node1_deformed = Q.rotateVector(deformed_node_i - C);
				V2D cf_node2_deformed = Q.rotateVector(deformed_node_j - C);

				// translational part of deformational local displacements
				UL(0) = cf_node1_deformed.x - rf_node1_position.x;
				UL(1) = cf_node1_deformed.y - rf_node1_position.y;
				UL(3) = cf_node2_deformed.x - rf_node2_position.x;
				UL(4) = cf_node2_deformed.y - rf_node2_position.y;

				// rotational part of deformational local rotation
				UL(2) = (Q * qn[0] * Q0.conjugate()).toRotationVector();
				UL(5) = (Q * qn[1] * Q0.conjugate()).toRotationVector();
			}
			else {
				dU = U;
				dU.addVector(1.0, U_commit, -1.0);
				UL = UL_commit; 
				UL.addMatrixVector(1.0, T, dU, 1.0);
			}

			computeBmatrix(node_i, node_j, L, B);  //rve nodes
			
			// section strain = B*U
			strain_s.addMatrixVector(0.0, B, UL, 1.0);

			//section compute strain - stress tangent
			int retval = section.compute(params, strain_s, do_implex, time_factor, stress_s, tangent_s);
			if (retval != 0) {
				asd_print_full("section !converged");
				return retval;
			}
			
			// RHS  
			RHSL.addMatrixTransposeVector(0.0, B, stress_s, L);
			RHS.addMatrixTransposeVector(0.0, T, RHSL, 1.0);

			// LHS
			BtC.addMatrixTransposeProduct(0.0, B, tangent_s, L);
			LHSL.addMatrixProduct(0.0, BtC, B, 1.0);
			TtK.addMatrixTransposeProduct(0.0, T, LHSL, 1.0);
			LHS.addMatrixProduct(0.0, TtK, T, 1.0);

			// SEMI-COMMIT
			if (params.implex && !do_implex) {
				UL_commit = UL;
			}
			
			return 0;
		} 
		inline int commitState() {
			qn_commit = qn;
			rn_commit = rn;
			return section.commitState();
		}
		inline void revertToLastCommit() {
			qn = qn_commit;
			rn = rn_commit;
			section.revertToLastCommit();  
		}
		inline void revertToStart() {
			qn = { Q2D::identity(), Q2D::identity() };
			qn_commit = { Q2D::identity(), Q2D::identity() };
			rn = V2D(0.0, 0.0);
			rn_commit = V2D(0.0, 0.0);
			UL_commit.Zero();
			section.revertToStart(); 
		}

	};

	template<int NFiber, int EType>
	inline int rve_process_element(bool flag, ElementComponent<NFiber, EType>& ele, const RVEStateVariables& rve, const ASDSteel1DMaterial::InputParameters& params, bool do_implex, double time_factor) {
		int retval = ele.compute(flag, rve, params, do_implex, time_factor);
		if (retval != 0) {
			asd_print_full("ele !converged");
			return retval;
		}
		auto& globals = Globals::instance();
		if (flag) {
			assembleRHS<EType>(flag, globals.element_RHS, globals.rve_R_el);
			assembleLHS<EType>(flag, globals.element_LHS, globals.rve_K_el);
		}
		else {

			assembleRHS<EType>(flag, globals.element_RHS, globals.rve_R);
			assembleLHS<EType>(flag, globals.element_LHS, globals.rve_K);
		}
		return 0;
	}

	class RVEModel
	{
	public:
		RVEStateVariables sv;
		ElementComponent<3, EType_Bot> e1;
		ElementComponent<1, EType_Mid> e2;
		ElementComponent<3, EType_Top> e3;
		ElementComponent<0, EType_El> e0;

		RVEModel() = default;
		
		inline int compute(bool flag, const ASDSteel1DMaterial::InputParameters& params, double U, bool do_implex, double time_factor, double& N, double& tangent) {

			// output
			int retval = -1;
			
			// start from previous commited state (not in implex commit)
			if (!(params.implex && !do_implex)) {
				revertToLastCommit();
			}
			
			// globals
			double tolR = params.tolR * params.sy;
			double tolU = params.tolU * params.length;
			//constexpr int max_iter = 200;
			auto& globals = Globals::instance();
			
			// utility for assembly
			// note: do it in reverse order, so that we can access the element RHS for element Bottom to get reaction
			auto lam_assemble = [this, flag, &globals, &params, &do_implex, &time_factor]() -> int {
				// zero R and K for element integration
				globals.rve_R.Zero();
				globals.rve_K.Zero();
				globals.rve_R_el.Zero();
				globals.rve_K_el.Zero();
				if (flag) {
					if (rve_process_element< 3, EType_Top>(flag, e3, sv, params, do_implex, time_factor) != 0) return -1;
					if (rve_process_element< 1, EType_Mid>(flag, e2, sv, params, do_implex, time_factor) != 0) return -1;
					if (rve_process_element< 3, EType_Bot>(flag, e1, sv, params, do_implex, time_factor) != 0) return -1;
					if (rve_process_element< 0, EType_El>(flag, e0, sv, params, do_implex, time_factor) != 0) return -1;
				}
				else {
					if (rve_process_element< 3, EType_Top>(flag, e3, sv, params, do_implex, time_factor) != 0) return -1;
					if (rve_process_element< 1, EType_Mid>(flag, e2, sv, params, do_implex, time_factor) != 0) return -1;
					if (rve_process_element< 3, EType_Bot>(flag, e1, sv, params, do_implex, time_factor) != 0) return -1;
					double I = M_PI * std::pow(params.radius, 4) / 4.0;
					double EI = params.E * I;
					double fact = 0.04 * params.length / (params.radius * 5.0);
					double penalty = 12.0 * EI * fact / std::pow(2.0*params.length, 3.0);
					double KReg = std::max(0.0, (params.length / params.lch_element - 1)) * penalty;
					globals.rve_K(6, 6) = globals.rve_K(6, 6) + KReg;
					globals.rve_R(6) = globals.rve_R(6)  + KReg * sv.UG(6);
				}                            
				return 0;
			};

			// impose BC
			if (flag) {
				sv.UG_el(10) = U;
			}
			else {
				sv.UG(7) = U;
			}
			
			// first residual
			if (lam_assemble() != 0) {
				asd_print_full("1lam assemble !=0");
				return -1;
			}

			double Rnorm = -1.0;
			double Unorm = -1.0;

			//for implicit correction convergence
			static std::array<Vector, 10> U_iterative = { Vector(7), Vector(7), Vector(7), Vector(7), Vector(7), Vector(7), Vector(7), Vector(7), Vector(7), Vector(7)};
			static std::array<Vector, 10> U_iterative_el = { Vector(10), Vector(10), Vector(10), Vector(10), Vector(10), Vector(10), Vector(10), Vector(10), Vector(10), Vector(10) };
			static std::array<double, 10> Rnorm_iterative = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	
			// newton loop for this step
			for (int iter = 0; iter < params.max_iter; ++iter) {

				// solve
				if (flag) {
					if (globals.rve_K_el.Solve(globals.rve_R_el, globals.rve_dU_el) != 0) {
						asd_print_full("Solve failed");
						break;
					}

					// update free DOFs
					for (int i = 0; i < 10; ++i) {
						sv.UG_el(i) -= globals.rve_dU_el(i);
					}

					// assemble
					if (lam_assemble() != 0) {
						asd_print_full("lam assemble !=0");
						break;
					}

					// convergence test

					Rnorm = globals.rve_R_el.Norm();
					Unorm = globals.rve_dU_el.Norm();
					asd_print("iter: " << iter << "R: " << Rnorm << ", U: " << Unorm);

					//for implicit correction convergence (max 10 iterations)
					if (params.implex && !do_implex) {
						U_iterative_el[iter] = sv.UG_el;
						Rnorm_iterative[iter] = Rnorm;
						if (iter >= 9) {
							auto minRiter = std::min_element(Rnorm_iterative.begin(), Rnorm_iterative.end());
							double minRnorm = *minRiter;
							if (!std::isfinite(minRnorm)) {
								Rnorm = std::numeric_limits<double>::max();
							}
							Rnorm = minRnorm;
							asd_print("minRnorm" << minRnorm);
							int pos_min = std::distance(Rnorm_iterative.begin(), minRiter);
							sv.UG_el = U_iterative_el[pos_min];
							lam_assemble();
							retval = 0;
							break;
						}
					}

					if (Rnorm < tolR || Unorm < tolU) {
						retval = 0;
						break;
					}
				}
				else {
					if (globals.rve_K.Solve(globals.rve_R, globals.rve_dU) != 0) {
						asd_print_full("Solve failed");
						break;
					}

					// update free DOFs
					for (int i = 0; i < 7; ++i) {
						sv.UG(i) -= globals.rve_dU(i);
					}

					// assemble
					if (lam_assemble() != 0) {
						asd_print_full("lam assemble !=0");
						break;
					}

					// convergence test

					Rnorm = globals.rve_R.Norm();
					Unorm = globals.rve_dU.Norm();
					asd_print("iter: " << iter << "R: " << Rnorm << ", U: " << Unorm);

					//for implicit correction convergence (max 10 iterations)
					if (params.implex && !do_implex) {
						U_iterative[iter] = sv.UG;
						Rnorm_iterative[iter] = Rnorm;
						if (iter >= 9) {
							auto minRiter = std::min_element(Rnorm_iterative.begin(), Rnorm_iterative.end());
							double minRnorm = *minRiter;
							if (!std::isfinite(minRnorm)) {
								Rnorm = std::numeric_limits<double>::max();
							}
							Rnorm = minRnorm;
							asd_print("minRnorm" << minRnorm);
							int pos_min = std::distance(Rnorm_iterative.begin(), minRiter);
							sv.UG = U_iterative[pos_min];
							lam_assemble();
							retval = 0;
							break;
						}
					}

					if (Rnorm < tolR || Unorm < tolU) {
						retval = 0;
						break;
					}

				}
				
				
			}
			
			if (retval != 0) {
				asd_print_full("rve !converged");
			}
			double Krc_dot_Kinv_Kcr = 0.0;
			// do it outside, element Bottom is the last one
			if (flag) {
				getRetainedComponents<EType_El>(flag, globals.element_LHS, globals.rve_Krr, globals.rve_Kcr_el, globals.rve_Krc_el);
				globals.rve_K_el.Solve(globals.rve_Kcr_el, globals.rve_Kinv_Kcr_el);
				Krc_dot_Kinv_Kcr = globals.rve_Krc_el ^ globals.rve_Kinv_Kcr_el;
			}
			else {

				getRetainedComponents<EType_Bot>(flag, globals.element_LHS, globals.rve_Krr, globals.rve_Kcr, globals.rve_Krc);
				globals.rve_K.Solve(globals.rve_Kcr, globals.rve_Kinv_Kcr);
				Krc_dot_Kinv_Kcr = globals.rve_Krc ^ globals.rve_Kinv_Kcr;
			}
			// get N from element_RHS (note: RHS = Fext-Fint)
			N = globals.element_RHS(1);
			// perform static condensation
			// Kred = Krr - (Krc *  inv(Kcc)*Kcr)
			
			
			globals.rve_Kred = globals.rve_Krr - Krc_dot_Kinv_Kcr;
			tangent = globals.rve_Kred;
			return retval;
		}
		inline int commitState() {
			
			int retval = 0;
			
			sv.UG_commit = sv.UG;
			sv.UG_el_commit = sv.UG_el;
			int el_retval;
			el_retval = e1.commitState();
			if (el_retval != 0) retval = el_retval;
			el_retval = e2.commitState();
			if (el_retval != 0) retval = el_retval;
			el_retval = e3.commitState();
			if (el_retval != 0) retval = el_retval;
			el_retval = e0.commitState();
			if (el_retval != 0) retval = el_retval;

			return retval;
		}
		inline void revertToLastCommit() {
			sv.UG = sv.UG_commit;
			sv.UG_el = sv.UG_el_commit;
			e0.revertToLastCommit();
			e1.revertToLastCommit();
			e2.revertToLastCommit();
			e3.revertToLastCommit();
		}
		inline void revertToStart() {
			sv.UG.Zero();
			sv.UG_commit.Zero();
			sv.UG_el.Zero();
			sv.UG_el_commit.Zero();
			e0.revertToStart();
			e1.revertToStart();
			e2.revertToStart();
			e3.revertToStart();
		}

	};
}

/**
	PIMPL pattern
	*/
class ASDSteel1DMaterialPIMPL
{
public:
	ASDSteel1DMaterialPIMPL() = default;
	ASDSteel1DMaterialPIMPL(const ASDSteel1DMaterialPIMPL&) = default;
public:
	SteelComponent steel_comp;
	RVEModel rve_m;	
};

void* OPS_ASDSteel1DMaterial()
{
	
	// some kudos
	static bool first_done = false;
	if (!first_done) {
		opserr << "Using ASDSteel1D - Developed by: Alessia Casalucci, Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
		first_done = true;
	} 
	static const char* msg = "uniaxialMaterial ASDSteel1D $tag $E $sy $su $eu $lch $r  <-implex $implex>  <-buckling $buckling> <-K_alpha $K_alpha> <-max_iter $max_iter> <-tolU $tolU> <-tolR $tolR>";

	// check arguments
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 7) {
		opserr <<
			"uniaxialMaterial ASDSteel1D Error: Few arguments (< 7).\n" << msg << "\n";
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
	double lch;
	double r;
	bool implex = false;
	bool buckling = false;
	double K_alpha = 0.5;
	double max_iter= 100;
	double tolU = 1.0e-6;
	double tolR = 1.0e-6;
	bool have_K_alpha = false;
	bool have_max_iter = false;
	bool have_tolU = false;
	bool have_tolR = false;


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
	if (!lam_get_dparam(&lch, "lch")) return nullptr;
	if (!lam_get_dparam(&r, "r")) return nullptr;
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
	
	// parse optional arguments
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* value = OPS_GetString();
		if (strcmp(value, "-implex") == 0) {
			implex = true;
		}
		if (strcmp(value, "-buckling") == 0) {
			buckling = true;
		}
		if (strcmp(value, "-K_alpha") == 0) {
			if (!lam_optional_double("K_alpha", K_alpha))
				return nullptr;
			have_K_alpha = true;
		}
		if (strcmp(value, "-max_iter") == 0) {
			if (!lam_optional_double("max_iter", max_iter))
				return nullptr;
			have_max_iter = true;
		}
		if (strcmp(value, "-tolU") == 0) {
			if (!lam_optional_double("tolU", tolU))
				return nullptr;
			have_tolU = true;
		}
		if (strcmp(value, "-tolR") == 0) {
			if (!lam_optional_double("tolR", tolR))
				return nullptr;
			have_tolR = true;
		}
	}

	// checks
	if (sy >= su) {
		opserr << "nDMaterial ASDSteel1D Error: invalid value for 'su' (" << su << "). It should be larger than sy.\n" << msg << "\n";
		return nullptr;
	}
	if (E == 0) {
		opserr << "nDMaterial ASDSteel1D Error: invalid value for 'E' (" << E << "). It should be non-zero.\n" << msg << "\n";
		return nullptr;
	}
	if (lch == 0) {
		opserr << "nDMaterial ASDSteel1D Error: invalid value for 'lch' (" << lch << "). It should be non-zero.\n" << msg << "\n";
		return nullptr;
	}
	if (r == 0) {
		opserr << "nDMaterial ASDSteel1D Error: invalid value for 'r' (" << r << "). It should be non-zero.\n" << msg << "\n";
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
	params.implex = implex;
	params.buckling = buckling;
	params.length = lch / 2.0; // consider half distance, the RVE uses symmetry
	params.radius = r;
	params.K_alpha = K_alpha;
	params.max_iter = max_iter;
	params.tolU = tolU;
	params.tolR = tolR;
	params.lch_element = params.length;
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

ASDSteel1DMaterial::ASDSteel1DMaterial(
	int _tag,
	const InputParameters& _params)
	: UniaxialMaterial(_tag, MAT_TAG_ASDSteel1DMaterial)
	, params(_params)
	, pdata(new ASDSteel1DMaterialPIMPL())
{
	// intialize C as C0
	C = getInitialTangent();
}

ASDSteel1DMaterial::ASDSteel1DMaterial()
	: UniaxialMaterial(0, MAT_TAG_ASDSteel1DMaterial)
	, pdata(new ASDSteel1DMaterialPIMPL())
{
}

ASDSteel1DMaterial::ASDSteel1DMaterial(const ASDSteel1DMaterial& other)
	: UniaxialMaterial(other.getTag(), MAT_TAG_ASDSteel1DMaterial)
	, params(other.params)
	, dtime_n(other.dtime_n)
	, dtime_n_commit(other.dtime_n_commit)
	, commit_done(other.commit_done)
	, strain(other.strain)
	, strain_commit(other.strain_commit)
	, stress(other.stress)
	, stress_commit(other.stress_commit)
	, C(other.C)
	, stress_rve(other.stress_rve)
	, stress_rve_commit(other.stress_rve_commit)
	, C_rve(other.C_rve)
	, energy(other.energy)
	, pdata(other.pdata ? new ASDSteel1DMaterialPIMPL(*other.pdata) : nullptr)
{
}

ASDSteel1DMaterial::~ASDSteel1DMaterial()
{
	if(pdata) delete pdata;
}

int ASDSteel1DMaterial::setTrialStrain(double v, double r)
{
	params.lch_element = ops_TheActiveElement ? ops_TheActiveElement->getCharacteristicLength()/2 : params.length;
	// save dT
	dtime_n = ops_Dt;
	if (!commit_done) {
		dtime_n_commit = dtime_n;
	}
	double time_factor = 1.0;
	if (params.implex && (dtime_n_commit > 0.0))
		time_factor = dtime_n / dtime_n_commit;
	
	// save macro strain
	strain = v;
	// homogenize micro response (stress/tangent)

	asd_print("MAT set trial (" << (int)params.implex << ")");
	int retval;
	// time factor for explicit extrapolation (todo: make a unique function)
	
	if (params.buckling) {
		//regularization for element length lower than physical length (weighted average)
		retval = homogenize(params.implex);
		C = C_rve;
		stress = stress_rve;
	}
	else {
		opserr << "HERE\n";		
		////computation of sigma and tangent of steel
		pdata->steel_comp.revertToLastCommit();
		pdata->steel_comp.compute(params, params.implex, time_factor, strain, stress, C);
	}
	
	// then apply damage on stress e C

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
		asd_print("MAT commit (" << (int)false << ")");
		int retval = homogenize(false);
		if (retval < 0) return retval;
		double sigma_macro = 0.0;
		double tangent_macro = 0.0;
		retval = pdata->steel_comp.compute(params, false, 1.0, strain, sigma_macro, tangent_macro);
	}

	// compute energy
	energy += 0.5 * (stress_commit + stress) * (strain - strain_commit);

	// forward to all components
	pdata->rve_m.commitState();
	pdata->steel_comp.commitState();

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
	// forward to all components
	pdata->rve_m.revertToLastCommit();
	pdata->steel_comp.revertToLastCommit();

	// state variables
	strain = strain_commit;
	stress = stress_commit;
	stress_rve = stress_rve_commit;

	// implex
	dtime_n = dtime_n_commit;

	// done
	return 0;
}

int ASDSteel1DMaterial::revertToStart(void)
{
	// forward to all components
	pdata->rve_m.revertToStart();
	pdata->steel_comp.revertToStart();

	// strain, stress and tangent
	strain = 0.0;
	strain_commit = 0.0;
	stress = 0.0;
	stress_commit = 0.0;
	stress_rve = 0.0;
	stress_rve_commit = 0.0;
	C = getInitialTangent();
	C_rve = getInitialTangent();

	// implex
	dtime_n = 0.0;
	dtime_n_commit = 0.0;
	commit_done = false;

	// output variables
	energy = 0.0;

	// done
	return 0;
}

UniaxialMaterial* ASDSteel1DMaterial::getCopy(void)
{
	// we can safely use the default copy-constructor
	return new ASDSteel1DMaterial(*this);
}

void ASDSteel1DMaterial::Print(OPS_Stream& s, int flag)
{
	s << "ASDSteel1D Material, tag: " << this->getTag() << "\n";
}

int ASDSteel1DMaterial::sendSelf(int commitTag, Channel &theChannel)
{
	//// aux
	//int counter;

	//// send DBL data
	//Vector ddata(InputParameters::NDATA + 1*StateVariablesSteel::NDATA + 10);
	//counter = 0;
	//ddata(counter++) = static_cast<double>(getTag());
	//ddata(counter++) = params.E;
	//ddata(counter++) = params.sy;
	//ddata(counter++) = params.H1;
	//ddata(counter++) = params.H2;
	//ddata(counter++) = params.gamma1;
	//ddata(counter++) = params.gamma2;
	//ddata(counter++) = static_cast<double>(params.implex);
	//steel.sendSelf(counter, ddata);
	//ddata(counter++) = dtime_n;
	//ddata(counter++) = dtime_n_commit;
	//ddata(counter++) = static_cast<double>(commit_done);
	//ddata(counter++) = strain;
	//ddata(counter++) = strain_commit;
	//ddata(counter++) = stress;
	//ddata(counter++) = stress_commit;
	//ddata(counter++) = C;
	//ddata(counter++) = energy;
	//if (theChannel.sendVector(getDbTag(), commitTag, ddata) < 0) {
	//	opserr << "ASDSteel1DMaterial::sendSelf() - failed to send DBL data\n";
	//	return -1;
	//}

	//// done
	return 0;
}

int ASDSteel1DMaterial::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	//// aux
	//int counter;

	//// recv DBL data
	//Vector ddata(InputParameters::NDATA + 1 * StateVariablesSteel::NDATA + 10);
	//if (theChannel.recvVector(getDbTag(), commitTag, ddata) < 0) {
	//	opserr << "ASDSteel1DMaterial::recvSelf() - failed to receive DBL data\n";
	//	return -1;
	//}
	//counter = 0;
	//setTag(ddata(counter++));
	//params.E = ddata(counter++);
	//params.sy = ddata(counter++);
	//params.H1 = ddata(counter++);
	//params.H2 = ddata(counter++);
	//params.gamma1 = ddata(counter++);
	//params.gamma2 = ddata(counter++);
	//params.implex = static_cast<bool>(ddata(counter++));
	//steel.recvSelf(counter, ddata);
	//dtime_n = ddata(counter++);
	//dtime_n_commit = ddata(counter++);
	//commit_done = static_cast<bool>(ddata(counter++));
	//strain = ddata(counter++);
	//strain_commit = ddata(counter++);
	//stress = ddata(counter++);
	//stress_commit = ddata(counter++);
	//C = ddata(counter++);
	//energy = ddata(counter++);

	//// done
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
	static std::vector<std::string> lb_buckling_ratio = { "BI" };

	// all outputs are 1D
	static Vector out1(1);

	// check specific responses
	if (argc > 0) {
		// 1000 - base steel output
		if (strcmp(argv[0], "BI") == 0 || strcmp(argv[0], "BucklingIndicator") == 0) {
			out1(0) = pdata->rve_m.sv.UG(6);
			return make_resp(1001, out1, &lb_buckling_ratio);
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
		out1(0) = pdata->rve_m.sv.UG(6);
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

int ASDSteel1DMaterial::homogenize(bool do_implex)
{	
	// return value
	int retval = 0;

	// check for elastic correction
	bool elastic_correction = params.lch_element > params.length;

	auto& globals = Globals::instance();
	if (elastic_correction) {
		globals.setRVENodes_el(params.length,params.lch_element);
	}
	else {
		globals.setRVENodes(params.length);
	}	

	// time factor for explicit extrapolation
	double time_factor = 1.0;
	if (params.implex && do_implex && (dtime_n_commit > 0.0))
		time_factor = dtime_n / dtime_n_commit;

	// from macro strain to micro strain
	double macro_strain = strain;
	double Uy;
	if (elastic_correction) {
		Uy = -macro_strain * params.lch_element;
	}
	else {
		Uy = -macro_strain * params.length;
	}
	double N = 0.0;
	double T = 0.0;
	double area = 0.0;
	retval = pdata->rve_m.compute(elastic_correction, params, Uy, do_implex, time_factor, N, T);

	if (retval != 0) {
		return retval;
	}

	// from micro stress/tangent to macro stress/tangent
	area = M_PI * params.radius * params.radius;

	stress_rve = -N/area;
	if (elastic_correction) {
		C_rve = T * params.lch_element / area;
	}
	else {
		C_rve = T * params.length / area;
	}
	

	// done
	return retval;
}




