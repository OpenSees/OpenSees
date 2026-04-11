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

// Alessia Casalucci, Massimo Petracca, Guido Camata - ASDEA Software, Italy
//
// A unified and efficient plastic-damage material model for steel bars including fracture, bond-slip, and buckling via multiscale homogenization
//

#include <ASDSteel1DMaterial.h>
#include <Channel.h>
#include <OPS_Globals.h>
#include <FEM_ObjectBroker.h>
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

//#define ASD_STEEL1D__VERBOSE

#ifdef ASD_STEEL1D__VERBOSE
#define asd_print_full(X) opserr  << __func__ << " (Line " << __LINE__ << "): " << X << "\n"
#define asd_print(X) opserr << X << "\n"
#else
#define asd_print_full(X)
#define asd_print(X)
#endif // ASD_STEEL1D__VERBOSE

// anonymous namespace for utilities
namespace {


	enum ErrorCodes {
		EC_Generic = -1,
		EC_IMPLEX_Error_Control = -10
	};

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
		int serializationDataSize() const;
		void serialize(Vector& data, int& pos);
		void deserialize(Vector& data, int& pos);
		inline int commitState() {
			// store the previously committed variables for next move from n to n - 1
			lambda_commit_old = lambda_commit;
			// state variables
			epl_commit = epl;
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
			epl = epl_commit;
			alpha1 = alpha1_commit;
			alpha2 = alpha2_commit;
			lambda = lambda_commit;
			strain = strain_commit;
			stress = stress_commit;
		}
		inline void revertToStart() {
			// state variables
			epl = 0.0;
			epl_commit = 0.0;
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
			// base steel response
			alpha1 = alpha1_commit;
			alpha2 = alpha2_commit;
			lambda = lambda_commit;
			epl = epl_commit;
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
				// extrapolate plastic flow direction
				sg = sg_commit;
				// extrapolate lambda
				//  xn + time_factor * (xn - xnn);
				double delta_lambda = time_factor * (lambda_commit - lambda_commit_old);
				lambda = lambda_commit + delta_lambda;
				// update plastic strain
				if (sg > 0.0)
					epl = epl_commit + sg * delta_lambda;
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
							sg = sign(lam_rel_stress());
							// update plastic multiplier
							lambda += delta_lambda;
							// update plastic strain
							if (sg > 0.0)
								epl += sg * delta_lambda;
							// compute tangent
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
		// positive plastic strain (for fracture)
		double epl = 0.0;
		double epl_commit = 0.0;
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
		static constexpr int NDATA = 14;
	};
	int SteelComponent::serializationDataSize() const
	{
		return NDATA; //= 14
	}

	void SteelComponent::serialize(Vector& data, int& pos)
	{
		data(pos++) = epl;
		data(pos++) = epl_commit;
		data(pos++) = alpha1;
		data(pos++) = alpha1_commit;
		data(pos++) = alpha2;
		data(pos++) = alpha2_commit;
		data(pos++) = lambda;
		data(pos++) = lambda_commit;
		data(pos++) = lambda_commit_old;
		data(pos++) = sg_commit;
		data(pos++) = strain;
		data(pos++) = strain_commit;
		data(pos++) = stress;
		data(pos++) = stress_commit;
	}

	void SteelComponent::deserialize(Vector& data, int& pos)
	{
		epl = data(pos++);
		epl_commit = data(pos++);
		alpha1 = data(pos++);
		alpha1_commit = data(pos++);
		alpha2 = data(pos++);
		alpha2_commit = data(pos++);
		lambda = data(pos++);
		lambda_commit = data(pos++);
		lambda_commit_old = data(pos++);
		sg_commit = data(pos++);
		strain = data(pos++);
		strain_commit = data(pos++);
		stress = data(pos++);
		stress_commit = data(pos++);
	}
	/*
	Series component.
	*/
	class SeriesComponent
	{
	public:
		double lch_anchor = 0.0;
		SteelComponent steel_material;
		UniaxialMaterial* slip_material = nullptr;
		int serializationDataSize() const;
		void serialize(Vector& data, int& pos, int commitTag, Channel& theChannel);
		void deserialize(Vector& data, int& pos, int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);
		void serialize_slip(int commitTag, Channel& theChannel);
		void deserialize_slip(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);
		
	public:
		using param_t = ASDSteel1DMaterial::InputParameters;
		inline void initialize_lch_anchor(const param_t& params) {
			if (lch_anchor == 0.0) { 
				lch_anchor = lch_anchor_compute(params);
				if (slip_material) {
					if (params.auto_regularization) {
						// if the slip material supports the lch_ref parameter, let's set it to the anchorage length
						Parameter lch_param(0, 0, 0, 0);
						const char* the_args[] = { "lch_ref" };
						if (slip_material->setParameter(the_args, 1, lch_param) != -1) {
							lch_param.update(lch_anchor);
						}
					}
				}
			}
		}
		SeriesComponent() = default;
		SeriesComponent(const SeriesComponent& c)
			: steel_material(c.steel_material)
			, slip_material(nullptr)
		{
			setSlipMaterial(c.slip_material);
		}
		~SeriesComponent() {
			if (slip_material)
				delete slip_material;
		}
	
		inline int commitState() {
			// store the previously committed variables for next move from n to n - 1

			// state variables

			int retval = 0;

			retval = steel_material.commitState();

			if (slip_material)
				retval = slip_material->commitState();
		
			// done 
			return retval;
		}
		inline void revertToLastCommit() {
			// state variables
			steel_material.revertToLastCommit();

			if (slip_material)
				slip_material->revertToLastCommit();
		
		}
		inline void revertToStart() {
			// state variables
			steel_material.revertToStart();

			if (slip_material)
				slip_material->revertToStart();
		
		}
		inline double lch_anchor_compute(const param_t& params) {
			//lch_anchor estimation 

			// compute the gradient
			double max_slip = params.radius * 100.0;
			double tol = params.sy * 1.0e-6;
			auto compute_gradient = [&](double slip) {
				double epsilon = max_slip * 1.0e-8;
				slip_material->setTrialStrain(slip);
				slip_material->commitState();
				double stress_0 = slip_material->getStress();
				slip_material->revertToStart();
				slip_material->setTrialStrain(slip + epsilon);
				slip_material->commitState();
				double stress_1 = slip_material->getStress();
				slip_material->revertToStart();
				double K = (stress_1 - stress_0) / epsilon;
				if (std::abs(K) < tol)
					K = 0.0;
				return K;
			};

			double X0 = 0.0;
			double Y0 = compute_gradient(X0);
			double S0 = sign(Y0);
			double X1 = max_slip;
			double Y1 = compute_gradient(X1);
			double S1 = sign(Y1);
			double XM = (X1 + X0)/2.0;
			bool found_peak = false;
			for (int i = 0; i < 20; i++) {
				double YM = compute_gradient(XM);				
				double SM = sign(YM);
				if (S0 != SM) {
					X1 = XM;
					Y1 = YM;
					S1 = SM;
					found_peak = true;
				}
				else if (SM != S1) {
					X0 = XM;
					Y0 = YM;
					S0 = SM;
					found_peak = true;
				}
				else {
					opserr << "Failed " << S0 << ", " << SM << ", " << S1 << "\n";
					break;
				}
				XM = (X1 + X0) / 2.0;
			}

			double tau_max = 0.0;
			double lch_anchor = 0.0;
			if (found_peak) {
				slip_material->setTrialStrain(XM);
				slip_material->commitState();
				tau_max = slip_material->getStress();
				slip_material->revertToStart();

				lch_anchor = (params.sy * params.radius) / (2.0 * tau_max);
			}
			else {
				opserr << "No peak found — using default lch_anchor\n";
				lch_anchor = 80.0 * params.radius;
			}
			return lch_anchor;
		}
		inline int compute(const param_t& params, bool do_implex, double time_factor, double _strain, double& sigma, double& tangent) {
			if (slip_material) {

				initialize_lch_anchor(params);

				// compute series response
				constexpr int MAX_ITER = 100;
				constexpr double  TOL = 1e-6;
				double scale_factor = 1.0;
				// if in RVE for buckling, only 1 element has the slip material (1/3 of the rve length)
				if (params.buckling) {
					if (params.lch_element > params.length) {
						scale_factor = params.length / (3.0 * params.lch_element);
					}
					else {
						scale_factor = 1.0 / 3.0;
					}
				}					
					
				_strain *= scale_factor;

				// initial guess for the slip strain
				double lch_ele = 2.0 * params.lch_element; // it was divided in constructor for RVE symmetry
				double strain_slip = slip_material->getStrain() / lch_anchor;
				double strain_steel = _strain - strain_slip;
				
				// iterative procedure to impose the iso-stress condition
				double sigma_steel;
				double tangent_steel;
				double Ktol = 1.0e-12 * params.E;
				double Rtol = TOL * params.sy;
				double Stol = TOL;
				for (int iter = 0; iter < MAX_ITER; ++iter) {
					int Tsteel = steel_material.compute(params, do_implex, time_factor, strain_steel, sigma_steel, tangent_steel);
					int Tslip = slip_material->setTrialStrain(strain_slip * lch_anchor);
					double sigma_slip = slip_material->getStress() * 2.0 * lch_anchor/ params.radius;
					double tangent_slip = slip_material->getTangent() * 2.0 * lch_anchor *lch_anchor / params.radius;
					double residual = sigma_slip - sigma_steel;
					double residual_derivative = tangent_steel + tangent_slip;
					if (std::abs(residual_derivative) < Ktol) {
						residual_derivative = params.E + slip_material->getInitialTangent() * 2.0 * lch_anchor * lch_anchor / params.radius;
					}
					double strain_increment = - residual / residual_derivative;
					strain_slip += strain_increment;
					strain_steel = _strain - strain_slip;

					if (std::abs(residual) < Rtol || std::abs(strain_increment) < Stol) {
						sigma = sigma_steel;
						if (std::abs(tangent_steel) < Ktol) tangent_steel = Ktol;
						if (std::abs(tangent_slip) < Ktol) tangent_slip = Ktol;
						tangent = scale_factor *  1.0 / (1.0 / tangent_steel + 1.0 / tangent_slip);
						asd_print("iter: " << iter);
						return 0;
					}				
				}
				asd_print_full("series failed");
				return -1;
			}
			else {
				return steel_material.compute(params, do_implex, time_factor, _strain, sigma, tangent);
			}
		}

		inline void setSlipMaterial(UniaxialMaterial* prototype) {
			if (slip_material) {
				delete slip_material;
				slip_material = nullptr;
			}
			if (prototype) {
				slip_material = prototype->getCopy();
				asd_print("setting slip material: " << slip_material->getClassType());
			}
		}

	};
	int SeriesComponent::serializationDataSize() const
	{
		return 3 + steel_material.serializationDataSize();
	}

	void SeriesComponent::serialize(Vector& data, int& pos, int commitTag, Channel& theChannel)
	{
		steel_material.serialize(data, pos);
		if (slip_material) {
			data(pos++) = static_cast<double>(slip_material->getClassTag());
			int mat_db_tag = slip_material->getDbTag();
			if (mat_db_tag == 0) {
				mat_db_tag = theChannel.getDbTag();
				if (mat_db_tag != 0)
					slip_material->setDbTag(mat_db_tag);
			}
			data(pos++) = static_cast<double>(mat_db_tag);
			data(pos++) = lch_anchor;
		}
		else {
			data(pos++) = -1.0;  // classTag
			data(pos++) = -1.0;  // dbTag
			data(pos++) = -1.0;  //lch_anchor
		}
	}
	
	void SeriesComponent::deserialize(Vector& data, int& pos, int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker) {
		steel_material.deserialize(data, pos);
		if (slip_material) {
			delete slip_material;
			slip_material = nullptr;
		}
		int classTag = static_cast<int>(data(pos++));
		int mat_db_tag = static_cast<int>(data(pos++));
		if (classTag >= 0) {
			UniaxialMaterial*  new_slip_material = theBroker.getNewUniaxialMaterial(classTag);
			if (!new_slip_material) {
				opserr << "SeriesComponent::deserialize - failed to get new UniaxialMaterial from broker\n";
				return;
			}
			new_slip_material->setDbTag(mat_db_tag);
			slip_material = new_slip_material;
		}
		lch_anchor = data(pos++);

	}

	void SeriesComponent::serialize_slip(int commitTag, Channel& theChannel)
	{
		if (slip_material) {
			int send_result = slip_material->sendSelf(commitTag, theChannel);
			if (send_result < 0) {
				opserr << "SeriesComponent:: serialize - failed to send slip material \n";
			}
		}
	}

	void SeriesComponent::deserialize_slip(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker) {

		if (slip_material) {
			if (slip_material->recvSelf(commitTag, theChannel, theBroker) < 0) {
				opserr << "SeriesComponent::deserialize - failed to receive slip material\n";
				delete slip_material;
				slip_material = nullptr;
				return;
			}
		}
	}

	/**
	Section component.
	*/
	template<int NFiber>
	class SectionComponent
	{
	private:
		SectionComponent() = delete;
	public:
		inline int serializationDataSize() const { return 0; }
		inline void serialize(Vector& data, int& pos) {}
		inline void deserialize(Vector& data, int& pos) {}
	};

	/**
	Section component with 1 fiber.
	*/
	template<>
	class SectionComponent<1>
	{
	public:
		SeriesComponent series;

	public:
		SectionComponent() = default;

		inline int commitState() {
			return series.commitState();
		}
		inline void revertToLastCommit() {
			series.revertToLastCommit();
		}
		inline void revertToStart() {
			series.revertToStart();
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
			int retval;
			retval = series.compute(params, do_implex, time_factor, strain(0), fiber_stress, fiber_tangent);

			double fiber_force = fiber_stress * area;
			stress(0) = fiber_force;
			tangent(0, 0) = fiber_tangent * area;
			// zero terms
			tangent(0, 1) = tangent(1, 0) = tangent(0, 2) = tangent(2, 0) = tangent(1, 2) = tangent(2, 1) = 0.0;
			// done
			return retval;
		}
		inline void setSlipMaterial(UniaxialMaterial* prototype)
		{
			series.setSlipMaterial(prototype);
		}

		inline int serializationDataSize() const { return series.steel_material.serializationDataSize(); }
		inline void serialize(Vector& data, int& pos) { series.steel_material.serialize(data, pos); }
		inline void deserialize(Vector& data, int& pos) { series.steel_material.deserialize(data, pos); }
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
		inline int serializationDataSize() const { 
			return fibers[0].serializationDataSize() +
				fibers[1].serializationDataSize() +
				fibers[2].serializationDataSize();
		}
		inline void serialize(Vector& data, int& pos) {
			for (int i = 0; i < fibers.size(); ++i) {
				fibers[i].serialize(data, pos);
			}
		}
		inline void deserialize(Vector& data, int& pos) {
			for (int i = 0; i < fibers.size(); ++i) {
				fibers[i].deserialize(data, pos);
			}
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
	inline void getElementDisplacementVector(const ASDSteel1DMaterial::InputParameters& params, bool elastic_correction, const Vector& G, Vector& L) {
		// G = global vector, size = 8 -> full = 12, semi-full 8 (7 free + 1 Uy imposed)
		// L = local vector , size = 6
		if (elastic_correction && params.auto_regularization) {
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
	inline void assembleRHS(const ASDSteel1DMaterial::InputParameters& params, bool elastic_correction, const Vector& L, Vector& G) {
		// G = global vector, size = 7 (-> full = 12, semi-full 8 (7 free + 1 Uy imposed) only free)
		// L = local vector , size = 6
		if (elastic_correction && params.auto_regularization) {
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
	inline void assembleLHS(const ASDSteel1DMaterial::InputParameters& params, bool elastic_correction, const Matrix& M, Matrix& G) {
		// G = global matrix, size = 7x7 -> full = 12x12, semi-full 8 (7 free + 1 Uy imposed)
		// L = local  vector , size = 6
		// M = local  matrix , size = 6x6
		if (elastic_correction && params.auto_regularization) {
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
	inline void getRetainedComponents(const ASDSteel1DMaterial::InputParameters& params,bool elastic_correction, const Matrix& M, double& Krr, Vector& Kcr, Vector& Krc) {
		if (elastic_correction && params.auto_regularization) {
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
	inline void computeLocalFrameLinearKin(const V2D& X1, const V2D& X2, Q2D& Qi, V2D& C, Matrix& T) {
		M2D Ri = computeLocalFrameInternal(X1, X2, Qi, C);
		//T matrix :  R
		for (int i = 0; i < 2;++i) {
			for (int j = 0; j < 2;++j) {
				T(i, j) = T(i + 3, j + 3) = Ri(i, j);
				T(i, j + 3) = T(i + 3, j) = 0.0;
			}
		}
		T(2, 2) = T(5, 5) = 1.0;
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

		int serializationDataSize() const;
		void serialize(Vector& data, int& pos);
		void deserialize(Vector& data, int& pos);
	};
	int RVEStateVariables::serializationDataSize() const
	{
		return 38;
	}

	void RVEStateVariables::serialize(Vector& data, int& pos)
	{
		for (int i = 0; i < 8; ++i) data(pos++) = UG[i];
		for (int i = 0; i < 8; ++i) data(pos++) = UG_commit[i];
		for (int i = 0; i < 11; ++i) data(pos++) = UG_el[i];
		for (int i = 0; i < 11; ++i) data(pos++) = UG_el_commit[i];
	}
	void RVEStateVariables::deserialize(Vector& data, int& pos)
	{
		for (int i = 0; i < 8; ++i) UG[i] = data(pos++);
		for (int i = 0; i < 8; ++i) UG_commit[i] = data(pos++);
		for (int i = 0; i < 11; ++i) UG_el[i] = data(pos++);
		for (int i = 0; i < 11; ++i) UG_el_commit[i] = data(pos++);
	}

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

		int serializationDataSize() const;
		void serialize(Vector& data, int& pos);
		void deserialize(Vector& data, int& pos);

		
		ElementComponent() = default;
		
		inline int compute(bool elastic_correction, const RVEStateVariables& rve, const ASDSteel1DMaterial::InputParameters& params, bool do_implex, double time_factor ) {

			// get global variables
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

			// determine node numbering (changes for elastic correction)
			int inode;
			int jnode;
			if (elastic_correction && params.auto_regularization) {
				inode = ETypeTraits<EType>::N1;
				jnode = ETypeTraits<EType>::N2;
			}
			else {
				inode = ETypeTraits<EType>::N1-1;
				jnode = ETypeTraits<EType>::N2-1;
			}
			
			// obtain global displacement vectors (size changes for elastic correction)
			if (elastic_correction && params.auto_regularization) {
				getElementDisplacementVector<EType>(params, elastic_correction, rve.UG_el, U);
				getElementDisplacementVector<EType>(params, elastic_correction, rve.UG_el_commit, U_commit);
			}										 
			else {									
				getElementDisplacementVector<EType>(params, elastic_correction, rve.UG, U);
				getElementDisplacementVector<EType>(params, elastic_correction, rve.UG_commit, U_commit);
			}
			

			// undeformed nodes
			V2D node_i;
			V2D node_j;
			if (elastic_correction && params.auto_regularization) {
				node_i = rve_nodes_el[inode];
				node_j = rve_nodes_el[jnode];
			}
			else {
				node_i = rve_nodes[inode];
				node_j = rve_nodes[jnode];
			}
			
			// current orientation and center
			Q2D Q;
			V2D C;
			if (EType == EType_El) {

				// use linear kinematics
				computeLocalFrameLinearKin(node_i, node_j, Q, C, T);  // undeformed nodes
				UL.addMatrixVector(0.0, T, U, 1.0);
			}
			else {
				
				// use nonlinear (corotational) kinematics

				if (!do_implex) {

					// computed deformed node position
					V2D deformed_node_i = node_i + V2D(U(0), U(1));
					V2D deformed_node_j = node_j + V2D(U(3), U(4));

					// compute local frame in current configuration
					computeLocalFrame(deformed_node_i, deformed_node_j, Q, C, T);

					// compute also local frame in the undeformed configuration
					Q2D Q0;
					V2D C0;
					computeLocalFrame(node_i, node_j, Q0, C0);

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

					// computed deformed node position with committed displacement
					V2D deformed_node_i_old = node_i + V2D(U_commit(0), U_commit(1));
					V2D deformed_node_j_old = node_j + V2D(U_commit(3), U_commit(4));

					// compute local frame in deformed committed configuration
					computeLocalFrame(deformed_node_i_old, deformed_node_j_old, Q, C, T);

					// transform only the incremental part with this explicit approximation
					dU = U;
					dU.addVector(1.0, U_commit, -1.0);
					UL = UL_commit;
					UL.addMatrixVector(1.0, T, dU, 1.0);
				}

			}
			
			// strain-displacement matrix in local frame
			double L = 0.0;
			computeBmatrix(node_i, node_j, L, B);  //rve nodes
			
			// section strain = B*U
			strain_s.addMatrixVector(0.0, B, UL, 1.0);
			int retval;
		
			retval = section.compute(params, strain_s, do_implex, time_factor, stress_s, tangent_s);
	
			
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
	inline int ElementComponent<NFiber, EType>::serializationDataSize() const{
		int size = 18 + section.serializationDataSize();		
		return size;
	}

	template<int NFiber, int EType>
	inline void ElementComponent<NFiber, EType>::serialize(Vector& data, int& pos) {
		for (int i = 0; i < 6; ++i) data(pos++) = UL_commit(i);
		for (int i = 0; i < 2; ++i) data(pos++) = qn[i].toRotationVector();
		for (int i = 0; i < 2; ++i) data(pos++) = qn_commit[i].toRotationVector();
		data(pos++) = rn.x;
		data(pos++) = rn.y;
		data(pos++) = rn_commit.x;
		data(pos++) = rn_commit.y;
		section.serialize(data, pos);
	}

	template<int NFiber, int EType>
	inline void ElementComponent<NFiber, EType>::deserialize(Vector& data, int& pos) {
		for (int i = 0; i < 6; ++i) UL_commit(i) = data(pos++);
		for (int i = 0; i < 2; ++i) qn[i] = Q2D::fromRotationVector(data(pos++));
		for (int i = 0; i < 2; ++i) qn_commit[i] = Q2D::fromRotationVector(data(pos++));
		rn.x = data(pos++);
		rn.y = data(pos++);
		rn_commit.x = data(pos++);
		rn_commit.y = data(pos++);
		section.deserialize(data, pos);
	}
	
	template<int NFiber, int EType>
	inline int rve_process_element(bool elastic_correction, ElementComponent<NFiber, EType>& ele, const RVEStateVariables& rve, const ASDSteel1DMaterial::InputParameters& params, bool do_implex, double time_factor) {
		int retval = ele.compute(elastic_correction, rve, params, do_implex, time_factor);
		if (retval != 0) {
			asd_print_full("ele !converged");
			return retval;
		}
		auto& globals = Globals::instance();
		if (elastic_correction && params.auto_regularization) {
			assembleRHS<EType>(params, elastic_correction, globals.element_RHS, globals.rve_R_el);
			assembleLHS<EType>(params, elastic_correction, globals.element_LHS, globals.rve_K_el);
		}
		else {

			assembleRHS<EType>(params, elastic_correction, globals.element_RHS, globals.rve_R);
			assembleLHS<EType>(params, elastic_correction, globals.element_LHS, globals.rve_K);
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
		ElementComponent<1, EType_El> e0; // nonlinear section (1 fiber), but linear kinematics
		int serializationDataSize() const;
		void serialize(Vector& data, int& pos);
		void deserialize(Vector& data, int& pos);
		RVEModel() = default;
		
		inline void setSlipMaterial(UniaxialMaterial* prototype) {
			e2.section.setSlipMaterial(prototype);
		}

		inline int compute(bool elastic_correction, const ASDSteel1DMaterial::InputParameters& params, double U, bool do_implex, double time_factor, double& N, double& tangent) {

			// output
			int retval = -1;
			
			// start from previous commited state (not in implex commit)
			if (!(params.implex && !do_implex)) {
				revertToLastCommit();
			}
			
			// globals
			double tolR = params.tolR * params.sy;
			double tolU = params.tolU * params.length;
			auto& globals = Globals::instance();
			
			// utility for assembly
			// note: do it in reverse order, so that we can access the element RHS for element Bottom to get reaction
			auto lam_assemble = [this, elastic_correction, &globals, &params, &do_implex, &time_factor]() -> int {
				// zero R and K for element integration
				globals.rve_R.Zero();
				globals.rve_K.Zero();
				globals.rve_R_el.Zero();
				globals.rve_K_el.Zero();
				if (elastic_correction && params.auto_regularization) {
					if (rve_process_element< 3, EType_Top>(elastic_correction, e3, sv, params, do_implex, time_factor) != 0) return -1;
					if (rve_process_element< 1, EType_Mid>(elastic_correction, e2, sv, params, do_implex, time_factor) != 0) return -1;
					if (rve_process_element< 3, EType_Bot>(elastic_correction, e1, sv, params, do_implex, time_factor) != 0) return -1;
					if (rve_process_element< 1, EType_El>(elastic_correction, e0, sv, params, do_implex, time_factor) != 0) return -1;
				}
				else {
					if (rve_process_element< 3, EType_Top>(elastic_correction, e3, sv, params, do_implex, time_factor) != 0) return -1;
					if (rve_process_element< 1, EType_Mid>(elastic_correction, e2, sv, params, do_implex, time_factor) != 0) return -1;
					if (rve_process_element< 3, EType_Bot>(elastic_correction, e1, sv, params, do_implex, time_factor) != 0) return -1;
					double I = M_PI * std::pow(params.radius, 4) / 4.0;
					double EI = params.E * I;
					double fact = 0.04 * params.length / (params.radius * 5.0);
					double penalty = 12.0 * EI * fact / std::pow(2.0*params.length, 3.0);
					double KReg;
					if (!params.auto_regularization) {
						KReg = 0.0;
					}
					else {
						KReg = std::max(0.0, (params.length / params.lch_element - 1)) * penalty;
					}
					globals.rve_K(6, 6) = globals.rve_K(6, 6) + KReg;
					globals.rve_R(6) = globals.rve_R(6)  + KReg * sv.UG(6);
				}                            
				return 0;
			};

			// impose BC
			if (elastic_correction && params.auto_regularization) {
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
				if (elastic_correction && params.auto_regularization) {
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
			if (elastic_correction && params.auto_regularization) {
				getRetainedComponents<EType_El>(params, elastic_correction, globals.element_LHS, globals.rve_Krr, globals.rve_Kcr_el, globals.rve_Krc_el);
				globals.rve_K_el.Solve(globals.rve_Kcr_el, globals.rve_Kinv_Kcr_el);
				Krc_dot_Kinv_Kcr = globals.rve_Krc_el ^ globals.rve_Kinv_Kcr_el;
			}
			else {

				getRetainedComponents<EType_Bot>(params, elastic_correction, globals.element_LHS, globals.rve_Krr, globals.rve_Kcr, globals.rve_Krc);
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
	int RVEModel::serializationDataSize() const {
		return sv.serializationDataSize()
			 + e1.serializationDataSize()
			 + e2.serializationDataSize()
			 + e3.serializationDataSize()
			 + e0.serializationDataSize();
	}
	void RVEModel::serialize(Vector& data, int& pos) {
		sv.serialize(data, pos);	
		e1.serialize(data, pos);
		e2.serialize(data, pos);
		e3.serialize(data, pos);
		e0.serialize(data, pos);
	}
	void RVEModel::deserialize(Vector& data, int& pos) {
		sv.deserialize(data, pos);
		e1.deserialize(data, pos);
		e2.deserialize(data, pos);
		e3.deserialize(data, pos);
		e0.deserialize(data, pos);
	}
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

}

/**
	PIMPL pattern
	*/
class ASDSteel1DMaterialPIMPL
{
public:
	ASDSteel1DMaterialPIMPL() = default;
	ASDSteel1DMaterialPIMPL(const ASDSteel1DMaterialPIMPL&) = default;
	inline void setSlipMaterial(UniaxialMaterial* prototype) {
		rve_m.setSlipMaterial(prototype);
	}
public:
	SeriesComponent steel_comp;
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
	static const char* msg = "uniaxialMaterial ASDSteel1D $tag $E $sy $su $eu  <-implex> <-implexControl $implexErrorTolerance $implexTimeReductionLimit> <-auto_regularization> <-buckling  $lch < $r>> <-fracture  <$r>> <-slip $matTag <$r>> <-K_alpha $K_alpha> <-max_iter $max_iter> <-tolU $tolU> <-tolR $tolR>";

	// check arguments
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 5) {
		opserr <<
			"uniaxialMaterial ASDSteel1D Error: Few arguments (< 5).\n" << msg << "\n";
		return nullptr;
	}

	// numData
	int numData = 1;

	// data
	int tag;
	double E = 0.0;
	double sy = 0.0;
	double su = 0.0;
	double eu = 0.0;
	double lch = 0.0;
	double r = 0.0; // default to 0.0, means not provided
	bool implex = false;
	bool implex_control = false;
	double implex_error_tolerance = 0.05;
	double implex_time_redution_limit = 0.01;
	bool auto_regularization = false;
	bool buckling = false;
	bool fracture = false;
	bool slip = false;
	double K_alpha = 0.5;
	double max_iter= 100;
	double tolU = 1.0e-6;
	double tolR = 1.0e-6;
	bool have_K_alpha = false;
	bool have_max_iter = false;
	bool have_tolU = false;
	bool have_tolR = false;
	UniaxialMaterial* matSlip = nullptr;


	// get tag
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "UniaxialMaterial ASDSteel1D Error: invalid 'tag'.\n";
		return nullptr;
	}

	// get steel base arguments
	auto lam_get_dparam = [&numData](double* val, const char* valname) -> bool {
		if (OPS_GetDouble(&numData, val) != 0) {
			opserr << "UniaxialMaterial ASDSteel1D Error: invalid '" << valname << "'.\n" << msg << "\n";
			return false;
		}
		if (*val <= 0.0) {
			opserr << "UniaxialMaterial ASDSteel1D Error: invalid value for '" << valname << "' (" << *val << "). It should be strictly positive.\n" << msg << "\n";
			return false;
		}
		return true;
	};
	if (!lam_get_dparam(&E, "E")) return nullptr;
	if (!lam_get_dparam(&sy, "sy")) return nullptr;
	if (!lam_get_dparam(&su, "su")) return nullptr;
	if (!lam_get_dparam(&eu, "eu")) return nullptr;
	auto lam_optional_double = [&numData](const char* variable, double& value) -> bool {
		if (OPS_GetNumRemainingInputArgs() > 0) {
			if (OPS_GetDouble(&numData, &value) < 0) {
				opserr << "UniaxialMaterial ASDSteel1D Error: failed to get '" << variable << "'.\n";
				return false;
			}
		}
		else {
			opserr << "UniaxialMaterial ASDSteel1D Error: '" << variable << "' requested but not provided.\n";
			return false;
		}
		return true;
	};
	
	// parse optional arguments
	int trials = 0;
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* value = OPS_GetString();
		if (strcmp(value, "-implex") == 0) {
			implex = true;
		}
		if (strcmp(value, "-implexControl") == 0) {
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
		if (strcmp(value, "-auto_regularization") == 0) {
			auto_regularization = true;
		}
		if (strcmp(value, "-buckling") == 0) {
			buckling = true;
			// lch is mandatory
			if (OPS_GetNumRemainingInputArgs() < 1) {
				opserr << "UniaxialMaterial ASDSteel1D: '-buckling' requires at least '$lch'\n";
				return nullptr;
			}
			if(!lam_optional_double("lch", lch)) return nullptr;
			// radius is optional (can be defined either here or in -slip)
			if (OPS_GetNumRemainingInputArgs() > 0) {
				double trial_radius;
				auto old_num_rem = OPS_GetNumRemainingInputArgs();
				if (OPS_GetDouble(&numData, &trial_radius) < 0) {
					// radius not provided, go back
					auto new_num_rem = OPS_GetNumRemainingInputArgs();
					if (new_num_rem < old_num_rem)
						OPS_ResetCurrentInputArg(-1);
				}
				else {
					// radius give, do checks
					if (trial_radius != r && r != 0.0) {
						opserr << "UniaxialMaterial ASDSteel1D: radius provied in -buckling (" << trial_radius << ") does not match the one provided in -slip or -fracture (" << r << "). Please use it only once.\n";
						return nullptr;
					}
					if (trial_radius <= 0.0) {
						opserr << "UniaxialMaterial ASDSteel1D: radius provied in -buckling should be strictly positive.\n";
						return nullptr;
					}
					r = trial_radius;
				}
			}
		}
		if (strcmp(value, "-fracture") == 0) {
			fracture = true;
			if (OPS_GetNumRemainingInputArgs() > 0) {
				double trial_radius;
				auto old_num_rem = OPS_GetNumRemainingInputArgs();
				if (OPS_GetDouble(&numData, &trial_radius) < 0) {
					// radius not provided, go back
					auto new_num_rem = OPS_GetNumRemainingInputArgs();
					if (new_num_rem < old_num_rem)
						OPS_ResetCurrentInputArg(-1);
				}
				else {
					// radius give, do checks
					if (trial_radius != r && r != 0.0) {
						opserr << "UniaxialMaterial ASDSteel1D: radius provied in -fracture (" << trial_radius << ") does not match the one provided in -slip or -buckling (" << r << "). Please use it only once.\n";
						return nullptr;
					}
					if (trial_radius <= 0.0) {
						opserr << "UniaxialMaterial ASDSteel1D: radius provied in -fracture should be strictly positive.\n";
						return nullptr;
					}
					r = trial_radius;
				}
			}
		}
		if (strcmp(value, "-slip") == 0) {
			slip = true;
			if (OPS_GetNumRemainingInputArgs() < 1) {
				opserr << "UniaxialMaterial ASDSteel1D: '-slip' requires '$matTag'.\n";
				return nullptr;
			}
			// get slip material
			int matTag;
			if (OPS_GetInt(&numData, &matTag) != 0) {
				opserr << "UniaxialMaterial ASDSteel1D: cannot get slip material tag\n";
				return nullptr;
			}
			// Retrieve the UniaxialMaterial*
			matSlip = OPS_getUniaxialMaterial(matTag);
			if (matSlip == nullptr) {
				opserr << "ASDSteel1D Error: No existing UniaxialMaterial with tag " << matTag << " for -slip option.\n";
				return nullptr;
			}
			
			// radius is optional (can be defined either here or in -buckling)
			if (OPS_GetNumRemainingInputArgs() > 0) {
				double trial_radius;
				auto old_num_rem = OPS_GetNumRemainingInputArgs();
				if (OPS_GetDouble(&numData, &trial_radius) < 0) {
					// radius not provided, go back
					auto new_num_rem = OPS_GetNumRemainingInputArgs();
					if (new_num_rem < old_num_rem)
						OPS_ResetCurrentInputArg(-1);
				}
				else {
					// radius give, do checks
					if (trial_radius != r && r != 0.0) {
						opserr << "UniaxialMaterial ASDSteel1D: radius provied in -slip (" << trial_radius << ") does not match the one provided in -buckling or -fracture (" << r << "). Please use it only once.\n";
						return nullptr;
					}
					if (trial_radius <= 0.0) {
						opserr << "UniaxialMaterial ASDSteel1D: radius provied in -slip should be strictly positive.\n";
						return nullptr;
					}
					r = trial_radius;
				}
			}
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
		opserr << "1DMaterial ASDSteel1D Error: invalid value for 'su' (" << su << "). It should be larger than sy  (" << sy << ").\n" << msg << "\n";
		return nullptr;
	}
	if (E == 0) {
		opserr << "1DMaterial ASDSteel1D Error: invalid value for 'E' (" << E << "). It should be non-zero.\n" << msg << "\n";
		return nullptr;
	}
	
	// obtain chaboche params from E, sy, su, eu
	// we want to use 2 hardening functions as per chaboche model.
	// so that the initial slope is close to E and the the stress apporaches su at eu
	//double dy = su - sy;
	//double H1 = E / 1000.0 / eu * dy / 40.0/1000.0;
	//double gamma1 = H1 / dy;
	//double H2 = H1 * 50;
	//double gamma2 = gamma1 * 50;
	//double alpha = 0.9;
	double sy_norm = sy / E;
	double su_norm = su / E;
	double dy_norm = su_norm - sy_norm;

	double n = 400.0;
	double m = 50.0;

	double H1_norm = (1.0 / n) / eu;
	double gamma1_norm = H1_norm / dy_norm;
	double H1 = H1_norm * E;
	double gamma1 = gamma1_norm ;

	double H2 = H1 * m;
	double gamma2 = gamma1 * m;
	double alpha = 0.9;
	ASDSteel1DMaterial::InputParameters params;
	params.E = E;
	params.sy = sy;
	params.eu = eu;
	params.H1 = H1 * alpha;
	params.gamma1 = gamma1;
	params.H2 = H2 * (1.0 - alpha);
	params.gamma2 = gamma2;
	params.implex = implex;
	params.implex_control = implex_control;
	params.implex_error_tolerance = implex_error_tolerance;
	params.implex_time_redution_limit = implex_time_redution_limit;
	params.auto_regularization = auto_regularization;
	params.buckling = buckling;
	params.fracture = fracture;
	params.slip = slip;
	params.length = buckling ? lch / 2.0 : 1.0; // consider half distance, the RVE uses symmetry
	params.radius = r;
	params.K_alpha = K_alpha;
	params.max_iter = max_iter;
	params.tolU = tolU;
	params.tolR = tolR;
	params.lch_element = params.length;
	// create the material
	UniaxialMaterial* instance = new ASDSteel1DMaterial(
		// tag
		tag, params, matSlip);
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
	const ASDSteel1DMaterial::InputParameters& _params,
	UniaxialMaterial* slip_material)
	: UniaxialMaterial(_tag, MAT_TAG_ASDSteel1DMaterial)
	, params(_params)
	, pdata(new ASDSteel1DMaterialPIMPL())
{
	// set the slip material if defined for the RVE...
	pdata->setSlipMaterial(slip_material);
	// ... and for the base material
	pdata->steel_comp.setSlipMaterial(slip_material);
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
	if (!dtime_is_user_defined) {
		dtime_n = ops_Dt;
		if (!commit_done) {
			dtime_0 = dtime_n;
			dtime_n_commit = dtime_n;
		}
	}
	// time factor for explicit extrapolation
	double time_factor = 1.0;
	if (params.implex && (dtime_n_commit > 0.0))
		time_factor = dtime_n / dtime_n_commit;
	
	// save macro strain
	strain = v;
	// homogenize micro response (stress/tangent)

	asd_print("MAT set trial (" << (int)params.implex << ")");
	int retval;
	
	double epl = 0.0;
	if (params.buckling) {  // se buckling farlo con l'energia ottenuta da rve (incremento spostamento * incremento forze (u - u commit)*(residuo-residuo commit) e normalizzarlo con una frazione della lunghezza* sigmay * Area
		//regularization for element length lower than physical length (weighted average)
		
		if (params.implex) {
			if (params.implex_control) {
				double area = M_PI * params.radius * params.radius;
				bool elastic_correction = params.lch_element > params.length;
				// implicit solution
				double retval = homogenize(false);             
				if (retval != 0) return retval;

				Vector U_impl = elastic_correction && params.auto_regularization ? pdata->rve_m.sv.UG_el : pdata->rve_m.sv.UG;
				double N_impl = N_rve_last;  
				double sigma_impl = N_impl / stress;

				pdata->rve_m.revertToLastCommit();

				// explicit solution
				retval = homogenize(true);
				if (retval != 0) return retval;

				Vector  U_expl = elastic_correction && params.auto_regularization ? pdata->rve_m.sv.UG_el : pdata->rve_m.sv.UG;
				double N_expl = N_rve_last;
				double sigma_expl = N_expl / stress;
				
				Vector dU_err = U_expl - U_impl;

				// L2 Norm 
				double norm_dU = dU_err.Norm();
				double implex_error_u = norm_dU / params.radius;
				double implex_error_N = (sigma_expl - sigma_impl) / params.sy;
				implex_error = std::max(implex_error_u, implex_error_N);
				if (implex_error > params.implex_error_tolerance) {
					if (dtime_n >= params.implex_time_redution_limit * dtime_0) {
						retval = EC_IMPLEX_Error_Control;
					}
				}				
			}
			else {
				retval = homogenize(true);
			}
		}
		else {
			retval = homogenize(false);
		}
		C = C_rve;
		stress = stress_rve;
		epl = pdata->rve_m.e2.section.series.steel_material.epl;


	}
	else {
		//computation of sigma and tangent of steel
		pdata->steel_comp.revertToLastCommit();
		
		if (params.implex) {
			if (params.implex_control) {
				double sigma_macro = 0.0;
				double tangent_macro = 0.0;
				// implicit solution
				int retval = pdata->steel_comp.compute(params, false, 1.0, strain, sigma_macro, tangent_macro);
				double sigma_impl = sigma_macro;
				pdata->steel_comp.revertToLastCommit();
				// explicit solution
				retval = pdata->steel_comp.compute(params, true, time_factor, strain, stress, C);
				double sigma_expl = stress;
				
				implex_error = std::abs(sigma_expl - sigma_impl) / params.sy; //semplifica
				if (implex_error > params.implex_error_tolerance) {
					if (dtime_n >= params.implex_time_redution_limit * dtime_0) {
						retval = EC_IMPLEX_Error_Control;
					}
				}
			}
			else {
				retval = pdata->steel_comp.compute(params, true, time_factor, strain, stress, C);
			}
		}

		else {
			retval = pdata->steel_comp.compute(params, params.implex, time_factor, strain, stress, C);
		}
		epl = pdata->steel_comp.steel_material.epl;
	}
	if (params.fracture) {
		double d = 0.0;
		double eupl = params.eu - params.sy / params.E;
		if (epl > eupl) {
			double epl_max;
			if (!params.auto_regularization) {
				epl_max = eupl + eupl;
			}
			else {
				epl_max = eupl + eupl * (16.0 * params.radius) / (2.0 * params.lch_element);
			}
			double G = (epl_max-eupl) * params.sy/4.0;
			double sigma_damaged = std::max(1.0e-4 * params.sy, params.sy * exp(-params.sy * (epl-eupl) / G));
			d = 1.0 - sigma_damaged / params.sy;
		}
		//apply damage on stress e C
		stress *= (1.0 - d);
		C *= (1.0 - d);  //tangent for implex, secant for implicit
	}
	
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
		if (params.buckling) {
			asd_print("MAT commit (" << (int)false << ")");
			double area = M_PI * params.radius * params.radius;
			bool elastic_correction = params.lch_element > params.length;
			Vector  U_expl = elastic_correction && params.auto_regularization ? pdata->rve_m.sv.UG_el : pdata->rve_m.sv.UG;
			double N_expl = N_rve_last;
			double sigma_expl = N_expl / stress;
			
			// implicit solution
			double retval = homogenize(false);

			Vector U_impl = elastic_correction && params.auto_regularization ? pdata->rve_m.sv.UG_el : pdata->rve_m.sv.UG;
			double N_impl = N_rve_last;
			double sigma_impl = N_impl / stress;

			Vector dU_err = U_expl - U_impl;

			// Norm 
			double norm_dU = dU_err.Norm();
			double implex_error_u = norm_dU / params.radius;
			double implex_error_N = (sigma_expl - sigma_impl) / params.sy;
			implex_error = std::max(implex_error_u, implex_error_N);
			GlobalParameters::instance().setMaxError(std::max(implex_error, GlobalParameters::instance().getMaxError()));
			GlobalParameters::instance().accumulateAverageError(implex_error);
			if (retval < 0) return retval;
		}
		else {
			double sigma_expl = stress;
			double sigma_macro = 0.0;
			double tangent_macro = 0.0;
			int retval = pdata->steel_comp.compute(params, false, 1.0, strain, sigma_macro, tangent_macro);
			double sigma_impl = sigma_macro;
			implex_error = std::abs(sigma_expl - sigma_impl) / params.sy;
			GlobalParameters::instance().setMaxError(std::max(implex_error, GlobalParameters::instance().getMaxError()));
			GlobalParameters::instance().accumulateAverageError(implex_error);
			if (retval < 0) return retval;
		}

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
	dtime_0 = 0.0;
	dtime_is_user_defined = false;
	// IMPL-EX error
	implex_error = 0.0;
	// Commit flag
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

void ASDSteel1DMaterial::Print(OPS_Stream& s, int  elastic_correction)
{
	s << "ASDSteel1D Material, tag: " << this->getTag() << "\n";
}

int ASDSteel1DMaterial::sendSelf(int commitTag, Channel &theChannel)
{
	//// aux
	int counter;

	//// send DBL data
	counter = 0;

	// variable DBL data size

	int nv_dbl = 13 + 
		params.NDATA +
		pdata->rve_m.serializationDataSize() +  
		pdata->steel_comp.serializationDataSize() +
		pdata->rve_m.e2.section.series.serializationDataSize();



	Vector ddata(nv_dbl);

	// this data
	ddata(counter++) = static_cast<double>(getTag());
	ddata(counter++) = dtime_n;
	ddata(counter++) = dtime_n_commit;
	ddata(counter++) = dtime_0;
	ddata(counter++) = static_cast<int>(dtime_is_user_defined);
	ddata(counter++) = implex_error;
	ddata(counter++) = static_cast<double>(commit_done);
	ddata(counter++) = strain;
	ddata(counter++) = strain_commit;
	ddata(counter++) = stress;
	ddata(counter++) = stress_commit;
	ddata(counter++) = C;
	ddata(counter++) = stress_rve;
	ddata(counter++) = stress_rve_commit;
	ddata(counter++) = C_rve;
	ddata(counter++) = energy;
	// params
	ddata(counter++) = params.E;
	ddata(counter++) = params.sy;
	ddata(counter++) = params.eu;
	ddata(counter++) = params.H1;
	ddata(counter++) = params.H2;
	ddata(counter++) = params.gamma1;
	ddata(counter++) = params.gamma2;
	ddata(counter++) = static_cast<double>(params.implex);
	ddata(counter++) = static_cast<int>(params.implex_control);
	ddata(counter++) = params.implex_error_tolerance;
	ddata(counter++) = params.implex_time_redution_limit;
	ddata(counter++) = static_cast<double>(params.auto_regularization);
	ddata(counter++) = static_cast<double>(params.buckling);
	ddata(counter++) = static_cast<double>(params.fracture);
	ddata(counter++) = static_cast<double>(params.slip);
	ddata(counter++) = params.radius;
	ddata(counter++) = params.length;
	ddata(counter++) = params.lch_element;
	ddata(counter++) = params.K_alpha;
	ddata(counter++) = params.max_iter;
	ddata(counter++) = params.tolU;
	ddata(counter++) = params.tolR;
	// rve
	pdata->rve_m.serialize(ddata, counter);
	// steel base (no buckling)
	pdata->steel_comp.serialize(ddata, counter, commitTag, theChannel);
	// slip 
	pdata->rve_m.e2.section.series.serialize(ddata, counter, commitTag, theChannel);
	
	
	
	if (theChannel.sendVector(getDbTag(), commitTag, ddata) < 0) {
		opserr << "ASDSteel1DMaterial::sendSelf() - failed to send DBL data\n";
		return -1;
	}

	pdata->steel_comp.serialize_slip(commitTag, theChannel);
	pdata->rve_m.e2.section.series.serialize_slip(commitTag, theChannel);


	// done
	return 0;
}

int ASDSteel1DMaterial::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	// aux
	int counter;

	// variable DBL data size
	int nv_dbl = 13 +
		params.NDATA +
		pdata->rve_m.serializationDataSize() +
		pdata->steel_comp.serializationDataSize() +
		pdata->rve_m.e2.section.series.serializationDataSize();

	//	pdata->rve_m.e0.section.series.serializationDataSize() + pdata->rve_m.e2.section.series.serializationDataSize() + pdata->steel_comp.serializationDataSize();
	Vector ddata(nv_dbl);

	// recv DBL data
	//Vector ddata(InputParameters::NDATA + 1 * StateVariablesSteel::NDATA + 10);
	if (theChannel.recvVector(getDbTag(), commitTag, ddata) < 0) {
		opserr << "ASDSteel1DMaterial::recvSelf() - failed to receive DBL data\n";
		return -1;
	}
	
	counter = 0;
	
	// this data
	setTag(ddata(counter++));
	dtime_n = ddata(counter++);
	dtime_n_commit = ddata(counter++);
	dtime_0 = ddata(counter++);
	dtime_is_user_defined = static_cast<bool>(ddata(counter++));
	implex_error = ddata(counter++);
	commit_done = static_cast<bool>(ddata(counter++));
	strain = ddata(counter++);
	strain_commit = ddata(counter++);
	stress = ddata(counter++);
	stress_commit = ddata(counter++);
	C = ddata(counter++);
	stress_rve_commit = ddata(counter++);
	stress_rve = ddata(counter++);
	C_rve = ddata(counter++);
	energy = ddata(counter++);

	//params
	params.E = ddata(counter++);
	params.sy = ddata(counter++);
	params.eu = ddata(counter++);
	params.H1 = ddata(counter++);
	params.H2 = ddata(counter++);
	params.gamma1 = ddata(counter++);
	params.gamma2 = ddata(counter++);
	params.implex = static_cast<bool>(ddata(counter++));
	params.implex_control = static_cast<bool>(ddata(counter++));
	params.implex_error_tolerance = ddata(counter++);
	params.implex_time_redution_limit = ddata(counter++);
	params.auto_regularization = static_cast<bool>(ddata(counter++));
	params.buckling = static_cast<bool>(ddata(counter++));
	params.fracture = static_cast<bool>(ddata(counter++));
	params.slip = static_cast<bool>(ddata(counter++));
	params.radius = ddata(counter++);
	params.length = ddata(counter++);
	params.lch_element = ddata(counter++);
	params.K_alpha = ddata(counter++);
	params.max_iter = ddata(counter++);
	params.tolU = ddata(counter++);
	params.tolR = ddata(counter++);
	// rve
	pdata->rve_m.deserialize(ddata, counter);
	// steel base (in case of no buckling)
	pdata->steel_comp.deserialize(ddata, counter, commitTag, theChannel, theBroker);
	pdata->steel_comp.deserialize_slip(commitTag, theChannel, theBroker);
	// slip 
	pdata->rve_m.e2.section.series.deserialize(ddata, counter, commitTag, theChannel, theBroker);
	pdata->rve_m.e2.section.series.deserialize_slip(commitTag, theChannel, theBroker);

	// done
	return 0;
}

int ASDSteel1DMaterial::setParameter(const char** argv, int argc, Parameter& param)
{

	// 1000 - time
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

	// default
	return -1;
}

int ASDSteel1DMaterial::updateParameter(int parameterID, Information& info)
{
	switch (parameterID) {

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
	static std::vector<std::string> lb_damage = { "D" };
	static std::vector<std::string> lb_eqpl_strain = { "PLE" };
	static std::vector<std::string> lb_slip_resp = { "Slip", "Tau" };
	static std::vector<std::string> lb_time = { "dTime", "dTimeCommit", "dTimeInitial" };
	static std::vector<std::string> lb_implex_error = { "Error" };



	// check specific responses
	if (argc > 0) {
		if (strcmp(argv[0], "BI") == 0 || strcmp(argv[0], "BucklingIndicator") == 0) {
			return make_resp(1001, getBucklingIndicator(), &lb_buckling_ratio);
		}
		if (strcmp(argv[0], "D") == 0 || strcmp(argv[0], "Damage") == 0) {
			return make_resp(1002, getDamage(), &lb_damage);
		}
		if (strcmp(argv[0], "PLE") == 0 || strcmp(argv[0], "EquivalentPlasticStrain") == 0) {
			return make_resp(1003, getEqPlStrain(), &lb_eqpl_strain);
		}
		if (strcmp(argv[0], "SlipResponse") == 0){
			return make_resp(1004, getSlipResponse(), &lb_slip_resp);
		}
		if (strcmp(argv[0], "time") == 0 || strcmp(argv[0], "Time") == 0) {
			return make_resp(1005, getTimeIncrements(), &lb_time);
		}
		if (strcmp(argv[0], "implexError") == 0 || strcmp(argv[0], "ImplexError") == 0) {
			return make_resp(1006, getImplexError(), &lb_implex_error);
		}
	}

	// otherwise return base-class response
	return UniaxialMaterial::setResponse(argv, argc, output);
}

int ASDSteel1DMaterial::getResponse(int responseID, Information& matInformation)
{
	// all outputs are 1D
	//static Vector out1(1);

	switch (responseID) {
		// 1000 - base steel output
	case 1001:		
		return matInformation.setVector(getBucklingIndicator());
		// 1002 - damage
	case 1002:
		return matInformation.setVector(getDamage());
	case 1003:
		return matInformation.setVector(getEqPlStrain());
	case 1004:
		return matInformation.setVector(getSlipResponse());
		// 1005 - internal time
	case 1005: return matInformation.setVector(getTimeIncrements());
	case 1006: return matInformation.setVector(getImplexError());
	default:
		break;
	}
	return UniaxialMaterial::getResponse(responseID, matInformation);
}

double ASDSteel1DMaterial::getEnergy(void)
{
	return energy;
}

const Vector& ASDSteel1DMaterial::getBucklingIndicator() const
{
	static Vector d(1);
	d.Zero();
	d(0) = pdata->rve_m.sv.UG(6);
	return d;
}

const Vector& ASDSteel1DMaterial::getDamage() const
{
	static Vector d(1);
	d.Zero();
	if (params.fracture) {
		double eupl = params.eu - params.sy / params.E;
		double epl = params.buckling ? pdata->rve_m.e2.section.series.steel_material.epl : pdata->steel_comp.steel_material.epl;
		if (epl > eupl) {
			double epl_max;
			if (!params.auto_regularization) {
				epl_max = eupl + eupl;
			}
			else {
				epl_max = eupl + eupl * (16.0 * params.radius) / (2.0 * params.lch_element);
			}
			double G = (epl_max - eupl) * params.sy / 4.0;
			double sigma_damaged = std::max(1.0e-4 * params.sy, params.sy * exp(-params.sy * (epl - eupl) / G));
			d(0) = 1.0 - sigma_damaged / params.sy;
		}
	}
	return d;
}

const Vector& ASDSteel1DMaterial::getEqPlStrain() const
{
	static Vector d(1);
	d.Zero();
	d(0) = params.buckling ? pdata->rve_m.e2.section.series.steel_material.epl : pdata->steel_comp.steel_material.epl;
	return d;
}

const Vector& ASDSteel1DMaterial::getSlipResponse() const
{
	static Vector d(2);
	d.Zero();
	if (params.slip) {
		if (params.buckling) {
			d(0) = pdata->rve_m.e2.section.series.slip_material->getStrain();
			d(1) = pdata->rve_m.e2.section.series.slip_material->getStress();
		}
		else {
			d(0) = pdata->steel_comp.slip_material->getStrain();
			d(1) = pdata->steel_comp.slip_material->getStress();
		}
	}
	return d;
}

const Vector& ASDSteel1DMaterial::getTimeIncrements() const
{
	static Vector d(3);
	d(0) = dtime_n;
	d(1) = dtime_n_commit;
	d(2) = dtime_0;
	return d;
}

const Vector& ASDSteel1DMaterial::getImplexError() const
{
	static Vector d(1);
	d(0) = implex_error;
	return d;
}

int ASDSteel1DMaterial::homogenize(bool do_implex)
{	
	// return value
	int retval = 0;

	// check for elastic correction
	bool elastic_correction = params.lch_element > params.length;

	auto& globals = Globals::instance();
	if (elastic_correction && params.auto_regularization) {
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
	if (elastic_correction && params.auto_regularization) {
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
	N_rve_last = N;
	// from micro stress/tangent to macro stress/tangent
	area = M_PI * params.radius * params.radius;

	stress_rve = -N/area;
	if (elastic_correction && params.auto_regularization) {
		C_rve = T * params.lch_element / area;
	}
	else {
		C_rve = T * params.length / area;
	}
	

	// done
	return retval;
}




