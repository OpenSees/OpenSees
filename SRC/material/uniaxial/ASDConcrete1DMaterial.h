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

#ifndef ASDConcrete1DMaterial_h
#define ASDConcrete1DMaterial_h

#include <UniaxialMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include <cmath>
#include <memory>
#include <vector>
#include <map>

class ASDConcrete1DMaterial : public UniaxialMaterial
{
public:
	// sub-classes

	// A point in the hardening law
	struct HardeningLawPoint {
		double x = 0.0; // total strain
		double y = 0.0; // final backbone stress (plasticity & damage)
		double d = 0.0; // damage variable
		double q = 0.0; // backbone stress in the effective-stress space (plasticity only)
		HardeningLawPoint() = default;
		HardeningLawPoint(double _x, double _y, double _d, double _q)
			: x(_x), y(_y), d(_d), q(_q) {}
		inline double totalStrain()const { return x; }
		inline double plasticStrain(double E)const { return x - q / E; }
		inline double stress()const { return y; }
		inline double effectiveStress()const { return q; }
		inline double elasticStress(double E) const { return E * x; }
		inline double crackingDamage()const { return d; }
	};

	// The Hardening Law Type
	enum class HardeningLawType {
		Tension = 0,
		Compression
	};

	// The Hardening Law Point Component (for output)
	enum class HardeningLawPointComponent {
		TotalStrain = 0,
		EffectiveStress,
		NominalStress
	};

	// The Hardening law
	class HardeningLaw {
	public:
		// default constructor
		HardeningLaw() = default;
		// full constructor
		HardeningLaw(
			int tag, HardeningLawType type,
			double E,
			const std::vector<double>& x,
			const std::vector<double>& y,
			const std::vector<double>& d);

		// regularizes the hardening curve according to the 'lch'
		// characteristic length of the parent element, and the 'lch_ref'
		// characteristic length of the input hardening curve
		void regularize(double lch, double lch_ref);
		// nullify previous regularization
		void deRegularize();
		// evaluate the hardening law at a certain strain
		HardeningLawPoint evaluateAt(double x) const;
		// get max stress value
		double computeMaxStress() const;
		// serialization
		int serializationDataSize() const;
		void serialize(Vector& data, int& pos);
		void deserialize(Vector& data, int& pos);
		// properties
		inline bool isValid() const { return m_valid; }
		inline int tag()const { return m_tag; }
		inline HardeningLawType type()const { return m_type; }
		inline const std::vector<HardeningLawPoint>& points()const { return m_points; }
		inline double strainTolerance()const { return m_xtolerance; }
		inline double stressTolerance()const { return m_ytolerance; }
		inline bool hasStrainSoftening()const { return m_fracture_energy_is_bounded; }
		inline double strainAtOnsetOfCrack()const {
			if (m_fracture_energy_is_bounded && m_softening_begin < m_points.size())
				return m_points[m_softening_begin].totalStrain();
			return 0.0;
		}

	private:
		// adjusts the points
		void adjust();
		// computes the fracture energy of the hardening curve
		void computeFractureEnergy();

	private:
		// tag (same as the parent material's tag)
		int m_tag = 0;
		// type
		HardeningLawType m_type = HardeningLawType::Tension;
		// the hardening points of the backbone curve in total-strain
		std::vector<HardeningLawPoint> m_points;
		// the fracture energy (computed only if there is softening)
		double m_fracture_energy = 0.0;
		// true if there is softening, false otherwise
		bool m_fracture_energy_is_bounded = false;
		// the location of the point in m_points where the softening begins
		std::size_t m_softening_begin = 0;
		// the location of the point in m_points where the softening ends
		std::size_t m_softening_end = 0;
		// false if the initialization (in constructor) fails, true otherwise
		bool m_valid = false;
		// tolerances
		double m_xtolerance = 1.0e-12;
		double m_ytolerance = 1.0e-12;
	};

	// The Hardening Law Storage. A static storage for original hardening laws
	class HardeningLawStorage {
	public:
		using PointerType = std::shared_ptr<HardeningLaw>;
		using MapType = std::map<int, PointerType>;

	public:
		HardeningLawStorage() = default;
		HardeningLawStorage(const HardeningLawStorage&) = delete;
		HardeningLawStorage& operator = (const HardeningLawStorage&) = delete;

	public:
		// Access to this singleton
		static HardeningLawStorage& instance();
		// store
		void store(const HardeningLaw& hl);
		// recover
		PointerType recover(int tag, HardeningLawType type);

	private:
		MapType m_tension;
		MapType m_compression;
	};

public:
	// life-cycle
	ASDConcrete1DMaterial(
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
		const HardeningLaw& _hc);
	ASDConcrete1DMaterial();
	~ASDConcrete1DMaterial();

	// info
	const char* getClassType(void) const { return "ASDConcrete1DMaterial"; }

	// set strain
	int setTrialStrain(double v, double r = 0.0);

	// get state
	double getStrain(void);
	double getStress(void);
	double getTangent(void);
	double getInitialTangent(void);

	// handle state
	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	// copy and others...
	UniaxialMaterial* getCopy(void);
	void Print(OPS_Stream& s, int flag = 0);

	// send/recv self
	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

	// parameters and responses
	int setParameter(const char** argv, int argc, Parameter& param);
	int updateParameter(int parameterID, Information& info);
	Response* setResponse(const char** argv, int argc, OPS_Stream& output);
	int getResponse(int responseID, Information& matInformation);
	double getEnergy(void);

private:
	// internal computation
	int compute(bool do_implex, bool do_tangent);
	Vector getHardeningLawVector(HardeningLawType ltype, HardeningLawPointComponent c) const;
	const Vector& getStrainMeasure() const;
	const Vector& getDamage() const;
	const Vector& getEquivalentPlasticStrain() const;
	const Vector& getCrackWidth() const;
	const Vector& getCrushWidth() const;
	const Vector& getImplexError() const;
	const Vector& getTimeIncrements() const;

 private:
	 // Young's modulus
	 double E = 0.0;
	 // Viscosity for rate-dependent damage
	 double eta = 0.0;
	 // True = use the IMPL-EX algorithm
	 bool implex = false;
	 // True = keep IMPL-EX error under control
	 bool implex_control = false;
	 // Maximum allowed IMPL-EX error (default = 5%)
	 double implex_error_tolerance = 0.05;
	 // Minimum allowed time step reduction factor under which IMPL-EX error is not controlled anymore (default = 1%)
	 double implex_time_redution_limit = 0.01;
	 // Scale factor for implex extrapolation
	 double implex_alpha = 1.0;
	 // True = use the tangent matrix, False (default) = use the secant matrix
	 bool tangent = false;
	 // True = automatically regularize the fracture energy using the element's characteristic length
	 bool auto_regularize = true;
	 bool regularization_done = false;
	 double lch = 1.0; // the parent-element's characteristic length
	 double lch_ref = 1.0; // the reference characteristic length (i.e. the size the specific-fracture-energy in the hardening-law is referred to)
	 // The hardening law for the tensile response
	 HardeningLaw ht;
	 // The hardening law for the compressive response
	 HardeningLaw hc;
	 // state variables - tension
	 double xt = 0.0;
	 double xt_commit = 0.0;
	 double xt_commit_old = 0.0;
	 // state variables - compression
	 double xc = 0.0;
	 double xc_commit = 0.0;
	 double xc_commit_old = 0.0;
	 // state variables - implex
	 double dtime_n = 0.0;
	 double dtime_n_commit = 0.0;
	 double dtime_0 = 0.0;
	 bool dtime_is_user_defined = false;
	 bool commit_done = false;
	 double implex_error = 0.0;
	 double PT_commit = 0.5;
	 // strain, stress and tangent
	 double strain = 0.0;
	 double strain_commit = 0.0;
	 double stress = 0.0;
	 double stress_commit = 0.0;
	 double stress_eff = 0.0;
	 double stress_eff_commit = 0.0;
	 double C = 0.0;
	 // other variables for output purposes
	 double dt_bar = 0.0;
	 double dc_bar = 0.0;
	 double energy = 0.0;
};

#endif
