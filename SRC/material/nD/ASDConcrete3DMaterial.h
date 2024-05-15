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

// Massimo Petracca - ASDEA Software, Italy
//
// A Generic Orthotropic Material Wrapper that can convert any
// nonlinear isotropic material into an orthotropic one by means of tensor
// mapping
//

#ifndef ASDConcrete3DMaterial_h
#define ASDConcrete3DMaterial_h

#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include <cmath>
#include <memory>
#include <vector>
#include <map>

class ASDConcrete3DMaterial : public NDMaterial 
{
public:
	// sub-classes

	// A simple 3D vector
	struct Vector3 {
		double x = 0.0;
		double y = 0.0;
		double z = 0.0;
		Vector3() = default;
		Vector3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
		inline double norm()const { return std::sqrt(x * x + y * y + z * z); }
		inline void normalize() {
			double n = norm();
			if (n > 0.0) {
				x /= n;
				y /= n;
				z /= n;
			}
		}
		inline double dot(const Vector3& b)const { return x * b.x + y * b.y + z * b.z; }
		inline Vector3 cross(const Vector3& b)const {
			return Vector3(b.z * y - b.y * z, b.x * z - b.z * x, b.y * x - b.x * y);
		}
		inline Vector3 operator-()const { return Vector3(-x, -y, -z); }
		inline Vector3 operator*(double s)const { return Vector3(x * s, y * s, z * s); }
		inline Vector3 operator+(const Vector3& b)const { return Vector3(x + b.x, y + b.y, z + b.z); }
		inline Vector3 operator-(const Vector3& b)const { return Vector3(x - b.x, y - b.y, z - b.z); }
	};

	// A class for performing stress decomposition
	struct StressDecomposition {
		StressDecomposition() = default;
		int compute(const Vector& S, double cdf);
		void recompose(const Vector& S, Vector& Sv) const;
		void recompose(Vector& Sv) const;
		Vector Si = Vector(3); // principal stressess 0>1>2
		Matrix V = Matrix(3, 3); // principal directions in columns
		Matrix PT = Matrix(6, 6); // tensile projector
		Matrix PC = Matrix(6, 6); // compressive projector
		Vector ST = Vector(6); // tensile stress
		Vector SC = Vector(6); // compressive stress
		double R = 0.0;
	};

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

	// The Crack Planes Storage. A static storage for crack normal vectors
	class CrackPlanesStorage {
	public:
		using Vector3List = std::vector<Vector3>;
		using Vector3ListPointer = std::shared_ptr<Vector3List>;
		using MapType = std::map<int, Vector3ListPointer>;

	private:
		CrackPlanesStorage() = default;
		CrackPlanesStorage(const CrackPlanesStorage&) = delete;
		CrackPlanesStorage& operator = (const CrackPlanesStorage&) = delete;

	public:
		// Access to this singleton
		static CrackPlanesStorage& instance();
		// Gets a pointer to a vector of normals with n90 subdivisions every pi/2.
		// It will return a pointer to an existing vector. If no normals have been generated
		// before for the current n90, it will be generated here.
		Vector3ListPointer get(int n90);

	private:
		MapType m_map;
	};

	// The Crack Planes. The instance that holds a pointer to a list of normal vectors
	class CrackPlanes {
	public:
		CrackPlanes() = default;
		CrackPlanes(int n90);
	public:
		// Sets the current normal. Finds the closest normal location in m_normals
		void setCurrentNormal(const Vector3& ni);
		// Returns the equivalent strain at the closest normal location
		double getCurrentEquivalentStrain()const;
		// Updates the equivalent strain at the closest normal location. The smoothing angle
		// should be provided in radians
		void updateCurrentEquivalentStrain(double x, double smooth_angle);
		// reset to zero
		void reset();
		// get count
		inline std::size_t count()const { return m_equivalent_strain.size(); }
		// get equivalent strain
		double getEquivalentStrainAtNormal(std::size_t i) const;
		// set equivalent strain
		void setEquivalentStrainAtNormal(std::size_t i, double x);
		// get normal
		const Vector3& getNormal(std::size_t i) const;
		// closest normal
		std::size_t getClosestNormal(const Vector3& N) const;
		// get at most 3 normals with maximum values
		std::vector<int> getMax3Normals(double smooth_angle) const;
		// serialization
		int serializationDataSize() const;
		void serialize(Vector& data, int& pos);
		void deserialize(Vector& data, int& pos);
	private:
		int m_n90 = 0;
		CrackPlanesStorage::Vector3ListPointer m_normals;
		std::vector<double> m_equivalent_strain = { 0.0 };
		Vector3 m_current_normal;
		std::size_t m_closest_normal_loc = 0;
	};

public:
	// life-cycle
	ASDConcrete3DMaterial(
		int _tag,
		double _E,
		double _v,
		double _rho,
		double _eta,
		double _Kc,
		bool _implex,
		bool _implex_control,
		double _implex_error_tolerance,
		double _implex_time_reduction_limit,
		double _implex_alpha,
		bool _tangent,
		bool _auto_regularize,
		double _lch_ref,
		const HardeningLaw& _ht,
		const HardeningLaw& _hc,
		double _cdf,
		int _nct,
		int _ncc,
		double _smoothing_angle);
	ASDConcrete3DMaterial();
	~ASDConcrete3DMaterial();

	// info
	const char* getClassType(void) const { return "ASDConcrete3DMaterial"; };

	// density
	double getRho(void);

	// set strain
	int setTrialStrain(const Vector &v);
	int setTrialStrain(const Vector& v, const Vector& r);
	int setTrialStrainIncr(const Vector& v);
	int setTrialStrainIncr(const Vector& v, const Vector& r);

	// get state
	const Vector& getStrain(void);
	const Vector& getStress(void);
	const Matrix& getTangent(void);
	const Matrix& getInitialTangent(void);

	// handle state
	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	// copy and others...
	NDMaterial* getCopy(void);
	NDMaterial* getCopy(const char* code);
	const char* getType(void) const;
	int getOrder(void) const;
	void Print(OPS_Stream &s, int flag=0);

	// send/recv self
	virtual int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

	// parameters and responses
	int setParameter(const char** argv, int argc, Parameter& param);
	int updateParameter(int parameterID, Information& info);
	Response* setResponse(const char** argv, int argc, OPS_Stream& output);
	int getResponse(int responseID, Information& matInformation);

private:
	// internal computation
	int compute(bool do_implex, bool do_tangent);
	double lublinerCriterion(double s1, double s2, double s3, double ft, double fc, double k1, double scale) const;
	double equivalentTensileStrainMeasure(double s1, double s2, double s3) const;
	double equivalentCompressiveStrainMeasure(double s1, double s2, double s3) const;
	Vector getHardeningLawVector(HardeningLawType ltype, HardeningLawPointComponent c) const;
	const Vector& getMaxStrainMeasure() const;
	const Vector& getAvgStrainMeasure() const;
	const Vector& getMaxDamage() const;
	const Vector& getAvgDamage() const;
	const Vector& getMaxEquivalentPlasticStrain() const;
	const Vector& getAvgEquivalentPlasticStrain() const;
	const Vector& getMaxCrackWidth() const;
	const Vector& getAvgCrackWidth() const;
	const Vector& getMaxCrushWidth() const;
	const Vector& getAvgCrushWidth() const;
	const Vector& getCrackPattern() const;
	const Vector& getCrushPattern() const;
	const Vector& getImplexError() const;
	const Vector& getTimeIncrements() const;

private:
	// Young's modulus
	double E = 0.0;
	// Poisson's ratio
	double v = 0.0;
	// Mass density
	double rho = 0.0;
	// Viscosity for rate-dependent damage
	double eta = 0.0;
	// Kc parameter for triaxial compression (default suggested by Lubliner et al.)
	double Kc = 2.0 / 3.0;
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
	// Compression/Tension peak ratio
	double fcft_ratio = 10;
	// cross-damage factor
	double cdf = 0.0;
	// number of normals and smoothing angle ( 0 means isotropic internal variables )
	int nct = 0;
	int ncc = 0;
	double smoothing_angle = 0.7854; // 45 deg
	// state variables - tension
	CrackPlanes svt;
	CrackPlanes svt_commit;
	CrackPlanes svt_commit_old;
	// state variables - compression
	CrackPlanes svc;
	CrackPlanes svc_commit;
	CrackPlanes svc_commit_old;
	// state variables - implex
	double dtime_n = 0.0;
	double dtime_n_commit = 0.0;
	double dtime_0 = 0.0;
	bool dtime_is_user_defined = false;
	bool commit_done = false;
	double implex_error = 0.0;
	Matrix PT_commit = Matrix(6, 6);
	double R_commit = 0.0;
	// strain, stress and tangent
	Vector strain = Vector(6);
	Vector strain_commit = Vector(6);
	Vector stress = Vector(6);
	Vector stress_eff = Vector(6);
	Vector stress_eff_commit = Vector(6);
	Matrix C = Matrix(6, 6);
	// other variables for output purposes
	double dt_bar = 0.0;
	double dc_bar = 0.0;
	Vector iso_crack_normal = Vector(3);
	Vector iso_crush_normal = Vector(3);
	double xt_max = 0.0;
	double xt_max_commit = 0.0;
	double xc_max = 0.0;
	double xc_max_commit = 0.0;
};

#endif
