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

// Massimo Petracca - ASDEA Software, Italy (2022)
//
// Series3DMaterial is an aggregation
// of NDMaterials (3D only) objects all considered acting in Series.
//

#include <Series3DMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <DummyStream.h>
#include <MaterialResponse.h>
#include <ID.h>
#include <cmath>
#include <string>

// interface to LAPACK
#ifdef _WIN32

extern "C" int  DGETRF(int* M, int* N, double* A, int* LDA,
	int* iPiv, int* INFO);

extern "C" int  DGETRS(char* TRANS,
	int* N, int* NRHS, double* A, int* LDA,
	int* iPiv, double* B, int* LDB, int* INFO);

#else

extern "C" int dgetrf_(int* M, int* N, double* A, int* LDA,
	int* iPiv, int* INFO);

extern "C" int dgetrs_(char* TRANS, int* N, int* NRHS, double* A, int* LDA,
	int* iPiv, double* B, int* LDB, int* INFO);

#endif

// namespace for utilities
namespace Series3DUtils {

	// Solver wrapper to factorize once, and solve multiple times
	class SolverWrapper {
	private:
		std::vector<double> A;
		std::vector<int> IPIV;
	public:
		bool factorize(const Matrix& M) {
			int n = M.noRows();
			if (n == 0 || n != M.noCols())
				return false;
			A.resize(static_cast<std::size_t>(n * n));
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					A[static_cast<std::size_t>(j * n + i)] = M(i, j);
			IPIV.resize(static_cast<std::size_t>(n));
			int info;
#ifdef _WIN32
			DGETRF(&n, &n, A.data(), &n, IPIV.data(), &info);
#else
			dgetrf_(&n, &n, A.data(), &n, IPIV.data(), &info);
#endif
			return (info == 0);
		}
		bool solve(const Vector& B, Vector& X) {
			int n = B.Size();
			if (n == 0 || n != static_cast<int>(IPIV.size()))
				return false;
			X = B;
			int nrhs = 1;
			int info;
#ifdef _WIN32
			DGETRS("N", &n, &nrhs, A.data(), &n, IPIV.data(), &X(0), &n, &info);
#else
			dgetrs_("N", &n, &nrhs, A.data(), &n, IPIV.data(), &X(0), &n, &info);
#endif
			return (info == 0);
		}
	};

	// a custom stream to store ResponseType
	class CustomStream : public DummyStream {
	public:
		std::vector<std::string> components;
		CustomStream() = default;
		int tag(const char* name, const char* value) {
			if (strcmp(name, "ResponseType") == 0)
				components.push_back(value);
			return 0;
		};
	};

	// Response wrapper
	class ResponseWrapper {
	public:
		std::string name;
		std::vector<std::string> components;
		std::vector<Response*> responses;
		ResponseWrapper() = default;
		ResponseWrapper(std::size_t n) : responses(n, nullptr) {}
		~ResponseWrapper() {
			for (Response* ires : responses) {
				if (ires)
					delete ires;
			}
		}
	};
}

void *OPS_Series3DMaterial(void)
{
	// info
	const char* info = "nDMaterial Series3D $tag    $tag1 $tag2 ... $tagN   <-weights $w1 $w2 ... $wN> <-maxIter $maxIter> <-relTol $relTol> <-absTol $absTol> <-verbose>";

	// utility for parsing
	static std::vector<char> my_buffer(1024);
	auto get_string_input = []() -> std::string {
		return OPS_GetStringFromAll(my_buffer.data(), static_cast<int>(my_buffer.size()));
	};

	// convert string to int
	auto to_int = [](const std::string& x, bool& ok) {
		std::size_t pos;
		int res = 0;
		try {
			res = std::stoi(x, &pos);
			ok = (pos == x.size());
		}
		catch (...) {
			ok = false;
		}
		return res;
	};

	// convert string to double
	auto to_double = [](const std::string& x, bool& ok) {
		std::size_t pos;
		double res = 0.0;
		try {
			res = std::stod(x, &pos);
			ok = (pos == x.size());
		}
		catch (...) {
			ok = false;
		}
		return res;
	};

	// check minimum required arguments
	int argc = OPS_GetNumRemainingInputArgs();
	if (argc < 2) {
		opserr <<
			"nDMaterial Series3D Error: Few arguments (< 2).\n"
			"Want: " << info << ".\n";
		return nullptr;
	}

	// the parsing stage
	// -1 = invalid, 0 = read -materials, 1 = read -weights
	int parsing_stage = 0;

	// begin parsing
	int counter = 0;
	bool ok;
	int tag = 0;
	std::vector<NDMaterial*> materials;
	std::vector<double> weights;
	int max_iter = 10;
	double rel_tol = 1.0e-4;
	double abs_tol = 1.0e-8;
	bool verbose = false;
	while (argc > 0) {
		// read next word
		std::string word = get_string_input();
		++counter;
		argc = OPS_GetNumRemainingInputArgs();
		// check tag first
		if (counter == 1) {
			tag = to_int(word, ok);
			if (!ok) {
				opserr << "nDMaterial Series3D Error: Cannot get the wrapper material tag (wrong word: \"" << word.c_str() << "\".\n";
				return nullptr;
			}
			continue;
		}
		// read material ids
		if (parsing_stage == 0) {
			int imaterial_tag = to_int(word, ok);
			if (ok) {
				NDMaterial* imaterial = OPS_getNDMaterial(imaterial_tag);
				if (imaterial == nullptr) {
					opserr <<
						"nDMaterial Series3D Error: No existing NDMaterial with tag " <<
						imaterial_tag << " for NDMaterial Series3D " << tag << ".\n";
					return nullptr;
				}
				materials.push_back(imaterial);
				continue;
			}
			else {
				parsing_stage = -1;
			}
		}
		// read weights
		if (parsing_stage == 1) {
			double iweight = to_double(word, ok);
			if (ok) {
				weights.push_back(iweight);
				continue;
			}
			else {
				parsing_stage = -1;
			}
		}
		// check others
		if (word == "-weights") {
			if (weights.size() > 0) {
				opserr << "nDMaterial Series3D Error: Cannot use the -weights keyword multiple times.\n";
				return nullptr;
			}
			parsing_stage = 1;
			continue;
		}
		if (word == "-maxIter") {
			if (argc == 0) {
				opserr << "nDMaterial Series3D Error: -maxIter keyword provided without a value.\n";
				return nullptr;
			}
			word = get_string_input();
			++counter;
			argc = OPS_GetNumRemainingInputArgs();
			max_iter = std::abs(to_int(word, ok));
			if (!ok) {
				opserr << "nDMaterial Series3D Error: -maxIter keyword provided with an invalid value(\"" << word.c_str() << "\").\n";
				return nullptr;
			}
		}
		if (word == "-relTol") {
			if (argc == 0) {
				opserr << "nDMaterial Series3D Error: -relTol keyword provided without a value.\n";
				return nullptr;
			}
			word = get_string_input();
			++counter;
			argc = OPS_GetNumRemainingInputArgs();
			rel_tol = std::abs(to_double(word, ok));
			if (!ok) {
				opserr << "nDMaterial Series3D Error: -relTol keyword provided with an invalid value(\"" << word.c_str() << "\").\n";
				return nullptr;
			}
		}
		if (word == "-absTol") {
			if (argc == 0) {
				opserr << "nDMaterial Series3D Error: -absTol keyword provided without a value.\n";
				return nullptr;
			}
			word = get_string_input();
			++counter;
			argc = OPS_GetNumRemainingInputArgs();
			abs_tol = std::abs(to_double(word, ok));
			if (!ok) {
				opserr << "nDMaterial Series3D Error: -absTol keyword provided with an invalid value(\"" << word.c_str() << "\").\n";
				return nullptr;
			}
		}
		if (word == "-verbose") {
			verbose = true;
		}
	}

	// some final check
	if (materials.size() == 0) {
		opserr << "nDMaterial Series3D Error: No material provided.\n";
		return nullptr;
	}
	if (weights.size() == 0) {
		// no user-defined weights... by default use unit weight for all materials
		weights.resize(materials.size(), 1.0);
	}
	else {
		if (weights.size() != materials.size()) {
			opserr << "nDMaterial Series3D Error: the number of materials (" << 
				static_cast<int>(materials.size()) << ") must be equal to the number of weights (" <<
				static_cast<int>(weights.size()) << ").\n";
			return nullptr;
		}
	}

	// create the series material wrapper
	NDMaterial* result = new Series3DMaterial(tag, materials, weights, max_iter, rel_tol, abs_tol, verbose);
	if (result == nullptr) {
		opserr << "WARNING could not create NDMaterial of type Series3D.\n";
		return nullptr;
	}
	return result;
}

Series3DMaterial::Series3DMaterial(
	int tag,
	const std::vector<NDMaterial*>& theMaterials,
	const std::vector<double>& theWeights,
	int theMaxNumberOfIterations,
	double theRelativeTolerance,
	double theAbsoluteTolerance,
	bool verbose)
	: NDMaterial(tag, ND_TAG_Series3DMaterial)
	, m_materials(theMaterials.size(), nullptr)
	, m_weights(theWeights)
	, m_max_iter(theMaxNumberOfIterations)
	, m_rel_tol(theRelativeTolerance)
	, m_abs_tol(theAbsoluteTolerance)
	, m_verbose(verbose)
	, m_mat_strain_commit(theMaterials.size(), Vector(6))
{
	// copy the materials
	for (std::size_t i = 0; i < theMaterials.size(); ++i) {
		NDMaterial* the_copy = theMaterials[i]->getCopy("ThreeDimensional");
		if (the_copy == 0) {
			opserr << 
				"nDMaterial Series3D Error: failed to get a (3D) copy of the material at location " << 
				static_cast<int>(i)+1 << " of " << static_cast<int>(theMaterials.size()) << "\n";
			exit(-1);
		}
		m_materials[i] = the_copy;
	}

	// compute initial tangent here and also the sabilization term
	static Matrix iCinv(6, 6);
	static Matrix Cinv(6, 6);
	Cinv.Zero();
	for (std::size_t i = 0; i < m_materials.size(); ++i) {
		const Matrix& iC = m_materials[i]->getInitialTangent();
		if (iC.Invert(iCinv) < 0) {
			opserr << 
				"nDMaterial Series3D Error: Cannot invert the initial tangent of material " << 
				static_cast<int>(i) + 1 << "\n";
			exit(-1);
		}
		Cinv.addMatrix(1.0, iCinv, m_weights[i]);
	}
	if (Cinv.Invert(m_initial_tangent) < 0) {
		opserr << "nDMaterial Series3D Error: Cannot invert the homogenized initial tangent.\n"
			"Make sure the materials are properly defined.\n";
		exit(-1);
	}
	m_tangent = m_initial_tangent;
}

Series3DMaterial::Series3DMaterial()
	: NDMaterial(0, ND_TAG_Series3DMaterial)
{
}

Series3DMaterial::~Series3DMaterial()
{ 
	for (NDMaterial* item : m_materials) {
		if (item) {
			delete item;
		}
	}
}

double Series3DMaterial::getRho(void)
{
	double rho = 0.0;
	for (std::size_t i = 0; i < m_materials.size(); ++i) {
		rho += m_materials[i]->getRho() * m_weights[i];
	}
	return rho;
}

int Series3DMaterial::setTrialStrain(const Vector& strain)
{
	// set homogenized strain
	m_strain = strain;
	// impose iso-stress condition
	IterativeTangentType ittype = IT_Tangent;
	bool converged = imposeIsoStressCondition(ittype);
	if (!converged) {
		ittype = IT_StabilizedTangent;
		converged = imposeIsoStressCondition(ittype);
		if (!converged) {
			ittype = IT_Initial;
			int aux = m_max_iter;
			m_max_iter = aux * 10;
			converged = imposeIsoStressCondition(ittype);
			m_max_iter = aux;
		}
		ittype = IT_Tangent;
	}
	// compute tangent
	computeHomogenizedStress();
	computeHomogenizedTangent(ittype);
	// done
	return converged ? 0 : -1;
}

const Vector &Series3DMaterial::getStrain(void)
{
	return m_strain;
}

const Vector &Series3DMaterial::getStress(void)
{
	return m_stress;
}

const Matrix &Series3DMaterial::getTangent(void)
{
	return m_tangent;
}

const Matrix &Series3DMaterial::getInitialTangent(void)
{
	return m_initial_tangent;
}

int Series3DMaterial::commitState(void)
{
	// commit each material
	int res = 0;
	for (std::size_t i = 0; i < m_materials.size(); ++i) {
		if (m_materials[i]->commitState() != 0) {
			res = -1;
		}
		m_mat_strain_commit[i] = m_materials[i]->getStrain();
	}
	// re-compute homogenized stress in case some sub-material
	// is using some special integration scheme (such as impl-ex)
	computeHomogenizedStress();
	// commit this
	m_lambda_commit = m_lambda;
	m_strain_commit = m_strain;
	m_stress_commit = m_stress;
	// done
	return res;
}

int Series3DMaterial::revertToLastCommit(void)
{
	// revert each material
	int res = 0;
	for (std::size_t i = 0; i < m_materials.size(); ++i) {
		if (m_materials[i]->revertToLastCommit() != 0) {
			res = -1;
		}
		m_materials[i]->setTrialStrain(m_mat_strain_commit[i]);
	}
	// revert this
	m_lambda = m_lambda_commit;
	m_strain = m_strain_commit;
	m_stress = m_stress_commit;
	// done
	return res;
}

int Series3DMaterial::revertToStart(void)
{
	int res = 0;
	// only if not called by the initial state analysis off command
	if (!ops_InitialStateAnalysis) { 
		m_lambda.Zero();
		m_lambda_commit.Zero();
		m_strain.Zero();
		m_strain_commit.Zero();
		m_stress.Zero();
		m_stress_commit.Zero();
		m_tangent.Zero();
		for (std::size_t i = 0; i < m_materials.size(); ++i) {
			if (m_materials[i]->revertToStart() != 0) {
				res = -1;
			}
			m_mat_strain_commit[i].Zero();
		}
	}
	return res;
}

NDMaterial* Series3DMaterial::getCopy(void)
{
	Series3DMaterial *theCopy = new Series3DMaterial();
	theCopy->setTag(getTag());
	
	theCopy->m_materials.resize(m_materials.size());
	theCopy->m_weights.resize(m_weights.size());
	theCopy->m_mat_strain_commit.resize(m_mat_strain_commit.size());
	for (std::size_t i = 0; i < m_materials.size(); ++i) {
		theCopy->m_materials[i] = m_materials[i]->getCopy("ThreeDimensional");
		theCopy->m_weights[i] = m_weights[i];
		theCopy->m_mat_strain_commit[i] = m_mat_strain_commit[i];
	}

	theCopy->m_max_iter = m_max_iter;
	theCopy->m_rel_tol = m_rel_tol;
	theCopy->m_abs_tol = m_abs_tol;
	theCopy->m_verbose = m_verbose;

	theCopy->m_lambda = m_lambda;
	theCopy->m_lambda_commit = m_lambda_commit;
	
	theCopy->m_strain = m_strain;
	theCopy->m_strain_commit = m_strain_commit;
	theCopy->m_stress = m_stress;
	theCopy->m_stress_commit = m_stress_commit;
	theCopy->m_tangent = m_tangent;
	theCopy->m_initial_tangent = m_initial_tangent;

	theCopy->m_stab = m_stab;

	return theCopy;
}

NDMaterial* Series3DMaterial::getCopy(const char* code)
{
	if (strcmp(code, "ThreeDimensional") == 0)
		return getCopy();
	return NDMaterial::getCopy(code);
}

const char* Series3DMaterial::getType(void) const
{
	return "ThreeDimensional";
}

int Series3DMaterial::getOrder(void) const
{
	return 6;
}

void Series3DMaterial::Print(OPS_Stream &s, int flag)
{
	s << "Series3D Material, tag: " << this->getTag() << "\n";
}

int Series3DMaterial::sendSelf(int commitTag, Channel &theChannel)
{
	int num_mat = static_cast<int>(m_materials.size());

	// send basic int data
	static ID I1(4);
	I1(0) = getTag();
	I1(1) = num_mat;
	I1(2) = m_max_iter;
	I1(3) = static_cast<int>(m_verbose);
	if (theChannel.sendID(getDbTag(), commitTag, I1) < 0) {
		opserr << "Series3DMaterial::sendSelf() - failed to send data (I1)\n";
		return -1;
	}

	// send double data
	static Vector D1;
	D1.resize(num_mat*2 /*ints for mats*/ + num_mat * 7 /*floats for mats*/ + 111 /*other fixed data*/);
	int counter = 0;
	/*ints for mats*/
	for (int i = 0; i < num_mat; ++i) {
		NDMaterial* imat = m_materials[static_cast<std::size_t>(i)];
		D1(counter++) = static_cast<double>(imat->getClassTag());
		int mat_db_tag = imat->getDbTag();
		if (mat_db_tag == 0) {
			mat_db_tag = theChannel.getDbTag();
			if (mat_db_tag != 0)
				imat->setDbTag(mat_db_tag);
		}
		D1(counter++) = static_cast<double>(mat_db_tag);
	}
	/*floats for mats*/
	for (int i = 0; i < num_mat; ++i) {
		D1(counter++) = m_weights[static_cast<std::size_t>(i)];
		const Vector& istrain_commit = m_mat_strain_commit[i];
		for (int j = 0; j < 6; ++j)
			D1(counter++) = istrain_commit(j);
	}
	/*other fixed data*/
	for (int i = 0; i < 6; ++i) 
		D1(counter++) = m_lambda(i);
	for (int i = 0; i < 6; ++i) 
		D1(counter++) = m_lambda_commit(i);
	for (int i = 0; i < 6; ++i)
		D1(counter++) = m_strain(i);
	for (int i = 0; i < 6; ++i)
		D1(counter++) = m_strain_commit(i);
	for (int i = 0; i < 6; ++i)
		D1(counter++) = m_stress(i);
	for (int i = 0; i < 6; ++i)
		D1(counter++) = m_stress_commit(i);
	for (int i = 0; i < 6; ++i)
		for(int j = 0; j < 6; ++j)
			D1(counter++) = m_tangent(i, j);
	for (int i = 0; i < 6; ++i)
		for (int j = 0; j < 6; ++j)
			D1(counter++) = m_initial_tangent(i, j);
	D1(counter++) = m_rel_tol;
	D1(counter++) = m_abs_tol;
	D1(counter++) = m_stab;
	if (theChannel.sendVector(getDbTag(), commitTag, D1) < 0) {
		opserr << "Series3DMaterial::sendSelf() - failed to send data (D1)\n";
		return -1;
	}

	// send all materials
	for (std::size_t i = 0; i < m_materials.size(); ++i) {
		if (m_materials[i]->sendSelf(commitTag, theChannel) < 0) {
			opserr << "Series3DMaterial::sendSelf() - failed to send material " << static_cast<int>(i) << "\n";
			return -1;
		}
	}

	// done
	return 0;
}

int Series3DMaterial::recvSelf(int commitTag, Channel & theChannel, FEM_ObjectBroker & theBroker)
{
	// receive basic int data,
	// reset and reallocate materials and weights vectors
	static ID I1(4);
	if (theChannel.recvID(getDbTag(), commitTag, I1) < 0) {
		opserr << "Series3DMaterial::recvSelf() - failed to receive data (I1)\n";
		return -1;
	}
	setTag(I1(0));
	int num_mat = I1(1);
	m_max_iter = I1(2);
	m_verbose = static_cast<bool>(I1(3));
	if (m_materials.size() > 0) {
		for (std::size_t i = 0; i < m_materials.size(); ++i) {
			if (m_materials[i])
				delete m_materials[i];
		}
		m_materials.clear();
		m_weights.clear();
	}
	m_materials.resize(static_cast<std::size_t>(num_mat), nullptr);
	m_weights.resize(static_cast<std::size_t>(num_mat), 0.0);
	m_mat_strain_commit.resize(static_cast<std::size_t>(num_mat), Vector(6));

	// receive double data
	static Vector D1;
	D1.resize(num_mat * 2 /*ints for mats*/ + num_mat * 7 /*floats for mats*/ + 111 /*other fixed data*/);
	int counter = 0;
	if (theChannel.recvVector(getDbTag(), commitTag, D1) < 0) {
		opserr << "Series3DMaterial::recvSelf() - failed to receive data (D1)\n";
		return -1;
	}
	/*ints for mats*/
	counter = 2 * num_mat; // left for later
	/*floats for mats*/
	for (int i = 0; i < num_mat; ++i) {
		m_weights[static_cast<std::size_t>(i)] = D1(counter++);
		Vector& istrain_commit = m_mat_strain_commit[i];
		for (int j = 0; j < 6; ++j)
			istrain_commit(j) = D1(counter++);
	}
	/*other fixed data*/
	for (int i = 0; i < 6; ++i)
		m_lambda(i) = D1(counter++);
	for (int i = 0; i < 6; ++i)
		m_lambda_commit(i) = D1(counter++);
	for (int i = 0; i < 6; ++i)
		m_strain(i) = D1(counter++);
	for (int i = 0; i < 6; ++i)
		m_strain_commit(i) = D1(counter++);
	for (int i = 0; i < 6; ++i)
		m_stress(i) = D1(counter++);
	for (int i = 0; i < 6; ++i)
		m_stress_commit(i) = D1(counter++);
	for (int i = 0; i < 6; ++i)
		for (int j = 0; j < 6; ++j)
			m_tangent(i, j) = D1(counter++);
	for (int i = 0; i < 6; ++i)
		for (int j = 0; j < 6; ++j)
			m_initial_tangent(i, j) = D1(counter++);
	m_rel_tol = D1(counter++);
	m_abs_tol = D1(counter++);
	m_stab = D1(counter++);

	// receive all materials
	counter = 0;
	for (int i = 0; i < num_mat; i++) {
		int mat_class_tag = D1(counter++);
		int mat_db_tag = D1(counter++);
		NDMaterial* the_material = theBroker.getNewNDMaterial(mat_class_tag);
		if (the_material == nullptr) {
			opserr << "Series3DMaterial::recvSelf() - could not get a new NDMaterial\n";
			return -1;
		}
		the_material->setDbTag(mat_db_tag);
		if (the_material->recvSelf(commitTag, theChannel, theBroker) < 0) {
			opserr << "Series3DMaterial::recvSelf() - failed to receive material\n";
			return -1;
		}
		m_materials[static_cast<std::size_t>(i)] = the_material;
	}

	// done
	return 0;
}

int Series3DMaterial::setParameter(const char** argv, int argc, Parameter& param)
{
	// forward to the sub-materials
	int res = -1;
	for (NDMaterial* item : m_materials) {
		if (item->setParameter(argv, argc, param) == 0)
			res = 0;
	}
	return res;
}

Response* Series3DMaterial::setResponse(const char** argv, int argc, OPS_Stream& output)
{
	if (argc > 0) {
		if (strcmp(argv[0], "material") == 0) {
			if (argc > 2) {
				int imat = atoi(argv[1]) - 1;
				if (imat >= 0 && imat < static_cast<int>(m_materials.size())) {
					return m_materials[static_cast<std::size_t>(imat)]->setResponse(&argv[2], argc - 2, output);
				}
			}
		}
		else if (strcmp(argv[0], "homogenized") == 0 && m_materials.size() > 0) {
			if (argc > 1) {
				std::string name = argv[1];
				std::shared_ptr<Series3DUtils::ResponseWrapper> wres;
				int wres_id = 0;
				// get from map ...
				for (auto& item : m_response_map) {
					if (item.second->name == name) {
						wres = item.second;
						wres_id = item.first;
						break;
					}
				}
				// ... or try make it and map it (generate id >= 1000)
				if (wres == nullptr) {
					wres = std::make_shared<Series3DUtils::ResponseWrapper>(m_materials.size());
					wres_id = m_response_map.size() > 0 ? m_response_map.rbegin()->first + 1 : 1000;
					wres->name = name;
					int response_size = 0;
					for (std::size_t i = 0; i < m_materials.size(); ++i) {
						NDMaterial* imat = m_materials[i];
						Series3DUtils::CustomStream ds;
						Response* ires = imat->setResponse(&argv[1], argc - 1, ds);
						if (ires) {
							int ires_size = ires->getInformation().getData().Size();
							if (response_size == 0)
								response_size = ires_size;
							if (ires_size == 0 || ires_size != response_size) {
								// validity check
								response_size = 0;
								wres = nullptr;
								break;
							}
							if (response_size > 0 && 
								wres->components.size() == 0 && 
								ds.components.size() == static_cast<std::size_t>(response_size)) {
								wres->components = ds.components;
							}
							wres->responses[i] = ires;
						}
					}
					// last validity check
					if (response_size == 0)
						wres = nullptr;
					// map it
					if (wres)
						m_response_map[wres_id] = wres;
				}
				// go on if valid
				if (wres) {
					output.tag("NdMaterialOutput");
					output.attr("matType", getClassType());
					output.attr("matTag", getTag());
					for (const auto& item : wres->components)
						output.tag("ResponseType", item.c_str());
					Vector data(static_cast<int>(wres->components.size()));
					MaterialResponse* resp = new MaterialResponse(this, wres_id, data);
					output.endTag();
					return resp;
				}
			}
		}
	}
	return NDMaterial::setResponse(argv, argc, output);
}

int Series3DMaterial::getResponse(int responseID, Information& matInformation)
{
	// check if it is a homogenized response
	auto iter = m_response_map.find(responseID);
	if (iter != m_response_map.end()) {
		auto wres = iter->second;
		if (matInformation.theVector) {
			matInformation.theVector->Zero();
			double wsum = 0.0;
			for (std::size_t i = 0; i < m_materials.size(); ++i) {
				if (wres->responses[i] == nullptr) continue;
				if (wres->responses[i]->getResponse() < 0) continue;
				const Vector& idata = wres->responses[i]->getInformation().getData();
				if (idata.Size() == matInformation.theVector->Size()) {
					matInformation.theVector->addVector(1.0, idata, m_weights[i]);
				}
				wsum += m_weights[i];
			}
			if (wsum > 0.0)
				(*matInformation.theVector) /= wsum;
			return 0;
		}
	}
	// default
	return NDMaterial::getResponse(responseID, matInformation);
}

bool Series3DMaterial::imposeIsoStressCondition(IterativeTangentType ittype)
{
	// declare the static solver
	static Series3DUtils::SolverWrapper solver;

	// restore committed lagrange multipliers
	m_lambda = m_lambda_commit;

	// restore committed state in all materials
	for (std::size_t i = 0; i < m_materials.size(); ++i)
		m_materials[i]->setTrialStrain(m_mat_strain_commit[i]);

	// compute initial residual norm
	double NR0 = computeResidualNorm();

	// Newton iteration to impose the iso-stress condition
	if (m_verbose) {
		opserr << "\n   Series3D (" << getTag() << ") - impose iso-stress condition\n";
	}
	bool converged = false;
	int iter;
	for (iter = 0; iter < m_max_iter; ++iter) {

		// compute denominator and factorize it
		const Matrix& D = computeDenominator(ittype);
		if (!solver.factorize(D)) {
			if (m_verbose) {
				opserr << "      singular D matrix\n";
			}
			break;
		}

		// compute ewg
		const Vector& ewg = computeWeightedStrainResidual();

		// solve for unknowns (strain vectors and lagrange multipliers)
		if (!solveForStrainVectors(ewg, ittype, solver, D)) {
			if (m_verbose) {
				opserr << "      cannot solve for dStrain and/or dLambda\n";
			}
			break;
		}

		// test convergence
		double NR = computeResidualNorm();
		double NRratio = NR0 > 0.0 ? NRratio = NR / NR0 : 1.0;
		if (m_verbose) {
			opserr << "      iter: " << iter + 1 << " - ratio: " << NRratio << " - norm: " << NR << "\n";
		}
		if (NRratio < m_rel_tol || NR < m_abs_tol) {
			converged = true;
			break;
		}
	}

	if (m_verbose) {
		if (converged)
			opserr << "      converged in " << iter << " iterations\n";
		else
			opserr << "      non-convergence in iso-stress constraint\n";
	}

	// done
	return converged;
}

double Series3DMaterial::computeResidualNorm()
{
	/*
	Computes the residual norm.
	The residual vector consists in N+1 6x1-blocks.
	The first N blocks are for the strain DOFs, while the last block
	is for the lagrange multipliers.
	It takes the following form:

	R:
	[     -w1 * (lambda + s1)     ]
	[     -w2 * (lambda + s2)     ]
	[             ...             ]
	[     -wN * (lambda + sN)     ]
	[ g - e1*w1 -e2*w2 ... -eN*wN ]

	where:
	- wi is the weight of the i-th material;
	- si is the stress vector of the i-th material;
	- ei is the strain vector of the i-th material;
	- lambda is the lagrange mutliplier vector;
	- g is the homogenized strain;
	*/

	static Vector NR2(6);
	static Vector NRtemp(6);

	double NR1 = 0.0;
	NR2 = m_strain;
	for (std::size_t i = 0; i < m_materials.size(); ++i) {
		double wi = m_weights[i];
		const Vector& si = m_materials[i]->getStress();
		const Vector& ei = m_materials[i]->getStrain();

		NRtemp = m_lambda;
		NRtemp.addVector(1.0, si, 1.0);
		NRtemp *= -wi;

		NR1 += NRtemp ^ NRtemp;

		NR2.addVector(1.0, ei, -wi);
	}

	return sqrt(NR1 + (NR2 ^ NR2));
}

const Matrix& Series3DMaterial::computeDenominator(IterativeTangentType ittype)
{
	/*
	Computes the denominator matrix D used to solve for the unknowns.
	It is a 6x6 matrix and takes the following form:

	D: SUM_i(wi * MULT_j(kj, j!=i))
	*/

	static Matrix D(6, 6);
	static Matrix Di(6, 6);
	static Matrix DiTemp(6, 6);

	D.Zero();
	for (std::size_t i = 0; i < m_materials.size(); ++i) {
		double wi = m_weights[i];
		Di.Zero();
		for (int q = 0; q < 6; q++)
			Di(q, q) = wi;
		for (std::size_t j = 0; j < m_materials.size(); ++j) {
			if (i != j) {
				const Matrix& Kj = getMaterialTangent(m_materials[j], ittype);
				DiTemp.addMatrixProduct(0.0, Di, Kj, 1.0);
				Di = DiTemp;
			}
		}
		D.addMatrix(1.0, Di, 1.0);
	}

	return D;
}

const Vector& Series3DMaterial::computeWeightedStrainResidual()
{
	/*
	Computes the 6x1 vector as the difference between
	the weighted-average of the materials' strain vectors and 
	the homogenized strain vector

	ewg: SUM_i(wi * ei) - g

	where:
	- wi is the weight of the i-th material;
	- ei is the strain vector of the i-th material;
	- g is the homogenized strain;
	*/

	static Vector ewg(6);

	ewg.addVector(0.0, m_strain, -1.0);
	for(std::size_t i = 0; i < m_materials.size(); ++i) {
		double wi = m_weights[i];
		const Vector& ei = m_materials[i]->getStrain();
		ewg.addVector(1.0, ei, wi);
	}

	return ewg;
}

bool Series3DMaterial::solveForLagrangeMultipliers(
	const Vector& ewg,
	IterativeTangentType ittype,
	Series3DUtils::SolverWrapper& solver)
{
	/*
	solves for the iterative correction of the lagrange multiplier 6x1 vector

	dLambda: D^-1 * (MULT_i(ki)*ewg - SUM_i( MULT_j(kj, j!=i) * wi*(lambda + si) ))

	where:
	- wi is the weight of the i-th material;
	- ei is the strain vector of the i-th material;
	- g is the homogenized strain;
	*/

	static Matrix A(6, 6);
	static Matrix ATemp(6, 6);
	static Vector B(6);
	static Matrix Bi(6, 6);
	static Matrix BiTemp(6, 6);
	static Vector lsi(6);
	static Vector Blsi(6);
	static Vector Aewg(6);
	static Vector dLambda(6);

	A.Zero();
	for (int q = 0; q < 6; ++q)
		A(q, q) = 1.0;
	B.Zero();
	for (std::size_t i = 0; i < m_materials.size(); ++i) {

		double wi = m_weights[i];
		const Matrix& ki = getMaterialTangent(m_materials[i], ittype);

		ATemp.addMatrixProduct(0.0, A, ki, 1.0);
		A = ATemp;

		Bi.Zero();
		for (int q = 0; q < 6; ++q)
			Bi(q, q) = 1.0;
		for (std::size_t j = 0; j < m_materials.size(); ++j) {
			if (i != j) {
				const Matrix& kj = getMaterialTangent(m_materials[j], ittype);
				BiTemp.addMatrixProduct(0.0, Bi, kj, 1.0);
				Bi = BiTemp;
			}
		}

		lsi = m_lambda;
		lsi.addVector(1.0, m_materials[i]->getStress(), 1.0);
		Blsi.addMatrixVector(0.0, Bi, lsi, wi);
		B.addVector(1.0, Blsi, -1.0);

	}

	Aewg.addMatrixVector(0.0, A, ewg, 1.0);
	B.addVector(1.0, Aewg, 1.0);
	if (!solver.solve(B, dLambda))
		return false;

	m_lambda.addVector(1.0, dLambda, 1.0);
	return true;
}

bool Series3DMaterial::solveForStrainVectors(
	const Vector& ewg, 
	IterativeTangentType ittype, 
	Series3DUtils::SolverWrapper& solver,
	const Matrix& D)
{
	/*
	solves for the iterative correction of the strain 6x1 vectors for each material

	dStrain_i: D^-1 * (
	     - MULT_j(kj, j!=i)*ewg
		 + SUM_j(  MULT_q(kq , q!=i and q!=j) * wj*(lambda + sj)  , j!=i)
		 - SUM_j(  MULT_q(kq , q!=i and q!=j) * wj  , j!=i) * (lambda + si)
		 )

	Notes:
		To unroll this loop we need to first store all the updates dStrain_i,
		and then we can call the setStrain method on the sub-materials,
		otherwise their stress and tangent will change after setting the new strain.
		To avoid repreted dynamic allocations, we create a static vector of Vector(6)
		of a reasonably large number, and only if necessary we resize it.
	*/

	// the storage for dStrain_i
	auto make_static_dStrain_vector = []() {
		constexpr std::size_t num = 100;
		std::vector<Vector> dStrain_storage(num);
		for (std::size_t i = 0; i < num; ++i)
			dStrain_storage[i].resize(6);
		return dStrain_storage;
	};
	static std::vector<Vector> dStrain_vector = make_static_dStrain_vector();
	if (m_materials.size() > dStrain_vector.size())
		dStrain_vector.resize(m_materials.size(), Vector(6));

	// other static variables
	static Matrix A(6, 6);
	static Vector B(6);
	static Matrix C(6, 6);
	static Matrix KKq(6, 6);
	static Vector lsi(6);
	static Vector lsj(6);
	static Matrix auxM(6, 6);

	// for each material, compute and store the strain correction dStrain
	for (std::size_t i = 0; i < m_materials.size(); ++i) {

		A.Zero();
		for (int q = 0; q < 6; ++q)
			A(q, q) = 1.0;
		B.Zero();
		C.Zero();

		for (std::size_t j = 0; j < m_materials.size(); ++j) {
			if (j != i) {

				const Matrix& kj = getMaterialTangent(m_materials[j], ittype);
				auxM.addMatrixProduct(0.0, A, kj, 1.0);
				A = auxM;

				KKq.Zero();
				for (int q = 0; q < 6; ++q)
					KKq(q, q) = 1.0;
				for (std::size_t q = 0; q < m_materials.size(); ++q) {
					if (q != i && q != j) {
						const Matrix& kq = getMaterialTangent(m_materials[q], ittype);
						auxM.addMatrixProduct(0.0, KKq, kq, 1.0);
						KKq = auxM;
					}
				}

				lsj = m_lambda;
				lsj.addVector(1.0, m_materials[j]->getStress(), 1.0);
				B.addMatrixVector(1.0, KKq, lsj, m_weights[j]);

				C.addMatrix(1.0, KKq, m_weights[j]);
			}
		}

		B.addMatrixVector(1.0, A, ewg, -1.0);
		lsi = m_lambda;
		lsi.addVector(1.0, m_materials[i]->getStress(), 1.0);
		B.addMatrixVector(1.0, C, lsi, -1.0);

		// solve for dStrain
		Vector& dStraini = dStrain_vector[i];
		if (!solver.solve(B, dStraini))
			return false;
	}

	// now, before updating the strain in each material,
	// we can solve for the lagrange multipliers.
	solveForLagrangeMultipliers(ewg, ittype, solver);

	// now we can update each material
	static Vector strain_new(6);
	for (std::size_t i = 0; i < m_materials.size(); ++i) {
		NDMaterial* imaterial = m_materials[i];
		strain_new = imaterial->getStrain();
		strain_new.addVector(1.0, dStrain_vector[i], 1.0);
		if (imaterial->setTrialStrain(strain_new) != 0)
			return false;
	}

	// done
	return true;
}

const Matrix& Series3DMaterial::getMaterialTangent(NDMaterial* mat, IterativeTangentType ittype) const
{
	static Matrix Kaux(6, 6);
	if (ittype == Series3DMaterial::IT_Initial) {
		return mat->getInitialTangent();
	}
	else if (ittype == Series3DMaterial::IT_StabilizedTangent) {
		Kaux.addMatrix(0.0, mat->getInitialTangent(), m_stab);
		Kaux.addMatrix(1.0, mat->getTangent(), 1.0 - m_stab);
		return Kaux;
	}
	else {
		return mat->getTangent();
	}
}

void Series3DMaterial::computeHomogenizedTangent(IterativeTangentType ittype)
{
	static Matrix iCinv(6, 6);
	static Matrix Cinv(6, 6);
	bool done;

	if (ittype == IT_Tangent) {
		done = true;
		Cinv.Zero();
		for (std::size_t i = 0; i < m_materials.size(); ++i) {
			NDMaterial* imaterial = m_materials[i];
			const Matrix& iC = getMaterialTangent(imaterial, ittype);
			if (iC.Invert(iCinv) != 0) {
				done = false;
				break;
			}
			Cinv.addMatrix(1.0, iCinv, m_weights[i]);
		}
		if (done) {
			done = (Cinv.Invert(m_tangent) == 0);
		}
		if (done) return; // ok
		ittype = IT_StabilizedTangent; // otherwise try the stabilized tangent
	}

	if (ittype == IT_StabilizedTangent) {
		done = true;
		Cinv.Zero();
		for (std::size_t i = 0; i < m_materials.size(); ++i) {
			NDMaterial* imaterial = m_materials[i];
			const Matrix& iC = getMaterialTangent(imaterial, ittype);
			if (iC.Invert(iCinv) != 0) {
				done = false;
				break;
			}
			Cinv.addMatrix(1.0, iCinv, m_weights[i]);
		}
		if (done) {
			done = (Cinv.Invert(m_tangent) == 0);
		}
		if (done) return; // ok
		ittype = IT_Initial; // otherwise go with initial tangent
	}

	// if none of the above worked, fallback to initial
	m_tangent = m_initial_tangent;
}

void Series3DMaterial::computeHomogenizedStress()
{
	// if the convergence is achieved, all materials should have the same stress.
	// we compute anyway the average, just because some material may use some non-standard
	// integration schemes...
	m_stress.Zero();
	double wsum = 0.0;
	for (std::size_t i = 0; i < m_materials.size(); ++i) {
		m_stress.addVector(1.0, m_materials[i]->getStress(), m_weights[i]);
		wsum += m_weights[i];
	}
	if (wsum > 0.0)
		m_stress /= wsum;
}
