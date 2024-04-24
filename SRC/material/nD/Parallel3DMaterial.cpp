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
// Parallel3DMaterial is an aggregation
// of NDMaterials (3D only) objects all considered acting in Parallel.
//

#include <Parallel3DMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <DummyStream.h>
#include <MaterialResponse.h>
#include <ID.h>
#include <cmath>
#include <string>

// namespace for utilities
namespace Parallel3DUtils {

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

void *OPS_Parallel3DMaterial(void)
{
	// info
	const char* info = "nDMaterial Parallel3D $tag    $tag1 $tag2 ... $tagN   <-weights $w1 $w2 ... $wN>";

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
			"nDMaterial Parallel3D Error: Few arguments (< 2).\n"
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
	while (argc > 0) {
		// read next word
		std::string word = get_string_input();
		++counter;
		argc = OPS_GetNumRemainingInputArgs();
		// check tag first
		if (counter == 1) {
			tag = to_int(word, ok);
			if (!ok) {
				opserr << "nDMaterial Parallel3D Error: Cannot get the wrapper material tag (wrong word: \"" << word.c_str() << "\".\n";
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
						"nDMaterial Parallel3D Error: No existing NDMaterial with tag " <<
						imaterial_tag << " for NDMaterial Parallel3D " << tag << ".\n";
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
				opserr << "nDMaterial Parallel3D Error: Cannot use the -weights keyword multiple times.\n";
				return nullptr;
			}
			parsing_stage = 1;
			continue;
		}
	}

	// some final check
	if (materials.size() == 0) {
		opserr << "nDMaterial Parallel3D Error: No material provided.\n";
		return nullptr;
	}
	if (weights.size() == 0) {
		// no user-defined weights... by default use unit weight for all materials
		weights.resize(materials.size(), 1.0);
	}
	else {
		if (weights.size() != materials.size()) {
			opserr << "nDMaterial Parallel3D Error: the number of materials (" << 
				static_cast<int>(materials.size()) << ") must be equal to the number of weights (" <<
				static_cast<int>(weights.size()) << ").\n";
			return nullptr;
		}
	}

	// create the series material wrapper
	NDMaterial* result = new Parallel3DMaterial(tag, materials, weights);
	if (result == nullptr) {
		opserr << "WARNING could not create NDMaterial of type Parallel3D.\n";
		return nullptr;
	}
	return result;
}

Parallel3DMaterial::Parallel3DMaterial(
	int tag,
	const std::vector<NDMaterial*>& theMaterials,
	const std::vector<double>& theWeights)
	: NDMaterial(tag, ND_TAG_Parallel3DMaterial)
	, m_materials(theMaterials.size(), nullptr)
	, m_weights(theWeights)
{
	// copy the materials
	for (std::size_t i = 0; i < theMaterials.size(); ++i) {
		NDMaterial* the_copy = theMaterials[i]->getCopy("ThreeDimensional");
		if (the_copy == 0) {
			opserr << 
				"nDMaterial Paralell3D Error: failed to get a (3D) copy of the material at location " << 
				static_cast<int>(i)+1 << " of " << static_cast<int>(theMaterials.size()) << "\n";
			exit(-1);
		}
		m_materials[i] = the_copy;
	}
	computeInitialTangent();
}

Parallel3DMaterial::Parallel3DMaterial()
	: NDMaterial(0, ND_TAG_Parallel3DMaterial)
{
}

Parallel3DMaterial::~Parallel3DMaterial()
{ 
	for (NDMaterial* item : m_materials) {
		if (item) {
			delete item;
		}
	}
}

double Parallel3DMaterial::getRho(void)
{
	double rho = 0.0;
	for (std::size_t i = 0; i < m_materials.size(); ++i) {
		rho += m_materials[i]->getRho() * m_weights[i];
	}
	return rho;
}

int Parallel3DMaterial::setTrialStrain(const Vector& v)
{
	m_strain = v;
	int result = 0;
	for (std::size_t i = 0; i < m_materials.size(); ++i) {
		if (!m_materials[i]->setTrialStrain(v))
			result = -1;
	}
	computeStress();
	computeTangent();
	return result;
}

int Parallel3DMaterial::setTrialStrain(const Vector& v, const Vector& r)
{
	m_strain = v;
	int result = 0;
	for (std::size_t i = 0; i < m_materials.size(); ++i) {
		if (!m_materials[i]->setTrialStrain(v, r))
			result = -1;
	}
	computeStress();
	computeTangent();
	return result;
}

const Vector &Parallel3DMaterial::getStrain(void)
{
	return m_strain;
}

const Vector &Parallel3DMaterial::getStress(void)
{
	return m_stress;
}

const Matrix &Parallel3DMaterial::getTangent(void)
{
	return m_tangent;
}

const Matrix &Parallel3DMaterial::getInitialTangent(void)
{
	return m_initial_tangent;
}

int Parallel3DMaterial::commitState(void)
{
	// commit each material
	int res = 0;
	for (std::size_t i = 0; i < m_materials.size(); ++i) {
		if (m_materials[i]->commitState() != 0)
			res = -1;
	}
	// commit this
	m_strain_commit = m_strain;
	// re-compute homogenized stress and tangent
	computeStress();
	computeTangent();
	// done
	return res;
}

int Parallel3DMaterial::revertToLastCommit(void)
{
	// revert each material
	int res = 0;
	for (std::size_t i = 0; i < m_materials.size(); ++i) {
		if (m_materials[i]->revertToLastCommit() != 0)
			res = -1;
	}
	// revert this
	m_strain = m_strain_commit;
	// re-compute homogenized stress and tangent
	computeStress();
	computeTangent();
	// done
	return res;
}

int Parallel3DMaterial::revertToStart(void)
{
	int res = 0;
	// only if not called by the initial state analysis off command
	if (!ops_InitialStateAnalysis) { 
		m_strain.Zero();
		m_strain_commit.Zero();
		m_stress.Zero();
		m_tangent.Zero();
		for (std::size_t i = 0; i < m_materials.size(); ++i) {
			if (m_materials[i]->revertToStart() != 0)
				res = -1;
		}
	}
	return res;
}

NDMaterial* Parallel3DMaterial::getCopy(void)
{
	return new Parallel3DMaterial(getTag(), m_materials, m_weights);
}

NDMaterial* Parallel3DMaterial::getCopy(const char* code)
{
	if (strcmp(code, "ThreeDimensional") == 0)
		return getCopy();
	return NDMaterial::getCopy(code);
}

const char* Parallel3DMaterial::getType(void) const
{
	return "ThreeDimensional";
}

int Parallel3DMaterial::getOrder(void) const
{
	return 6;
}

void Parallel3DMaterial::Print(OPS_Stream &s, int flag)
{
	s << "Parallel3D Material, tag: " << this->getTag() << "\n";
}

int Parallel3DMaterial::sendSelf(int commitTag, Channel &theChannel)
{
	int num_mat = static_cast<int>(m_materials.size());

	// send basic int data
	static ID I1(2);
	I1(0) = getTag();
	I1(1) = num_mat;
	if (theChannel.sendID(getDbTag(), commitTag, I1) < 0) {
		opserr << "Parallel3DMaterial::sendSelf() - failed to send data (I1)\n";
		return -1;
	}

	// send double data
	static Vector D1;
	D1.resize(num_mat*2 /*ints for mats*/ + num_mat /*floats for mats*/ + 90 /*other fixed data*/);
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
	}
	/*other fixed data*/
	for (int i = 0; i < 6; ++i) 
		D1(counter++) = m_strain(i);
	for (int i = 0; i < 6; ++i) 
		D1(counter++) = m_strain_commit(i);
	for (int i = 0; i < 6; ++i)
		D1(counter++) = m_stress(i);
	for (int i = 0; i < 6; ++i)
		for(int j = 0; j < 6; ++j)
			D1(counter++) = m_tangent(i, j);
	for (int i = 0; i < 6; ++i)
		for (int j = 0; j < 6; ++j)
			D1(counter++) = m_initial_tangent(i, j);
	if (theChannel.sendVector(getDbTag(), commitTag, D1) < 0) {
		opserr << "Parallel3DMaterial::sendSelf() - failed to send data (D1)\n";
		return -1;
	}

	// send all materials
	for (std::size_t i = 0; i < m_materials.size(); ++i) {
		if (m_materials[i]->sendSelf(commitTag, theChannel) < 0) {
			opserr << "Parallel3DMaterial::sendSelf() - failed to send material " << static_cast<int>(i) << "\n";
			return -1;
		}
	}

	// done
	return 0;
}

int Parallel3DMaterial::recvSelf(int commitTag, Channel & theChannel, FEM_ObjectBroker & theBroker)
{
	// receive basic int data,
	// reset and reallocate materials and weights vectors
	static ID I1(2);
	if (theChannel.recvID(getDbTag(), commitTag, I1) < 0) {
		opserr << "Parallel3DMaterial::recvSelf() - failed to receive data (I1)\n";
		return -1;
	}
	setTag(I1(0));
	int num_mat = I1(1);
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

	// receive double data
	static Vector D1;
	D1.resize(num_mat * 2 /*ints for mats*/ + num_mat /*floats for mats*/ + 90 /*other fixed data*/);
	int counter = 0;
	if (theChannel.recvVector(getDbTag(), commitTag, D1) < 0) {
		opserr << "Parallel3DMaterial::recvSelf() - failed to receive data (D1)\n";
		return -1;
	}
	/*ints for mats*/
	counter = 2 * num_mat; // left for later
	/*floats for mats*/
	for (int i = 0; i < num_mat; ++i) {
		m_weights[static_cast<std::size_t>(i)] = D1(counter++);
	}
	/*other fixed data*/
	for (int i = 0; i < 6; ++i)
		m_strain(i) = D1(counter++);
	for (int i = 0; i < 6; ++i)
		m_strain_commit(i) = D1(counter++);
	for (int i = 0; i < 6; ++i)
		m_stress(i) = D1(counter++);
	for (int i = 0; i < 6; ++i)
		for (int j = 0; j < 6; ++j)
			m_tangent(i, j) = D1(counter++);
	for (int i = 0; i < 6; ++i)
		for (int j = 0; j < 6; ++j)
			m_initial_tangent(i, j) = D1(counter++);

	// receive all materials
	counter = 0;
	for (int i = 0; i < num_mat; i++) {
		int mat_class_tag = D1(counter++);
		int mat_db_tag = D1(counter++);
		NDMaterial* the_material = theBroker.getNewNDMaterial(mat_class_tag);
		if (the_material == nullptr) {
			opserr << "Parallel3DMaterial::recvSelf() - could not get a new NDMaterial\n";
			return -1;
		}
		the_material->setDbTag(mat_db_tag);
		if (the_material->recvSelf(commitTag, theChannel, theBroker) < 0) {
			opserr << "Parallel3DMaterial::recvSelf() - failed to receive material\n";
			return -1;
		}
		m_materials[static_cast<std::size_t>(i)] = the_material;
	}

	// done
	return 0;
}

int Parallel3DMaterial::setParameter(const char** argv, int argc, Parameter& param)
{
	// forward to the sub-materials
	int res = -1;
	for (NDMaterial* item : m_materials) {
		if (item->setParameter(argv, argc, param) == 0)
			res = 0;
	}
	return res;
}

Response* Parallel3DMaterial::setResponse(const char** argv, int argc, OPS_Stream& output)
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
				std::shared_ptr<Parallel3DUtils::ResponseWrapper> wres;
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
					wres = std::make_shared<Parallel3DUtils::ResponseWrapper>(m_materials.size());
					wres_id = m_response_map.size() > 0 ? m_response_map.rbegin()->first + 1 : 1000;
					wres->name = name;
					int response_size = 0;
					for (std::size_t i = 0; i < m_materials.size(); ++i) {
						NDMaterial* imat = m_materials[i];
						Parallel3DUtils::CustomStream ds;
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

int Parallel3DMaterial::getResponse(int responseID, Information& matInformation)
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

void Parallel3DMaterial::computeInitialTangent()
{
	m_initial_tangent.Zero();
	for (std::size_t i = 0; i < m_materials.size(); ++i) {
		m_initial_tangent.addMatrix(1.0, m_materials[i]->getInitialTangent(), m_weights[i]);
	}
}

void Parallel3DMaterial::computeTangent()
{
	m_tangent.Zero();
	for (std::size_t i = 0; i < m_materials.size(); ++i) {
		m_tangent.addMatrix(1.0, m_materials[i]->getTangent(), m_weights[i]);
	}
}

void Parallel3DMaterial::computeStress()
{
	m_stress.Zero();
	for (std::size_t i = 0; i < m_materials.size(); ++i) {
		m_stress.addVector(1.0, m_materials[i]->getStress(), m_weights[i]);
	}
}

