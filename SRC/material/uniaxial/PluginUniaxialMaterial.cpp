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

// Written: Massimo Petracca 
// Created: 02/2020
// Revision: A

#include <PluginUniaxialMaterial.h>
#include <PluginFramework.h>
#include <elementAPI.h>
#include <Element.h>
#include <ID.h>
#include <Channel.h>
#include <MaterialResponse.h>
#include <actor/message/Message.h>
#include <string>
#include <exception>
#include <limits>

namespace details {

	inline PluginMaterialDescriptor* getDescriptor(const std::string& plugin_library, const std::string& plugin_function)
	{
		PluginMaterialDescriptor* descriptor = 0;
		try {
			descriptor = PluginFramework::instance().getMaterialDescriptor(plugin_library, plugin_function);
		}
		catch (const std::exception & ex) {
			opserr << "PluginUniaxialMaterial internal error while getting the requested procedure:\n"
				<< ex.what() << "\n";
			return 0;
		}
		catch (...) {
			opserr << "PluginUniaxialMaterial unknown exception while getting the requested procedure:\n";
			return 0;
		}
		return descriptor;
	}

}

void*
OPS_PluginUniaxialMaterial(void)
{
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 3) { // check for mandatory arguments (at least 3)
		opserr << "Want: uniaxialMaterial plugin $tag   $pluginLibrary $pluginFunction    <+ plugin dependent arguments>\n";
		return 0;
	}

	int idata[1];
	int num_data;

	// get tag
	num_data = 1;
	if (OPS_GetInt(&num_data, idata) != 0) {
		opserr << "WARNING invalid integer tag: uniaxialMaterial plugin \n";
		return 0;
	}
	int tag = idata[0];

	// get library and function names
	std::string plugin_library = OPS_GetString();
	std::string plugin_function = OPS_GetString();
	
	// get plugin material descriptor
	PluginMaterialDescriptor* descriptor = details::getDescriptor(plugin_library, plugin_function);
	if (descriptor == 0)
		return 0;

	// allocate the PluginMaterialData. now material data has the tag, the pointer to
	// the procedure, and all data properly inizialized to zero
	PluginMaterialData* data = PluginFramework::instance().makeMaterialData(descriptor->procedure, tag);
	if (data == 0) {
		opserr << "PluginUniaxialMaterial Error: Failed to initialize material data\n";
		return 0;
	}
	data->mat_type = PF_MAT_TYPE_UNIAXIAL;

	// call to get initialization info. Here the function must set initialization info for proper allocation
	// (message)
	int rcode = 0;
	data->proc(data, PF_MAT_GET_INIT_INFO, &rcode);
	if (rcode != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_GET_INIT_INFO\n";
		delete data;
		return 0;
	}
	if (data->mat_type != PF_MAT_TYPE_UNIAXIAL) {
		opserr << 
			"PluginUniaxialMaterial Error: Failed in job = PF_MAT_GET_INIT_INFO.\n"
			"mat_type (" << data->mat_type << ") not set correctly to PF_MAT_TYPE_UNIAXIAL (" << PF_MAT_TYPE_UNIAXIAL << ").\n"
			"Make sure the plugin material you selected is a valid uniaxial material plugin";
		delete data;
		return 0;
	}

	// now we can parse the message and get information about this material plugin.
	// information such as the input arguments, the available results/parameters/variables.
	// note: do it only once for each plugin material type
	if (data->message == 0) {
		opserr << 
			"PluginUniaxialMaterial Error: message field is null.\n"
			"Please you need to fill the message field to parse input arguments\n";
		delete data;
		return 0;
	}
	if (descriptor->parseMessage(data->message) != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed to parse message.\n";
		delete data;
		return 0;
	}
	data->message = 0; // don't need it anymore, the plugin is free to do whatever it likes

	// get number of arguments for next allocation
	data->n_param = static_cast<int>(descriptor->arguments.size()); 

	// now we can allocate data for input arguments
	if (PluginFramework::instance().allocateData(data) != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed to allocate material data\n";
		delete data;
		return 0;
	}

	// we correctly parsed the message. now we can parse the tcl command
	if (PluginFramework::instance().parseTclCommand(descriptor, data) != 0) {
		opserr << "PluginUniaxialMaterial Error: cannot parse the TCL command.\n";
		PluginFramework::instance().releaseData(data);
		delete data;
		return 0;
	}

	// call for initialization. Here the function must initialize the sate of the plugin material
	rcode = 0;
	data->proc(data, PF_MAT_INITIALIZE, &rcode);
	if (rcode != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_INITIALIZE\n";
		delete data;
		return 0;
	}

	// create the new material with the current data, it will complete 
	// the setup of pointers in data
	PluginUniaxialMaterial* material = new PluginUniaxialMaterial(descriptor, data);

	// done
	return material;
}

PluginUniaxialMaterial::PluginUniaxialMaterial()
	: UniaxialMaterial(0, MAT_TAG_PluginUniaxialMaterial)
	, m_descriptor(0)
	, m_data(0)
	, m_lch(1.0)
	, m_lch_calculated(false)
{
}

PluginUniaxialMaterial::PluginUniaxialMaterial(PluginMaterialDescriptor* descr, PluginMaterialData* d)
	: UniaxialMaterial(d->tag, MAT_TAG_PluginUniaxialMaterial)
	, m_descriptor(descr)
	, m_data(d)
	, m_lch(1.0)
	, m_lch_calculated(false)
{
	// save a pointer to lch in m_data
	m_data->lch = &m_lch;
}

PluginUniaxialMaterial::~PluginUniaxialMaterial()
{
	int rcode = 0;
	m_data->proc(m_data, PF_MAT_FINALIZE, &rcode);
	if (rcode != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_FINALIZE\n";
	}
	PluginFramework::instance().releaseData(m_data);
	delete m_data;
}

void PluginUniaxialMaterial::Print(OPS_Stream& s, int flag)
{
	s << "PluginUniaxialMaterial - Tag: " << getTag() << "\n";
}

const char* PluginUniaxialMaterial::getClassType() const
{
	return "PluginUniaxialMaterial";
}

int PluginUniaxialMaterial::sendSelf(int commitTag, Channel& theChannel)
{
	// serialize to string
	int rcode = 0;
	m_data->proc(m_data, PF_MAT_SERIALIZE, &rcode);
	if (rcode != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_SERIALIZE\n";
		return -1;
	}
	if (m_data->message == 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_SERIALIZE\n";
		return -1;
	}

	// copy message in a local string and reset the original pointer to 0
	std::string msg_copy(m_data->message);
	m_data->message = 0;

	// get message size
	int msg_data_size = static_cast<int>(msg_copy.size() + 1);

	// send message size
	ID idata(1);
	idata(0) = msg_data_size;
	if (theChannel.sendID(0, commitTag, idata) < 0) {
		opserr << "PluginUniaxialMaterial::sendSelf() - failed to send message size\n";
		return -1;
	}

	// send message
	Message msg((char*)msg_copy.c_str(), msg_data_size);
	if (theChannel.sendMsg(0, commitTag, msg) < 0) {
		opserr << "PluginUniaxialMaterial::sendSelf() - failed to send message\n";
		return -1;
	}

	// done
	return 0;
}

int PluginUniaxialMaterial::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	// recv message size
	ID idata(1);
	if (theChannel.recvID(0, commitTag, idata) < 0) {
		opserr << "PluginUniaxialMaterial::recvSelf() - failed to recv message size\n";
		return -1;
	}
	int msg_data_size = idata(0);

	// recv message
	std::vector<char> msg_data_chars(msg_data_size);
	Message msg(msg_data_chars.data(), msg_data_size);
	if (theChannel.recvMsg(0, commitTag, msg) < 0) {
		opserr << "PluginUniaxialMaterial::recvSelf() - failed to recv message\n";
		return -1;
	}

	// set the message pointer
	m_data->message = msg_data_chars.data();

	// deserialize from string
	int rcode = 0;
	m_data->proc(m_data, PF_MAT_DESERIALIZE, &rcode);
	if (rcode != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_DESERIALIZE\n";
		return -1;
	}
	
	// reset the original pointer to 0
	m_data->message = 0;

	// done
	return 0;
}

int PluginUniaxialMaterial::setParameter(const char** argv, int argc, Parameter& param)
{
	// PF_TODO
	return -1;
}

int PluginUniaxialMaterial::updateParameter(int parameterID, Information& info)
{
	// PF_TODO
	return 0;
}

int PluginUniaxialMaterial::activateParameter(int parameterID)
{
	// PF_TODO
	return -0;
}

int PluginUniaxialMaterial::setVariable(const char* variable, Information&)
{
	// PF_TODO
	return -1;
}

int PluginUniaxialMaterial::getVariable(const char* variable, Information&)
{
	// PF_TODO
	return -1;
}

int PluginUniaxialMaterial::setTrialStrain(double strain, double strainRate)
{
	// compute lch only once
	if (!m_lch_calculated && (ops_TheActiveElement != 0)) {
		m_lch = ops_TheActiveElement->getCharacteristicLength();
		m_lch_calculated = true;
	}
	// set
	m_data->strain = &strain;
	m_data->strain_rate = &strainRate;
	// call
	int rcode = 0;
	m_data->proc(m_data, PF_MAT_COMPUTE, &rcode);
	// unset
	m_data->strain = 0;
	m_data->strain_rate = 0;
	// check
	if (rcode != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_COMPUTE\n";
	}
	// done
	return rcode;
}

int PluginUniaxialMaterial::setTrialStrain(double strain, double temperature, double strainRate)
{
	// set
	m_data->temperature = &temperature;
	// call
	int res = setTrialStrain(strain, strainRate);
	// unset
	m_data->temperature = 0;
	// done
	return res;
}

double PluginUniaxialMaterial::getStrain()
{
	// set
	double res;
	m_data->response = &res;
	// call
	int rcode = 0;
	m_data->proc(m_data, PF_MAT_GET_STRAIN, &rcode);
	// unset
	m_data->response = 0;
	// check
	if (rcode != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_GET_STRAIN\n";
		exit(-1);
	}
	// done
	return res;
}

double PluginUniaxialMaterial::getStrainRate()
{
	// set
	double res;
	m_data->response = &res;
	// call
	int rcode = 0;
	m_data->proc(m_data, PF_MAT_GET_STRAIN_RATE, &rcode);
	// unset
	m_data->response = 0;
	// check
	if (rcode != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_GET_STRAIN_RATE\n";
		exit(-1);
	}
	// done
	return res;
}

double PluginUniaxialMaterial::getStress()
{
	// set
	double res;
	m_data->response = &res;
	// call
	int rcode = 0;
	m_data->proc(m_data, PF_MAT_GET_STRESS, &rcode);
	// unset
	m_data->response = 0;
	// check
	if (rcode != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_GET_STRESS\n";
		exit(-1);
	}
	// done
	return res;
}

double PluginUniaxialMaterial::getTangent()
{
	// set
	double res;
	m_data->response = &res;
	// call
	int rcode = 0;
	m_data->proc(m_data, PF_MAT_GET_TANGENT, &rcode);
	// unset
	m_data->response = 0;
	// check
	if (rcode != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_GET_TANGENT\n";
		exit(-1);
	}
	// done
	return res;
}

double PluginUniaxialMaterial::getInitialTangent()
{
	// set
	double res;
	m_data->response = &res;
	// call
	int rcode = 0;
	m_data->proc(m_data, PF_MAT_GET_INITIAL_TANGENT, &rcode);
	// unset
	m_data->response = 0;
	// check
	if (rcode != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_GET_INITIAL_TANGENT\n";
		exit(-1);
	}
	// done
	return res;
}

double PluginUniaxialMaterial::getDampTangent()
{
	// set
	double res;
	m_data->response = &res;
	// call
	int rcode = 0;
	m_data->proc(m_data, PF_MAT_GET_DAMP_TANGENT, &rcode);
	// unset
	m_data->response = 0;
	// check
	if (rcode != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_GET_DAMP_TANGENT\n";
		exit(-1);
	}
	// done
	return res;
}

double PluginUniaxialMaterial::getRho()
{
	// set
	double res;
	m_data->response = &res;
	// call
	int rcode = 0;
	m_data->proc(m_data, PF_MAT_GET_RHO, &rcode);
	// unset
	m_data->response = 0;
	// check
	if (rcode != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_GET_RHO\n";
		exit(-1);
	}
	// done
	return res;
}

int PluginUniaxialMaterial::commitState()
{
	// call
	int rcode = 0;
	m_data->proc(m_data, PF_MAT_COMMIT, &rcode);
	if (rcode != 0) 
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_COMMIT\n";
	return rcode;
}

int PluginUniaxialMaterial::revertToLastCommit()
{
	// call
	int rcode = 0;
	m_data->proc(m_data, PF_MAT_REVERT, &rcode);
	if (rcode != 0)
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_REVERT\n";
	return rcode;
}

int PluginUniaxialMaterial::revertToStart()
{
	// call
	int rcode = 0;
	m_data->proc(m_data, PF_MAT_REVERT_TO_START, &rcode);
	if (rcode != 0)
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_REVERT_TO_START\n";
	return rcode;
}

UniaxialMaterial* PluginUniaxialMaterial::getCopy()
{
	// allocate the PluginMaterialData. now material data has the tag, the pointer to
	// the procedure, and all data properly inizialized to zero
	PluginMaterialData* data = PluginFramework::instance().makeMaterialData(m_descriptor->procedure, getTag());
	if (data == 0) {
		opserr << "PluginUniaxialMaterial Error: Failed to initialize material data (in getCopy)\n";
		return 0;
	}

	// call for initialization info (PF_MAT_GET_INIT_INFO) avoided here. Can copy those values
	data->mat_type = m_data->mat_type;
	data->n_param = m_data->n_param;

	// now we can allocate data for input arguments and state variables
	if (PluginFramework::instance().allocateData(data) != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed to allocate material data (in getCopy)\n";
		delete data;
		return 0;
	}

	// copy data param and state. parseTclCommand avoided here. Can copy those values
	for (int i = 0; i < m_data->n_param; i++)
		data->param[i] = m_data->param[i];

	// call for initialization. Here the function must initialize the sate of the plugin material
	int rcode = 0;
	data->proc(data, PF_MAT_INITIALIZE, &rcode);
	if (rcode != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_INITIALIZE\n";
		delete data;
		return 0;
	}

	// create the new material with the current data, it will complete 
	// the setup of pointers in data
	PluginUniaxialMaterial* material = new PluginUniaxialMaterial(m_descriptor, data);

	// done
	return material;
}

Response* PluginUniaxialMaterial::setResponse(const char** argv, int argc, OPS_Stream& theOutputStream)
{
	if (argc < 1)
		return 0;

	std::string iarg(argv[0]);
	for (std::size_t i = 0; i < m_descriptor->responses.size(); i++) {
		const PluginResponseDescriptor& d = m_descriptor->responses[i];
		if ((d.name == iarg) && (d.components.size() > 0)) {

			theOutputStream.tag("NdMaterialOutput");
			theOutputStream.attr("matType", this->getClassType());
			theOutputStream.attr("matTag", this->getTag());

			for(std::size_t j = 0; j < d.components.size(); j++)
				theOutputStream.tag("ResponseType", d.components[j].c_str());

			if (d.components.size() == 1) {
				return new MaterialResponse(this, d.id, 0.0);
			}
			else {
				return new MaterialResponse(this, d.id, Vector(static_cast<int>(d.components.size())));
			}

			theOutputStream.endTag();
		}
	}
	
	// call the base class implementation
	return UniaxialMaterial::setResponse(argv, argc, theOutputStream);
}

int PluginUniaxialMaterial::getResponse(int responseID, Information& matInformation)
{
	for (std::size_t i = 0; i < m_descriptor->responses.size(); i++) {
		const PluginResponseDescriptor& d = m_descriptor->responses[i];
		if (d.id == responseID) {

			// set
			m_data->response_id = responseID;
			double value_scalar = 0.0;
			Vector value_vector;
			if (d.components.size() == 1) {
				m_data->response = &value_scalar;
			}
			else {
				value_vector.resize(d.components.size());
				m_data->response = &value_vector(0);
			}

			// call
			int rcode = 0;
			m_data->proc(m_data, PF_MAT_GET_RESPONSE, &rcode);

			// unset
			m_data->response_id = 0;
			m_data->response = 0;

			// check
			if (rcode != 0) {
				opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_GET_RESPONSE\n";
				exit(-1);
			}

			// done
			if (d.components.size() == 1) 
				return matInformation.setDouble(value_scalar);
			else 
				return matInformation.setVector(value_vector);
		}
	}

	// call the base class implementation
	return UniaxialMaterial::getResponse(responseID, matInformation);
}

int PluginUniaxialMaterial::getResponseSensitivity(int responseID, int gradIndex, Information& info)
{
	// PF_TODO
	return -1;
}

bool PluginUniaxialMaterial::hasFailed()
{
	// set
	double res;
	m_data->response = &res;
	// call
	int rcode = 0;
	m_data->proc(m_data, PF_MAT_GET_IS_FAILED, &rcode);
	// unset
	m_data->response = 0;
	// check
	if (rcode != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_GET_IS_FAILED\n";
		exit(-1);
	}
	// done
	return std::abs(res) > std::numeric_limits<double>::epsilon();
}

double PluginUniaxialMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
	// PF_TODO
	return 0.0;
}

double PluginUniaxialMaterial::getStrainSensitivity(int gradIndex)
{
	// PF_TODO
	return 0.0;
}

double PluginUniaxialMaterial::getTangentSensitivity(int gradIndex)
{
	// PF_TODO
	return 0.0;
}

double PluginUniaxialMaterial::getInitialTangentSensitivity(int gradIndex)
{
	// PF_TODO
	return 0.0;
}

double PluginUniaxialMaterial::getDampTangentSensitivity(int gradIndex)
{
	// PF_TODO
	return 0.0;
}

double PluginUniaxialMaterial::getRhoSensitivity(int gradIndex)
{
	// PF_TODO
	return 0.0;
}

int PluginUniaxialMaterial::commitSensitivity(double strainGradient, int gradIndex, int numGrads)
{
	// PF_TODO
	return 0;
}

double PluginUniaxialMaterial::getEnergy()
{
	// set
	double res;
	m_data->response = &res;
	// call
	int rcode = 0;
	m_data->proc(m_data, PF_MAT_GET_ENERGY, &rcode);
	// unset
	m_data->response = 0;
	// check
	if (rcode != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_GET_ENERGY\n";
		exit(-1);
	}
	// done
	return res;
}
