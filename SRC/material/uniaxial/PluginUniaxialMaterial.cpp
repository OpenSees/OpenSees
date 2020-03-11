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
#include <string>
#include <exception>
#include <limits>

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

	// allocate the PluginMaterialData. now material data has the tag, the pointer to
	// the procedure, and all data properly inizialized to zero
	PluginMaterialData* data = PluginFramework::instance().makeMaterialData(descriptor->procedure, tag);
	if (data == 0) {
		opserr << "PluginUniaxialMaterial Error: Failed to initialize material data\n";
		return 0;
	}

	// call for initialization. Here the function must set initialization info for proper allocation
	// (mat_type, n_param, n_state,  message)
	int rcode = 0;
	data->proc(data, PF_MAT_INITIALIZE, &rcode);
	if (rcode != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in job = PF_MAT_INITIALIZE\n";
		delete data;
		return 0;
	}
	if (data->mat_type != PF_MAT_TYPE_UNIAXIAL) {
		opserr << 
			"PluginUniaxialMaterial Error: Failed in job = PF_MAT_INITIALIZE.\n"
			"mat_type (" << data->mat_type << ") not set correctly to PF_MAT_TYPE_UNIAXIAL (" << PF_MAT_TYPE_UNIAXIAL << ").\n"
			"Make sure the plugin material you selected is a valid uniaxial material plugin";
		delete data;
		return 0;
	}

	// now we can allocate data for input arguments and state variables
	if (PluginFramework::instance().allocateData(data) != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed to allocate material data\n";
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
		PluginFramework::instance().releaseData(data);
		delete data;
		return 0;
	}
	if (descriptor->parseMessage(data->message) != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed to parse message.\n";
		PluginFramework::instance().releaseData(data);
		delete data;
		return 0;
	}
	data->message = 0; // don't need it anymore, the plugin is free to do whatever it likes
	if (descriptor->arguments.size() != static_cast<std::size_t>(data->n_param)) {
		opserr << "PluginUniaxialMaterial Error: n_param does not match the number of arguments (A records) in message.\n";
		PluginFramework::instance().releaseData(data);
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
{
}

PluginUniaxialMaterial::PluginUniaxialMaterial(PluginMaterialDescriptor* descr, PluginMaterialData* d)
	: UniaxialMaterial(d->tag, MAT_TAG_PluginUniaxialMaterial)
	, m_descriptor(descr)
	, m_data(d)
{
	// call revert to start to inizialize variable now that are properly allocated
	if (revertToStart() != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed in revertToStart for first initialization\n";
		exit(-1);
	}
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
	// PF_TODO
	return -1;
}

int PluginUniaxialMaterial::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	// PF_TODO
	return -1;
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

	// call for initialization (PF_MAT_INITIALIZE) avoided here. Can copy those values
	data->mat_type = m_data->mat_type;
	data->n_param = m_data->n_param;

	// now we can allocate data for input arguments and state variables
	if (PluginFramework::instance().allocateData(data) != 0) {
		opserr << "PluginUniaxialMaterial Error: Failed to allocate material data (in getCopy)\n";
		delete data;
		return 0;
	}

	// copy data param and state
	for (int i = 0; i < m_data->n_param; i++)
		data->param[i] = m_data->param[i];

	// create the new material with the current data, it will complete 
	// the setup of pointers in data
	return new PluginUniaxialMaterial(m_descriptor, data);
}

Response* PluginUniaxialMaterial::setResponse(const char** argv, int argc, OPS_Stream& theOutputStream)
{
	// PF_TODO
	return 0;
}

int PluginUniaxialMaterial::getResponse(int responseID, Information& matInformation)
{
	// PF_TODO
	return -1;
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
