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

#include <PluginFramework.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <string>
#include <sstream>
#include <exception>
#include <iostream>

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <dlfcn.h>
#endif // defined(_WIN32)
#include <stdlib.h>
#include <sys/stat.h>

#define PG_VERBOSE

#ifdef PG_VERBOSE
#define pg_print(X) opserr << X
#else
#define pg_print(X)
#endif

#define TRY_ALLOC(BLOCK) \
	try { \
		BLOCK; \
	} \
	catch(const std::exception &ex) { \
		opserr << "Failed to allocate (" << ex.what() << ")\n"; \
		return -1; \
	} \
	catch(...) { \
		opserr << "Failed to allocate (...)\n"; \
		return -1; \
	}

#define ENUM_STR(X) #X
#define CASE_ENUM_STR(X) case X : return ENUM_STR(X)

namespace details_plugin_framework {

	/**
	common utilities
	*/

	class formatter
	{
		std::stringstream ss;
		formatter(const formatter&) {}
		formatter& operator = (const formatter&) { return *this;  }
	public:
		formatter() {}
		template<typename T>
		inline formatter& operator << (const T& x) {
			ss << x;
			return *this;
		}
		operator std::string() const { return ss.str(); }
	};

	namespace mem {

		template<typename T>
		inline int allocate(T*& x, std::size_t N) {
			if (N == 1) {
				TRY_ALLOC(x = new T);
				if (x == 0) {
					opserr << "Failed to allocate (...)\n";
					return -1;
				}
				(*x) = T(0);
			}
			else if (N > 1) {
				TRY_ALLOC(x = new T[N]);
				if (x == 0) {
					opserr << "Failed to allocate (...)\n";
					return -1;
				}
				for (std::size_t i = 0; i < N; i++)
					x[i] = T(0);
			}
			return 0;
		}

		template<typename T>
		inline void release(T*& x, std::size_t N) {
			if (N == 1) {
				delete x;
			}
			else if (N > 1) {
				delete[] x;
			}
			x = 0;
		}

	}

	namespace strings {

		inline std::string trim(const std::string& x, const std::string c)
		{
			size_t p0 = x.find_first_not_of(c);
			if (p0 == std::string::npos)
				return "";
			size_t p1 = x.find_last_not_of(c);
			return x.substr(p0, p1 + 1 - p0);
		}

		inline std::string trim(const std::string& x)
		{
			static const std::string c = " \t\r\n";
			return trim(x, c);
		}

		inline void split(const std::string& text, char sep, std::vector<std::string>& tokens, bool skip_empty = true, bool do_trim = true)
		{
			if (tokens.size() > 0) tokens.clear();
			std::size_t start = 0, end = 0;
			while (true)
			{
				end = text.find(sep, start);
				if (end == std::string::npos)
				{
					if (start < text.size())
						tokens.push_back(do_trim ? trim(text.substr(start)) : text.substr(start));
					else if (!skip_empty)
						tokens.push_back("");
					break;
				}
				std::string subs = text.substr(start, end - start);
				if (do_trim)
					subs = trim(subs);
				if (skip_empty)
				{
					if (subs.size() > 0)
						tokens.push_back(subs);
				}
				else
				{
					tokens.push_back(subs);
				}
				start = end + 1;
			}
		}

	}

	/**
	All utilities for loading libraries and symbols at runtime
	*/

	namespace libload {

		inline std::string prependLibPrefix(const std::string& lib_name) {
			return formatter() << "lib" << lib_name;
		}

		std::string appendSharedLibPostfix(const std::string& lib_name) {
#if defined(_WIN32)
			return formatter() << lib_name << ".dll";
#elif defined(_MACOSX)
			return formatter() << lib_name << ".dylib";
#else
			return formatter() << lib_name << ".so";
#endif
		}

		int unloadLib(void** libHandle) {

#if defined(_WIN32)

			HINSTANCE hLib = (HINSTANCE)(*libHandle);
			if (hLib != NULL) {
				if (!FreeLibrary((HMODULE)hLib)) {
					opserr << "cannot unload library\n";
					return -3;
				}
				(*libHandle) = NULL;
			}

#else

			if ((*libHandle) != NULL) {
				dlclose(*libHandle);
				(*libHandle) = NULL;
			}

#endif

			return 0;
		}

		int loadLib(const std::string& lib_name, void** libHandle) {

			int result = 0;
			*libHandle = NULL;

#if defined(_WIN32)

			HINSTANCE hLib = LoadLibraryA(lib_name.c_str());
			if (hLib != NULL) {
				*libHandle = (void*)hLib;
			}
			else {
				opserr << "cannot load library \"" << lib_name.c_str() << "\"\n";
				result = -1;
			}

#else

			* libHandle = dlopen(lib_name.c_str(), RTLD_NOW);
			char* error = dlerror();
			if (*libHandle == NULL) {
				opserr << "cannot load library \"" << lib_name.c_str() << "\"" << "\n";
				opserr << "internal error : " << error << "\n";
				result = -1; // no lib exists
			}

#endif

			return result;
		}

		int loadSym(void* libHandle, const std::string& sym_name, void** funcHandle) {

			int result = 0;
			*funcHandle = NULL;

#if defined(_WIN32)

			HINSTANCE hLib = (HINSTANCE)libHandle;
			if (hLib != NULL) {
				(*funcHandle) = (void*)GetProcAddress((HMODULE)hLib, sym_name.c_str());
				if (*funcHandle == NULL) {
					std::stringstream ss;
					ss << sym_name << "_";
					std::string sym_name_ = ss.str();
					(*funcHandle) = (void*)GetProcAddress((HMODULE)hLib, sym_name_.c_str());
				}
				if (*funcHandle == NULL) {
					opserr << "cannot load symbol \"" << sym_name.c_str() << "\"\n";
					result = -2;
				}
			}
			else {
				opserr << "cannot load symbol \"" << sym_name.c_str() << "\" libHandle is null\n";
				result = -1;
			}

#else

			if (libHandle != NULL) {
				(*funcHandle) = dlsym(libHandle, sym_name.c_str());
				char* error = dlerror();
				if (*funcHandle == NULL) {
					std::stringstream ss;
					ss << sym_name << "_";
					std::string sym_name_ = ss.str();
					(*funcHandle) = dlsym(libHandle, sym_name_.c_str());
					error = dlerror();
				}
				if (*funcHandle == NULL) {
					opserr << "cannot load symbol \"" << sym_name.c_str() << "\"" << "\n";
					opserr << "internal error: " << error << "\n";
					return -2;
				}
			}
			else {
				opserr << "cannot load symbol \"" << sym_name.c_str() << "\" libHandle is null" << "\n";
				return -1;
			}

#endif

			return result;
		}

	}

}

namespace details = details_plugin_framework;

/**
PluginArgumentDescriptor implementation
*/

PluginArgumentDescriptor::PluginArgumentDescriptor()
	: name("")
	, is_optional(false)
	, default_value(0.0)
{
}

PluginArgumentDescriptor::PluginArgumentDescriptor(const std::string& the_name)
	: name(the_name)
	, is_optional(false)
	, default_value(0.0)
{
}

PluginArgumentDescriptor::PluginArgumentDescriptor(const std::string& the_name, double def)
	: name(the_name)
	, is_optional(true)
	, default_value(def)
{
}

/**
PluginMaterialDescriptor implementation
*/

PluginMaterialDescriptor::PluginMaterialDescriptor(PluginLibrary* lib, const std::string& name, PluginMaterialProc p)
	: library(lib)
	, procedure_name(name)
	, procedure(p)
	, message_parsed(false)
{
}

int PluginMaterialDescriptor::parseMessage(const char* m)
{
	// do this only once
	if (message_parsed)
		return 0;
	message_parsed = true;

	// get records
	std::string mess(m);
	std::vector<std::string> records;
	details::strings::split(mess, ';', records);

	// processed data, avoid duplicates
	std::set<std::string> processed_arguments;
	std::set<std::string> processed_responses;
	std::set<std::string> processed_parameters;
	std::set<std::string> processed_variables;

	// process each record
	for (std::size_t i = 0; i < records.size(); i++) {
		const std::string irec = records[i];

		// get fields
		std::vector<std::string> fields;
		details::strings::split(irec, '|', fields);
		if (fields.size() == 0)
			continue;

		// process each field
		const std::string f0 = fields[0];
		if (f0 == "A") {

			// A|<name>|<default>
			if (fields.size() < 2) {
				opserr << "PluginMaterialDescriptor Error: the A record must have at least 1 field";
				return -1;
			}

			const std::string& name = fields[1];
			double val = 0.0;
			bool has_def_val = false;
			if (fields.size() > 2) {
				std::stringstream conv(fields[2]);
				if (conv >> val)
					has_def_val = true;
			}

			size_t nold = processed_arguments.size();
			processed_arguments.insert(name);
			if (nold == processed_arguments.size()) {
				opserr << "PluginMaterialDescriptor Error: the A records must have unique names";
				return -1;
			}

			if (has_def_val)
				arguments.push_back(PluginArgumentDescriptor(name, val));
			else
				arguments.push_back(PluginArgumentDescriptor(name));
			
		}
		else if (f0 == "R") {

			// R|<id>|<name>|<comp_1_name>|<comp_2_name>|...|<comp_N_name>
			if (fields.size() < 4) {
				opserr << "PluginMaterialDescriptor Error: the R record must have at least 3 fields (id, name, and at least 1 component)";
				return -1;
			}

			int id;
			if (!(std::stringstream(fields[1]) >> id)) {
				opserr << "PluginMaterialDescriptor Error: Cannot get the response id from the R record";
				return -1;
			}

			const std::string& name = fields[2];

			std::size_t nold = processed_responses.size();
			processed_responses.insert(name);
			if (nold == processed_responses.size()) {
				opserr << "PluginMaterialDescriptor Error: the R records must have unique names";
				return -1;
			}

			std::size_t ncomp = fields.size() - 3;

			PluginResponseDescriptor resp;
			resp.id = id;
			resp.name = name;
			resp.components.resize(ncomp);
			for (std::size_t j = 3; j < fields.size(); j++)
				resp.components[j-3] = fields[j];

			responses.push_back(resp);

		}
		else if (f0 == "P") {

			// P|<id>|<name>
			if (fields.size() != 3) {
				opserr << "PluginMaterialDescriptor Error: the P record must have 2 fields (id and name)";
				return -1;
			}

			int id;
			if (!(std::stringstream(fields[1]) >> id)) {
				opserr << "PluginMaterialDescriptor Error: Cannot get the response id from the P record";
				return -1;
			}

			const std::string& name = fields[2];

			std::size_t nold = processed_parameters.size();
			processed_parameters.insert(name);
			if (nold == processed_parameters.size()) {
				opserr << "PluginMaterialDescriptor Error: the P records must have unique names";
				return -1;
			}

			PluginResponseDescriptor resp;
			resp.id = id;
			resp.name = name;

			parameters.push_back(resp);

		}
		else if (f0 == "V") {

			// V|<id>|<name>|<n_comp>
			if (fields.size() != 4) {
				opserr << "PluginMaterialDescriptor Error: the V record must have 3 fields (id, name and number of components)";
				return -1;
			}

			int id;
			if (!(std::stringstream(fields[1]) >> id)) {
				opserr << "PluginMaterialDescriptor Error: Cannot get the response id from the V record";
				return -1;
			}

			const std::string& name = fields[2];

			std::size_t nold = processed_variables.size();
			processed_variables.insert(name);
			if (nold == processed_variables.size()) {
				opserr << "PluginMaterialDescriptor Error: the V records must have unique names";
				return -1;
			}

			std::size_t ncomp;
			if (!(std::stringstream(fields[3]) >> ncomp)) {
				opserr << "PluginMaterialDescriptor Error: Cannot get the response number of components from the V record";
				return -1;
			}
			if (ncomp < 1) {
				opserr << "PluginMaterialDescriptor Error: number of components for the V record must be at least 1";
				return -1;
			}

			PluginResponseDescriptor resp;
			resp.id = id;
			resp.name = name;
			resp.components.resize(ncomp);

			variables.push_back(resp);

		}
	}

	// done
	return 0;
}

/**
PluginLibrary implementation
*/

PluginLibrary::PluginLibrary(const std::string& name, void* lh)
	: library_name(name)
	, library_handle(lh)
{
	pg_print("Plugin Library \"" << library_name.c_str() << "\" correctly loaded\n");
}

PluginLibrary::~PluginLibrary()
{
	// release pointers to material descriptors
	for (material_map_t::iterator it = materials.begin(); it != materials.end(); ++it) 
		if (it->second) 
			delete it->second;
	materials.clear();

	// unload library
	if (library_handle) {
		if (details::libload::unloadLib((void**)&library_handle) == 0) {
			pg_print("Plugin Library \"" << library_name.c_str() << "\" correctly un-loaded\n");
		}
		else {
			pg_print("Plugin Library \"" << library_name.c_str() << "\" could not be un-loaded!\n");
		}
	}
	else {
		pg_print("Plugin Library \"" << library_name.c_str() << "\" deleted, no lib handle to un-load...\n");
	}
}

/**
PluginFramework implementation
*/

PluginFramework::PluginFramework()
{
	// it is a singleton so it will be called only once
	opserr <<
		"PluginFramework - Written by ASDEA Software Technology: M.Petracca, G.Camata\n"
		"ASDEA Software Technology: https://asdeasoft.net \n"
		"STKO (Scientific ToolKit for OpenSees): https://asdeasoft.net/stko/ \n"
		"If you use this tool, please cite us:\n"
		"Petracca, M., Candeloro, F., & Camata, G. (2017). \"STKO user manual\". ASDEA Software Technology, Pescara Italy.\n";
}

PluginFramework& PluginFramework::operator=(const PluginFramework&)
{
	return *this;
}

PluginFramework::~PluginFramework()
{
	// release all the loaded plugins
	for (plugin_map_t::iterator it = plugins.begin(); it != plugins.end(); ++it) {
		if (it->second) {
			delete it->second;
		}
	}
	plugins.clear();
}

PluginFramework& PluginFramework::instance()
{
	static PluginFramework _instance;
	return _instance;
}

const char* PluginFramework::materialJobTypeToString(PluginMaterialJobType x)
{
	switch (x)
	{
		CASE_ENUM_STR(PF_MAT_GET_INIT_INFO);
		CASE_ENUM_STR(PF_MAT_INITIALIZE);
		CASE_ENUM_STR(PF_MAT_FINALIZE);

		CASE_ENUM_STR(PF_MAT_COMMIT);
		CASE_ENUM_STR(PF_MAT_REVERT);
		CASE_ENUM_STR(PF_MAT_REVERT_TO_START);

		CASE_ENUM_STR(PF_MAT_SERIALIZE);
		CASE_ENUM_STR(PF_MAT_DESERIALIZE);

		CASE_ENUM_STR(PF_MAT_COMPUTE);

		CASE_ENUM_STR(PF_MAT_GET_STRAIN);
		CASE_ENUM_STR(PF_MAT_GET_STRAIN_RATE);
		CASE_ENUM_STR(PF_MAT_GET_STRESS);
		CASE_ENUM_STR(PF_MAT_GET_TANGENT);
		CASE_ENUM_STR(PF_MAT_GET_INITIAL_TANGENT);
		CASE_ENUM_STR(PF_MAT_GET_DAMP_TANGENT);
		CASE_ENUM_STR(PF_MAT_GET_RHO);
		CASE_ENUM_STR(PF_MAT_GET_ENERGY);
		CASE_ENUM_STR(PF_MAT_GET_IS_FAILED);
		CASE_ENUM_STR(PF_MAT_GET_RESPONSE);

		CASE_ENUM_STR(PF_MAT_GET_VARIABLE);
		CASE_ENUM_STR(PF_MAT_SET_VARIABLE);

		CASE_ENUM_STR(PF_MAT_UPDATE_PARAMETER);
		CASE_ENUM_STR(PF_MAT_ACTIVATE_PARAMETER);

		CASE_ENUM_STR(PF_MAT_GET_RESPONSE_SENSITIVITY);
		CASE_ENUM_STR(PF_MAT_GET_STRAIN_SENSITIVITY);
		CASE_ENUM_STR(PF_MAT_GET_STRESS_SENSITIVITY);
		CASE_ENUM_STR(PF_MAT_GET_TANGENT_SENSITIVITY);
		CASE_ENUM_STR(PF_MAT_GET_INITIAL_TANGENT_SENSITIVITY);
		CASE_ENUM_STR(PF_MAT_GET_DAMP_TANGENT_SENSITIVITY);
		CASE_ENUM_STR(PF_MAT_GET_RHO_SENSITIVITY);
		CASE_ENUM_STR(PF_MAT_COMMIT_SENSITIVITY);
	default:
		return "Unknown";
	}
}

PluginMaterialDescriptor* PluginFramework::getMaterialDescriptor(const std::string& library_name, const std::string& function_name)
{
	// find the library or laod it
	PluginLibrary* plugin_library = 0;
	{
		plugin_map_t::const_iterator it = plugins.find(library_name);
		if (it == plugins.end()) {

			// library not yet loaded
			void* library_handle = 0;
			std::string trial_name = details::libload::appendSharedLibPostfix(library_name);
			if (details::libload::loadLib(trial_name, &library_handle) != 0) {
				std::string trial_name_2 = details::libload::prependLibPrefix(trial_name);
				if (details::libload::loadLib(trial_name_2, &library_handle) != 0) {
					throw std::runtime_error(details::formatter()
						<< "PluginFramework Error: cannot load the specified library, tried both \""
						<< trial_name << "\" and \"" << trial_name_2 << "\"");
				}
			}
			if (library_handle == 0) {
				throw std::runtime_error("PluginFramework Error: library handle is null. This should never happen!");
			}

			plugin_library = new PluginLibrary(library_name, library_handle);
			plugins[library_name] = plugin_library;
			opserr << "loaded plugin library : " << library_name.c_str() << "\n";
		}
		else {

			// we have the library
			plugin_library = it->second;

		}
		if (plugin_library == 0) {
			throw std::runtime_error("PluginFramework Error: plugin library is a null pointer. This should never happen!");
		}
	}

	// find the material descriptor or create it and load the procedure
	PluginMaterialDescriptor* descriptor = 0;
	{
		PluginLibrary::material_map_t::const_iterator it = plugin_library->materials.find(function_name);
		if (it == plugin_library->materials.end()) {

			// function not yet loaded
			PluginMaterialProc plugin_function;
			if (details::libload::loadSym(plugin_library->library_handle, function_name, (void**)&plugin_function) != 0) {
				throw std::runtime_error(details::formatter()
					<< "PluginFramework Error: cannot load the specified function \""
					<< function_name << "\"");
			}
			if (plugin_function == 0) {
				throw std::runtime_error("PluginFramework Error: function handle is null. This should never happen!");
			}
			descriptor = new PluginMaterialDescriptor(plugin_library, function_name, plugin_function);

			plugin_library->materials[function_name] = descriptor;
		}
		else {

			// we have the function
			descriptor = it->second;

		}
	}

	// done
	return descriptor;
}

PluginMaterialData* PluginFramework::makeMaterialData(PluginMaterialProc p, int tag)
{
	PluginMaterialData* d = new PluginMaterialData();
	d->proc = p;

	// opaque pointer to custom data
	d->custom = 0;

	// description
	d->message = 0;
	d->tag = tag;
	d->mat_type = 0;
	d->n_param = 0;
	
	// output parameter/response/variable
	d->response_id = 0;
	d->response_size = 0;
	d->response = 0;

	// for sensitivity
	d->sens_strain_gradient = 0;
	d->sens_grad_index = 0;
	d->sens_conditonal = 0;
	d->sens_num_grad = 0;
	d->sens_unused_placeholder = 0;

	// parameters and state variables
	d->param = 0;

	// input
	d->strain = 0;
	d->strain_rate = 0;
	d->temperature = 0;
	d->lch = 1.0;
	d->dT = 0.0;

	return d;
}

int PluginFramework::allocateData(PluginMaterialData* d)
{
	if (details::mem::allocate(d->param, d->n_param) != 0) return -1;
	return 0;
}

int PluginFramework::releaseData(PluginMaterialData* d)
{
	details::mem::release(d->param, d->n_param);
	return 0;
}

int PluginFramework::parseTclCommand(const PluginMaterialDescriptor* descriptor, PluginMaterialData* d)
{
	// initial checks
	if (descriptor == 0) {
		opserr << "PluginFramework - parseTclCommand Error: material descriptor is null\n";
		return -1;
	}
	if (d == 0) {
		opserr << "PluginFramework - parseTclCommand Error: material data is null\n";
		return -1;
	}

	// make sure the number of parameters is even
	// because we need name-value pairs
	int numArgs = OPS_GetNumRemainingInputArgs();
	if ((numArgs % 2) != 0) {
		opserr << "PluginFramework - parseTclCommand Error: the number of arguments must be even (i.e. N name-value pairs)\n";
		return -1;
	}

	// parse input arguments
	std::map<std::string, double> input_args;
	for (int i = 0; i < numArgs / 2; i++) {

		std::string name = details::strings::trim(OPS_GetString(), "-");
		if (name.empty()) {
			opserr << "PluginFramework - parseTclCommand Error: the name of an argument cannot be empty\n";
			return -1;
		}

		int n_data = 1;
		double value = 0.0;
		if (OPS_GetDouble(&n_data, &value) != 0) {
			opserr 
				<< "PluginFramework - parseTclCommand Error: Cannot obtain a double-precision number for the argument -"
				<< name.c_str() << "\n";
			return -1;
		}

		size_t n = input_args.size();
		input_args[name] = value;
		if (n == input_args.size()) {
			opserr
				<< "PluginFramework - parseTclCommand Error: The argument -"
				<< name.c_str() << " was specified multiple times.\n";
			return -1;
		}
	}

	// now fill the param vector in data according to the values provided by the user
	// and the optional ones for which we have default values.
	// here we check that all mandatory data have been correctly specified.
	for (std::size_t i = 0; i < descriptor->arguments.size(); i++) {
		const PluginArgumentDescriptor& iarg = descriptor->arguments[i];
		std::map<std::string, double>::const_iterator it = input_args.find(iarg.name);
		if (it != input_args.end()) {
			// the user specified this value
			d->param[i] = it->second;
		}
		else {
			// the user did not specify this value. let's see if it was optional
			if (iarg.is_optional) {
				d->param[i] = iarg.default_value;
			}
			else {
				opserr
					<< "PluginFramework - parseTclCommand Error: The argument -"
					<< iarg.name.c_str() << " is mandatory, please define this argument.\n";
				return -1;
			}
		}
	}

	// done
	return 0;
}

