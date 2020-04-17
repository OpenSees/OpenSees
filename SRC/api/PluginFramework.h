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

#ifndef PluginFramework_h
#define PluginFramework_h

// Written: Massimo Petracca 
// Created: 02/2020
// Revision: A
//
// Description: This file contains the PluginFramework manager
// class, which is a singleton that handles the loading of
// plugins

#include <PluginFrameworkAPI.h>
#include <vector>
#include <set>
#include <map>
#include <sstream>
#include <limits>
#include <iomanip>


/**
The serializer utility used to serialize heterogeneous data in a string
used to send/receive copies of plugin wrappers to/from processes 
*/
#define fmt_double std::setprecision(std::numeric_limits<double>::digits10 + 1)
class PluginSerializer
{
private:
	std::stringstream ss;

public:
	struct darray_wrapper {
		double* x;
		size_t n;
		darray_wrapper(double* _x, std::size_t _n)
			: x(_x), n(_n)
		{}
	};

public:
	PluginSerializer() {}
	PluginSerializer(const std::string& x)
		: ss(x) {}
	inline std::string str()const {
		return ss.str();
	}
	inline std::string remaining() {
		auto pos = ss.tellg();
		ss.seekg(0, ss.end);
		auto n = ss.tellg() - pos - 1; // exclude last endline
		ss.seekg(pos);
		return ss.str().substr(ss.tellg(), n);
	}
	explicit operator bool() const {
		return !ss.fail();
	}
	bool operator!() const {
		return ss.fail();
	}

public:
	inline PluginSerializer& operator << (const char* x) {
		ss << x << '\n';
		return *this;
	}
	inline PluginSerializer& operator << (const std::string& x) {
		ss << x << '\n';
		return *this;
	}
	inline PluginSerializer& operator << (bool x) {
		ss << x << '\n';
		return *this;
	}
	inline PluginSerializer& operator << (int x) {
		ss << x << '\n';
		return *this;
	}
	inline PluginSerializer& operator << (std::size_t x) {
		ss << x << '\n';
		return *this;
	}
	inline PluginSerializer& operator << (double x) {
		ss << fmt_double << x << '\n';
		return *this;
	}
	inline PluginSerializer& operator << (const darray_wrapper& x) {
		if (x.n > 0) {
			ss << fmt_double << x.x[0];
			for (std::size_t i = 1; i < x.n; ++i)
				ss << " " << fmt_double << x.x[i];
		}
		ss << '\n';
		return *this;
	}

public:
	inline PluginSerializer& operator >> (std::string& x) {
		std::getline(ss, x, '\n');
		return *this;
	}
	inline PluginSerializer& operator >> (bool& x) {
		ss >> x;
		return *this;
	}
	inline PluginSerializer& operator >> (int& x) {
		ss >> x;
		return *this;
	}
	inline PluginSerializer& operator >> (std::size_t& x) {
		ss >> x;
		return *this;
	}
	inline PluginSerializer& operator >> (double& x) {
		ss >> x;
		return *this;
	}
	inline PluginSerializer& operator >> (darray_wrapper& x) {
		for (std::size_t i = 0; i < x.n; ++i) {
			ss >> x.x[i];
		}
		static std::string dummy;
		std::getline(ss, dummy, '\n');
		return *this;
	}
};

/**
Descriptor for an input argument
*/
class PluginArgumentDescriptor
{
public:
	PluginArgumentDescriptor();
	PluginArgumentDescriptor(const std::string& the_name);
	PluginArgumentDescriptor(const std::string &the_name, double def);

public:
	/// argument name
	std::string name;
	/// data for optional values
	bool is_optional;
	double default_value;
};

/**
Descriptor for a response
*/
class PluginResponseDescriptor
{
public:
	/// response id
	int id;
	/// response name
	std::string name;
	/// all components, its length defined the size of the Vector of results
	std::vector<std::string> components;
};

class PluginLibrary; // forward declaration

/**
Descriptor for a material plugin
*/
class PluginMaterialDescriptor
{
public:
	PluginMaterialDescriptor(PluginLibrary* lib, const std::string &name, PluginMaterialProc p);
	int parseMessage(const char* m);

public:
	/// pointer to the parent library
	PluginLibrary* library;
	/// procedure name
	std::string procedure_name;
	/// pointer to the plugin material procedure
	PluginMaterialProc procedure;
	/// input arguments
	std::vector<PluginArgumentDescriptor> arguments;
	/// responses
	std::vector<PluginResponseDescriptor> responses;
	/// parameters
	std::vector<PluginResponseDescriptor> parameters;
	/// variables
	std::vector<PluginResponseDescriptor> variables;
private:
	bool message_parsed;
};

/**
A simple class that contains the plugin library handle, and all the functions
loaded from the library mapped to their name
*/
class PluginLibrary
{
public:
	typedef std::map<std::string, PluginMaterialDescriptor*> material_map_t;

public:
	PluginLibrary(const std::string& name, void* lh);
	~PluginLibrary();

public:
	/// the library name
	std::string library_name;
	/// the library handle
	void* library_handle;
	/// all materials loaded from this library
	material_map_t materials;
};

/**
The PluginFramework is the manager class for loading plugin libraries and procedures.
It is implemented as a singleton. Everythin is handled internally.
The only public function is getMaterialDescriptor. It may throw an exception if it cannot load
the library or the function
*/
class PluginFramework
{
private:
	typedef std::map<std::string, PluginLibrary*> plugin_map_t;

private:
	PluginFramework();
	PluginFramework(const PluginFramework&);
	PluginFramework& operator = (const PluginFramework&);

public:
	~PluginFramework();
	static PluginFramework& instance();
	static const char* materialJobTypeToString(PluginMaterialJobType x);

public:
	///
	/// \brief getMaterialDescriptor returns a pointer to the requested PluginMaterialDescriptor
	/// \param library_name the plug-in library name
	/// \param function_name the plug-in function name name in the plug-in library
	/// \note throws a runtime_error exception in case something goes wrong while loading the library
	///
	PluginMaterialDescriptor* getMaterialDescriptor(const std::string& library_name, const std::string& function_name);
	///
	/// \brief makeMaterialData creates a new material data
	/// and initializes its fields to default values.
	/// \param p the plug-in function pointer
	/// \param tag the material tag
	///
	PluginMaterialData* makeMaterialData(PluginMaterialProc p, int tag);
	///
	/// \brief allocateData allocates the fields in d based on the material type specified in d
	/// during the first call to the material procedure with job = PF_MAT_INITIALIZE
	/// \param d the data
	///
	int allocateData(PluginMaterialData* d);
	///
	/// \brief releaseData destroys the fields in d previously allocated by a call to
	/// \ref allocateData. This method should be called after the call to the material
	/// procedure with job = PF_MAT_FINALIZE
	/// \param d the data
	///
	int releaseData(PluginMaterialData* d);
	///
	/// \brief parseTclCommand parses the TCL command according to the arguments in
	/// the descriptor, and places the obtained values in the data param vector
	///
	int parseTclCommand(const PluginMaterialDescriptor* descriptor, PluginMaterialData* d);

private:
	/// all the loaded plugins
	plugin_map_t plugins;
};

#endif // PluginFramework_h

