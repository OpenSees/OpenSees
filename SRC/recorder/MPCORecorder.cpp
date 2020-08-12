/*************************************************************************************

notes, todo-list and workarounds

note 1:
some beams define a variable number of integration points. this is accounted for by
this recorder. The only requirement is that the element should provide the
response 'integrationPoints'.
Now DispBeamColumn(2d or 3d)WithSensitivity do not implement it, furthermore they
define their class tag in their own cpp file...
For those elements the variable number of integration points is not supported.
<update this node>

note 2:
gauss point / section index in the output from elem->setResponse are 1 based, while
indices for fibers are 0-based. In the mpco result file we make everything 0-based.

note 3:
the most expensive part is the writeSection method. However MPCORecored calls this
function only after a domain change. Some optimizations can be done....

note 4:
see $MP(2017/04/20)

$WO:SHELL_SEC_KEYWORD
workaround for shells, they used the "material" keyword for their section!
it would be better to use "section" as in beams, but we cannot modify those files...
\todo check the implementation of shell elements, if in future versions those elements
will be fixed and will handle the section keyword, remove this workaround

note 5:
NDMaterial in setResponse gives unknownStress for stress components if type == PlabeFiber!
todo: add auto-component naming in case of duplicated components!
in STKO components are assumed all different!

**************************************************************************************/

// some definitions

/* 
loads hdf5 shared library at runtime. if uncommented, hdf5 will be linked dynamically.
If uncommeted, then you need to define the macro H5_BUILT_AS_DYNAMIC_LIB to tell hdf5 that
we want to dynamic link (shared library).
Then we need to set the HDF5 include directory.
And finally for the linker: hdf5 or libhdf5
and path to HDF5 lib dir
*/
#define MPCO_HDF5_LOADED_AT_RUNTIME 

/* if hdf5 is loaded at runtime, this macro makes the process of loading hdf5 verbose */
#define MPCO_LIBLOADER_VERBOSE

/* max number of iterations to guess the number of cross sections in elements */
#define MPCO_MAX_TRIAL_NSEC 1000

/* max number of iterations to guess the number of fibers in cross sections */
#define MPCO_MAX_TRIAL_NFIB 1000000

//#define MPCO_WRITE_SECTION_IS_VERBOSE
//#define MPCO_WRITE_LOC_AX_IS_VERBOSE
//#define MPCO_TIMING

/*
enables SWMR (Single Writer - Multiple Readers) to allow reading this database from multiple processes
while opensees is writing data. Warning: this is a new feature in hdf5 version 1.10.0.
*/
#define MPCO_USE_SWMR

// opensees
#include "MPCORecorder.h"
#include "Channel.h"
#include "Domain.h"
#include "OPS_Globals.h"
#include "elementAPI.h"
#include "ID.h"
#include "Pressure_ConstraintIter.h"
#include "Pressure_Constraint.h"
#include "Node.h"
#include "NodeIter.h"
#include "Element.h"
#include "ElementIter.h"
#include "Vector.h"
#include "Matrix.h"
#include "Response.h"
#include "Message.h"
#include "CompositeResponse.h"
#include "section/SectionForceDeformation.h"
#include "MeshRegion.h"
// hdf5
#if !defined(MPCO_HDF5_LOADED_AT_RUNTIME)
// include hdf5 headers
#include "hdf5.h"
#include "hdf5_hl.h"
#else
// include system headers for dynamic loading of libraries
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
#endif // MPCO_HDF5_LOADED_AT_RUNTIME
// std
#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <limits>
#include <stdint.h>

/*************************************************************************************

macros

**************************************************************************************/

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// element tags not defined in an accessible .h file... \todo check them...
#ifndef ELE_TAG_FourNodeQuadWithSensitivity
#define ELE_TAG_FourNodeQuadWithSensitivity 100000011
#endif // !ELE_TAG_FourNodeQuadWithSensitivity
#ifndef ELE_TAG_DispBeamColumn2dWithSensitivity
#define ELE_TAG_DispBeamColumn2dWithSensitivity 102030  
#endif // !ELE_TAG_DispBeamColumn2dWithSensitivity
#ifndef ELE_TAG_DispBeamColumn3dWithSensitivity
#define ELE_TAG_DispBeamColumn3dWithSensitivity 1110000  
#endif // !ELE_TAG_DispBeamColumn3dWithSensitivity
#ifndef ELE_TAG_BbarBrickWithSensitivity
# define ELE_TAG_BbarBrickWithSensitivity 1984587234  
#endif // !ELE_TAG_BbarBrickWithSensitivity
#ifndef ELE_TAG_AC3D8HexWithSensitivity
#define ELE_TAG_AC3D8HexWithSensitivity 100001  
#endif // !ELE_TAG_AC3D8HexWithSensitivity

/*warning: the following 2 tags from ./ulBeamColumn are set to same value!!! this will not work properly!*/
#ifndef Ele_TAG_Elastic2dGNL
#define Ele_TAG_Elastic2dGNL -1  
#endif // !Ele_TAG_Elastic2dGNL
#ifndef TAG_InelasticYS2DGNL
#define TAG_InelasticYS2DGNL -1  
#endif // !TAG_InelasticYS2DGNL

#define OPS_STREAM_TAGS_MPCORecorder_ElementOutputDescriptorStream 1001 /** \todo should we move it to classTags file ? */

#ifdef MPCO_HDF5_LOADED_AT_RUNTIME

/*
some basic types used in HDF5
*/

typedef int64_t hid_t;
typedef int herr_t;
typedef unsigned long long 	hsize_t;
typedef int H5T_str_t; // enum (int) in hdf5
typedef unsigned int H5F_libver_t; // enum (uint) in hdf5
typedef unsigned int H5F_scope_t; // enum (uint) in hdf5

/*
HDF5 version info
*/

#define H5_VERS_MAJOR	1	/* For major interface/format changes  	     */
#define H5_VERS_MINOR	10	/* For minor interface/format changes  	     */
#define H5_VERS_RELEASE	1	/* For tweaks, bug-fixes, or development     */
#define H5_VERS_SUBRELEASE ""	/* For pre-releases like snap0       */
                              /* Empty string for real releases.           */
#define H5_VERS_INFO    "HDF5 library version: 1.10.1"      /* Full version string */

/*
cout wrapper for library loader verbosity
*/

#ifdef MPCO_LIBLOADER_VERBOSE
#define MPCO_LIBLOADER_COUT(X) std::cout << X
#else
#define MPCO_LIBLOADER_COUT(X)
#endif

/*
utilities for run-time loading of shared libraries
*/

namespace libload {

	std::string appendSharedLibPostfix(const std::string &lib_name) {
		std::stringstream ss;
		ss << lib_name;
#if defined(_WIN32)
		ss << ".dll";
#elif defined(_MACOSX)
		ss << ".dylib";
#else
		ss << ".so";
#endif
		return ss.str();
	}

	int unloadLib(void **libHandle) {

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

	int loadLib(const std::string &lib_name, void **libHandle) {

		int result = 0;
		*libHandle = NULL;

#if defined(_WIN32)

		HINSTANCE hLib = LoadLibraryA(lib_name.c_str());
		if (hLib != NULL) {
			*libHandle = (void *)hLib;
		}
		else {
			opserr << "cannot load library \"" << lib_name.c_str() << "\"\n";
			result = -1;
		}

#else

		*libHandle = dlopen(lib_name.c_str(), RTLD_NOW);
		char *error = dlerror();
		if (*libHandle == NULL) {
			opserr << "cannot load library \"" << lib_name.c_str() << "\"" << "\n";
			opserr << "internal error : " << error << "\n";
			result = -1; // no lib exists
		}

#endif

		return result;
	}

	int loadSym(void *libHandle, const std::string &sym_name, void **funcHandle) {

		int result = 0;
		*funcHandle = NULL;

#if defined(_WIN32)

		HINSTANCE hLib = (HINSTANCE)libHandle;
		if (hLib != NULL) {
			(*funcHandle) = (void *)GetProcAddress((HMODULE)hLib, sym_name.c_str());
			if (*funcHandle == NULL) {
				std::stringstream ss;
				ss << sym_name << "_";
				std::string sym_name_ = ss.str();
				(*funcHandle) = (void *)GetProcAddress((HMODULE)hLib, sym_name_.c_str());
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
			char *error = dlerror();
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

/*
macro for symbol loading
*/

#define MPCO_LIBLOADER_LOAD_SYM(X) \
if (libload::loadSym(LibraryLoader::lib_handle, #X, (void **)&ptr_##X) != 0) return; \
MPCO_LIBLOADER_COUT("loaded symbol \"" << #X << "\"." << "\n")

/*
Library Loader singleton.
Lazy loading of the HDF5 library and symbols.
loading is delayed to the first time LibraryLoader::instance() method is called.
the loaded library is then freed upon destruction of this singleton instance
*/

class LibraryLoader {
private:
	LibraryLoader()
		: loaded(false)
		, lib_handle(NULL)
	{
		std::string lib_name = libload::appendSharedLibPostfix("hdf5");
		MPCO_LIBLOADER_COUT("loading \"" << lib_name.c_str() << "\" ..." << "\n");

		if (libload::loadLib(lib_name, &lib_handle) != 0) {
			MPCO_LIBLOADER_COUT("\"" << lib_name.c_str() << "\" not found, trying with \"lib" << lib_name.c_str() << "\"\n");
			std::stringstream ss_for_prefix;
			ss_for_prefix << "lib" << lib_name;
			std::string lib_name_2 = ss_for_prefix.str();
			MPCO_LIBLOADER_COUT("loading \"" << lib_name_2.c_str() << "\" ..." << "\n");
			if (libload::loadLib(lib_name_2, &lib_handle) != 0) {
				opserr << "MPCORecorder error:\n"
					"Cannot load HDF5 (tried both \"" << lib_name.c_str() << "\" and \"" << lib_name_2.c_str() << "\" names for the shared library).\n"
					"HDF5 version " << H5_VERS_MAJOR << "." << H5_VERS_MINOR << " must be installed on your machine.\n"
					"pre-built binaries can be downloaded from: https://www.hdfgroup.org/downloads/hdf5/" << "\n";
				return;
			}
		}
		std::cout << "MPCORecorder: done loading HDF5\n";
		loaded = true;

		// functions
		MPCO_LIBLOADER_LOAD_SYM(H5check_version);
		MPCO_LIBLOADER_LOAD_SYM(H5open);
		MPCO_LIBLOADER_LOAD_SYM(H5Screate_simple);
		MPCO_LIBLOADER_LOAD_SYM(H5Sclose);
		MPCO_LIBLOADER_LOAD_SYM(H5Acreate2);
		MPCO_LIBLOADER_LOAD_SYM(H5Awrite);
		MPCO_LIBLOADER_LOAD_SYM(H5Aclose);
		MPCO_LIBLOADER_LOAD_SYM(H5Tcopy);
		MPCO_LIBLOADER_LOAD_SYM(H5Tset_size);
		MPCO_LIBLOADER_LOAD_SYM(H5Tset_strpad);
		MPCO_LIBLOADER_LOAD_SYM(H5Gcreate2);
		MPCO_LIBLOADER_LOAD_SYM(H5Gclose);
		MPCO_LIBLOADER_LOAD_SYM(H5Dcreate2);
		MPCO_LIBLOADER_LOAD_SYM(H5Dclose);
		MPCO_LIBLOADER_LOAD_SYM(H5Dwrite);
		MPCO_LIBLOADER_LOAD_SYM(H5Pcreate);
		MPCO_LIBLOADER_LOAD_SYM(H5Pclose);
		MPCO_LIBLOADER_LOAD_SYM(H5Pset_link_creation_order);
		MPCO_LIBLOADER_LOAD_SYM(H5Pset_libver_bounds);
		MPCO_LIBLOADER_LOAD_SYM(H5Fcreate);
		MPCO_LIBLOADER_LOAD_SYM(H5Fflush);
		MPCO_LIBLOADER_LOAD_SYM(H5Fclose);
#ifdef MPCO_USE_SWMR
		MPCO_LIBLOADER_LOAD_SYM(H5Fstart_swmr_write);
#endif // MPCO_USE_SWMR

		// static variables
		MPCO_LIBLOADER_LOAD_SYM(H5T_STD_I32LE_g);
		MPCO_LIBLOADER_LOAD_SYM(H5T_NATIVE_INT_g);
		MPCO_LIBLOADER_LOAD_SYM(H5T_IEEE_F64LE_g);
		MPCO_LIBLOADER_LOAD_SYM(H5T_NATIVE_DOUBLE_g);
		MPCO_LIBLOADER_LOAD_SYM(H5T_C_S1_g);
		MPCO_LIBLOADER_LOAD_SYM(H5P_CLS_FILE_CREATE_ID_g);
		MPCO_LIBLOADER_LOAD_SYM(H5P_CLS_FILE_ACCESS_ID_g);
		MPCO_LIBLOADER_LOAD_SYM(H5P_CLS_GROUP_CREATE_ID_g);
	}
	~LibraryLoader() {
		if (loaded) {
			MPCO_LIBLOADER_COUT("destroying lib handle: " << lib_handle << " ... ");
			libload::unloadLib(&lib_handle);
			MPCO_LIBLOADER_COUT(lib_handle << "\n");
		}
	}
	LibraryLoader(const LibraryLoader &other);
	LibraryLoader &operator = (const LibraryLoader &other);

public:
	static LibraryLoader &instance() {
		static LibraryLoader _instance;
		return _instance;
	}

public:
	// misc
	bool loaded;
	// lib handle
	void *lib_handle;
	// loaded function pointers
	herr_t (*ptr_H5check_version)(unsigned majnum, unsigned minnum, unsigned relnum);
	herr_t (*ptr_H5open)(void);
	hid_t  (*ptr_H5Screate_simple)(int rank, const hsize_t dims[], const hsize_t maxdims[]);
	herr_t (*ptr_H5Sclose)(hid_t space_id);
	hid_t  (*ptr_H5Acreate2)(hid_t loc_id, const char *attr_name, hid_t type_id, hid_t space_id, hid_t acpl_id, hid_t aapl_id);
	herr_t (*ptr_H5Awrite)(hid_t attr_id, hid_t type_id, const void *buf);
	herr_t (*ptr_H5Aclose)(hid_t attr_id);
	hid_t  (*ptr_H5Tcopy)(hid_t type_id);
	herr_t (*ptr_H5Tset_size)(hid_t type_id, size_t size);
	herr_t (*ptr_H5Tset_strpad)(hid_t type_id, H5T_str_t strpad);
	hid_t  (*ptr_H5Gcreate2)(hid_t loc_id, const char *name, hid_t lcpl_id, hid_t gcpl_id, hid_t gapl_id);
	herr_t (*ptr_H5Gclose)(hid_t group_id);
	hid_t  (*ptr_H5Dcreate2)(hid_t loc_id, const char *name, hid_t type_id, hid_t space_id, hid_t lcpl_id, hid_t dcpl_id, hid_t dapl_id);
	herr_t (*ptr_H5Dclose)(hid_t dset_id);
	herr_t (*ptr_H5Dwrite)(hid_t dset_id, hid_t mem_type_id, hid_t mem_space_id, hid_t file_space_id, hid_t plist_id, const void *buf);
	hid_t  (*ptr_H5Pcreate)(hid_t cls_id);
	herr_t (*ptr_H5Pclose)(hid_t plist_id);
	herr_t (*ptr_H5Pset_link_creation_order)(hid_t plist_id, unsigned crt_order_flags);
	herr_t (*ptr_H5Pset_libver_bounds)(hid_t plist_id, H5F_libver_t low, H5F_libver_t high);
	hid_t  (*ptr_H5Fcreate)(const char *filename, unsigned flags, hid_t create_plist, hid_t access_plist);
	herr_t (*ptr_H5Fflush)(hid_t object_id, H5F_scope_t scope);
	herr_t (*ptr_H5Fclose)(hid_t file_id);
#ifdef MPCO_USE_SWMR
	herr_t (*ptr_H5Fstart_swmr_write)(hid_t file_id);
#endif // MPCO_USE_SWMR

	// loaded global variable pointers
	hid_t *ptr_H5T_STD_I32LE_g;
	hid_t *ptr_H5T_NATIVE_INT_g;
	hid_t *ptr_H5T_IEEE_F64LE_g;
	hid_t *ptr_H5T_NATIVE_DOUBLE_g;
	hid_t *ptr_H5T_C_S1_g;
	hid_t *ptr_H5P_CLS_FILE_CREATE_ID_g;
	hid_t *ptr_H5P_CLS_FILE_ACCESS_ID_g;
	hid_t *ptr_H5P_CLS_GROUP_CREATE_ID_g;
};

/*
macros for loaded HDF5 functions
*/

#define H5check_version (*LibraryLoader::instance().ptr_H5check_version)
#define H5check()	H5check_version(H5_VERS_MAJOR,H5_VERS_MINOR, H5_VERS_RELEASE)

#define H5open (*LibraryLoader::instance().ptr_H5open)
#define H5OPEN  H5open(),
#define H5CHECK  H5check(),

#define H5Screate_simple (*LibraryLoader::instance().ptr_H5Screate_simple)
#define H5Sclose (*LibraryLoader::instance().ptr_H5Sclose)

#define H5Acreate2 (*LibraryLoader::instance().ptr_H5Acreate2)
#define H5Acreate H5Acreate2
#define H5Awrite (*LibraryLoader::instance().ptr_H5Awrite)
#define H5Aclose (*LibraryLoader::instance().ptr_H5Aclose)

#define H5Tcopy (*LibraryLoader::instance().ptr_H5Tcopy)
#define H5Tset_size (*LibraryLoader::instance().ptr_H5Tset_size)
#define H5Tset_strpad (*LibraryLoader::instance().ptr_H5Tset_strpad)

#define H5Gcreate2 (*LibraryLoader::instance().ptr_H5Gcreate2)
#define H5Gclose (*LibraryLoader::instance().ptr_H5Gclose)
#define H5Gcreate H5Gcreate2

#define H5Dcreate2 (*LibraryLoader::instance().ptr_H5Dcreate2)
#define H5Dclose (*LibraryLoader::instance().ptr_H5Dclose)
#define H5Dwrite (*LibraryLoader::instance().ptr_H5Dwrite)
#define H5Dcreate H5Dcreate2

#define H5Pcreate (*LibraryLoader::instance().ptr_H5Pcreate)
#define H5Pclose (*LibraryLoader::instance().ptr_H5Pclose)
#define H5Pset_link_creation_order (*LibraryLoader::instance().ptr_H5Pset_link_creation_order)
#define H5Pset_libver_bounds (*LibraryLoader::instance().ptr_H5Pset_libver_bounds)

#define H5Fcreate (*LibraryLoader::instance().ptr_H5Fcreate)
#define H5Fflush (*LibraryLoader::instance().ptr_H5Fflush)
#define H5Fclose (*LibraryLoader::instance().ptr_H5Fclose)
#ifdef MPCO_USE_SWMR
#define H5Fstart_swmr_write (*LibraryLoader::instance().ptr_H5Fstart_swmr_write)
#endif // MPCO_USE_SWMR

/*
macros for loaded HDF5 global variables
*/

#define H5T_STD_I32LE_g (*LibraryLoader::instance().ptr_H5T_STD_I32LE_g)
#define H5T_STD_I32LE (H5OPEN H5T_STD_I32LE_g)
#define H5T_NATIVE_INT_g (*LibraryLoader::instance().ptr_H5T_NATIVE_INT_g)
#define H5T_NATIVE_INT (H5OPEN H5T_NATIVE_INT_g)
#define H5T_IEEE_F64LE_g (*LibraryLoader::instance().ptr_H5T_IEEE_F64LE_g)
#define H5T_IEEE_F64LE (H5OPEN H5T_IEEE_F64LE_g)
#define H5T_NATIVE_DOUBLE_g (*LibraryLoader::instance().ptr_H5T_NATIVE_DOUBLE_g)
#define H5T_NATIVE_DOUBLE (H5OPEN H5T_NATIVE_DOUBLE_g)
#define H5T_C_S1_g (*LibraryLoader::instance().ptr_H5T_C_S1_g)
#define H5T_C_S1 (H5OPEN H5T_C_S1_g)
#define H5P_CLS_FILE_CREATE_ID_g (*LibraryLoader::instance().ptr_H5P_CLS_FILE_CREATE_ID_g)
#define H5P_FILE_CREATE (H5OPEN H5P_CLS_FILE_CREATE_ID_g)
#define H5P_CLS_FILE_ACCESS_ID_g (*LibraryLoader::instance().ptr_H5P_CLS_FILE_ACCESS_ID_g)
#define H5P_FILE_ACCESS (H5OPEN H5P_CLS_FILE_ACCESS_ID_g)
#define H5P_CLS_GROUP_CREATE_ID_g (*LibraryLoader::instance().ptr_H5P_CLS_GROUP_CREATE_ID_g)
#define H5P_GROUP_CREATE (H5OPEN H5P_CLS_GROUP_CREATE_ID_g)

/*
some other useful things defined in HDF5 headers
*/

#define H5F_ACC_TRUNC (H5CHECK H5OPEN 0x0002u)	/*overwrite existing files   */

#define H5S_ALL (hid_t)0

#define H5P_DEFAULT (hid_t)0 

#define H5P_CRT_ORDER_TRACKED           0x0001
#define H5P_CRT_ORDER_INDEXED           0x0002

// this is an enum in hdf5: H5F_libver_t
#define H5F_LIBVER_LATEST 1

// this is an enum in hdf5: H5T_str_t
#define H5T_STR_NULLTERM 0

// this is an enum in hdf5: H5F_scope_t
#define H5F_SCOPE_LOCAL 0

#endif // MPCO_HDF5_LOADED_AT_RUNTIME

#define HID_INVALID -1

/*namespace for utilities*/
namespace utils {

	/*
	utilities for parsing
	*/
	namespace parsing {

		enum option_type {
			opt_none,
			opt_result_on_nodes,
			opt_result_on_nodes_sens,
			opt_result_on_elements,
			opt_time,
			opt_region
		};

	}

	/*
	utilities for string
	*/
	namespace strings {

		template<class T>
		inline std::string to_string(T x) {
			std::stringstream ss;
			ss << x;
			return ss.str();
		}

		inline void split(const std::string &text, char sep, std::vector<std::string> &tokens, bool skip_empty = false)
		{
			if (tokens.size() > 0) tokens.clear();
			std::size_t start = 0, end = 0;
			while (true)
			{
				end = text.find(sep, start);
				if (end == std::string::npos)
				{
					if (start < text.size())
						tokens.push_back(text.substr(start));
					else if (!skip_empty)
						tokens.push_back("");
					break;
				}
				std::string subs = text.substr(start, end - start);
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

		struct concatenator {
			template<class T>
			inline concatenator &operator << (const T &x) {
				ss << x;
				return *this;
			}
			std::stringstream ss;
		};

	}

	/*
	utilities for shells
	*/
	namespace shell {

		inline bool isShellElementTag(int ele_tag) {
			/**
			used to check if an element is a shell.
			used in workaournd for shell sections.
			note that the keyword "section" doesn't work with shells (only in beams)
			shells use "material" keyword, and then the "fiber" keyword
			*/
			return (
				ele_tag == ELE_TAG_ShellMITC4 ||
				ele_tag == ELE_TAG_ShellMITC4Thermal ||
				ele_tag == ELE_TAG_ShellMITC9 ||
				ele_tag == ELE_TAG_ShellDKGQ ||
				ele_tag == ELE_TAG_ShellNLDKGQ ||
				ele_tag == ELE_TAG_ShellNLDKGQThermal ||
				ele_tag == ELE_TAG_ShellDKGT ||
				ele_tag == ELE_TAG_ShellNLDKGT ||
				ele_tag == ELE_TAG_ShellANDeS ||
				ele_tag == ELE_TAG_ASDShellQ4 ||
				ele_tag == ELE_TAG_ASDShellT3
				);
		}

	}

	/*
	utilities for local axes
	*/
	namespace locax {

		struct quaternion {
			double w;
			double x;
			double y;
			double z;
			inline const double *data()const { return &w; }
		};

		inline quaternion quatFromMat(const Vector &vx, const Vector &vy, const Vector &vz) {
			// vx = first column of rot matrix
			// vy = second column
			// vz = third column
			double xx = vx(0);
			double yy = vy(1);
			double zz = vz(2);
			double tr = xx + yy + zz;
			quaternion Q;
			if ((tr > xx) && (tr > yy) && (tr > zz))
			{
				double S = std::sqrt(tr + 1.0) * 2.0;
				Q.w = 0.25 * S;
				Q.x = (vz(1) - vy(2)) / S;
				Q.y = (vx(2) - vz(0)) / S;
				Q.z = (vy(0) - vx(1)) / S;
			}
			else if ((xx > yy) && (xx > zz))
			{
				double S = std::sqrt(1.0 + xx - yy - zz) * 2.0;
				Q.w = (vz(1) - vy(2)) / S;
				Q.x = 0.25 * S;
				Q.y = (vx(1) + vy(0)) / S;
				Q.z = (vx(2) + vz(0)) / S;
			}
			else if (yy > zz)
			{
				double S = std::sqrt(1.0 + yy - xx - zz) * 2.0;
				Q.w = (vx(2) - vz(0)) / S;
				Q.x = (vx(1) + vy(0)) / S;
				Q.y = 0.25 * S;
				Q.z = (vy(2) + vz(1)) / S;
			}
			else
			{
				double S = std::sqrt(1.0 + zz - xx - yy) * 2.0;
				Q.w = (vy(0) - vx(1)) / S;
				Q.x = (vx(2) + vz(0)) / S;
				Q.y = (vy(2) + vz(1)) / S;
				Q.z = 0.25 * S;
			}
			// safe normalization
			double squared_norm = Q.x*Q.x + Q.y*Q.y + Q.z*Q.z + Q.w*Q.w;
			if (squared_norm > 0.0 && squared_norm != 1.0) {
				double n = std::sqrt(squared_norm);
				Q.x /= n;
				Q.y /= n;
				Q.z /= n;
				Q.w /= n;
			}
			return Q;
		}
	}

	/*
	misc utilities
	*/
	namespace misc {

		template<class T>
		inline bool areVectorsEqual(const std::vector<T> &a, const std::vector<T> &b) {
			if (a.size() != b.size())
				return false;
			for (size_t i = 0; i < a.size(); i++)
				if (a[i] != b[i])
					return false;
			return true;
		}

		template<>
		inline bool areVectorsEqual<std::string>(const std::vector<std::string> &a, const std::vector<std::string> &b) {
			if (a.size() != b.size())
				return false;
			for (size_t i = 0; i < a.size(); i++) {
				const std::string &sa = a[i];
				const std::string &sb = b[i];
				if (sa.size() != sb.size())
					return false;
				if (sa != sb)
					return false;
			}
			return true;
		}

		inline bool lessThanWithTol(double a, double b, double tol) {
			return (a - b) < -tol;
		}

		inline bool greaterThanWithTol(double a, double b, double tol) {
			return (a - b) > tol;
		}

		void bufferNodeResponseVec3u(size_t node_counter, int ndim, const Vector &iresponse, std::vector<double> &buffer)
		{
			size_t j = node_counter * ndim;
			buffer[j] = iresponse[0];
			if (ndim > 1) {
				buffer[j + 1] = iresponse[1];
				if (ndim > 2) {
					buffer[j + 2] = iresponse[2];
				}
			}
		}

		void bufferNodeResponseVec3r(size_t node_counter, int ndim, const Vector &iresponse, std::vector<double> &buffer)
		{
			if (ndim == 2) {
				size_t j = node_counter; // node_counter * 1 (1 rotation component in 2D)
				if (iresponse.Size() > 2)
					buffer[j] = iresponse[2];
				else
					buffer[j] = 0.0;
			}
			else if (ndim == 3) {
				size_t j = node_counter * 3; // node_counter * 3 (3 rotation components in 3D)
				if (iresponse.Size() > 5) {
					buffer[j] = iresponse[3];
					buffer[j + 1] = iresponse[4];
					buffer[j + 2] = iresponse[5];
				}
				else {
					buffer[j] = 0.0;
					buffer[j + 1] = 0.0;
					buffer[j + 2] = 0.0;
				}
				
			}
		}
	}

}

#define MPCO_MAKE_STRING(X) (utils::strings::concatenator() << X).ss.str()

/*utilities for element results. need to stay here before the h5 namespace*/
namespace mpco {

	namespace element {

		struct FiberData
		{
			FiberData()
				: y(0.0), z(0.0), a(0.0) {}
			FiberData(double _y, double _z, double _a)
				: y(_y), z(_z), a(_a) {}
			inline bool operator < (const FiberData &other) const {
				const double rel_tol = 1.0e-5;
				double tol = std::max(std::abs(y), std::abs(other.y))*rel_tol;
				if (utils::misc::lessThanWithTol(y, other.y, tol)) return true;
				if (utils::misc::greaterThanWithTol(y, other.y, tol)) return false;
				tol = std::max(std::abs(z), std::abs(other.z))*rel_tol;
				if (utils::misc::lessThanWithTol(z, other.z, tol)) return true;
				if (utils::misc::greaterThanWithTol(z, other.z, tol)) return false;
				tol = std::max(std::abs(a), std::abs(other.a))*rel_tol;
				if (utils::misc::lessThanWithTol(a, other.a, tol)) return true;
				if (utils::misc::greaterThanWithTol(a, other.a, tol)) return false;
				return false; // equal
			}
			inline bool operator > (const FiberData &other) const {
				const double rel_tol = 1.0e-5;
				double tol = std::max(std::abs(y), std::abs(other.y))*rel_tol;
				if (utils::misc::greaterThanWithTol(y, other.y, tol)) return true;
				if (utils::misc::lessThanWithTol(y, other.y, tol)) return false;
				tol = std::max(std::abs(z), std::abs(other.z))*rel_tol;
				if (utils::misc::greaterThanWithTol(z, other.z, tol)) return true;
				if (utils::misc::lessThanWithTol(z, other.z, tol)) return false;
				tol = std::max(std::abs(a), std::abs(other.a))*rel_tol;
				if (utils::misc::greaterThanWithTol(a, other.a, tol)) return true;
				if (utils::misc::lessThanWithTol(a, other.a, tol)) return false;
				return false; // equal
			}
			inline const double *data()const { return &y; }
			double y;
			double z;
			double a;
		};

		typedef std::vector<FiberData> FiberDataCollection;

		typedef std::vector<int> FiberDataMaterialCollection;

		struct FiberSectionData
		{
			FiberSectionData() : fibers() {}
			inline bool operator < (const FiberSectionData &other) const {
				if (fibers.size() < other.fibers.size()) return true;
				if (fibers.size() > other.fibers.size()) return false;
				for (size_t i = 0; i < fibers.size(); i++) {
					const mpco::element::FiberData &a = fibers[i];
					const mpco::element::FiberData &b = other.fibers[i];
					if (a < b) return true;
					if (a > b) return false;
				}
				if (materials.size() < other.materials.size()) return true;
				if (materials.size() > other.materials.size()) return false;
				for (size_t i = 0; i < materials.size(); i++) {
					int a = materials[i];
					int b = other.materials[i];
					if (a < b) return true;
					if (a > b) return false;
				}
				return false; // equal
			}
			inline bool operator > (const FiberSectionData &other) const {
				if (fibers.size() > other.fibers.size()) return true;
				if (fibers.size() < other.fibers.size()) return false;
				for (size_t i = 0; i < fibers.size(); i++) {
					const mpco::element::FiberData &a = fibers[i];
					const mpco::element::FiberData &b = other.fibers[i];
					if (a > b) return true;
					if (a < b) return false;
				}
				if (materials.size() > other.materials.size()) return true;
				if (materials.size() < other.materials.size()) return false;
				for (size_t i = 0; i < materials.size(); i++) {
					int a = materials[i];
					int b = other.materials[i];
					if (a > b) return true;
					if (a < b) return false;
				}
				return false; // equal
			}
			FiberDataCollection fibers;
			FiberDataMaterialCollection materials;
		};

		struct ElemGaussPair
		{
			ElemGaussPair()
				: elem_id(0), gauss_id(0) {}
			ElemGaussPair(int _elem_id, int _gauss_id)
				: elem_id(_elem_id), gauss_id(_gauss_id) {}
			inline const int *data() const { return &elem_id; }
			int elem_id;
			int gauss_id;
		};

		struct SectionAssignment
		{
			SectionAssignment()
				: is_new(true), name("UnkownClassType"), fiber_section_data(), assignments() {}
			bool is_new;
			std::string name;
			FiberSectionData fiber_section_data;
			std::vector<ElemGaussPair> assignments;
		};

	}

}

/*mpco - hdf5 interface utilities*/
namespace h5 {

	namespace attribute {

		// low level functions for c interface

		herr_t writei(hid_t obj, const char *attr_name, const int *attr_data, hsize_t data_size)
		{
			herr_t status;
			hsize_t dim[1] = { data_size };
			hid_t space = H5Screate_simple(1, dim, NULL);
			hid_t attr = H5Acreate(obj, attr_name, H5T_STD_I32LE, space, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Awrite(attr, H5T_NATIVE_INT, attr_data);
			status = H5Aclose(attr);
			status = H5Sclose(space);
			return status;
		}
		herr_t writed(hid_t obj, const char *attr_name, const double *attr_data, hsize_t data_size)
		{
			herr_t status;
			hsize_t dim[1] = { data_size };
			hid_t space = H5Screate_simple(1, dim, NULL);
			hid_t attr = H5Acreate(obj, attr_name, H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Awrite(attr, H5T_NATIVE_DOUBLE, attr_data);
			status = H5Aclose(attr);
			status = H5Sclose(space);
			return status;
		}
		herr_t writes(hid_t obj, const char *attr_name, const char *attr_data, hsize_t data_size)
		{
			herr_t status;
			hsize_t dim[1] = { 1 };
			hid_t space = H5Screate_simple(1, dim, NULL);
			hid_t atype = H5Tcopy(H5T_C_S1);
			status = H5Tset_size(atype, data_size);
			status = H5Tset_strpad(atype, H5T_STR_NULLTERM);
			hid_t attr = H5Acreate(obj, attr_name, atype, space, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Awrite(attr, atype, attr_data);
			status = H5Aclose(attr);
			status = H5Sclose(space);
			return status;
		}

		// higher level c++ utils

		herr_t write(hid_t obj, const char *attr_name, int attr_data)
		{
			return writei(obj, attr_name, &attr_data, 1);
		}
		herr_t write(hid_t obj, const char *attr_name, const std::vector<int> &attr_data)
		{
			if (attr_data.size() > 0) {
				return writei(obj, attr_name, &attr_data[0], attr_data.size());
			}
			return 0;
		}
		herr_t write(hid_t obj, const char *attr_name, double attr_data)
		{
			return writed(obj, attr_name, &attr_data, 1);
		}
		herr_t write(hid_t obj, const char *attr_name, std::vector<double> &attr_data)
		{
			if (attr_data.size() > 0) {
				return writed(obj, attr_name, &attr_data[0], attr_data.size());
			}
			return 0;
		}
		herr_t write(hid_t obj, const char *attr_name, const std::string &attr_data)
		{
			if (!attr_data.empty()) {
				return writes(obj, attr_name, attr_data.c_str(), attr_data.size());
			}
			return 0;
		}

	}

	namespace group {

		// low level functions for c interface

		hid_t create(hid_t loc_id, const char *name, hid_t lcpl_id, hid_t gcpl_id, hid_t gapl_id) {
			return H5Gcreate(loc_id, name, lcpl_id, gcpl_id, gapl_id);
		}
		herr_t close(hid_t group_id) {
			if (group_id == HID_INVALID) return -1;
			return H5Gclose(group_id);
		}

		// higher level c++ utils

		hid_t createResultGroup(hid_t obj, hid_t gplist, const std::string &result_name,
			const std::string &disp_name, const std::string &components, int num_components,
			const std::string &dimension, const std::string &description,
			int result_type, int data_type)
		{
			herr_t status;
			hid_t h_gp_result = h5::group::create(obj, result_name.c_str(), H5P_DEFAULT, gplist, H5P_DEFAULT);
			status = h5::attribute::write(h_gp_result, "DISPLAY_NAME", disp_name);
			status = h5::attribute::write(h_gp_result, "COMPONENTS", components);
			status = h5::attribute::write(h_gp_result, "DIMENSION", dimension);
			status = h5::attribute::write(h_gp_result, "DESCRIPTION", description);
			status = h5::attribute::write(h_gp_result, "TYPE", result_type);
			status = h5::attribute::write(h_gp_result, "DATA_TYPE", data_type);
			return h_gp_result;
		}

	}

	namespace dataset {

		// low level functions for c interface

		herr_t close(hid_t dataset_id) {
			if (dataset_id == HID_INVALID) return -1;
			return H5Dclose(dataset_id);
		}
		hid_t createAndWrited1(hid_t obj, const char *name, const double *data, hsize_t data_size) {
			// error flags
			herr_t status;
			// create the dataspace
			hsize_t dim[1] = { data_size };
			hid_t space = H5Screate_simple(1, dim, NULL);
			// create the dataset and write data to it.
			hid_t dset = H5Dcreate(obj, name, H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
			// close and release resources
			status = H5Sclose(space);
			return dset;
		}
		hid_t createAndWrited2(hid_t obj, const char *name, const double *data, hsize_t rows, hsize_t cols)
		{
			// error flags
			herr_t status;
			// create the dataspace
			hsize_t dim[2] = { rows, cols };
			hid_t space = H5Screate_simple(2, dim, NULL);
			// create the dataset and write data to it.
			hid_t dset = H5Dcreate(obj, name, H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
			// close and release resources
			status = H5Sclose(space);
			return dset;
		}
		hid_t createAndWritei1(hid_t obj, const char *name, const int *data, hsize_t data_size)
		{
			// error flags
			herr_t status;
			// create the dataspace
			hsize_t dim[1] = { data_size };
			hid_t space = H5Screate_simple(1, dim, NULL);
			// create the dataset and write data to it.
			hid_t dset = H5Dcreate(obj, name, H5T_STD_I32LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
			// close and release resources
			status = H5Sclose(space);
			return dset;
		}
		hid_t createAndWritei2(hid_t obj, const char *name, const int *data, hsize_t rows, hsize_t cols)
		{
			// error flags
			herr_t status;
			// create the dataspace
			hsize_t dim[2] = { rows, cols };
			hid_t space = H5Screate_simple(2, dim, NULL);
			// create the dataset and write data to it.
			hid_t dset = H5Dcreate(obj, name, H5T_STD_I32LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
			// close and release resources
			status = H5Sclose(space);
			return dset;
		}
		hid_t createAndWrites(hid_t obj, const char *name, const char *data, hsize_t data_size)
		{
			// error flags
			herr_t status;
			// create the dataspace
			hsize_t dim[1] = { 1 };
			hid_t space = H5Screate_simple(1, dim, NULL);
			hid_t atype = H5Tcopy(H5T_C_S1);
			status = H5Tset_size(atype, data_size);
			status = H5Tset_strpad(atype, H5T_STR_NULLTERM);
			// create the dataset and write data to it.
			hid_t dset = H5Dcreate(obj, name, atype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			status = H5Dwrite(dset, atype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
			// close and release resources
			status = H5Sclose(space);
			return dset;
		}

		// higher level c++ utils

		hid_t createAndWrite(hid_t obj, const char *name, double data)
		{
			return createAndWrited1(obj, name, &data, 1);
		}
		hid_t createAndWrite(hid_t obj, const char *name, const std::vector<double> &data)
		{
			if (data.size() > 0) {
				return createAndWrited1(obj, name, &data[0], data.size());
			}
			return HID_INVALID;
		}
		hid_t createAndWrite(hid_t obj, const char *name, const std::vector<double> &data, size_t rows, size_t cols)
		{
			if (data.size() > 0 && data.size() == rows*cols) {
				return createAndWrited2(obj, name, &data[0], rows, cols);
			}
			return HID_INVALID;
		}
		hid_t createAndWrite(hid_t obj, const char *name, int data)
		{
			return createAndWritei1(obj, name, &data, 1);
		}
		hid_t createAndWrite(hid_t obj, const char *name, const std::vector<int> &data)
		{
			if (data.size() > 0) {
				return createAndWritei1(obj, name, &data[0], data.size());
			}
			return HID_INVALID;
		}
		hid_t createAndWrite(hid_t obj, const char *name, const std::vector<int> &data, size_t rows, size_t cols)
		{
			if (data.size() > 0 && data.size() == rows*cols) {
				return createAndWritei2(obj, name, &data[0], rows, cols);
			}
			return HID_INVALID;
		}
		hid_t createAndWrite(hid_t obj, const char *name, const std::string &data)
		{
			if (!data.empty()) {
				return createAndWrites(obj, name, data.c_str(), data.size());
			}
			return HID_INVALID;
		}
		hid_t createAndWrite(hid_t obj, const char *name, const std::vector<utils::locax::quaternion> &data)
		{
			if (data.size() > 0) {
				return createAndWrited2(obj, name, data[0].data(), data.size(), 4);
			}
			return HID_INVALID;
		}
		hid_t createAndWrite(hid_t obj, const char *name, const std::vector<mpco::element::FiberData> &data)
		{
			if (data.size() > 0) {
				return createAndWrited2(obj, name, data[0].data(), data.size(), 3);
			}
			return HID_INVALID;
		}
		hid_t createAndWrite(hid_t obj, const char *name, const std::vector<mpco::element::ElemGaussPair> &data)
		{
			if (data.size() > 0) {
				return createAndWritei2(obj, name, data[0].data(), data.size(), 2);
			}
			return HID_INVALID;
		}

	}

	namespace file {

		// low level functions for c interface

		herr_t close(hid_t file_id) {
			return H5Fclose(file_id);
		}
		herr_t flush(hid_t file_id) {
			return H5Fflush(file_id, H5F_SCOPE_LOCAL);
		}
		hid_t create(const char *filename, hid_t create_plist, hid_t acc_plist) {
			return H5Fcreate(filename, H5F_ACC_TRUNC, create_plist, acc_plist);
		}
#ifdef MPCO_USE_SWMR
		herr_t startSWMR(hid_t file_id) {
			return H5Fstart_swmr_write(file_id);
		}
#endif // MPCO_USE_SWMR
	}

	namespace plist {

		enum CreateOptions {
			FileCreate,
			FileAccess,
			GroupCreate
		};

		// low level functions for c interface

		herr_t close(hid_t plist_id) {
			if (plist_id == HID_INVALID) return -1;
			return H5Pclose(plist_id);
		}
		hid_t crate(int opt) {
			switch (opt) {
			case FileCreate:
				return H5Pcreate(H5P_FILE_CREATE);
			case FileAccess:
				return H5Pcreate(H5P_FILE_ACCESS);
			case GroupCreate:
				return H5Pcreate(H5P_GROUP_CREATE);
			default:
				return HID_INVALID;
			}
		}
		herr_t setLinkCreationOrder(hid_t plist_id, unsigned int crt_order_flags) {
			return H5Pset_link_creation_order(plist_id, crt_order_flags);
		}
		herr_t setLibVerBounds(hid_t plist_id, unsigned int minor, unsigned int major) {
			return H5Pset_libver_bounds(plist_id, (H5F_libver_t)minor, (H5F_libver_t)major);
		}
	}

}

/*mpco enums*/
namespace mpco {

	struct ResultType {
		enum Enum {
			Generic = 0,
			Modal
		};
	};

	struct ResultDataType {
		enum Enum {
			Scalar = 0,
			Vectorial,
			Tensorial,
			BeamForceDeformation,
			ShellForceDeformation
		};
	};

	struct NodalResultType {
		enum Enum {
			Displacement,
			Rotation,
			Velocity,
			AngularVelocity,
			Acceleration,
			AngularAcceleration,
			ReactionForce,
			ReactionMoment,
			ReactionForceIncludingInertia,
			ReactionMomentIncludingInertia,
			RayleighForce,
			RayleighMoment,
			Pressure,
			ModesOfVibration,
			ModesOfVibrationRotational,
			// sensitivity
			DisplacementSensitivity,
			RotationSensitivity,
			VelocitySensitivity,
			AngularVelocitySensitivity,
			AccelerationSensitivity,
			AngularAccelerationSensitivity
		};
	};

	struct ElementGeometryType {
		enum Enum {
			Custom = 0,
			Line_2N = 1,
			Line_3N,
			Triangle_3N = 100,
			Triangle_6N,
			Quadrilateral_4N = 200,
			Quadrilateral_8N,
			Quadrilateral_9N,
			Tetrahedron_4N = 300,
			Tetrahedron_10N,
			Hexahedron_8N = 400,
			Hexahedron_20N,
			Hexahedron_27N
		};
	};

	struct ElementIntegrationRuleType {
		enum Enum {
			NoIntegrationRule = 0,
			Line_GaussLegendre_1 = 1,
			Line_GaussLegendre_2,
			Line_GaussLegendre_3,
			Triangle_GaussLegendre_1 = 100,
			Triangle_GaussLegendre_2,
			Triangle_GaussLegendre_2B,
			Triangle_GaussLegendre_2C,
			Quadrilateral_GaussLegendre_1 = 200,
			Quadrilateral_GaussLegendre_2,
			Quadrilateral_GaussLegendre_3,
			Tetrahedron_GaussLegendre_1 = 300,
			Tetrahedron_GaussLegendre_2,
			Hexahedron_GaussLegendre_1 = 400,
			Hexahedron_GaussLegendre_2,
			Hexahedron_GaussLegendre_3,
			CustomIntegrationRule = 1000
		};
	};

	struct ElementOutputDescriptorType {
		/**
		these paths are intepreted as
		results on elements (either constant over element (# metadata = 1) or per element node (# metadata = # element nodes)
		[Element]
		these paths are interpreted as
		results on gauss points (# metadata = # element gauss points)
		[Element.Gauss]
		[Element.Gauss.Material]
		[Element.Gauss.Section]
		[Element.Gauss.Section.Material]
		[Element.Material]
		[Element.Section]
		[Element.Section.Material]
		these paths are interpreted as
		results on fibers (sub-integration points) (# metadata = # element gauss points AND (for each metadata) multiplicity = # fibers)
		[Element.Gauss.Section.Fiber]
		[Element.Gauss.Section.Fiber.Material]
		[Element.Section.Fiber]
		[Element.Section.Fiber.Material]
		*/
		enum Enum {
			Element = 0,
			Gauss,
			Section,
			Fiber,
			Material
		};

		static const char* toString(mpco::ElementOutputDescriptorType::Enum _type) {
			switch (_type) {
			case mpco::ElementOutputDescriptorType::Element: return "ELEMENT";
			case mpco::ElementOutputDescriptorType::Gauss: return "GAUSS";
			case mpco::ElementOutputDescriptorType::Section: return "SECTION";
			case mpco::ElementOutputDescriptorType::Fiber: return "FIBER";
			case mpco::ElementOutputDescriptorType::Material: return "MATERIAL";
			default:
				return "UnknownOutput";
			}
		}
	};
}

/*mpco utilities*/
namespace mpco {

	/*
	holds information for output frequency
	*/
	struct OutputFrequency {
		enum IncrementType {
			DeltaTime,
			NumberOfSteps
		};

		IncrementType type;
		double dt;
		int nsteps;
		double last_time;
		int last_step;

		OutputFrequency() : type(DeltaTime), dt(0.0), nsteps(1), last_time(0.0), last_step(0) {}
		void reset() {
			type = DeltaTime;
			dt = 0.0;
			nsteps = 1;
		}
	};

	/*
	a simple timer
	*/
	class Timer
	{
	public:
		Timer(const std::string &task_name)
			: m_task_name(task_name), m_t0(0), m_t1(0) {}
		inline void start() {
			m_t0 = clock();
			m_t1 = m_t0;
		}
		inline void stop() {
			m_t1 = clock();
			std::cout << "          MPCORecorder. Task = \"" << m_task_name << "\". Elapsed time: " << double(m_t1 - m_t0) / double(CLOCKS_PER_SEC) << "\n";
		}
	private:
		std::string m_task_name;
		clock_t m_t0;
		clock_t m_t1;
	};

	/*
	holds current informations
	*/
	class ProcessInfo
	{
	public:
		ProcessInfo()
			// domain and model information
			: domain(0)
			, num_dimensions(0)
			, current_model_stage_id(-1)
			// some handles in the hdf5 file
			, h_file_id(HID_INVALID)
			, h_file_proplist(HID_INVALID)
#ifdef MPCO_USE_SWMR
			, h_file_acc_proplist(HID_INVALID)
#endif // MPCO_USE_SWMR
			, h_group_proplist(HID_INVALID)
			// time step info
			, current_time_step_id(0)
			, current_time_step(0.0)
			// misc
			, eigen_first_initialization_done(false)
			, record_eigen_on_this_step(false)
			, eigen_last_time_set(0.0)
			, eigen_last_values()
		{}
	public:
		// domain and model information
		Domain *domain;
		int num_dimensions;
		int current_model_stage_id;
		// some handles in the hdf5 file
		hid_t h_file_id;
		hid_t h_file_proplist;
#ifdef MPCO_USE_SWMR
		hid_t h_file_acc_proplist;
#endif // MPCO_USE_SWMR
		hid_t h_group_proplist;
		// time step info
		int current_time_step_id;
		double current_time_step;
		// misc
		bool eigen_first_initialization_done;
		bool record_eigen_on_this_step;
		double eigen_last_time_set;
		Vector eigen_last_values;
	};

}

/*utilities for node results*/
namespace mpco {

	namespace node {

		class ResultRecorder
		{
		public:
			ResultRecorder(const mpco::ProcessInfo &info)
				: m_initialized(false)
				, m_ndim(info.num_dimensions)
				, m_result_name("")
				, m_result_display_name("")
				, m_components_name("")
				, m_num_components(0)
				, m_dimension("")
				, m_description("")
				, m_result_type(mpco::ResultType::Generic)
				, m_result_data_type(mpco::ResultDataType::Scalar)
			{}
			virtual ~ResultRecorder() {}
			virtual int getReactionFlag()const { return -1; }
			virtual int record(mpco::ProcessInfo &info, std::vector<Node*> &nodes)
			{
				/*
				error flags
				*/
				int retval = 0;
				herr_t status = 0;
				/*
				quick return
				*/
				if (m_num_components < 1)
					return retval;
				/*
				operations performed only once
				*/
				if (!m_initialized) {
					/*
					create result group
					*/
					hid_t h_gp_result = h5::group::createResultGroup(info.h_file_id, info.h_group_proplist,
						m_result_name, m_result_display_name,
						m_components_name, m_num_components, m_dimension, m_description,
						m_result_type, m_result_data_type);
					/*
					create the id dataset
					*/
					std::vector<int> buffer_id(nodes.size());
					for (size_t i = 0; i < nodes.size(); i++)
						buffer_id[i] = nodes[i]->getTag();
					hid_t h_dset_id = h5::dataset::createAndWrite(h_gp_result, "ID", buffer_id, buffer_id.size(), 1);
					/*
					create the data group
					*/
					hid_t h_gp_data = h5::group::create(h_gp_result, "DATA", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
					/*
					done
					*/
					status = h5::group::close(h_gp_data);
					status = h5::dataset::close(h_dset_id);
					status = h5::group::close(h_gp_result);
					m_initialized = true;
				}
				/*
				create the dataset for this timestep
				*/
				std::vector<double> buffer_data(nodes.size() * m_num_components);
				bufferResponse(info, nodes, buffer_data);
				std::stringstream ss_dset_name;
				ss_dset_name << m_result_name << "/DATA/STEP_" << info.current_time_step_id;
				std::string dset_name = ss_dset_name.str();
				hid_t h_dset_data = h5::dataset::createAndWrite(info.h_file_id, dset_name.c_str(), buffer_data, nodes.size(), m_num_components);
				status = h5::attribute::write(h_dset_data, "STEP", info.current_time_step_id);
				status = h5::attribute::write(h_dset_data, "TIME", info.current_time_step);
				status = h5::dataset::close(h_dset_data);
				/*
				return
				*/
				return retval;
			}
		protected:
			virtual void bufferResponse(mpco::ProcessInfo &info, std::vector<Node*> &nodes, std::vector<double> &buffer)const = 0;
		protected:
			bool m_initialized;
			int m_ndim;
			std::string m_result_name;
			std::string m_result_display_name;
			std::string m_components_name;
			int m_num_components;
			std::string m_dimension;
			std::string m_description;
			mpco::ResultType::Enum m_result_type;
			mpco::ResultDataType::Enum m_result_data_type;
		};

		class ResultRecorderDisplacement : public ResultRecorder
		{
		public:
			ResultRecorderDisplacement(const mpco::ProcessInfo &info)
				: ResultRecorder(info)
			{
				std::stringstream ss_buffer;
				ss_buffer << "MODEL_STAGE[" << info.current_model_stage_id << "]/RESULTS/ON_NODES/DISPLACEMENT";
				m_result_name = ss_buffer.str();
				m_result_display_name = "Displacement";
				m_num_components = 0;
				if (m_ndim == 1) {
					m_components_name = "Ux";
					m_num_components = 1;
					m_result_data_type = mpco::ResultDataType::Scalar;
				}
				else if (m_ndim == 2) {
					m_components_name = "Ux,Uy";
					m_num_components = 2;
					m_result_data_type = mpco::ResultDataType::Vectorial;
				}
				else if (m_ndim == 3) {
					m_components_name = "Ux,Uy,Uz";
					m_num_components = 3;
					m_result_data_type = mpco::ResultDataType::Vectorial;
				}
				m_dimension = "L";
				m_description = "Nodal displacement field";
				m_result_type = mpco::ResultType::Generic;
			}
		protected:
			virtual void bufferResponse(mpco::ProcessInfo &info, std::vector<Node*> &nodes, std::vector<double> &buffer)const {
				for (size_t i = 0; i < nodes.size(); i++)
					utils::misc::bufferNodeResponseVec3u(i, m_ndim, nodes[i]->getTrialDisp(), buffer);
			}
		};

		class ResultRecorderRotation : public ResultRecorder
		{
		public:
			ResultRecorderRotation(const mpco::ProcessInfo &info)
				: ResultRecorder(info)
			{
				std::stringstream ss_buffer;
				ss_buffer << "MODEL_STAGE[" << info.current_model_stage_id << "]/RESULTS/ON_NODES/ROTATION";
				m_result_name = ss_buffer.str();
				m_result_display_name = "Rotation";
				m_num_components = 0;
				if (m_ndim == 2) {
					m_components_name = "Rz";
					m_num_components = 1;
					m_result_data_type = mpco::ResultDataType::Scalar;
				}
				else {
					m_components_name = "Rx,Ry,Rz";
					m_num_components = 3;
					m_result_data_type = mpco::ResultDataType::Vectorial;
				}
				m_dimension = "A";
				m_description = "Nodal rotation field";
				m_result_type = mpco::ResultType::Generic;
			}
		protected:
			virtual void bufferResponse(mpco::ProcessInfo &info, std::vector<Node*> &nodes, std::vector<double> &buffer)const {
				for (size_t i = 0; i < nodes.size(); i++)
					utils::misc::bufferNodeResponseVec3r(i, m_ndim, nodes[i]->getTrialDisp(), buffer);
			}
		};

		class ResultRecorderPressure : public ResultRecorder
		{
		public:
			ResultRecorderPressure(const mpco::ProcessInfo &info)
				: ResultRecorder(info)
			{
				std::stringstream ss_buffer;
				ss_buffer << "MODEL_STAGE[" << info.current_model_stage_id << "]/RESULTS/ON_NODES/PRESSURE";
				m_result_name = ss_buffer.str();
				m_result_display_name = "Pressure";
				m_components_name = "p";
				m_num_components = 1;
				m_result_data_type = mpco::ResultDataType::Scalar;
				m_dimension = "F/L^2";
				m_description = "Nodal pressure field";
				m_result_type = mpco::ResultType::Generic;
			}
		protected:
			virtual void bufferResponse(mpco::ProcessInfo &info, std::vector<Node*> &nodes, std::vector<double> &buffer)const {
				for (size_t i = 0; i < nodes.size(); i++) {
					double pressure = 0.0;
					Node *inode = nodes[i];
					const Vector &vel = inode->getTrialVel();
					if (inode->getCrds().Size() == 2) {
						if (vel.Size() == 3) {
							pressure = vel[2];
						}
					}
					else if (inode->getCrds().Size() == 3) {
						if (vel.Size() == 4) {
							pressure = vel[3];
						}
					}
					buffer[i] = pressure;
				}
			}
		};

		class ResultRecorderReactionForce : public ResultRecorder
		{
		public:
			ResultRecorderReactionForce(const mpco::ProcessInfo &info)
				: ResultRecorder(info)
			{
				std::stringstream ss_buffer;
				ss_buffer << "MODEL_STAGE[" << info.current_model_stage_id << "]/RESULTS/ON_NODES/REACTION_FORCE";
				m_result_name = ss_buffer.str();
				m_result_display_name = "Reaction Force";
				m_num_components = 0;
				if (m_ndim == 1) {
					m_components_name = "RFx";
					m_num_components = 1;
					m_result_data_type = mpco::ResultDataType::Scalar;
				}
				else if (m_ndim == 2) {
					m_components_name = "RFx,RFy";
					m_num_components = 2;
					m_result_data_type = mpco::ResultDataType::Vectorial;
				}
				else if (m_ndim == 3) {
					m_components_name = "RFx,RFy,RFz";
					m_num_components = 3;
					m_result_data_type = mpco::ResultDataType::Vectorial;
				}
				m_dimension = "F";
				m_description = "Nodal reaction force field";
				m_result_type = mpco::ResultType::Generic;
			}
		protected:
			virtual int getReactionFlag()const { return 0; }
			virtual void bufferResponse(mpco::ProcessInfo &info, std::vector<Node*> &nodes, std::vector<double> &buffer)const {
				for (size_t i = 0; i < nodes.size(); i++)
					utils::misc::bufferNodeResponseVec3u(i, m_ndim, nodes[i]->getReaction(), buffer);
			}
		};

		class ResultRecorderReactionMoment : public ResultRecorder
		{
		public:
			ResultRecorderReactionMoment(const mpco::ProcessInfo &info)
				: ResultRecorder(info)
			{
				std::stringstream ss_buffer;
				ss_buffer << "MODEL_STAGE[" << info.current_model_stage_id << "]/RESULTS/ON_NODES/REACTION_MOMENT";
				m_result_name = ss_buffer.str();
				m_result_display_name = "Reaction Moment";
				m_num_components = 0;
				if (m_ndim == 2) {
					m_components_name = "RMz";
					m_num_components = 1;
					m_result_data_type = mpco::ResultDataType::Scalar;
				}
				else {
					m_components_name = "RMx,RMy,RMz";
					m_num_components = 3;
					m_result_data_type = mpco::ResultDataType::Vectorial;
				}
				m_dimension = "F*L";
				m_description = "Nodal reaction moment field";
				m_result_type = mpco::ResultType::Generic;
			}
		protected:
			virtual int getReactionFlag()const { return 0; }
			virtual void bufferResponse(mpco::ProcessInfo &info, std::vector<Node*> &nodes, std::vector<double> &buffer)const {
				for (size_t i = 0; i < nodes.size(); i++)
					utils::misc::bufferNodeResponseVec3r(i, m_ndim, nodes[i]->getReaction(), buffer);
			}
		};

		class ResultRecorderReactionForceIncIntertia : public ResultRecorderReactionForce
		{
		public:
			ResultRecorderReactionForceIncIntertia(const mpco::ProcessInfo &info)
				: ResultRecorderReactionForce(info)
			{
				std::stringstream ss_buffer;
				ss_buffer << "MODEL_STAGE[" << info.current_model_stage_id << "]/RESULTS/ON_NODES/REACTION_FORCE_INCLUDING_INERTIA";
				m_result_name = ss_buffer.str();
				m_result_display_name = "Reaction Force Including Inertia";
				m_description = "Nodal reaction force field including inertia";
			}
		protected:
			virtual int getReactionFlag()const { return 1; }
		};

		class ResultRecorderReactionMomentIncIntertia : public ResultRecorderReactionMoment
		{
		public:
			ResultRecorderReactionMomentIncIntertia(const mpco::ProcessInfo &info)
				: ResultRecorderReactionMoment(info)
			{
				std::stringstream ss_buffer;
				ss_buffer << "MODEL_STAGE[" << info.current_model_stage_id << "]/RESULTS/ON_NODES/REACTION_MOMENT_INCLUDING_INERTIA";
				m_result_name = ss_buffer.str();
				m_result_display_name = "Reaction Moment Including Inertia";
				m_description = "Nodal reaction moment field including inertia";
			}
		protected:
			virtual int getReactionFlag()const { return 1; }
		};

		class ResultRecorderReactionForceRayleigh : public ResultRecorderReactionForce
		{
		public:
			ResultRecorderReactionForceRayleigh(const mpco::ProcessInfo &info)
				: ResultRecorderReactionForce(info)
			{
				std::stringstream ss_buffer;
				ss_buffer << "MODEL_STAGE[" << info.current_model_stage_id << "]/RESULTS/ON_NODES/RAYLEIGH_FORCE";
				m_result_name = ss_buffer.str();
				m_result_display_name = "Rayleigh Force";
				m_description = "Nodal Rayleigh force field";
			}
		protected:
			virtual int getReactionFlag()const { return 2; }
		};

		class ResultRecorderReactionMomentRayleigh : public ResultRecorderReactionMoment
		{
		public:
			ResultRecorderReactionMomentRayleigh(const mpco::ProcessInfo &info)
				: ResultRecorderReactionMoment(info)
			{
				std::stringstream ss_buffer;
				ss_buffer << "MODEL_STAGE[" << info.current_model_stage_id << "]/RESULTS/ON_NODES/RAYLEIGH_MOMENT";
				m_result_name = ss_buffer.str();
				m_result_display_name = "Rayleigh Moment";
				m_description = "Nodal Rayleigh moment field";
			}
		protected:
			virtual int getReactionFlag()const { return 2; }
		};

		class ResultRecorderVelocity : public ResultRecorder
		{
		public:
			ResultRecorderVelocity(const mpco::ProcessInfo &info)
				: ResultRecorder(info)
			{
				std::stringstream ss_buffer;
				ss_buffer << "MODEL_STAGE[" << info.current_model_stage_id << "]/RESULTS/ON_NODES/VELOCITY";
				m_result_name = ss_buffer.str();
				m_result_display_name = "Velocity";
				m_num_components = 0;
				if (m_ndim == 1) {
					m_components_name = "Vx";
					m_num_components = 1;
					m_result_data_type = mpco::ResultDataType::Scalar;
				}
				else if (m_ndim == 2) {
					m_components_name = "Vx,Vy";
					m_num_components = 2;
					m_result_data_type = mpco::ResultDataType::Vectorial;
				}
				else if (m_ndim == 3) {
					m_components_name = "Vx,Vy,Vz";
					m_num_components = 3;
					m_result_data_type = mpco::ResultDataType::Vectorial;
				}
				m_dimension = "L/t";
				m_description = "Nodal velocity field";
				m_result_type = mpco::ResultType::Generic;
			}
		protected:
			virtual void bufferResponse(mpco::ProcessInfo &info, std::vector<Node*> &nodes, std::vector<double> &buffer)const {
				for (size_t i = 0; i < nodes.size(); i++)
					utils::misc::bufferNodeResponseVec3u(i, m_ndim, nodes[i]->getTrialVel(), buffer);
			}
		};

		class ResultRecorderAngularVelocity : public ResultRecorder
		{
		public:
			ResultRecorderAngularVelocity(const mpco::ProcessInfo &info)
				: ResultRecorder(info)
			{
				std::stringstream ss_buffer;
				ss_buffer << "MODEL_STAGE[" << info.current_model_stage_id << "]/RESULTS/ON_NODES/ANGULAR_VELOCITY";
				m_result_name = ss_buffer.str();
				m_result_display_name = "Angular Velocity";
				m_num_components = 0;
				if (m_ndim == 2) {
					m_components_name = "RVz";
					m_num_components = 1;
					m_result_data_type = mpco::ResultDataType::Scalar;
				}
				else {
					m_components_name = "RVx,RVy,RVz";
					m_num_components = 3;
					m_result_data_type = mpco::ResultDataType::Vectorial;
				}
				m_dimension = "L/t*A";
				m_description = "Nodal angular velocity field";
				m_result_type = mpco::ResultType::Generic;
			}
		protected:
			virtual void bufferResponse(mpco::ProcessInfo &info, std::vector<Node*> &nodes, std::vector<double> &buffer)const {
				for (size_t i = 0; i < nodes.size(); i++)
					utils::misc::bufferNodeResponseVec3r(i, m_ndim, nodes[i]->getTrialVel(), buffer);
			}
		};

		class ResultRecorderAcceleration : public ResultRecorder
		{
		public:
			ResultRecorderAcceleration(const mpco::ProcessInfo &info)
				: ResultRecorder(info)
			{
				std::stringstream ss_buffer;
				ss_buffer << "MODEL_STAGE[" << info.current_model_stage_id << "]/RESULTS/ON_NODES/ACCELERATION";
				m_result_name = ss_buffer.str();
				m_result_display_name = "Acceleration";
				m_num_components = 0;
				if (m_ndim == 1) {
					m_components_name = "Ax";
					m_num_components = 1;
					m_result_data_type = mpco::ResultDataType::Scalar;
				}
				else if (m_ndim == 2) {
					m_components_name = "Ax,Ay";
					m_num_components = 2;
					m_result_data_type = mpco::ResultDataType::Vectorial;
				}
				else if (m_ndim == 3) {
					m_components_name = "Ax,Ay,Az";
					m_num_components = 3;
					m_result_data_type = mpco::ResultDataType::Vectorial;
				}
				m_dimension = "L/t^2";
				m_description = "Nodal acceleration field";
				m_result_type = mpco::ResultType::Generic;
			}
		protected:
			virtual void bufferResponse(mpco::ProcessInfo &info, std::vector<Node*> &nodes, std::vector<double> &buffer)const {
				for (size_t i = 0; i < nodes.size(); i++)
					utils::misc::bufferNodeResponseVec3u(i, m_ndim, nodes[i]->getTrialAccel(), buffer);
			}
		};

		class ResultRecorderAngularAcceleration : public ResultRecorder
		{
		public:
			ResultRecorderAngularAcceleration(const mpco::ProcessInfo &info)
				: ResultRecorder(info)
			{
				std::stringstream ss_buffer;
				ss_buffer << "MODEL_STAGE[" << info.current_model_stage_id << "]/RESULTS/ON_NODES/ANGULAR_ACCELERATION";
				m_result_name = ss_buffer.str();
				m_result_display_name = "Angular Acceleration";
				m_num_components = 0;
				if (m_ndim == 2) {
					m_components_name = "RAz";
					m_num_components = 1;
					m_result_data_type = mpco::ResultDataType::Scalar;
				}
				else {
					m_components_name = "RAx,RAy,RAz";
					m_num_components = 3;
					m_result_data_type = mpco::ResultDataType::Vectorial;
				}
				m_dimension = "L/t^2*A";
				m_description = "Nodal angular acceleration field";
				m_result_type = mpco::ResultType::Generic;
			}
		protected:
			virtual void bufferResponse(mpco::ProcessInfo &info, std::vector<Node*> &nodes, std::vector<double> &buffer)const {
				for (size_t i = 0; i < nodes.size(); i++)
					utils::misc::bufferNodeResponseVec3r(i, m_ndim, nodes[i]->getTrialAccel(), buffer);
			}
		};

		class ResultRecorderModesOfVibration : public ResultRecorder
		{
		public:
			ResultRecorderModesOfVibration(const mpco::ProcessInfo &info)
				: ResultRecorder(info)
				, m_current_mode_to_buffer(0)
			{
				std::stringstream ss_buffer;
				ss_buffer << "MODEL_STAGE[" << info.current_model_stage_id << "]/RESULTS/ON_NODES/MODES_OF_VIBRATION(U)";
				m_result_name = ss_buffer.str();
				m_result_display_name = "Modes of Vibration (U)";
				m_num_components = 0;
				if (m_ndim == 1) {
					m_components_name = "Ux";
					m_num_components = 1;
					m_result_data_type = mpco::ResultDataType::Scalar;
				}
				else if (m_ndim == 2) {
					m_components_name = "Ux,Uy";
					m_num_components = 2;
					m_result_data_type = mpco::ResultDataType::Vectorial;
				}
				else if (m_ndim == 3) {
					m_components_name = "Ux,Uy,Uz";
					m_num_components = 3;
					m_result_data_type = mpco::ResultDataType::Vectorial;
				}
				m_dimension = "";
				m_description = "Eigenvector (translational components)";
				m_result_type = mpco::ResultType::Modal;
			}
			virtual int record(mpco::ProcessInfo &info, std::vector<Node*> &nodes)
			{
				/*
				error flags
				*/
				int retval = 0;
				herr_t status = 0;
				/*
				quick return
				*/
				if (m_num_components < 1)
					return retval;
				/*
				quick return
				*/
				if (!info.record_eigen_on_this_step)
					return retval;
				/*
				operations performed only once
				*/
				if (!m_initialized) {
					/*
					create result group
					*/
					hid_t h_gp_result = h5::group::createResultGroup(info.h_file_id, info.h_group_proplist,
						m_result_name, m_result_display_name,
						m_components_name, m_num_components, m_dimension, m_description,
						m_result_type, m_result_data_type);
					/*
					create the id dataset
					*/
					std::vector<int> buffer_id(nodes.size());
					for (size_t i = 0; i < nodes.size(); i++)
						buffer_id[i] = nodes[i]->getTag();
					hid_t h_dset_id = h5::dataset::createAndWrite(h_gp_result, "ID", buffer_id, buffer_id.size(), 1);
					/*
					create the data group
					*/
					hid_t h_gp_data = h5::group::create(h_gp_result, "DATA", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
					/*
					done
					*/
					status = h5::group::close(h_gp_data);
					status = h5::dataset::close(h_dset_id);
					status = h5::group::close(h_gp_result);
					m_initialized = true;
				}
				/*
				create the timestep group
				*/
				std::stringstream ss_gp_step_name;
				ss_gp_step_name << m_result_name << "/DATA/STEP_" << info.current_time_step_id;
				std::string gp_step_name = ss_gp_step_name.str();
				hid_t h_gp_step = h5::group::create(info.h_file_id, gp_step_name.c_str(), H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
				status = h5::attribute::write(h_gp_step, "STEP", info.current_time_step_id);
				status = h5::attribute::write(h_gp_step, "TIME", info.current_time_step);
				/*
				allocate buffer response
				*/
				std::vector<double> buffer(nodes.size() * m_num_components);
				/*
				for each mode
				*/
				int num_eigen = *OPS_GetNumEigen();
				for (int k = 0; k < num_eigen; k++) {
					/*
					buffer response
					*/
					m_current_mode_to_buffer = k;
					bufferResponse(info, nodes, buffer);
					/*
					compute additional infos
					*/
					double lambda = k < info.eigen_last_values.Size() ? info.eigen_last_values[k] : 0.0;
					double omega = std::sqrt(lambda);
					double freq = omega / (2.0 * M_PI);
					double period = 1.0 / freq;
					/*
					create modal group
					*/
					std::stringstream ss_dset_name;
					ss_dset_name << "MODE_" << k;
					std::string dset_name = ss_dset_name.str();
					hid_t h_dset_data = h5::dataset::createAndWrite(h_gp_step, dset_name.c_str(), buffer, nodes.size(), m_num_components);
					status = h5::attribute::write(h_dset_data, "MODE", k);
					status = h5::attribute::write(h_dset_data, "LAMBDA", lambda);
					status = h5::attribute::write(h_dset_data, "OMEGA", omega);
					status = h5::attribute::write(h_dset_data, "FREQUENCY", freq);
					status = h5::attribute::write(h_dset_data, "PERIOD", period);
					status = h5::dataset::close(h_dset_data);
				}
				/*
				close the timestep group
				*/
				status = h5::group::close(h_gp_step);
				/*
				return
				*/
				return retval;
			}
		protected:
			virtual void bufferResponse(mpco::ProcessInfo &info, std::vector<Node*> &nodes, std::vector<double> &buffer)const {
				for (size_t i = 0; i < nodes.size(); i++) {
					size_t j = i * m_ndim;
					const Matrix &inode_eigenvec = nodes[i]->getEigenvectors();
					buffer[j] = inode_eigenvec(0, m_current_mode_to_buffer);
					if (m_ndim > 1) {
						buffer[j + 1] = inode_eigenvec(1, m_current_mode_to_buffer);
						if (m_ndim > 2) {
							buffer[j + 2] = inode_eigenvec(2, m_current_mode_to_buffer);
						}
					}
				}
			}
		protected:
			int m_current_mode_to_buffer;
		};

		class ResultRecorderModesOfVibrationRotational : public ResultRecorderModesOfVibration
		{
		public:
			ResultRecorderModesOfVibrationRotational(const mpco::ProcessInfo &info)
				: ResultRecorderModesOfVibration(info)
			{
				std::stringstream ss_buffer;
				ss_buffer << "MODEL_STAGE[" << info.current_model_stage_id << "]/RESULTS/ON_NODES/MODES_OF_VIBRATION(R)";
				m_result_name = ss_buffer.str();
				m_result_display_name = "Modes of Vibration (R)";
				m_num_components = 0;
				if (m_ndim == 2) {
					m_components_name = "Rz";
					m_num_components = 1;
					m_result_data_type = mpco::ResultDataType::Scalar;
				}
				else {
					m_components_name = "Rx,Ry,Rz";
					m_num_components = 3;
					m_result_data_type = mpco::ResultDataType::Vectorial;
				}
				m_dimension = "";
				m_description = "Eigenvector (rotational components)";
				m_result_type = mpco::ResultType::Modal;
			}
		protected:
			virtual void bufferResponse(mpco::ProcessInfo &info, std::vector<Node*> &nodes, std::vector<double> &buffer)const {
				for (size_t i = 0; i < nodes.size(); i++) {
					const Matrix &inode_eigenvec = nodes[i]->getEigenvectors();
					if (m_ndim == 2 && inode_eigenvec.noRows() > 2) {
						size_t j = i;
						buffer[j] = inode_eigenvec(2, m_current_mode_to_buffer);
					}
					else if (m_ndim == 3 && inode_eigenvec.noRows() > 5) {
						size_t j = i * 3;
						buffer[j] = inode_eigenvec(3, m_current_mode_to_buffer);
						buffer[j + 1] = inode_eigenvec(4, m_current_mode_to_buffer);
						buffer[j + 2] = inode_eigenvec(5, m_current_mode_to_buffer);
					}
				}
			}
		};

		class ResultRecorderDisplacementSensitivity : public ResultRecorder
		{
		public:
			ResultRecorderDisplacementSensitivity(const mpco::ProcessInfo &info, int grad)
				: ResultRecorder(info)
				, m_grad(grad)
			{
				std::stringstream ss_buffer;
				ss_buffer << "MODEL_STAGE[" << info.current_model_stage_id << "]/RESULTS/ON_NODES/DISPLACEMENT_SENSITIVITY_" << m_grad;
				m_result_name = ss_buffer.str();
				m_result_display_name = MPCO_MAKE_STRING("Displacement Sensitivity d_U/d_q" << m_grad);
				m_num_components = 0;
				if (m_ndim == 1) {
					m_components_name = MPCO_MAKE_STRING("d_Ux/d_q" << m_grad);
					m_num_components = 1;
					m_result_data_type = mpco::ResultDataType::Scalar;
				}
				else if (m_ndim == 2) {
					m_components_name = MPCO_MAKE_STRING("d_Ux/d_q" << m_grad << ",d_Uy/d_q" << m_grad);
					m_num_components = 2;
					m_result_data_type = mpco::ResultDataType::Vectorial;
				}
				else if (m_ndim == 3) {
					m_components_name = MPCO_MAKE_STRING("d_Ux/d_q" << m_grad << ",d_Uy/d_q" << m_grad << ",d_Uz/d_q" << m_grad);
					m_num_components = 3;
					m_result_data_type = mpco::ResultDataType::Vectorial;
				}
				m_description = MPCO_MAKE_STRING("Nodal displacement sensitivity field d_U/d_q" << m_grad);
				m_result_type = mpco::ResultType::Generic;
			}
		protected:
			virtual void bufferResponse(mpco::ProcessInfo &info, std::vector<Node*> &nodes, std::vector<double> &buffer)const {
				if (m_grad < 1 || m_grad > info.domain->getNumParameters()) {
					opserr << "MPCORecorder sensitivity parameter is out of range: grad index = " << m_grad << ", domain->getNumParameters() = " << info.domain->getNumParameters() << "\n";
					return;
				}
				for (size_t i = 0; i < nodes.size(); i++) {
					Node *inode = nodes[i];
					size_t j = i * m_ndim;
					buffer[j] = inode->getDispSensitivity(1, m_grad - 1);
					if (m_ndim > 1) {
						buffer[j + 1] = inode->getDispSensitivity(2, m_grad - 1);
						if (m_ndim > 2) {
							buffer[j + 2] = inode->getDispSensitivity(3, m_grad - 1);
						}
					}
				}
			}
		private:
			int m_grad;
		};

		class ResultRecorderRotationSensitivity : public ResultRecorder
		{
		public:
			ResultRecorderRotationSensitivity(const mpco::ProcessInfo &info, int grad)
				: ResultRecorder(info)
				, m_grad(grad)
			{
				std::stringstream ss_buffer;
				ss_buffer << "MODEL_STAGE[" << info.current_model_stage_id << "]/RESULTS/ON_NODES/ROTATION_SENSITIVITY_" << m_grad;;
				m_result_name = ss_buffer.str();
				m_result_display_name = MPCO_MAKE_STRING("Rotation Sensitivity d_R/d_q" << m_grad);
				m_num_components = 0;
				if (m_ndim == 2) {
					m_components_name = MPCO_MAKE_STRING("d_Rz/d_q" << m_grad);
					m_num_components = 1;
					m_result_data_type = mpco::ResultDataType::Scalar;
				}
				else {
					m_components_name = MPCO_MAKE_STRING("d_Rx/d_q" << m_grad << ",d_Ry/d_q" << m_grad << ",d_Rz/d_q" << m_grad);
					m_num_components = 3;
					m_result_data_type = mpco::ResultDataType::Vectorial;
				}
				m_description = MPCO_MAKE_STRING("Nodal rotation sensitivity field d_R/d_q" << m_grad);
				m_result_type = mpco::ResultType::Generic;
			}
		protected:
			virtual void bufferResponse(mpco::ProcessInfo &info, std::vector<Node*> &nodes, std::vector<double> &buffer)const {
				if (m_grad < 1 || m_grad > info.domain->getNumParameters()) {
					opserr << "MPCORecorder sensitivity parameter is out of range: grad index = " << m_grad << ", domain->getNumParameters() = " << info.domain->getNumParameters() << "\n";
					return;
				}
				for (size_t i = 0; i < nodes.size(); i++) {
					Node *inode = nodes[i];
					const Vector &iresponse = inode->getTrialDisp(); // aux, just to check for rotational components
					if (m_ndim == 2) {
						size_t j = i; // node_counter * 1 (1 rotation component in 2D)
						if (iresponse.Size() > 2)
							buffer[j] = inode->getDispSensitivity(3, m_grad - 1);
						else
							buffer[j] = 0.0;
					}
					else if (m_ndim == 3) {
						size_t j = i * 3; // node_counter * 3 (3 rotation components in 3D)
						if (iresponse.Size() > 5) {
							buffer[j] = inode->getDispSensitivity(4, m_grad - 1);
							buffer[j + 1] = inode->getDispSensitivity(5, m_grad - 1);
							buffer[j + 2] = inode->getDispSensitivity(6, m_grad - 1);
						}
						else {
							buffer[j] = 0.0;
							buffer[j + 1] = 0.0;
							buffer[j + 2] = 0.0;
						}
					}
				}
			}
		private:
			int m_grad;
		};

		class ResultRecorderVelocitySensitivity : public ResultRecorder
		{
		public:
			ResultRecorderVelocitySensitivity(const mpco::ProcessInfo &info, int grad)
				: ResultRecorder(info)
				, m_grad(grad)
			{
				std::stringstream ss_buffer;
				ss_buffer << "MODEL_STAGE[" << info.current_model_stage_id << "]/RESULTS/ON_NODES/VELOCITY_SENSITIVITY_" << m_grad;
				m_result_name = ss_buffer.str();
				m_result_display_name = MPCO_MAKE_STRING("Velocity Sensitivity d_V/d_q" << m_grad);
				m_num_components = 0;
				if (m_ndim == 1) {
					m_components_name = MPCO_MAKE_STRING("d_Vx/d_q" << m_grad);
					m_num_components = 1;
					m_result_data_type = mpco::ResultDataType::Scalar;
				}
				else if (m_ndim == 2) {
					m_components_name = MPCO_MAKE_STRING("d_Vx/d_q" << m_grad << ",d_Vy/d_q" << m_grad);
					m_num_components = 2;
					m_result_data_type = mpco::ResultDataType::Vectorial;
				}
				else if (m_ndim == 3) {
					m_components_name = MPCO_MAKE_STRING("d_Vx/d_q" << m_grad << ",d_Vy/d_q" << m_grad << ",d_Vz/d_q" << m_grad);
					m_num_components = 3;
					m_result_data_type = mpco::ResultDataType::Vectorial;
				}
				m_description = MPCO_MAKE_STRING("Nodal velocity sensitivity field d_V/d_q" << m_grad);
				m_result_type = mpco::ResultType::Generic;
			}
		protected:
			virtual void bufferResponse(mpco::ProcessInfo &info, std::vector<Node*> &nodes, std::vector<double> &buffer)const {
				if (m_grad < 1 || m_grad > info.domain->getNumParameters()) {
					opserr << "MPCORecorder sensitivity parameter is out of range: grad index = " << m_grad << ", domain->getNumParameters() = " << info.domain->getNumParameters() << "\n";
					return;
				}
				for (size_t i = 0; i < nodes.size(); i++) {
					Node *inode = nodes[i];
					size_t j = i * m_ndim;
					buffer[j] = inode->getVelSensitivity(1, m_grad - 1);
					if (m_ndim > 1) {
						buffer[j + 1] = inode->getVelSensitivity(2, m_grad - 1);
						if (m_ndim > 2) {
							buffer[j + 2] = inode->getVelSensitivity(3, m_grad - 1);
						}
					}
				}
			}
		private:
			int m_grad;
		};

		class ResultRecorderAngularVelocitySensitivity : public ResultRecorder
		{
		public:
			ResultRecorderAngularVelocitySensitivity(const mpco::ProcessInfo &info, int grad)
				: ResultRecorder(info)
				, m_grad(grad)
			{
				std::stringstream ss_buffer;
				ss_buffer << "MODEL_STAGE[" << info.current_model_stage_id << "]/RESULTS/ON_NODES/ANGULAR_VELOCITY_SENSITIVITY_" << m_grad;;
				m_result_name = ss_buffer.str();
				m_result_display_name = MPCO_MAKE_STRING("Angular Velocity Sensitivity d_RV/d_q" << m_grad);
				m_num_components = 0;
				if (m_ndim == 2) {
					m_components_name = MPCO_MAKE_STRING("d_RVz/d_q" << m_grad);
					m_num_components = 1;
					m_result_data_type = mpco::ResultDataType::Scalar;
				}
				else {
					m_components_name = MPCO_MAKE_STRING("d_RVx/d_q" << m_grad << ",d_RVy/d_q" << m_grad << ",d_RVz/d_q" << m_grad);
					m_num_components = 3;
					m_result_data_type = mpco::ResultDataType::Vectorial;
				}
				m_description = MPCO_MAKE_STRING("Nodal angular velocity sensitivity field d_RV/d_q" << m_grad);
				m_result_type = mpco::ResultType::Generic;
			}
		protected:
			virtual void bufferResponse(mpco::ProcessInfo &info, std::vector<Node*> &nodes, std::vector<double> &buffer)const {
				if (m_grad < 1 || m_grad > info.domain->getNumParameters()) {
					opserr << "MPCORecorder sensitivity parameter is out of range: grad index = " << m_grad << ", domain->getNumParameters() = " << info.domain->getNumParameters() << "\n";
					return;
				}
				for (size_t i = 0; i < nodes.size(); i++) {
					Node *inode = nodes[i];
					const Vector &iresponse = inode->getTrialVel(); // aux, just to check for angular velocityal components
					if (m_ndim == 2) {
						size_t j = i; // node_counter * 1 (1 angular velocity component in 2D)
						if (iresponse.Size() > 2)
							buffer[j] = inode->getVelSensitivity(3, m_grad - 1);
						else
							buffer[j] = 0.0;
					}
					else if (m_ndim == 3) {
						size_t j = i * 3; // node_counter * 3 (3 angular velocity components in 3D)
						if (iresponse.Size() > 5) {
							buffer[j] = inode->getVelSensitivity(4, m_grad - 1);
							buffer[j + 1] = inode->getVelSensitivity(5, m_grad - 1);
							buffer[j + 2] = inode->getVelSensitivity(6, m_grad - 1);
						}
						else {
							buffer[j] = 0.0;
							buffer[j + 1] = 0.0;
							buffer[j + 2] = 0.0;
						}
					}
				}
			}
		private:
			int m_grad;
		};

		class ResultRecorderAccelerationSensitivity : public ResultRecorder
		{
		public:
			ResultRecorderAccelerationSensitivity(const mpco::ProcessInfo &info, int grad)
				: ResultRecorder(info)
				, m_grad(grad)
			{
				std::stringstream ss_buffer;
				ss_buffer << "MODEL_STAGE[" << info.current_model_stage_id << "]/RESULTS/ON_NODES/ACCELERATION_SENSITIVITY_" << m_grad;
				m_result_name = ss_buffer.str();
				m_result_display_name = MPCO_MAKE_STRING("Acceleration Sensitivity d_A/d_q" << m_grad);
				m_num_components = 0;
				if (m_ndim == 1) {
					m_components_name = MPCO_MAKE_STRING("d_Ax/d_q" << m_grad);
					m_num_components = 1;
					m_result_data_type = mpco::ResultDataType::Scalar;
				}
				else if (m_ndim == 2) {
					m_components_name = MPCO_MAKE_STRING("d_Ax/d_q" << m_grad << ",d_Ay/d_q" << m_grad);
					m_num_components = 2;
					m_result_data_type = mpco::ResultDataType::Vectorial;
				}
				else if (m_ndim == 3) {
					m_components_name = MPCO_MAKE_STRING("d_Ax/d_q" << m_grad << ",d_Ay/d_q" << m_grad << ",d_Az/d_q" << m_grad);
					m_num_components = 3;
					m_result_data_type = mpco::ResultDataType::Vectorial;
				}
				m_description = MPCO_MAKE_STRING("Nodal acceleration sensitivity field d_A/d_q" << m_grad);
				m_result_type = mpco::ResultType::Generic;
			}
		protected:
			virtual void bufferResponse(mpco::ProcessInfo &info, std::vector<Node*> &nodes, std::vector<double> &buffer)const {
				if (m_grad < 1 || m_grad > info.domain->getNumParameters()) {
					opserr << "MPCORecorder sensitivity parameter is out of range: grad index = " << m_grad << ", domain->getNumParameters() = " << info.domain->getNumParameters() << "\n";
					return;
				}
				for (size_t i = 0; i < nodes.size(); i++) {
					Node *inode = nodes[i];
					size_t j = i * m_ndim;
					buffer[j] = inode->getVelSensitivity(1, m_grad - 1);
					if (m_ndim > 1) {
						buffer[j + 1] = inode->getVelSensitivity(2, m_grad - 1);
						if (m_ndim > 2) {
							buffer[j + 2] = inode->getVelSensitivity(3, m_grad - 1);
						}
					}
				}
			}
		private:
			int m_grad;
		};

		class ResultRecorderAngularAccelerationSensitivity : public ResultRecorder
		{
		public:
			ResultRecorderAngularAccelerationSensitivity(const mpco::ProcessInfo &info, int grad)
				: ResultRecorder(info)
				, m_grad(grad)
			{
				std::stringstream ss_buffer;
				ss_buffer << "MODEL_STAGE[" << info.current_model_stage_id << "]/RESULTS/ON_NODES/ANGULAR_ACCELERATION_SENSITIVITY_" << m_grad;;
				m_result_name = ss_buffer.str();
				m_result_display_name = MPCO_MAKE_STRING("Angular Acceleration Sensitivity d_RA/d_q" << m_grad);
				m_num_components = 0;
				if (m_ndim == 2) {
					m_components_name = MPCO_MAKE_STRING("d_RAz/d_q" << m_grad);
					m_num_components = 1;
					m_result_data_type = mpco::ResultDataType::Scalar;
				}
				else {
					m_components_name = MPCO_MAKE_STRING("d_RAx/d_q" << m_grad << ",d_RAy/d_q" << m_grad << ",d_RAz/d_q" << m_grad);
					m_num_components = 3;
					m_result_data_type = mpco::ResultDataType::Vectorial;
				}
				m_description = MPCO_MAKE_STRING("Nodal angular acceleration sensitivity field dRA/d_q" << m_grad);
				m_result_type = mpco::ResultType::Generic;
			}
		protected:
			virtual void bufferResponse(mpco::ProcessInfo &info, std::vector<Node*> &nodes, std::vector<double> &buffer)const {
				if (m_grad < 1 || m_grad > info.domain->getNumParameters()) {
					opserr << "MPCORecorder sensitivity parameter is out of range: grad index = " << m_grad << ", domain->getNumParameters() = " << info.domain->getNumParameters() << "\n";
					return;
				}
				for (size_t i = 0; i < nodes.size(); i++) {
					Node *inode = nodes[i];
					const Vector &iresponse = inode->getTrialVel(); // aux, just to check for angular accelerational components
					if (m_ndim == 2) {
						size_t j = i; // node_counter * 1 (1 angular acceleration component in 2D)
						if (iresponse.Size() > 2)
							buffer[j] = inode->getVelSensitivity(3, m_grad - 1);
						else
							buffer[j] = 0.0;
					}
					else if (m_ndim == 3) {
						size_t j = i * 3; // node_counter * 3 (3 angular acceleration components in 3D)
						if (iresponse.Size() > 5) {
							buffer[j] = inode->getVelSensitivity(4, m_grad - 1);
							buffer[j + 1] = inode->getVelSensitivity(5, m_grad - 1);
							buffer[j + 2] = inode->getVelSensitivity(6, m_grad - 1);
						}
						else {
							buffer[j] = 0.0;
							buffer[j + 1] = 0.0;
							buffer[j + 2] = 0.0;
						}
					}
				}
			}
		private:
			int m_grad;
		};

		typedef std::map<mpco::NodalResultType::Enum, ResultRecorder*> ResultRecorderMap;

	}

}

/*utilities for element results*/
namespace mpco {

	namespace element {

		struct OutputDescriptorHeader
		{
			OutputDescriptorHeader()
				: num_columns(0)
				, num_components()
				, gauss_id()
				, components_path()
				, components()
				, multiplicity()
			{}

			inline bool operator == (const OutputDescriptorHeader &other) const {
				if (this != &other) {
					if (*this < other) return false;
					if (*this > other) return false;
					return true;
				}
				return true;
			}

			inline bool operator < (const OutputDescriptorHeader &other) const {
				if (num_columns < other.num_columns) return true;
				if (num_columns > other.num_columns) return false;
				if (multiplicity.size() < other.multiplicity.size()) return true;
				if (multiplicity.size() > other.multiplicity.size()) return false;
				for (size_t i = 0; i < multiplicity.size(); i++) {
					if (multiplicity[i] < other.multiplicity[i]) return true;
					if (multiplicity[i] > other.multiplicity[i]) return false;
				}
				if (num_components.size() < other.num_components.size()) return true;
				if (num_components.size() > other.num_components.size()) return false;
				for (size_t i = 0; i < num_components.size(); i++) {
					if (num_components[i] < other.num_components[i]) return true;
					if (num_components[i] > other.num_components[i]) return false;
				}
				if (gauss_id.size() < other.gauss_id.size()) return true;
				if (gauss_id.size() > other.gauss_id.size()) return false;
				for (size_t i = 0; i < gauss_id.size(); i++) {
					if (gauss_id[i] < other.gauss_id[i]) return true;
					if (gauss_id[i] > other.gauss_id[i]) return false;
				}
				if (components_path.size() < other.components_path.size()) return true;
				if (components_path.size() > other.components_path.size()) return false;
				for (size_t i = 0; i < components_path.size(); i++) {
					if (components_path[i] < other.components_path[i]) return true;
					if (components_path[i] > other.components_path[i]) return false;
				}
				if (components.size() < other.components.size()) return true;
				if (components.size() > other.components.size()) return false;
				for (size_t i = 0; i < components.size(); i++) {
					if (components[i] < other.components[i]) return true;
					if (components[i] > other.components[i]) return false;
				}
				return false;
			}

			inline bool operator > (const OutputDescriptorHeader &other) const {
				if (num_columns > other.num_columns) return true;
				if (num_columns < other.num_columns) return false;
				if (multiplicity.size() > other.multiplicity.size()) return true;
				if (multiplicity.size() < other.multiplicity.size()) return false;
				for (size_t i = 0; i > multiplicity.size(); i++) {
					if (multiplicity[i] > other.multiplicity[i]) return true;
					if (multiplicity[i] < other.multiplicity[i]) return false;
				}
				if (num_components.size() > other.num_components.size()) return true;
				if (num_components.size() < other.num_components.size()) return false;
				for (size_t i = 0; i > num_components.size(); i++) {
					if (num_components[i] > other.num_components[i]) return true;
					if (num_components[i] < other.num_components[i]) return false;
				}
				if (gauss_id.size() > other.gauss_id.size()) return true;
				if (gauss_id.size() < other.gauss_id.size()) return false;
				for (size_t i = 0; i > gauss_id.size(); i++) {
					if (gauss_id[i] > other.gauss_id[i]) return true;
					if (gauss_id[i] < other.gauss_id[i]) return false;
				}
				if (components_path.size() > other.components_path.size()) return true;
				if (components_path.size() < other.components_path.size()) return false;
				for (size_t i = 0; i > components_path.size(); i++) {
					if (components_path[i] > other.components_path[i]) return true;
					if (components_path[i] < other.components_path[i]) return false;
				}
				if (components.size() > other.components.size()) return true;
				if (components.size() < other.components.size()) return false;
				for (size_t i = 0; i > components.size(); i++) {
					if (components[i] > other.components[i]) return true;
					if (components[i] < other.components[i]) return false;
				}
				return false;
			}

			inline std::string toString()const {
				std::stringstream ss;
				ss << "MPCORecorder Element Output Descriptor Header:\n";
				ss << "Columns: " << num_columns << "\n";
				ss << "N. Components:\n";
				for (size_t i = 0; i < num_components.size(); i++) {
					ss << num_components[i] << "   ";
				}
				ss << "\n\nPaths:\n";
				for (std::vector<std::vector<int> >::const_iterator
					it1 = components_path.begin(); it1 != components_path.end(); ++it1) {
					for (std::vector<int>::const_iterator
						it2 = it1->begin(); it2 != it1->end(); ++it2) {
						ss << mpco::ElementOutputDescriptorType::toString((mpco::ElementOutputDescriptorType::Enum)*it2) << ".";
					}
					ss << "\n";
				}
				ss << "\nComponents:\n";
				for (std::vector<std::vector<std::string> >::const_iterator
					it1 = components.begin(); it1 != components.end(); ++it1) {
					for (std::vector<std::string>::const_iterator
						it2 = it1->begin(); it2 != it1->end(); ++it2) {
						ss << *it2 << "   ";
					}
					ss << "\n";
				}
				return ss.str();
			}

			void workaroundForSizeInconsistency(int data_size) {
				/** WORKAROUND
				ouch! in some cases (see some material models), the material in the setResponse method, gives a vector of size N
				but without specifying the response type, that we use to get the components!. In those cases there will be
				a mismatch between response->getInformation().getData().size() and header.num_components !!!
				note: we can apply this workaround if and only if all items have zero components
				*/
				bool all_zero = true;
				for (size_t i = 0; i < num_components.size(); i++) {
					if (num_components[i] != 0) {
						all_zero = false;
						break;
					}
				}
				if (all_zero && data_size > 0) {
					if (num_components.size() == 0) { // empty
						num_components.resize(1);
						num_components[0] = data_size;
						gauss_id.resize(1);
						gauss_id[0] = -1;
						components_path.resize(1);
						components_path[0].push_back(mpco::ElementOutputDescriptorType::Element);
						components.resize(1);
						components[0].resize((size_t)data_size);
						for (size_t i = 0; i < (size_t)data_size; i++) {
							std::stringstream aux; aux << "C" << i + 1; components[0][i] = aux.str();
						}
						multiplicity.resize(1);
						multiplicity[0] = 1;
						num_columns = data_size;
					}
					else {
						int num_items = (int)num_components.size();
						/* only if we can equally divide the data size by the number of items,
						i.e. we assume that each gauss point have the same components
						*/
						if (data_size % num_items == 0) {
							int n_comp = data_size / num_items;
							for (size_t i = 0; i < num_components.size(); i++) {
								num_components[i] = n_comp;
								multiplicity[i] = 1;
								components[i].resize((size_t)n_comp);
								for (int j = 0; j < n_comp; j++) {
									std::stringstream aux; aux << "C" << j + 1; components[i][j] = aux.str();
								}
							}
							// we leave gauss and comp paths as they are
							num_columns = data_size;
						}
					}
				}
			}

			void workaroundForDuplicatedComponents()
			{
				for (size_t i = 0; i < components.size(); i++) {
					std::vector<std::string> &i_components = components[i];
					if (i_components.size() > 0) {
						std::map<std::string, int> aux;
						for (size_t j = 0; j < i_components.size(); j++) {
							int &dupl_counter = aux[i_components[j]];
							if (dupl_counter > 0) {
								std::stringstream ss;
								ss << i_components[j] << "(" << dupl_counter << ")";
								i_components[j] = ss.str();
							}
							dupl_counter++;
						}
					}
				}
			}

			int num_columns;
			std::vector<int> num_components;
			std::vector<int> gauss_id; // note: 0-based
			std::vector<std::vector<int> > components_path;
			std::vector<std::vector<std::string> > components;
			std::vector<int> multiplicity;
		};

		class OutputDescriptor
		{
		public:
			OutputDescriptor()
				: type(mpco::ElementOutputDescriptorType::Element)
				, tag(0)
				, dummy_section_flag(false)
				, gp_number(0)
				, gp_eta(0.0)
				, fib_y(0.0)
				, fib_z(0.0)
				, fib_a(0.0)
				, components()
				, items()
			{}

			OutputDescriptor(const OutputDescriptor &other)
				: type(other.type)
				, tag(other.tag)
				, dummy_section_flag(other.dummy_section_flag)
				, gp_number(other.gp_number)
				, gp_eta(other.gp_eta)
				, fib_y(other.fib_y)
				, fib_z(other.fib_z)
				, fib_a(other.fib_a)
				, components(other.components)
				, items()
			{
				// make a deep copy of subitems
				items.resize(other.items.size());
				for (size_t i = 0; i < other.items.size(); i++)
					items[i] = new OutputDescriptor(*other.items[i]);
			}

			~OutputDescriptor() {
				for (size_t i = 0; i < items.size(); i++)
					if (items[i])
						delete items[i];
			}

			OutputDescriptor &operator = (const OutputDescriptor &other) {
				if (this != &other) {
					type = other.type;
					tag = other.tag;
					dummy_section_flag = other.dummy_section_flag;
					gp_number = other.gp_number;
					gp_eta = other.gp_eta;
					fib_y = other.fib_y;
					fib_z = other.fib_z;
					fib_a = other.fib_a;
					components = other.components;
					// make a deep copy of subitems
					items.resize(other.items.size());
					for (size_t i = 0; i < other.items.size(); i++)
						items[i] = new OutputDescriptor(*other.items[i]);
				}
				return *this;
			}

		public:
			// generic
			mpco::ElementOutputDescriptorType::Enum type;
			// for material or section
			int tag;
			bool dummy_section_flag;
			// for gauss point
			int gp_number;
			double gp_eta; // use only eta (for beams)
						   // for fibers
			double fib_y;
			double fib_z;
			double fib_a;
			// components (if items.size() == 0)
			std::vector<std::string> components;
			// subitems
			std::vector<OutputDescriptor*> items;

		public:
			std::string toString()const {
				std::stringstream ss;
				/*
				this prints a xml-like structure, just for debugging purposes
				*/
				ss << "----------------------------------------------------\n";
				ss << "OutputDescriptor info\n";
				ss << "----------------------------------------------------\n";
				ss << "XML-like structure\n";
				printInfo(0, ss);
				ss << "----------------------------------------------------\n";
				return ss.str();
			}

			void getGaussLocations(std::vector<double> &x) const {
				x.clear();
				appendGaussLocation(x);
				if (x.size() == 1) {
					x[0] = 0.0;
				}
				else {
					double xmin = std::numeric_limits<double>::max();
					double xmax = -xmin;
					for (size_t i = 0; i < x.size(); i++) {
						double ieta = x[i];
						if (ieta < xmin)
							xmin = ieta;
						else if (ieta > xmax)
							xmax = ieta;
					}
					double span = xmax - xmin;
					if (span == 0.0) {
						for (size_t i = 0; i < x.size(); i++)
							x[i] = 0.0;
					}
					else {
						for (size_t i = 0; i < x.size(); i++)
							x[i] = 2.0*(x[i] - xmin) / span - 1.0;
					}
				}
			}

			void getFiberData(std::vector<mpco::element::FiberData> &data,
				std::vector<int> &data_mat_id,
				std::vector<int> &sec_id,
				std::vector<int> &gp_id,
				std::vector<bool> &dummy_sec_flags)const {
				data.clear();
				data_mat_id.clear();
				sec_id.clear();
				gp_id.clear();
				int *temp_gp = 0;
				int *temp_sec = 0;
				bool *temp_dummy = 0;
				appendFiberData(data, data_mat_id,
					sec_id, gp_id, dummy_sec_flags, temp_sec, temp_gp, temp_dummy);
				if (temp_gp)
					delete temp_gp;
				if (temp_sec)
					delete temp_sec;
				if (temp_dummy)
					delete temp_dummy;
			}

			mpco::element::OutputDescriptorHeader makeHeader()const {
				mpco::element::OutputDescriptorHeader header;
				std::list<int> temp_path;
				int temp_gp_id(-1);
				makeHeaderInternal(header, temp_path, temp_gp_id);
				return header;
			}

			void fixFloatingFiberOutput() {
				/*
				for some reason, some fiber-based cross section do not write the SectionOutput-tag before the FiberOutput-tag.
				this make things complicated and not robust. this is a workaround to fix this problem.
				*/
				fixFloatingFiberOutputInternal();
			}

			int getNextGpTag() {
				int next_gp_tag = -1;
				getNextGpTagInternal(next_gp_tag);
				return next_gp_tag + 1;
			}

			void purge() {
				mergeGaussInternal();
				mergeSecInternal();
			}

		private:
			void printInfo(int level, std::stringstream &ss) const {
				std::stringstream ss_indent;
				for (int i = 0; i < level; i++) ss_indent << "\t";
				std::string indent = ss_indent.str();
				ss << indent << "<" << mpco::ElementOutputDescriptorType::toString(this->type);
				if (this->type == mpco::ElementOutputDescriptorType::Gauss) {
					ss << " number=\"" << this->gp_number << "\" eta=\"" << this->gp_eta << "\"";
				}
				else if (this->type == mpco::ElementOutputDescriptorType::Section) {
					ss << " tag=\"" << this->tag << "\"";
				}
				else if (this->type == mpco::ElementOutputDescriptorType::Material) {
					ss << " tag=\"" << this->tag << "\"";
				}
				ss << ">\n";
				for (int i = 0; i < this->components.size(); i++) {
					ss << indent << "\t" << components[i] << "\n";
				}
				for (int i = 0; i < items.size(); i++) {
					items[i]->printInfo(level + 1, ss);
				}
				ss << indent << "</" << mpco::ElementOutputDescriptorType::toString(this->type) << ">\n";
			}

			void makeHeaderInternal(mpco::element::OutputDescriptorHeader &header, std::list<int> &temp_path, int &temp_gp_id) const {
				if (type == mpco::ElementOutputDescriptorType::Gauss)
					temp_gp_id = gp_number;
				temp_path.push_back((int)type);
				if (components.size() > 0 || items.size() == 0) {
					header.num_columns += (int)components.size();
					int next_num_components = (int)components.size();
					int next_gauss_id = temp_gp_id;
					std::vector<int> next_components_path(temp_path.begin(), temp_path.end());
					if (header.multiplicity.size() == 0) {
						header.num_components.push_back(next_num_components);
						header.gauss_id.push_back(next_gauss_id);
						header.components_path.push_back(next_components_path);
						header.components.push_back(components);
						header.multiplicity.push_back(1);
					}
					else {
						bool equal_to_previous = false;
						size_t last_index = header.multiplicity.size() - 1;
						if (header.num_components[last_index] == next_num_components) {
							if (header.gauss_id[last_index] == next_gauss_id) {
								if (utils::misc::areVectorsEqual(header.components_path[last_index], next_components_path)) {
									if (utils::misc::areVectorsEqual(header.components[last_index], components)) {
										equal_to_previous = true;
									}
								}
							}
						}
						if (equal_to_previous) {
							header.multiplicity[last_index]++;
						}
						else {
							header.num_components.push_back(next_num_components);
							header.gauss_id.push_back(next_gauss_id);
							header.components_path.push_back(next_components_path);
							header.components.push_back(components);
							header.multiplicity.push_back(1);
						}
					}
				}
				for (size_t i = 0; i < items.size(); i++)
					items[i]->makeHeaderInternal(header, temp_path, temp_gp_id);
				temp_path.pop_back();
			}

			void appendGaussLocation(std::vector<double> &x) const {
				if (type == mpco::ElementOutputDescriptorType::Gauss)
					x.push_back(gp_eta);
				for (size_t i = 0; i < items.size(); i++)
					items[i]->appendGaussLocation(x);
			}

			void appendFiberData(std::vector<mpco::element::FiberData> &data, std::vector<int> &data_mat_id,
				std::vector<int> &sec_id, std::vector<int> &gp_id, std::vector<bool> &dummy_sec_flags,
				int* &temp_sec, int* &temp_gp, bool* &temp_dummy) const {
				if (type == mpco::ElementOutputDescriptorType::Gauss) {
					if (temp_gp == 0)
						temp_gp = new int();
					*temp_gp = gp_number;
				}
				else if (type == mpco::ElementOutputDescriptorType::Section) {
					if (temp_sec == 0)
						temp_sec = new int();
					if (temp_dummy == 0)
						temp_dummy = new bool();
					*temp_sec = tag;
					*temp_dummy = dummy_section_flag;
				}
				else if (type == mpco::ElementOutputDescriptorType::Fiber) {
					data.push_back(mpco::element::FiberData(fib_y, fib_z, fib_a));
					int fiber_mat_tag = -1;
					if (items.size() == 1) {
						if (items[0]->type == mpco::ElementOutputDescriptorType::Material) {
							fiber_mat_tag = items[0]->tag;
						}
					}
					data_mat_id.push_back(fiber_mat_tag);
				}
				if (type == mpco::ElementOutputDescriptorType::Fiber || items.size() == 0) {
					if (temp_sec)
						sec_id.push_back(*temp_sec);
					if (temp_gp)
						gp_id.push_back(*temp_gp);
					if (temp_dummy)
						dummy_sec_flags.push_back(*temp_dummy);
					/** \bug-fixed note here we MUST exit, avoiding asking to sub-items.
					for example if this item is a fiber, it may have sub-times like material!
					*/
					return;
				}
				// ask subitems
				for (size_t i = 0; i < items.size(); i++) {
					items[i]->appendFiberData(data, data_mat_id,
						sec_id, gp_id, dummy_sec_flags, temp_sec, temp_gp, temp_dummy);
				}
			}

			void fixFloatingFiberOutputInternal() {
				if (items.size() > 0) {
					if (type != mpco::ElementOutputDescriptorType::Section) {
						if (items[0]->type == mpco::ElementOutputDescriptorType::Fiber) {
							/* check only the first one. items are of the same type... */
							OutputDescriptor *dummy_section_level = new OutputDescriptor();
							dummy_section_level->type = mpco::ElementOutputDescriptorType::Section;
							dummy_section_level->tag = -123456;
							dummy_section_level->dummy_section_flag = true;
							dummy_section_level->items = items;
							items.clear();
							items.push_back(dummy_section_level);
						}
					}
					for (size_t i = 0; i < items.size(); i++)
						items[i]->fixFloatingFiberOutputInternal();
				}
			}

			void getNextGpTagInternal(int &next_gp_tag) {
				if (type == mpco::ElementOutputDescriptorType::Gauss) {
					if (next_gp_tag < gp_number)
						next_gp_tag = gp_number;
				}
				else {
					for (size_t i = 0; i < items.size(); i++)
						items[i]->getNextGpTagInternal(next_gp_tag);
				}
			}

			void mergeGaussInternal() {
				/* if multiple gauss items with same id exist, merge their contents
				\todo: check if the order of fibers is preserved..
				note mandatory now, because when asking for multiple sections, duplicate gauss points are produced
				but only the first one is filled
				*/
				if (items.size() > 0) {
					if (items[0]->type == mpco::ElementOutputDescriptorType::Gauss) {
						// all sub items are gauss descriptors, let's merge them
						std::map<int, OutputDescriptor*> aux;
						for (size_t i = 0; i < items.size(); i++) {
							OutputDescriptor* curr_item = items[i];
							std::map<int, OutputDescriptor*>::iterator it = aux.find(curr_item->gp_number);
							if (it == aux.end()) { // non existing
								aux[curr_item->gp_number] = curr_item;
							}
							else { // existing, move curr_item->items into existing_item->items
								OutputDescriptor* existing_item = it->second;
								for (size_t j = 0; j < curr_item->items.size(); j++)
									existing_item->items.push_back(curr_item->items[j]);
								curr_item->items.clear();
							}
						}
						items.clear();
						for (std::map<int, OutputDescriptor*>::iterator it = aux.begin(); it != aux.end(); ++it) {
							items.push_back(it->second);
						}
					}
					else {
						for (size_t i = 0; i < items.size(); i++)
							items[i]->mergeGaussInternal();
					}
				}
			}
			void mergeSecInternal() {
				if (items.size() > 0) {
					if (items[0]->type == mpco::ElementOutputDescriptorType::Section) {
						// all sub items are section descriptors, let's merge them
						std::map<int, OutputDescriptor*> aux;
						for (size_t i = 0; i < items.size(); i++) {
							OutputDescriptor* curr_item = items[i];
							std::map<int, OutputDescriptor*>::iterator it = aux.find(curr_item->tag);
							if (it == aux.end()) { // non existing
								aux[curr_item->tag] = curr_item;
							}
							else { // existing, move curr_item->items into existing_item->items
								OutputDescriptor* existing_item = it->second;
								for (size_t j = 0; j < curr_item->items.size(); j++)
									existing_item->items.push_back(curr_item->items[j]);
								curr_item->items.clear();
							}
						}
						items.clear();
						for (std::map<int, OutputDescriptor*>::iterator it = aux.begin(); it != aux.end(); ++it) {
							items.push_back(it->second);
						}
					}
					else {
						for (size_t i = 0; i < items.size(); i++)
							items[i]->mergeSecInternal();
					}
				}
			}
		};

		class OutputDescriptorStream : public OPS_Stream
		{
		public:
			OutputDescriptorStream(mpco::element::OutputDescriptor * _d)
				: OPS_Stream(OPS_STREAM_TAGS_MPCORecorder_ElementOutputDescriptorStream)
				, descr(_d)
				, current_level(0)
				, pending_close_tag(false)
			{}
			~OutputDescriptorStream() {}

			int tag(const char *name) {
				// get the element output descriptor at current level and id
				mpco::element::OutputDescriptor *eo_curr_lev = descr;
				for (int i = 1; i <= current_level; i++) {
					if (eo_curr_lev->items.size() == 0) {
						opserr << "MPCORecorder Error: cannot set attribute(name, int), empty item list.\n";
						exit(-1);
					}
					eo_curr_lev = eo_curr_lev->items[eo_curr_lev->items.size() - 1];
				}
				if (current_level == 0) {
					if (strcmp(name, "ElementOutput") == 0) {
						/** nothing to do. this is the root of the result tree*/
					}
					/* gauss output is the first entry */
					else if (strcmp(name, "GaussPoint") == 0 || strcmp(name, "GaussPointOutput") == 0) {
						mpco::element::OutputDescriptor *eo_new_curr_lev = new mpco::element::OutputDescriptor();
						eo_new_curr_lev->type = mpco::ElementOutputDescriptorType::Gauss;
						ensureItemsOfUniformType(eo_curr_lev, eo_new_curr_lev);
						eo_curr_lev->items.push_back(eo_new_curr_lev);
						current_level++;
					}
					else {
						/*opserr <<
						"MPCORecorder Error: invalid tag at level 0:\n"
						"expected \"GaussPoint\" or \"GaussPointOutput\", given \"" << name << "\"\n";
						exit(-1);*/
						/*
						let's try a last workaround for an inconsistency problem found in SSPbrick:
						in that element, in setResponse, the material's setResponse method is called without opening a GaussPoint tag...
						*/
						//simulate: tag("GaussPoint");
						{
							mpco::element::OutputDescriptor *eo_new_curr_lev = new mpco::element::OutputDescriptor();
							eo_new_curr_lev->type = mpco::ElementOutputDescriptorType::Gauss;
							ensureItemsOfUniformType(eo_curr_lev, eo_new_curr_lev);
							eo_curr_lev->items.push_back(eo_new_curr_lev);
							current_level++;
						}
						// create the attribute 
						attr("number", descr->getNextGpTag());
						/*
						close the gauss tag. we cannot do it here. Just set the pending_close_tag to true, so that
						when the 'real' tag gets closed, we automatically close the 'manual' gauss tag
						*/
						pending_close_tag = true;
						// recursion: recal this request now inside a gauss point tag
						tag(name);
					}
				}
				else if (current_level > 0) {
					if (strcmp(name, "NdMaterialOutput") == 0 || strcmp(name, "UniaxialMaterialOutput") == 0) {
						// its parent can be anything but ElementOutput, ok if current_level > 1
						mpco::element::OutputDescriptor *eo_new_curr_lev = new mpco::element::OutputDescriptor();
						eo_new_curr_lev->type = mpco::ElementOutputDescriptorType::Material;
						ensureItemsOfUniformType(eo_curr_lev, eo_new_curr_lev);
						eo_curr_lev->items.push_back(eo_new_curr_lev);
						current_level++;
					}
					else if (strcmp(name, "SectionOutput") == 0 || strcmp(name, "SectionForceDeformation") == 0) {
						// its parent can be GaussOutput or another SectionOutput
						if (!(eo_curr_lev->type == mpco::ElementOutputDescriptorType::Gauss || eo_curr_lev->type == mpco::ElementOutputDescriptorType::Section)) {
							opserr <<
								"MPCORecorder Error: invalid parent for \"" << name << "\" tag:\n"
								"expected \"GaussOutput\" or \"GaussPointOutput\""
								" or \"SectionOutput\" or \"SectionForceDeformation\", parent tag = \""
								<< mpco::ElementOutputDescriptorType::toString(eo_curr_lev->type) << "\"\n";
							exit(-1);
						}
						mpco::element::OutputDescriptor *eo_new_curr_lev = new mpco::element::OutputDescriptor();
						eo_new_curr_lev->type = mpco::ElementOutputDescriptorType::Section;
						ensureItemsOfUniformType(eo_curr_lev, eo_new_curr_lev);
						eo_curr_lev->items.push_back(eo_new_curr_lev);
						current_level++;
					}
					else if (strcmp(name, "FiberOutput") == 0) {
						// its parent can be only a SectionOutput
						if (!(eo_curr_lev->type == mpco::ElementOutputDescriptorType::Section || eo_curr_lev->type == mpco::ElementOutputDescriptorType::Gauss)) {
							opserr <<
								"MPCORecorder Error: invalid parent for \"" << name << "\" tag:\n"
								"expected \"GaussOutput\" or \"GaussPointOutput\""
								" or \"SectionOutput\" or \"SectionForceDeformation\", parent tag = \""
								<< mpco::ElementOutputDescriptorType::toString(eo_curr_lev->type) << "\"\n";
							exit(-1);
						}
						mpco::element::OutputDescriptor *eo_new_curr_lev = new mpco::element::OutputDescriptor();
						eo_new_curr_lev->type = mpco::ElementOutputDescriptorType::Fiber;
						ensureItemsOfUniformType(eo_curr_lev, eo_new_curr_lev);
						eo_curr_lev->items.push_back(eo_new_curr_lev);
						current_level++;
					}
					else {
						/*
						let's try a last workaround for an inconsistency problem found in ForceBeamColumn3d:
						in that element, in setResponse, if multiple sections are requested the GaussOutput tag is not
						closed...
						*/
						if ((strcmp(name, "GaussPoint") == 0 || strcmp(name, "GaussPointOutput") == 0) && current_level == 1) {
							endTag();
							tag(name);
						}
						else {
							opserr <<
								"MPCORecorder Error: invalid tag at level " << current_level << ":\n"
								"expected \"NdMaterialOutput\" or \"UniaxialMaterialOutput\""
								" or \"SectionOutput\" or \"SectionForceDeformation\""
								" or \"FiberOutput\""
								", given \"" << name << "\"\n";
							exit(-1);
						}
					}
				}
				return 0;
			};

			int tag(const char *name, const char *value) {
				mpco::element::OutputDescriptor *eo_curr_lev = descr;
				for (int i = 1; i <= current_level; i++) {
					if (eo_curr_lev->items.size() == 0) {
						opserr << "MPCORecorder Error: cannot set attribute(name, int), empty item list.\n";
						exit(-1);
					}
					eo_curr_lev = eo_curr_lev->items[eo_curr_lev->items.size() - 1];
				}
				if (strcmp(name, "ResponseType") == 0)
					eo_curr_lev->components.push_back(value);
				return 0;
			};

			int endTag() {
				if (current_level > 0) {
					// decrement the current level
					current_level--;
				}
				if (pending_close_tag) {
					// once more...
					if (current_level > 0) {
						// decrement the current level
						current_level--;
					}
					pending_close_tag = false;
				}
				return 0;
			};

			int attr(const char *name, int value) {
				if (current_level > 0) {
					mpco::element::OutputDescriptor *eo_curr_lev = descr;
					for (int i = 1; i <= current_level; i++) {
						if (eo_curr_lev->items.size() == 0) {
							opserr << "MPCORecorder Error: cannot set attribute(name, int), empty item list.\n";
							exit(-1);
						}
						eo_curr_lev = eo_curr_lev->items[eo_curr_lev->items.size() - 1];
					}
					if (eo_curr_lev->type == mpco::ElementOutputDescriptorType::Gauss) {
						if (strcmp(name, "number") == 0)
							eo_curr_lev->gp_number = value - 1; // note: make it 0-based
					}
					else if (eo_curr_lev->type == mpco::ElementOutputDescriptorType::Material) {
						if (strcmp(name, "tag") == 0 || strcmp(name, "matTag") == 0)
							eo_curr_lev->tag = value;
					}
					else if (eo_curr_lev->type == mpco::ElementOutputDescriptorType::Section) {
						if (strcmp(name, "tag") == 0 || strcmp(name, "secTag") == 0)
							eo_curr_lev->tag = value;
					}
				}
				return 0;
			};

			int attr(const char *name, double value) {
				if (current_level > 0) {
					mpco::element::OutputDescriptor *eo_curr_lev = descr;
					for (int i = 1; i <= current_level; i++) {
						if (eo_curr_lev->items.size() == 0) {
							opserr << "MPCORecorder Error: cannot set attribute(name, int), empty item list.\n";
							exit(-1);
						}
						eo_curr_lev = eo_curr_lev->items[eo_curr_lev->items.size() - 1];
					}
					if (eo_curr_lev->type == mpco::ElementOutputDescriptorType::Gauss) {
						if (strcmp(name, "eta") == 0)
							eo_curr_lev->gp_eta = value;
					}
					else if (eo_curr_lev->type == mpco::ElementOutputDescriptorType::Fiber) {
						if (strcmp(name, "yLoc") == 0)
							eo_curr_lev->fib_y = value;
						else if (strcmp(name, "zLoc") == 0)
							eo_curr_lev->fib_z = value;
						else if ((strcmp(name, "area") == 0) || (strcmp(name, "thickness") == 0))
							eo_curr_lev->fib_a = value;
					}
				}
				return 0;
			};

			int attr(const char *name, const char *value) { return 0; };
			int write(Vector &data) { return 0; };

			OPS_Stream& write(const char *s, int n) { return *this; };
			OPS_Stream& write(const unsigned char *s, int n) { return *this; };
			OPS_Stream& write(const signed char *s, int n) { return *this; };
			OPS_Stream& write(const void *s, int n) { return *this; };
			OPS_Stream& operator<<(char c) { return *this; };
			OPS_Stream& operator<<(unsigned char c) { return *this; };
			OPS_Stream& operator<<(signed char c) { return *this; };
			OPS_Stream& operator<<(const char *s) { return *this; };
			OPS_Stream& operator<<(const unsigned char *s) { return *this; };
			OPS_Stream& operator<<(const signed char *s) { return *this; };
			OPS_Stream& operator<<(const void *p) { return *this; };
			OPS_Stream& operator<<(int n) { return *this; };
			OPS_Stream& operator<<(unsigned int n) { return *this; };
			OPS_Stream& operator<<(long n) { return *this; };
			OPS_Stream& operator<<(unsigned long n) { return *this; };
			OPS_Stream& operator<<(short n) { return *this; };
			OPS_Stream& operator<<(unsigned short n) { return *this; };
			OPS_Stream& operator<<(bool b) { return *this; };
			OPS_Stream& operator<<(double n) { return *this; };
			OPS_Stream& operator<<(float n) { return *this; };

			int sendSelf(int commitTag, Channel &theChannel) { return 0; };
			int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) { return 0; };

		public:
			void finalizeSetResponse() {
				while (current_level > 0) {
					endTag();
				}
				/*
				this method was originally automatically called when a call to endTag determines
				the end of an ElementOutput. Unfortunately we found out that some elements (ForceBeamColumn3d)
				do not call endTag after a tag("GaussPointOutput") in case of output from multiple sections.
				So we need to call it manually.
				*/
				/*
				do some workaround...
				*/
				descr->fixFloatingFiberOutput();
			}

		private:
			void ensureItemsOfUniformType(mpco::element::OutputDescriptor *parent, mpco::element::OutputDescriptor *child) {
				if (parent->items.size() > 0) {
					if (child->type != parent->items.back()->type) {
						opserr << "MPCORecorder Error: (mpco::element::OutputDescriptor) "
							"Responses at the same level of the response tree must be of the same type.\n"
							"Expected: " << mpco::ElementOutputDescriptorType::toString(parent->items.back()->type)
							<< ", given: " << mpco::ElementOutputDescriptorType::toString(child->type) << "\n";
						exit(-1);
					}
				}
			}

		public:
			mpco::element::OutputDescriptor *descr;
			int current_level;
			bool pending_close_tag;
		};

		class OutputResponse
		{
		public:
			OutputResponse()
				: response(0)
				, element(0)
			{}
			OutputResponse(Element *_elem, Response *_resp)
				: response(_resp)
				, element(_elem)
			{}
			Response *response;
			Element *element;
		};

		struct OutputResponseCollection
		{
			OutputResponseCollection()
				: is_new(true)
				, dir_name("")
				, initialized(false)
				, items() {}
			bool is_new;
			std::string dir_name;
			bool initialized;
			std::vector<OutputResponse> items;
		};

		struct OutputWithSameCustomIntRuleCollection
		{
			typedef std::map<OutputDescriptorHeader, OutputResponseCollection> collection_type;
			collection_type items;
		};

		struct OutputWithSameIntRuleCollection
		{
			typedef std::map<int, OutputWithSameCustomIntRuleCollection> collection_type;
			collection_type items;
		};

		struct OutputWithSameClassTagCollection
		{
			typedef std::map<ElementIntegrationRuleType::Enum, OutputWithSameIntRuleCollection> collection_type;
			collection_type items;
		};

		struct ResultRecorder
		{
			typedef std::map<int, OutputWithSameClassTagCollection> collection_type;
			ResultRecorder()
				: initialized(false)
				, result_request()
				, response_map()
			{}
			bool initialized;
			std::vector<std::string> result_request;
			collection_type response_map;
		};

		typedef std::vector<ResultRecorder> ResultRecorderCollection;

		/*************************************************************************************

		utilities for element mapping based on:
		class tag
		integration rule
		default or custom integration rule
		<element group>

		In this way all elements in <element group> share the same:
		1) number of nodes, 2) number and location of integration points

		**************************************************************************************/

		struct ElementIntegrationRule
		{
			ElementIntegrationRule()
				: int_rule_type(ElementIntegrationRuleType::CustomIntegrationRule), x()
			{}
			ElementIntegrationRule(ElementIntegrationRuleType::Enum _int_rule_type)
				: int_rule_type(_int_rule_type), x()
			{}
			inline bool operator < (const ElementIntegrationRule &other) const {
				const double rel_tol = 1.0e-5;
				if (int_rule_type < other.int_rule_type) return true;
				if (int_rule_type > other.int_rule_type) return false;
				if (x.size() < other.x.size()) return true;
				if (x.size() > other.x.size()) return false;
				for (size_t i = 0; i < x.size(); i++) {
					double tol = std::max(std::abs(x[i]), std::abs(other.x[i]))*rel_tol;
					if (utils::misc::lessThanWithTol(x[i], other.x[i], rel_tol)) return true;
					if (utils::misc::greaterThanWithTol(x[i], other.x[i], rel_tol)) return false;
					// continue with the loop
				}
				return false; // everything is equal
			}
			inline bool operator > (const ElementIntegrationRule &other) const {
				const double rel_tol = 1.0e-5;
				if (int_rule_type > other.int_rule_type) return true;
				if (int_rule_type < other.int_rule_type) return false;
				if (x.size() > other.x.size()) return true;
				if (x.size() < other.x.size()) return false;
				for (size_t i = 0; i < x.size(); i++) {
					double tol = std::max(std::abs(x[i]), std::abs(other.x[i]))*rel_tol;
					if (utils::misc::greaterThanWithTol(x[i], other.x[i], rel_tol)) return true;
					if (utils::misc::lessThanWithTol(x[i], other.x[i], rel_tol)) return false;
					// continue with the loop
				}
				return false; // everything is equal
			}
			ElementIntegrationRuleType::Enum int_rule_type;
			std::vector<double> x;
		};

		struct ElementWithSameCustomIntRuleCollection
		{
			typedef std::vector<Element*> collection_type;

			ElementWithSameCustomIntRuleCollection()
				: is_new(true), custom_int_rule_index(0), name(""), items()
			{}

			bool is_new;
			int custom_int_rule_index;
			std::string name;
			std::vector<Element*> items;
		};

		struct ElementWithSameIntRuleCollection
		{
			typedef std::map<int, ElementWithSameCustomIntRuleCollection> submap_type;

			ElementWithSameIntRuleCollection()
				: is_new(true), int_rule_type(ElementIntegrationRuleType::CustomIntegrationRule), items()
			{}

			bool is_new;
			ElementIntegrationRuleType::Enum int_rule_type;
			std::map<int, ElementWithSameCustomIntRuleCollection> items;
		};

		struct ElementWithSameClassTagCollection
		{
			typedef std::map<ElementIntegrationRuleType::Enum, ElementWithSameIntRuleCollection> submap_type;

			ElementWithSameClassTagCollection()
				: is_new(true), class_tag(0), class_name("unknown"), num_nodes(0), geom_type(ElementGeometryType::Custom), items() {}

			bool is_new;
			int class_tag;
			std::string class_name;
			int num_nodes;
			ElementGeometryType::Enum geom_type;
			std::map<mpco::ElementIntegrationRuleType::Enum, ElementWithSameIntRuleCollection> items;
		};

		struct ElementCollection
		{
			typedef std::map<int, ElementWithSameClassTagCollection> submap_type;

			ElementCollection()
				: registered_custom_rules(), items()
			{}

			void mapElements(Domain *d, bool has_region, const std::vector<int> &subset) {
				/*
				clear previous mappings
				*/
				registered_custom_rules.clear();
				items.clear();
				/*
				quick return
				*/
				if (d == 0) return;
				/*
				auxiliary map to register custom rules
				*/
				std::map<ElementIntegrationRule, int> aux_map_custom_rules;
				/*
				loop over all elements in the domain
				*/
				size_t subset_elem_counter(0);
				ElementIter* element_iter = &(d->getElements());
				Element* current_element = 0;
				//while ((current_element = (*element_iter)()) != 0) {
				while (true) {
					/*
					get next element
					*/
					if (has_region) {
						if (subset_elem_counter == subset.size())
							break;
						current_element = d->getElement(subset[subset_elem_counter++]);
						if (current_element == 0)
							continue; // skip null and go to next iteration
					}
					else {
						current_element = (*element_iter)();
						if (current_element == 0) 
							break; 
					}
					/*
					get class tag, geometry and integration rule type
					*/
					int elem_type = current_element->getClassTag();
					ElementGeometryType::Enum geom_type;
					ElementIntegrationRuleType::Enum int_rule_type;
					getGeometryAndIntRuleByClassTag(elem_type, geom_type, int_rule_type);
					/*
					map by class tag
					*/
					ElementWithSameClassTagCollection &elem_coll_by_tag = items[elem_type];
					if (elem_coll_by_tag.is_new) {
						elem_coll_by_tag.class_tag = elem_type;
						elem_coll_by_tag.class_name = current_element->getClassType();
						elem_coll_by_tag.num_nodes = current_element->getNumExternalNodes();
						elem_coll_by_tag.geom_type = geom_type;
						elem_coll_by_tag.is_new = false;
					}
					/*
					make sure that every element with the same tag have the same number of nodes
					*/
					if (current_element->getNumExternalNodes() != elem_coll_by_tag.num_nodes) {
						opserr << "MPCORecorder Error while mapping elements: elements with different number of nodes "
							"exist withing the same class tag. This is not supported\n";
						exit(-1);
					}
					/*
					create the integration rule
					*/
					ElementIntegrationRule int_rule(int_rule_type);
					if (int_rule_type == ElementIntegrationRuleType::CustomIntegrationRule)
						getCustomGaussPointLocations(current_element, int_rule);
					/*
					if this is a custom rule, register it
					*/
					int custom_int_rule_index = 0;
					if (int_rule_type == ElementIntegrationRuleType::CustomIntegrationRule) {
						std::map<ElementIntegrationRule, int>::iterator reg_int_rule_iter = aux_map_custom_rules.find(int_rule);
						if (reg_int_rule_iter == aux_map_custom_rules.end()) {
							// new rule, define a new index starting from 1
							custom_int_rule_index = (int)aux_map_custom_rules.size() + 1;
							aux_map_custom_rules[int_rule] = custom_int_rule_index;
						}
						else {
							custom_int_rule_index = reg_int_rule_iter->second;
						}
					}
					/*
					map by integration rule
					*/
					ElementWithSameIntRuleCollection &elem_coll_by_rule = elem_coll_by_tag.items[int_rule_type];
					if (elem_coll_by_rule.is_new) {
						elem_coll_by_rule.int_rule_type = int_rule_type;
						elem_coll_by_rule.is_new = false;
					}
					/*
					map by custom integration rule index
					*/
					ElementWithSameCustomIntRuleCollection &elem_coll_by_custom_rule = elem_coll_by_rule.items[custom_int_rule_index];
					if (elem_coll_by_custom_rule.is_new) {
						elem_coll_by_custom_rule.custom_int_rule_index = custom_int_rule_index;
						elem_coll_by_custom_rule.is_new = false;
					}
					/*
					finally add this element
					*/
					elem_coll_by_custom_rule.items.push_back(current_element);
				}
				/*
				now fill the custom integration rule map
				*/
				for (std::map<ElementIntegrationRule, int>::iterator it = aux_map_custom_rules.begin();
					it != aux_map_custom_rules.end(); ++it) {
					registered_custom_rules[it->second] = it->first;
				}
			}

			void getGeometryAndIntRuleByClassTag(
				int elem_class_tag,
				ElementGeometryType::Enum &geom_type,
				ElementIntegrationRuleType::Enum &int_type) {
				/*
				set default values. custom geometry (i.e. point cloud)
				and no integration rule
				*/
				geom_type = ElementGeometryType::Custom;
				int_type = ElementIntegrationRuleType::NoIntegrationRule;
				/*
				2-node line with 1 gp
				*/
				if (
					// ./adpter actuators
					elem_class_tag == ELE_TAG_Actuator ||
					elem_class_tag == ELE_TAG_ActuatorCorot ||
					// ./truss
					elem_class_tag == ELE_TAG_Truss ||
					elem_class_tag == ELE_TAG_Truss2 ||
					elem_class_tag == ELE_TAG_TrussSection ||
					elem_class_tag == ELE_TAG_CorotTruss ||
					elem_class_tag == ELE_TAG_CorotTruss2 ||
					elem_class_tag == ELE_TAG_CorotTrussSection ||
					// ./zeroLength
					elem_class_tag == ELE_TAG_ZeroLength ||
					elem_class_tag == ELE_TAG_ZeroLengthSection ||
					elem_class_tag == ELE_TAG_ZeroLengthND ||
					elem_class_tag == ELE_TAG_CoupledZeroLength ||
					elem_class_tag == ELE_TAG_ZeroLengthRocking ||
					elem_class_tag == ELE_TAG_ZeroLengthContact2D ||
					elem_class_tag == ELE_TAG_ZeroLengthContact3D ||
					elem_class_tag == ELE_Tag_ZeroLengthImpact3D ||
					// ./elasticBeamColumn
					elem_class_tag == ELE_TAG_ElasticBeam2d ||
					elem_class_tag == ELE_TAG_ElasticBeam3d ||
					elem_class_tag == ELE_TAG_ElasticTimoshenkoBeam2d ||
					elem_class_tag == ELE_TAG_ElasticTimoshenkoBeam3d ||
					elem_class_tag == ELE_TAG_ModElasticBeam2d ||
					// .elastomericBearing
					elem_class_tag == ELE_TAG_ElastomericBearingBoucWen2d ||
					elem_class_tag == ELE_TAG_ElastomericBearingBoucWen3d ||
					elem_class_tag == ELE_TAG_ElastomericBearingBoucWenMod3d ||
					elem_class_tag == ELE_TAG_ElastomericBearingPlasticity2d ||
					elem_class_tag == ELE_TAG_ElastomericBearingPlasticity3d ||
					elem_class_tag == ELE_TAG_ElastomericBearingUFRP2d ||
					elem_class_tag == ELE_TAG_ElastomericBearingUFRP3d ||
					elem_class_tag == ELE_TAG_ElastomericX ||
					elem_class_tag == ELE_TAG_HDR ||
					elem_class_tag == ELE_TAG_LeadRubberX ||
					// ./ulBeamColumn
					/*warning: these two could go to the beam with custom integration rule, but they do not define everything properly! check in future versions */
					elem_class_tag == Ele_TAG_Elastic2dGNL || /*warning: no integrationPoints, no resp for all sections, no section response*/
					elem_class_tag == TAG_InelasticYS2DGNL /*warning: no integrationPoints, no resp for all sections, no section response*/
														   // todo... others
					) {
					geom_type = ElementGeometryType::Line_2N;
					int_type = ElementIntegrationRuleType::Line_GaussLegendre_1;
				}
				/*
				2-node beams with custom number of gp
				*/
				else if (
					// ./dispBeamColumn
					elem_class_tag == ELE_TAG_DispBeamColumn2d || /*warning: no integrationPoints*/
					elem_class_tag == ELE_TAG_DispBeamColumn3d || /*warning: no integrationPoints*/
					elem_class_tag == ELE_TAG_DispBeamColumn2dWithSensitivity || /*warning: no integrationPoints, no resp for all sections*/
					elem_class_tag == ELE_TAG_DispBeamColumn3dWithSensitivity || /*warning: no integrationPoints, no resp for all sections*/
																				 // ./dispBeamColumnInt
					elem_class_tag == ELE_TAG_DispBeamColumn2dInt || /*warning: no integrationPoints, no resp for all sections*/
					elem_class_tag == ELE_TAG_DispBeamColumn2dThermal || /*warning: no integrationPoints, no resp for all sections*/
																		 // ./forceBeamColumn
					elem_class_tag == ELE_TAG_ElasticForceBeamColumn2d || /*warning: no resp for all sections*/
					elem_class_tag == ELE_TAG_ElasticForceBeamColumn3d || /*warning: no resp for all sections*/
					elem_class_tag == ELE_TAG_ElasticForceBeamColumnWarping2d || /*warning: no resp for all sections*/
					elem_class_tag == ELE_TAG_ForceBeamColumn2d || /* <- OK! this one defines everything ! good job*/
					elem_class_tag == ELE_TAG_ForceBeamColumn3d || /* <- OK! this one defines everything ! good job*/
					elem_class_tag == ELE_TAG_ForceBeamColumnCBDI2d || /* <- OK! this one defines everything ! good job*/
					elem_class_tag == ELE_TAG_ForceBeamColumnWarping2d /* <- OK! this one defines everything ! good job*/
					) {
					geom_type = ElementGeometryType::Line_2N;
					int_type = ElementIntegrationRuleType::CustomIntegrationRule;
				}
				/*
				3-node triangle with 1 gp
				*/
				else if (
					// ./triangle
					elem_class_tag == ELE_TAG_Tri31
					) {
					geom_type = ElementGeometryType::Triangle_3N;
					int_type = ElementIntegrationRuleType::Triangle_GaussLegendre_1;
				}
				/*
				3-node triangle with 3 gp
				*/
				else if (
					// ./shell
					elem_class_tag == ELE_TAG_ShellANDeS
					) {
					geom_type = ElementGeometryType::Triangle_3N;
					int_type = ElementIntegrationRuleType::Triangle_GaussLegendre_2B;
				}
				/*
				3-node triangle with 4 gp
				*/
				else if (
					// ./shell
					elem_class_tag == ELE_TAG_ShellDKGT ||
					elem_class_tag == ELE_TAG_ShellNLDKGT ||
					elem_class_tag == ELE_TAG_ASDShellT3
					) {
					geom_type = ElementGeometryType::Triangle_3N;
					int_type = ElementIntegrationRuleType::Triangle_GaussLegendre_2C;
				}
				/*
				4-node quadrilateral with 1 gp
				*/
				else if (
					// ./UWelements
					elem_class_tag == ELE_TAG_SSPquad ||
					elem_class_tag == ELE_TAG_SSPquadUP
					)
				{
					geom_type = ElementGeometryType::Quadrilateral_4N;
					int_type = ElementIntegrationRuleType::Quadrilateral_GaussLegendre_1;
				}
				/*
				4-node quadrilateral with 2x2 gp
				*/
				else if (
					// ./quad
					elem_class_tag == ELE_TAG_ConstantPressureVolumeQuad ||
					elem_class_tag == ELE_TAG_EnhancedQuad ||
					elem_class_tag == ELE_TAG_FourNodeQuad ||
					elem_class_tag == ELE_TAG_FourNodeQuad3d ||
					elem_class_tag == ELE_TAG_FourNodeQuadWithSensitivity ||
					// ./shell
					elem_class_tag == ELE_TAG_ShellDKGQ ||
					elem_class_tag == ELE_TAG_ShellNLDKGQ ||
					elem_class_tag == ELE_TAG_ShellMITC4 ||
					elem_class_tag == ELE_TAG_ShellMITC4Thermal ||
					elem_class_tag == ELE_TAG_ASDShellQ4 ||
					// ./up
					elem_class_tag == ELE_TAG_BBarFourNodeQuadUP ||
					elem_class_tag == ELE_TAG_FourNodeQuadUP
					) {
					geom_type = ElementGeometryType::Quadrilateral_4N;
					int_type = ElementIntegrationRuleType::Quadrilateral_GaussLegendre_2;
				}
				/*
				9-node quadrilateral with 3x3 gp
				*/
				else if (
					// ./quad
					elem_class_tag == ELE_TAG_NineNodeMixedQuad ||
					// ./shell
					elem_class_tag == ELE_TAG_ShellMITC9 ||
					// ./up
					elem_class_tag == ELE_TAG_Nine_Four_Node_QuadUP
					) {
					geom_type = ElementGeometryType::Quadrilateral_9N;
					int_type = ElementIntegrationRuleType::Quadrilateral_GaussLegendre_3;
				}
				/*
				4-node tetrahedron with 1x1x1 gp
				*/
				else if (
					// ./tetrahedron
					elem_class_tag == ELE_TAG_FourNodeTetrahedron
					)
				{
					geom_type = ElementGeometryType::Tetrahedron_4N;
					int_type = ElementIntegrationRuleType::Tetrahedron_GaussLegendre_1;
				}
				/*
				8-node hexahedron with 1x1x1 gp
				*/
				else if (
					// ./UWelements
					elem_class_tag == ELE_TAG_SSPbrick ||
					elem_class_tag == ELE_TAG_SSPbrickUP
					)
				{
					geom_type = ElementGeometryType::Hexahedron_8N;
					int_type = ElementIntegrationRuleType::Hexahedron_GaussLegendre_1;
				}
				/*
				8-node hexahedron with 2x2x2 gp
				*/
				else if (
					// ./brick
					elem_class_tag == ELE_TAG_BbarBrick ||
					elem_class_tag == ELE_TAG_BbarBrickWithSensitivity ||
					elem_class_tag == ELE_TAG_Brick ||
					// ./up
					elem_class_tag == ELE_TAG_BBarBrickUP ||
					elem_class_tag == ELE_TAG_BrickUP ||
					// ./XMUelements
					elem_class_tag == ELE_TAG_AC3D8HexWithSensitivity
					)
				{
					geom_type = ElementGeometryType::Hexahedron_8N;
					int_type = ElementIntegrationRuleType::Hexahedron_GaussLegendre_2;
				}
				/*
				20-node hexahedron with 3x3x3 gp
				*/
				else if (
					// ./brick
					elem_class_tag == ELE_TAG_Twenty_Node_Brick ||
					// ./up
					elem_class_tag == ELE_TAG_Twenty_Eight_Node_BrickUP
					)
				{
					geom_type = ElementGeometryType::Hexahedron_20N;
					int_type = ElementIntegrationRuleType::Hexahedron_GaussLegendre_3;
				}
			}

			void getCustomGaussPointLocations(Element *elem, ElementIntegrationRule &rule) {
				/*
				clear any existing locations
				*/
				rule.x.clear();
				/*
				ask for integrationPoints ...
				*/
				{
					//std::cout << "get custom gp: trying with \"integrationPoints\"...\n";
					bool done = false;
					std::string request = "integrationPoints";
					int argc = 1;
					const char **argv = new const char*[argc];
					argv[0] = request.c_str();
					OutputDescriptor eo_descriptor;
					OutputDescriptorStream eo_stream(&eo_descriptor);
					Response *eo_response = elem->setResponse(argv, argc, eo_stream);
					eo_stream.finalizeSetResponse();
					if (eo_response) {
						eo_response->getResponse();
						const Vector &data = eo_response->getInformation().getData();
						if (data.Size() > 0) {
							if (data.Size() == 1) {
								rule.x.resize(1);
								rule.x[0] = 0.0;
							}
							else {
								double x_min = std::numeric_limits<double>::max();
								double x_max = -x_min;
								for (int i = 0; i < data.Size(); i++) {
									double ieta = data[i];
									if (ieta < x_min)
										x_min = ieta;
									else if (ieta > x_max)
										x_max = ieta;
								}
								double span = x_max - x_min;
								//elem->getNodePtrs(); // HERE
								if (span == 0.0) {
									rule.x.resize((size_t)data.Size(), 0.0);
								}
								else {
									rule.x.resize((size_t)data.Size());
									for (int i = 0; i < data.Size(); i++)
										rule.x[(size_t)i] = 2.0*(data[i] - x_min) / span - 1.0;
								}
							}
							done = true;
						}
						delete eo_response;
					}
					delete[] argv;
					if (done)
						return;
				}
				/*
				..., otherwise, ask for a dummy response on all sections ...
				*/
				{
					//std::cout << "get custom gp: trying with \"section(all)\"...\n";
					bool done = false;
					std::string request1 = "section";
					std::string request2 = "dummy";
					int argc = 2;
					const char **argv = new const char*[argc];
					argv[0] = request1.c_str();
					argv[1] = request2.c_str();
					OutputDescriptor eo_descriptor;
					OutputDescriptorStream eo_stream(&eo_descriptor);
					Response *eo_response = elem->setResponse(argv, argc, eo_stream);
					eo_stream.finalizeSetResponse();
					if (eo_response)
						delete eo_response; // we dont need it now
					eo_descriptor.getGaussLocations(rule.x);
					if (rule.x.size() > 0)
						done = true;
					delete[] argv;
					if (done)
						return;
				}
				/*
				..., otherwise, ask for a dummy response on all sections, findind out what is the number of gauss points
				*/
				{
					//std::cout << "get custom gp: trying with \"section(1,2,..,N)\"...\n";
					bool done = false;
					std::string request1 = "section";
					std::string request3 = "dummy";
					int argc = 3;
					const char **argv = new const char*[argc];
					argv[0] = request1.c_str();
					argv[2] = request3.c_str();

					int trial_num = 0;
					while (true) {
						trial_num++;
						if (trial_num > MPCO_MAX_TRIAL_NSEC) {
							// we should never get here, or at least we hope, anyway we need a limit!
							//opserr << "MPCORecorder warning: iterative guess of ngp: reached maximum number of iteration, giving up...\n";
							break;
						}
						std::stringstream ss_trial_num; ss_trial_num << trial_num;
						std::string s_trial_num = ss_trial_num.str();
						argv[1] = s_trial_num.c_str();
						OutputDescriptor eo_descriptor;
						OutputDescriptorStream eo_stream(&eo_descriptor);
						Response *eo_response = elem->setResponse(argv, argc, eo_stream);
						eo_stream.finalizeSetResponse();
						if (eo_response)
							delete eo_response; // we dont need it now
						std::vector<double> trial_x;
						eo_descriptor.getGaussLocations(trial_x);
						if (trial_x.size() > 0) {
							if (trial_x.size() > 1) {
								// we should never get here!
								opserr << "MPCORecorder warning: iterative guess of ngp: expected 1 trial gauss location, given = "
									<< (int)trial_x.size()
									<< "\nonly the first one will be considered\n";
							}
							rule.x.push_back(trial_x[0]);
						}
						else {
							// we reached the maximum number of gauss points for this element
							break;
						}
					}
					if (rule.x.size() > 0) {
						if (rule.x.size() == 1) {
							rule.x[0] = 0.0;
						}
						else {
							double x_min = std::numeric_limits<double>::max();
							double x_max = -x_min;
							for (size_t i = 0; i < rule.x.size(); i++) {
								double ieta = rule.x[i];
								if (ieta < x_min)
									x_min = ieta;
								else if (ieta > x_max)
									x_max = ieta;
							}
							double span = x_max - x_min;
							if (span == 0.0) {
								for (size_t i = 0; i < rule.x.size(); i++)
									rule.x[i] = 0.0;
							}
							else {
								for (size_t i = 0; i < rule.x.size(); i++)
									rule.x[i] = 2.0*(rule.x[i] - x_min) / span - 1.0;
							}
						}
						done = true;
					}

					delete[] argv;
					if (done)
						return;
				}
				/*
				..., otherwise, sorry we did our best... try to implement all kinds of response in this element
				to make things work smoothly
				*/
				opserr << "MPCORecorder warning: cannot get custom integration rule from element "
					<< elem->getTag() << "[class = " << elem->getClassType() << "]\n";
			}

			std::map<int, ElementIntegrationRule> registered_custom_rules;
			std::map<int, ElementWithSameClassTagCollection> items;
		};
	}

}

/*************************************************************************************

private_data class.
private storage class for MPCORecorder

**************************************************************************************/

class MPCORecorder::private_data
{
public:
	private_data()
		: filename()
		, initialized(false)
		, first_domain_changed_done(false)
		, info()
		, output_freq()
		, has_region(false)
		, node_set()
		, elem_set()
		, nodes()
		, elements()
		, nodal_results_requests()
		, sens_grad_indices()
		, nodal_recorders()
		, elemental_results_requests()
		, elemental_recorders()
		, elemental_responses()
		, elem_ngauss_nfiber_info()
		, send_self_count(0)
		, p_id(0)
	{}
public:
	// misc
	std::string filename;
	bool initialized;

	// domain changed stuff...
	bool first_domain_changed_done;
	int domain_changed_stamp;

	// info
	mpco::ProcessInfo info;

	// output frequency
	mpco::OutputFrequency output_freq;

	// nodes and elements
	bool has_region;
	std::vector<int> node_set;
	std::vector<int> elem_set;
	std::vector<Node*> nodes;
	mpco::element::ElementCollection elements;

	// nodal recorders
	std::vector<mpco::NodalResultType::Enum> nodal_results_requests;
	std::vector<int> sens_grad_indices;
	mpco::node::ResultRecorderMap nodal_recorders;

	// elemental recorders
	std::vector< std::vector< std::string > > elemental_results_requests;
	mpco::element::ResultRecorderCollection elemental_recorders;
	std::vector<Response*> elemental_responses;

	/*
	an auxiliary map to store, for each element the number of gauss points (with sections) and the number of fibers.
	key = element tag, value = vector with size equal to number of gauss points, and each item is
	a pair with the initial id of the fiber (typically 0 but some sections use 1-based fiber indexing) and
	the number of fibers for that gauss point
	*/
	std::map<int, std::vector<std::pair<int, int> > > elem_ngauss_nfiber_info;

	// parallel stuff
	int send_self_count;
	int p_id;
};

/*************************************************************************************

MPCORecorder class implementation

**************************************************************************************/

MPCORecorder::MPCORecorder()
	: Recorder(RECORDER_TAGS_MPCORecorder)
	, m_data(new private_data())
{
}

MPCORecorder::~MPCORecorder()
{
	if (m_data->initialized) {
		/*
		error flags
		*/
		herr_t status;
		/*
		close file
		*/
		status = h5::file::close(m_data->info.h_file_id);
		if (status < 0) {
			opserr << "MPCORecorder Error: cannot close file on destructor\n";
		}
		/*
		close property lists
		*/
		status = h5::plist::close(m_data->info.h_group_proplist);
		status = h5::plist::close(m_data->info.h_file_proplist);
#ifdef MPCO_USE_SWMR
		status = h5::plist::close(m_data->info.h_file_acc_proplist);
#endif // MPCO_USE_SWMR
		/*
		delete nodal recorders
		*/
		clearNodeRecorders();
		/*
		delete elemental recorders-responses
		*/
		clearElementRecorders();
	}
	delete m_data;
}

int MPCORecorder::record(int commitTag, double timeStamp)
{
	/*
	error flags
	*/
	int retval = 0;
	herr_t status = 0;
	/*
	save current time step info
	*/
	m_data->info.current_time_step = timeStamp;
	m_data->info.current_time_step_id = commitTag;
	/*
	check whether to record or not this timestep
	*/
	bool do_record = false;
	if (!m_data->initialized) {
		// always record the first time
		do_record = true;
	}
	else {
		if (m_data->output_freq.type == mpco::OutputFrequency::NumberOfSteps) {
			if ((m_data->info.current_time_step_id - m_data->output_freq.last_step) >= m_data->output_freq.nsteps) {
				do_record = true;
			}
		}
		else if (m_data->output_freq.type == mpco::OutputFrequency::DeltaTime) {
			if (std::abs(m_data->info.current_time_step - m_data->output_freq.last_time) >= m_data->output_freq.dt) {
				do_record = true;
			}
		}
	}
	if (do_record) {
		m_data->output_freq.last_step = m_data->info.current_time_step_id;
		m_data->output_freq.last_time = m_data->info.current_time_step;
	}
	else {
		return retval;
	}
	/*
	perform initilization on first call
	*/
	if (!m_data->initialized) {
		retval = initialize();
		if (retval) {
			opserr << "MPCRecorder Error: cannot initialize the recorder\n";
			return retval;
		}
	}
	/*
	check domain changed and perform related initializations
	*/
	bool rebuild_model = false;
	if (!m_data->first_domain_changed_done) {
		m_data->info.current_model_stage_id = m_data->info.domain->hasDomainChanged();
		m_data->first_domain_changed_done = true;
		rebuild_model = true;
	}
	else {
		int new_domain_stamp = m_data->info.domain->hasDomainChanged();
		if (new_domain_stamp != m_data->info.current_model_stage_id) {
			m_data->info.current_model_stage_id = new_domain_stamp;
			rebuild_model = true;
		}
	}
	if (rebuild_model) {
		retval = writeModel();
		if (retval) {
			opserr << "MPCRecorder Error: cannot write model\n";
			return retval;
		}
	}
	/*
	write node results
	*/
	retval = recordResultsOnNodes();
	if (retval) {
		opserr << "MPCRecorder Error: cannot record results on nodes\n";
		return retval;
	}
	/*
	write element results
	*/
	retval = recordResultsOnElements();
	if (retval) {
		opserr << "MPCRecorder Error: cannot record results on elements\n";
		return retval;
	}
	/*
	flush file
	*/ 
	status = h5::file::flush(m_data->info.h_file_id);
	if (status < 0) {
		opserr << "MPCORecorder Error: cannot flush file on record()\n";
		retval = -1;
		return retval;
	}
	return retval;
}

int MPCORecorder::restart(void)
{
	return 0;
}

int MPCORecorder::domainChanged(void)
{
	return 0;
}

int MPCORecorder::setDomain(Domain &theDomain)
{
	if (m_data->info.domain) {
		if (m_data->info.domain != &theDomain) {
			opserr << "MPCORecorder Error: setDomain with an existing domain is not implemented yet...\n";
			exit(-1);
		}
	}
	m_data->info.domain = &theDomain;
	return 0;
}

int MPCORecorder::sendSelf(int commitTag, Channel &theChannel)
{
	if (theChannel.isDatastore() == 1) {
		opserr << "MPCORecorder::sendSelf() - does not send data to a datastore\n";
		return -1;
	}

	m_data->send_self_count++;

	std::stringstream _info;
	_info << "MPCORecorder sendSelf from: " << m_data->p_id << ", send self count = " << m_data->send_self_count << "\n";
	std::cout << _info.str();

	// element results requests
	std::string elem_res_merged_string;
	{
		std::stringstream ss;
		for (size_t i = 0; i < m_data->elemental_results_requests.size(); i++) {
			if (i > 0) {
				ss << ':'; // char separator for results
			}
			const std::vector<std::string> &i_request = m_data->elemental_results_requests[i];
			for (size_t j = 0; j < i_request.size(); j++) {
				if (j > 0) {
					ss << '.'; // char separator for jth token
				}
			}
		}
		elem_res_merged_string = ss.str();
	}
	
	// send misc info
	{
		ID aux(8);
		aux(0) = getTag();
		aux(1) = m_data->send_self_count; // use the send self counter as p_id in the receiver
		aux(2) = static_cast<int>(m_data->filename.size());
		aux(3) = static_cast<int>(m_data->nodal_results_requests.size());
		aux(4) = static_cast<int>(elem_res_merged_string.size());
		aux(5) = static_cast<int>(m_data->has_region);
		aux(6) = static_cast<int>(m_data->node_set.size());
		aux(7) = static_cast<int>(m_data->elem_set.size());
		if (theChannel.sendID(0, commitTag, aux) < 0) {
			opserr << "MPCORecorder::sendSelf() - failed to send misc info\n";
			return -1;
		}
	}

	// send filename
	if (m_data->filename.size() > 0) {
		std::vector<char> aux(m_data->filename.begin(), m_data->filename.end());
		Message msg(aux.data(), static_cast<int>(m_data->filename.size()));
		if (theChannel.sendMsg(0, commitTag, msg) < 0) {
			opserr << "MPCORecorder::sendSelf() - failed to send filename\n";
			return -1;
		}
	}

	// send output frequency
	{
		Vector aux(3);
		aux(0) = static_cast<double>(m_data->output_freq.type);
		aux(1) = static_cast<double>(m_data->output_freq.dt);
		aux(2) = static_cast<double>(m_data->output_freq.nsteps);
		if (theChannel.sendVector(0, commitTag, aux) < 0) {
			opserr << "MPCORecorder::sendSelf() - failed to send output frequency\n";
			return -1;
		}
	}

	// send node result requests
	if (m_data->nodal_results_requests.size() > 0) {
		ID aux(static_cast<int>(m_data->nodal_results_requests.size()));
		for (size_t i = 0; i < m_data->nodal_results_requests.size(); i++)
			aux(static_cast<int>(i)) = static_cast<int>(m_data->nodal_results_requests[i]);
		if (theChannel.sendID(0, commitTag, aux) < 0) {
			opserr << "MPCORecorder::sendSelf() - failed to send node result requests\n";
			return -1;
		}
	}

	// send node result requests (sens grad indices)
	if (m_data->sens_grad_indices.size() > 0) {
		ID aux(static_cast<int>(m_data->sens_grad_indices.size()));
		for (size_t i = 0; i < m_data->sens_grad_indices.size(); i++)
			aux(static_cast<int>(i)) = static_cast<int>(m_data->sens_grad_indices[i]);
		if (theChannel.sendID(0, commitTag, aux) < 0) {
			opserr << "MPCORecorder::sendSelf() - failed to send node result requests (sensitivity parameter indices)\n";
			return -1;
		}
	}

	// send element result requests
	if (elem_res_merged_string.size() > 0) {
		std::vector<char> aux(elem_res_merged_string.begin(), elem_res_merged_string.end());
		Message msg(aux.data(), static_cast<int>(elem_res_merged_string.size()));
		if (theChannel.sendMsg(0, commitTag, msg) < 0) {
			opserr << "MPCORecorder::sendSelf() - failed to send element result requests\n";
			return -1;
		}
	}

	if (m_data->has_region) {
		// send node set
		if (m_data->node_set.size()) {
			ID aux(static_cast<int>(m_data->node_set.size()));
			for (size_t i = 0; i < m_data->node_set.size(); i++)
				aux(static_cast<int>(i)) = m_data->node_set[i];
			if (theChannel.sendID(0, commitTag, aux) < 0) {
				opserr << "MPCORecorder::sendSelf() - failed to send node set\n";
				return -1;
			}
		}

		// send elem set
		if (m_data->elem_set.size()) {
			ID aux(static_cast<int>(m_data->elem_set.size()));
			for (size_t i = 0; i < m_data->elem_set.size(); i++)
				aux(static_cast<int>(i)) = m_data->elem_set[i];
			if (theChannel.sendID(0, commitTag, aux) < 0) {
				opserr << "MPCORecorder::sendSelf() - failed to send element set\n";
				return -1;
			}
		}
	}

	return 0;
}

int MPCORecorder::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	if (theChannel.isDatastore() == 1) {
		opserr << "MPCORecorder::sendSelf() - does not send data to a datastore\n";
		return -1;
	}

	if (m_data) {
		delete m_data;
		m_data = 0;
	}

	m_data = new private_data();
	m_data->send_self_count = -1;

	// recv misc info
	size_t aux_filename_size(0);
	size_t aux_node_res_size(0);
	size_t aux_res_merged_string_size(0);
	size_t aux_node_set_size(0);
	size_t aux_elem_set_size(0);
	{
		ID aux(8);
		if (theChannel.recvID(0, commitTag, aux) < 0) {
			opserr << "MPCORecorder::sendSelf() - failed to recv misc info\n";
			return -1;
		}
		setTag(aux(0));
		m_data->p_id = aux(1); // use the send self counter as p_id in the receiver
		aux_filename_size = aux(2);
		aux_node_res_size = aux(3);
		aux_res_merged_string_size = aux(4);
		m_data->has_region = aux(5) != 0;
		aux_node_set_size = static_cast<size_t>(aux(6));
		aux_elem_set_size = static_cast<size_t>(aux(7));
	}

	// recv filename
	if (aux_filename_size > 0) {
		std::vector<char> aux(aux_filename_size);
		Message msg(aux.data(), static_cast<int>(aux_filename_size));
		if (theChannel.recvMsg(0, commitTag, msg) < 0) {
			opserr << "MPCORecorder::sendSelf() - failed to recv filename\n";
			return -1;
		}
		m_data->filename = std::string(aux.begin(), aux.end());
	}

	// recv output frequency
	{
		Vector aux(3);
		if (theChannel.recvVector(0, commitTag, aux) < 0) {
			opserr << "MPCORecorder::sendSelf() - failed to recv output frequency\n";
			return -1;
		}
		m_data->output_freq.type = static_cast<mpco::OutputFrequency::IncrementType>(static_cast<int>(aux(0)));
		m_data->output_freq.dt = aux(1);
		m_data->output_freq.nsteps = static_cast<int>(aux(2));
	}

	// recv node result requests
	if (aux_node_res_size > 0) {
		ID aux(static_cast<int>(aux_node_res_size));
		if (theChannel.recvID(0, commitTag, aux) < 0) {
			opserr << "MPCORecorder::sendSelf() - failed to recv node result requests\n";
			return -1;
		}
		m_data->nodal_results_requests.resize(aux_node_res_size);
		for (size_t i = 0; i < aux_node_res_size; i++)
			m_data->nodal_results_requests[i] = static_cast<mpco::NodalResultType::Enum>(aux(static_cast<int>(i)));
	}

	// recv node result requests (sens grad indices)
	if (m_data->sens_grad_indices.size() > 0) {
		ID aux(static_cast<int>(aux_node_res_size));
		if (theChannel.recvID(0, commitTag, aux) < 0) {
			opserr << "MPCORecorder::sendSelf() - failed to recv node result requests (sensitivity parameter indices)\n";
			return -1;
		}
		m_data->sens_grad_indices.resize(aux_node_res_size);
		for (size_t i = 0; i < aux_node_res_size; i++)
			m_data->sens_grad_indices[i] = aux(static_cast<int>(i));
	}

	// recv element result requests
	if (aux_res_merged_string_size > 0) {
		std::vector<char> aux(aux_res_merged_string_size);
		Message msg(aux.data(), static_cast<int>(aux_res_merged_string_size));
		if (theChannel.sendMsg(0, commitTag, msg) < 0) {
			opserr << "MPCORecorder::sendSelf() - failed to recv element result requests\n";
			return -1;
		}
		std::vector<std::string> aux_1;
		utils::strings::split(std::string(aux.begin(), aux.end()), ':', aux_1, true);
		for (size_t i = 0; i < aux_1.size(); i++) {
			std::vector<std::string> aux_2;
			utils::strings::split(aux_1[i], '.', aux_2, true);
			m_data->elemental_results_requests.push_back(aux_2);
		}
	}

	if (m_data->has_region) {
		// recv node set
		if (aux_node_set_size > 0) {
			ID aux(static_cast<int>(aux_node_set_size));
			if (theChannel.recvID(0, commitTag, aux) < 0) {
				opserr << "MPCORecorder::sendSelf() - failed to recv node set\n";
				return -1;
			}
			m_data->node_set.resize(aux_node_set_size);
			for (size_t i = 0; i < aux_node_set_size; i++)
				m_data->node_set[i] = aux(static_cast<int>(i));
		}

		// recv elem set
		if (aux_elem_set_size > 0) {
			ID aux(static_cast<int>(aux_elem_set_size));
			if (theChannel.recvID(0, commitTag, aux) < 0) {
				opserr << "MPCORecorder::sendSelf() - failed to recv elem set\n";
				return -1;
			}
			m_data->elem_set.resize(aux_elem_set_size);
			for (size_t i = 0; i < aux_elem_set_size; i++)
				m_data->elem_set[i] = aux(static_cast<int>(i));
		}
	}

	std::stringstream _info;
	_info << "MPCORecorder recvSelf from: " << m_data->p_id << ", send self count = " << m_data->send_self_count << "\n";
	_info << "filename = " << m_data->filename << "\n";
	_info << "freq: ";
	if (m_data->output_freq.type == mpco::OutputFrequency::DeltaTime) 
		_info << "dt = " << m_data->output_freq.dt << "\n";
	else 
		_info << "nsteps = " << m_data->output_freq.nsteps << "\n";
	_info << "nodal results [" << m_data->nodal_results_requests.size() << "]" << "\n";
	for (size_t i = 0; i < m_data->nodal_results_requests.size(); i++) {
		_info << "   " << m_data->nodal_results_requests[i] << "\n";
	}
	_info << "elemental results [" << m_data->elemental_results_requests.size() << "]" << "\n";
	for (size_t i = 0; i < m_data->elemental_results_requests.size(); i++) {
		const std::vector<std::string> &ireq = m_data->elemental_results_requests[i];
		for (size_t j = 0; j < ireq.size(); j++) {
			if (j > 0)
				_info << ".";
			_info << ireq[j];
		}
		_info << "\n";
	}
	_info << "node set:" << "\n" << "[";
	{
		int n_counter(0);
		for (size_t i = 0; i < m_data->node_set.size(); i++) {
			_info << m_data->node_set[i] << " ";
			n_counter++;
			if (n_counter >= 5) {
				_info << "\n";
				n_counter = 0;
			}
		}
	}
	_info << "elem set:" << "\n" << "[";
	{
		int n_counter(0);
		for (size_t i = 0; i < m_data->elem_set.size(); i++) {
			_info << m_data->elem_set[i] << " ";
			n_counter++;
			if (n_counter >= 5) {
				_info << "\n";
				n_counter = 0;
			}
		}
	}
	_info << "]" << "\n";
	std::cout << _info.str();

	return 0;
}

int MPCORecorder::initialize()
{
	/*
	error flags
	*/
	int retval = 0;
	herr_t status = 0;
	/*
	initialize only once
	*/
	if (m_data->initialized)
		return retval;
	/*
	get the dimension of the problem
	*/
	if (m_data->info.domain == 0) {
		opserr << "MPCORecorder error: null domain\n";
		exit(-1);
	}
	if (m_data->info.num_dimensions < 1) {
		NodeIter& node_iter = m_data->info.domain->getNodes();
		Node* current_node = 0;
		while ((current_node = node_iter()) != 0) {
			int trial_dim = current_node->getCrds().Size();
			if (m_data->info.num_dimensions < 1) {
				m_data->info.num_dimensions = trial_dim;
			}
			else {
				if (m_data->info.num_dimensions != trial_dim) {
					opserr << "MPCORecorder error: nodes with different dimensions, this is not supported\n";
					exit(-1);
				}
			}
		}
		if (m_data->info.num_dimensions < 1 || m_data->info.num_dimensions > 3) {
			opserr << "MPCORecorder error: invalid spatial dimension. expected 1 2 or 3, given: " << m_data->info.num_dimensions << "\n";
			exit(-1);
		}
	}
	/*
	create property lists (enable link creation order tracking)
	*/
	m_data->info.h_file_proplist = h5::plist::crate(h5::plist::FileCreate);
	status = h5::plist::setLinkCreationOrder(m_data->info.h_file_proplist, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED);
	m_data->info.h_group_proplist = h5::plist::crate(h5::plist::GroupCreate);
	status = h5::plist::setLinkCreationOrder(m_data->info.h_group_proplist, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED);
	/*
	open or create the output file
	*/
	std::string the_filename;
	{
		std::vector<std::string> filename_words;
		utils::strings::split(m_data->filename, '.', filename_words, true);
		if (filename_words.size() == 0) 
			filename_words.push_back("DefaultName");
		if (filename_words.back() != "mpco")
			filename_words.push_back("mpco");
		std::stringstream ss_filename;
		for (size_t i = 0; i < filename_words.size() - 1; i++)
			ss_filename << filename_words[i] << '.';
		if (m_data->send_self_count != 0) { // > 0 -> we are in p0, < 0 -> we are in secondary procs, = 0 -> not in parallel
			ss_filename << "p" << m_data->p_id << '.';
		}
		ss_filename << filename_words.back();
		the_filename = ss_filename.str();
	}
#ifdef MPCO_USE_SWMR
	m_data->info.h_file_acc_proplist = h5::plist::crate(h5::plist::FileAccess);
	status = h5::plist::setLibVerBounds(m_data->info.h_file_acc_proplist, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
	m_data->info.h_file_id = h5::file::create(the_filename.c_str(), m_data->info.h_file_proplist, m_data->info.h_file_acc_proplist);
	status = h5::file::startSWMR(m_data->info.h_file_id);
#else
	m_data->info.h_file_id = h5::file::create(the_filename.c_str(), m_data->info.h_file_proplist, H5P_DEFAULT);
#endif // MPCO_USE_SWMR
	if (m_data->info.h_file_id == HID_INVALID) {
		opserr << "MPCORecorder Error: cannot create or open file: \"" << the_filename.c_str() << "\"";
		exit(-1);
	}
	/*
	create info group and metadata
	*/
	hid_t h_gp_info = h5::group::create(m_data->info.h_file_id, "INFO", H5P_DEFAULT, m_data->info.h_group_proplist, H5P_DEFAULT);
	hid_t h_dset_dim = h5::dataset::createAndWrite(h_gp_info, "SPATIAL_DIM", m_data->info.num_dimensions);
	status = h5::dataset::close(h_dset_dim);
	status = h5::group::close(h_gp_info);	
	/*
	mark as initialized and return
	*/
	m_data->initialized = true;
	return retval;
}

int MPCORecorder::writeModel()
{
	/*
	error flags
	*/
	int retval = 0;
	herr_t status;
	/*
	create model and results groups
	*/
	std::stringstream ss_current_model_stage_name;
	ss_current_model_stage_name << "MODEL_STAGE[" << m_data->info.current_model_stage_id << "]";
	std::string model_stage_dir = ss_current_model_stage_name.str();
	hid_t h_gp_model_stage = h5::group::create(m_data->info.h_file_id, model_stage_dir.c_str(), H5P_DEFAULT, m_data->info.h_group_proplist, H5P_DEFAULT);
	hid_t h_gp_model = h5::group::create(h_gp_model_stage, "MODEL", H5P_DEFAULT, m_data->info.h_group_proplist, H5P_DEFAULT);
	hid_t h_gp_results = h5::group::create(h_gp_model_stage, "RESULTS", H5P_DEFAULT, m_data->info.h_group_proplist, H5P_DEFAULT);
	hid_t h_gp_results_on_nodes = h5::group::create(h_gp_results, "ON_NODES", H5P_DEFAULT, m_data->info.h_group_proplist, H5P_DEFAULT);
	hid_t h_gp_results_on_elements = h5::group::create(h_gp_results, "ON_ELEMENTS", H5P_DEFAULT, m_data->info.h_group_proplist, H5P_DEFAULT);
	/*
	attributes
	*/
	status = h5::attribute::write(h_gp_model_stage, "STEP", m_data->info.current_time_step_id);
	status = h5::attribute::write(h_gp_model_stage, "TIME", m_data->info.current_time_step);
	/*
	close groups
	*/
	status = h5::group::close(h_gp_results_on_elements);
	status = h5::group::close(h_gp_results_on_nodes);
	status = h5::group::close(h_gp_results);
	status = h5::group::close(h_gp_model);
	status = h5::group::close(h_gp_model_stage);
	/*
	write info.nodes
	*/
	retval = writeModelNodes();
	if (retval)
		return retval;
	/*
	write elements
	*/
	retval = writeModelElements();
	if (retval)
		return retval;
	/*
	write local axes
	*/
	retval = writeModelLocalAxes();
	if (retval)
		return retval;
	/*
	write cross sections
	*/
	retval = writeSections();
	if (retval)
		return retval;
	/*
	write sets
	*/
	retval = writeSets();
	if (retval)
		return retval;
	/*
	initialize node recorders
	*/
	retval = initNodeRecorders();
	if (retval)
		return retval;
	/*
	initialize element recorders
	*/
	retval = initElementRecorders();
	if (retval)
		return retval;
	/*
	return
	*/
	return retval;
}

int MPCORecorder::writeModelNodes()
{
#ifdef MPCO_TIMING
	mpco::Timer timer("writeModelNodes"); timer.start();
#endif // MPCO_TIMING
	/*
	error flags
	*/
	int retval = 0;
	herr_t status = 0;
	/*
	get all nodes except pressure nodes
	*/
	m_data->nodes.clear();
	ID pressure_node_tags(0, m_data->info.domain->getNumPCs());
	Pressure_ConstraintIter& pressure_constraint_iter = m_data->info.domain->getPCs();
	Pressure_Constraint* pressure_constraint = 0;
	while ((pressure_constraint = pressure_constraint_iter()) != 0) {
		Node* pressure_node = pressure_constraint->getPressureNode();
		if (pressure_node) {
			pressure_node_tags.insert(pressure_node->getTag());
		}
	}
	if (m_data->has_region) { // only nodes in the requested regions
		for (std::vector<int>::const_iterator id_iter = m_data->node_set.begin(); id_iter != m_data->node_set.end(); ++id_iter) {
			int node_id = *id_iter;
			if (pressure_node_tags.getLocationOrdered(node_id) < 0) {
				Node *current_node = m_data->info.domain->getNode(node_id);
				if (current_node != 0) {
					m_data->nodes.push_back(current_node);
				}
			}
		}
	}
	else { // all nodes
		NodeIter& node_iter = m_data->info.domain->getNodes();
		Node* current_node = 0;
		while ((current_node = node_iter()) != 0) {
			int nd = current_node->getTag();
			if (pressure_node_tags.getLocationOrdered(nd) < 0) {
				m_data->nodes.push_back(current_node);
			}
		}
	}

	/*
	quick return
	*/
	size_t num_nodes = m_data->nodes.size();
	if (num_nodes == 0) {
		opserr << "MPCORecorder Error: no nodes to write\n";
		return -1;
	}
	/*
	allocate data for id and coordinates datasets
	*/
	std::vector<int> buffer_id(num_nodes);
	std::vector<double> buffer_coord(num_nodes*m_data->info.num_dimensions);
	for (size_t i = 0; i < num_nodes; i++) {
		Node *inode = m_data->nodes[i];
		buffer_id[i] = inode->getTag();
		size_t j = i * m_data->info.num_dimensions;
		const Vector &coord = inode->getCrds();
		buffer_coord[j] = coord[0];
		if (m_data->info.num_dimensions > 1) {
			buffer_coord[j + 1] = coord[1];
			if (m_data->info.num_dimensions > 2) {
				buffer_coord[j + 2] = coord[2];
			}
		}
	}
	/*
	create the nodes group
	*/
	std::stringstream ss_gp_nodes_dir;
	ss_gp_nodes_dir << "MODEL_STAGE[" << m_data->info.current_model_stage_id << "]/MODEL/NODES";
	std::string gp_nodes_dir = ss_gp_nodes_dir.str();
	hid_t h_gp_nodes = h5::group::create(m_data->info.h_file_id, gp_nodes_dir.c_str(), H5P_DEFAULT, m_data->info.h_group_proplist, H5P_DEFAULT);
	/*
	create the ID dataset
	*/
	hid_t dset_id = h5::dataset::createAndWrite(h_gp_nodes, "ID", buffer_id);
	status = h5::dataset::close(dset_id);
	/*
	create the COORDINATES dataspace
	*/
	hid_t dset_coord = h5::dataset::createAndWrite(h_gp_nodes, "COORDINATES", buffer_coord, num_nodes, (hsize_t)m_data->info.num_dimensions);
	status = h5::dataset::close(dset_coord);
	/*
	close nodes group
	*/
	status = h5::group::close(h_gp_nodes);
	/*
	return
	*/
#ifdef MPCO_TIMING
	timer.stop();
#endif // MPCO_TIMING
	return retval;
}

int MPCORecorder::writeModelElements()
{
#ifdef MPCO_TIMING
	mpco::Timer timer("writeModelElements"); timer.start();
#endif // MPCO_TIMING
	/*
	error flags
	*/
	int retval = 0;
	herr_t status = 0;
	/*
	create element group
	*/
	std::stringstream ss_gp_elements_dir;
	ss_gp_elements_dir << "MODEL_STAGE[" << m_data->info.current_model_stage_id << "]/MODEL/ELEMENTS";
	std::string gp_elements_dir = ss_gp_elements_dir.str();
	hid_t h_gp_elements = h5::group::create(m_data->info.h_file_id, gp_elements_dir.c_str(), H5P_DEFAULT, m_data->info.h_group_proplist, H5P_DEFAULT);
	/*
	map all elements
	*/
	m_data->elements.mapElements(m_data->info.domain, m_data->has_region, m_data->elem_set);
	/*
	loop over all element group-by-tag
	*/
	for (mpco::element::ElementCollection::submap_type::iterator it1 = m_data->elements.items.begin();
		it1 != m_data->elements.items.end(); ++it1) {
		mpco::element::ElementWithSameClassTagCollection &elem_by_tag = it1->second;
		/*
		loop over all element group-by-integration rule
		*/
		for (mpco::element::ElementWithSameClassTagCollection::submap_type::iterator it2 = elem_by_tag.items.begin();
			it2 != elem_by_tag.items.end(); ++it2) {
			mpco::element::ElementWithSameIntRuleCollection &elem_by_rule = it2->second;
			/*
			loop over all element group-by-custom(or default) integration rule
			*/
			for (mpco::element::ElementWithSameIntRuleCollection::submap_type::iterator it3 = elem_by_rule.items.begin();
				it3 != elem_by_rule.items.end(); ++it3) {
				mpco::element::ElementWithSameCustomIntRuleCollection &elem_by_custom_rule = it3->second;
				/*
				create a dataset with all elements in this group
				*/
				if (elem_by_custom_rule.items.size() > 0) {
					/*
					create a name for this dataset using the following format <class_tag>-<class_name>[<integration_rule>:<custom_rule_index>]
					*/
					std::stringstream ss_dset_name;
					ss_dset_name << elem_by_tag.class_tag << "-" << elem_by_tag.class_name
						<< "[" << elem_by_rule.int_rule_type << ":"
						<< elem_by_custom_rule.custom_int_rule_index << "]";
					elem_by_custom_rule.name = ss_dset_name.str();
					/*
					create the dataset and write data to it.
					the first column is the element id, the remaining columns are nodal connectivity
					*/
					std::vector<int> buffer(elem_by_custom_rule.items.size()*(1 + elem_by_tag.num_nodes));
					size_t offset = 0;
					for (mpco::element::ElementWithSameCustomIntRuleCollection::collection_type::iterator it4 = elem_by_custom_rule.items.begin();
						it4 != elem_by_custom_rule.items.end(); ++it4) {
						Element *elem = *it4;
						buffer[offset] = elem->getTag();
						for (size_t j = 0; j < elem_by_tag.num_nodes; j++)
							buffer[j + 1 + offset] = elem->getExternalNodes()((int)j);
						offset += (1 + elem_by_tag.num_nodes);
					}
					hid_t dset_id = h5::dataset::createAndWrite(h_gp_elements, elem_by_custom_rule.name.c_str(), buffer, elem_by_custom_rule.items.size(), (hsize_t)(1 + elem_by_tag.num_nodes));
					/*
					add attributes
					*/
					status = h5::attribute::write(dset_id, "GEOMETRY", (int)elem_by_tag.geom_type);
					status = h5::attribute::write(dset_id, "INTEGRATION_RULE", (int)elem_by_rule.int_rule_type);
					if (elem_by_rule.int_rule_type == mpco::ElementIntegrationRuleType::CustomIntegrationRule) {
						status = h5::attribute::write(dset_id, "CUSTOM_INTEGRATION_RULE", elem_by_custom_rule.custom_int_rule_index);
						if (elem_by_custom_rule.custom_int_rule_index != 0) {
							mpco::element::ElementIntegrationRule &custom_rule = m_data->elements.registered_custom_rules[elem_by_custom_rule.custom_int_rule_index];
							h5::attribute::write(dset_id, "GP_X", custom_rule.x);
						}
					}
					/*
					close dataset
					*/
					status = h5::dataset::close(dset_id);
				}
			}
		}
	}
	/*
	close element group
	*/
	status = h5::group::close(h_gp_elements);
	/*
	return
	*/
#ifdef MPCO_TIMING
	timer.stop();
#endif // MPCO_TIMING
	return retval;
}

int MPCORecorder::writeModelLocalAxes()
{
#ifdef MPCO_TIMING
	mpco::Timer timer("writeModelLocalAxes"); timer.start();
#endif // MPCO_TIMING
#ifdef MPCO_WRITE_LOC_AX_IS_VERBOSE
	std::stringstream ss;
	ss << "local axes:\n";
#endif // MPCO_WRITE_LOC_AX_IS_VERBOSE
	/*
	error flags
	*/
	int retval = 0;
	herr_t status = 0;
	/*
	allocate vectors for ID and local axes, only for elements that provide a response for local axes
	*/
	std::vector<int> buffer_el_id;
	std::vector<utils::locax::quaternion> buffer_el_data;
	/*
	element request
	*/
	int argc = 1;
	const char **argv = new const char*[argc];
	argv[0] = "localAxes";
	Vector vx(3);
	Vector vy(3);
	Vector vz(3);
	/*
	loop over all element group-by-tag
	*/
	for (mpco::element::ElementCollection::submap_type::iterator it1 = m_data->elements.items.begin();
		it1 != m_data->elements.items.end(); ++it1) {
		mpco::element::ElementWithSameClassTagCollection &elem_by_tag = it1->second;
		/*
		loop over all element group-by-integration rule
		*/
		for (mpco::element::ElementWithSameClassTagCollection::submap_type::iterator it2 = elem_by_tag.items.begin();
			it2 != elem_by_tag.items.end(); ++it2) {
			mpco::element::ElementWithSameIntRuleCollection &elem_by_rule = it2->second;
			/*
			loop over all element group-by-custom(or default) integration rule
			*/
			for (mpco::element::ElementWithSameIntRuleCollection::submap_type::iterator it3 = elem_by_rule.items.begin();
				it3 != elem_by_rule.items.end(); ++it3) {
				mpco::element::ElementWithSameCustomIntRuleCollection &elem_by_custom_rule = it3->second;
				/*
				loop over all alements
				*/
				for (std::vector<Element*>::iterator it4 = elem_by_custom_rule.items.begin();
					it4 != elem_by_custom_rule.items.end(); ++it4) {
					Element *elem = *it4;
					/*
					get the element response for the ith section and the ith fiber
					*/
					mpco::element::OutputDescriptor eo_descriptor;
					mpco::element::OutputDescriptorStream eo_stream(&eo_descriptor);
					Response *eo_response = elem->setResponse(argv, argc, eo_stream);
					if (eo_response) {
						eo_response->getResponse();
						const Vector &local_axes_packed_9 = eo_response->getInformation().getData();
						if (local_axes_packed_9.Size() != 9) {
							opserr << "MPCORecorder Error: invalid response data size for local axes\n";
							continue;
						}
						for (int i = 0; i < 3; i++) {
							vx(i) = local_axes_packed_9(i);
							vy(i) = local_axes_packed_9(i + 3);
							vz(i) = local_axes_packed_9(i + 6);
						}
						buffer_el_id.push_back(elem->getTag());
						buffer_el_data.push_back(utils::locax::quatFromMat(vx, vy, vz));
#ifdef MPCO_WRITE_LOC_AX_IS_VERBOSE
						const utils::locax::quaternion &qq = buffer_el_data.back();
						ss << "E: " << elem->getTag() << " "
							<< "[" << vx(0) << ", " << vx(1) << ", " << vx(2) << "]"
							<< "[" << vy(0) << ", " << vy(1) << ", " << vy(2) << "]"
							<< "[" << vz(0) << ", " << vz(1) << ", " << vz(2) << "]"
							<< "[" << qq.w << ", " << qq.x << ", " << qq.y << ", " << qq.z << "]"
							<< "\n";
#endif // MPCO_WRITE_LOC_AX_IS_VERBOSE
					}
				}
			}
		}
	}
#ifdef MPCO_WRITE_LOC_AX_IS_VERBOSE
	std::cout << ss.str();
#endif // MPCO_WRITE_LOC_AX_IS_VERBOSE
	/*
	free argv
	*/
	delete[] argv;
	/*
	quick return
	*/
	if (buffer_el_id.size() != buffer_el_data.size()) {
		opserr << "MPCORecorder Error: buffer_el_id.size() != buffer_el_data.size()\n";
		return -1;
	}
	if (buffer_el_id.size() == 0)
		return retval;
	/*
	create local axes group
	*/
	std::stringstream ss_gp_local_axes_dir;
	ss_gp_local_axes_dir << "MODEL_STAGE[" << m_data->info.current_model_stage_id << "]/MODEL/LOCAL_AXES";
	std::string gp_local_axes_dir = ss_gp_local_axes_dir.str();
	hid_t h_gp_local_axes = h5::group::create(m_data->info.h_file_id, gp_local_axes_dir.c_str(), H5P_DEFAULT, m_data->info.h_group_proplist, H5P_DEFAULT);
	/*
	ID dataset
	*/
	hid_t dset_id = h5::dataset::createAndWrite(h_gp_local_axes, "ID", buffer_el_id);
	status = h5::dataset::close(dset_id);
	/*
	local axes dataset
	*/
	hid_t dset_axes = h5::dataset::createAndWrite(h_gp_local_axes, "LOCAL_AXES", buffer_el_data);
	status = h5::dataset::close(dset_axes);
	/*
	close element group
	*/
	status = h5::group::close(h_gp_local_axes);
	/*
	return
	*/
#ifdef MPCO_TIMING
	timer.stop();
#endif // MPCO_TIMING
	return retval;
}

int MPCORecorder::writeSections()
{
#ifdef MPCO_TIMING
	mpco::Timer timer("writeSections"); timer.start();
#endif // MPCO_TIMING
	/*
	error flags
	*/
	int retval = 0;
	herr_t status = 0;
	/*
	auxiliary variables
	*/
	std::map<int, mpco::element::SectionAssignment> sec_assignments;
	std::map<mpco::element::FiberSectionData, mpco::element::SectionAssignment> aux_sec_assignments;
	std::map<int, std::string> sec_id_and_names;
	std::map<int, std::string> mat_id_and_names;
	/*
	an auxiliary map to store, for each element the number of gauss points (with sections) and the number of fibers.
	key = element tag, value = vector with size equal to the number of gauss points, and each item is the number of fibers
	for that gauss point
	*/
	m_data->elem_ngauss_nfiber_info.clear();
	/*
	create section assignments group
	*/
	std::stringstream ss_gp_section_assignments_dir;
	ss_gp_section_assignments_dir << "MODEL_STAGE[" << m_data->info.current_model_stage_id << "]/MODEL/SECTION_ASSIGNMENTS";
	std::string gp_section_assignments_dir = ss_gp_section_assignments_dir.str();
	hid_t h_gp_section_assignments = h5::group::create(m_data->info.h_file_id, gp_section_assignments_dir.c_str(), H5P_DEFAULT, m_data->info.h_group_proplist, H5P_DEFAULT);
	/*
	loop over all element group-by-tag
	*/
	for (mpco::element::ElementCollection::submap_type::iterator it1 = m_data->elements.items.begin();
		it1 != m_data->elements.items.end(); ++it1) {
		mpco::element::ElementWithSameClassTagCollection &elem_by_tag = it1->second;
		/*
		loop over all element group-by-integration rule
		*/
		for (mpco::element::ElementWithSameClassTagCollection::submap_type::iterator it2 = elem_by_tag.items.begin();
			it2 != elem_by_tag.items.end(); ++it2) {
			mpco::element::ElementWithSameIntRuleCollection &elem_by_rule = it2->second;
			/*
			loop over all element group-by-custom(or default) integration rule
			*/
			for (mpco::element::ElementWithSameIntRuleCollection::submap_type::iterator it3 = elem_by_rule.items.begin();
				it3 != elem_by_rule.items.end(); ++it3) {
				mpco::element::ElementWithSameCustomIntRuleCollection &elem_by_custom_rule = it3->second;
				/*
				loop over all alements 
				*/
				for (std::vector<Element*>::iterator it4 = elem_by_custom_rule.items.begin(); 
					it4 != elem_by_custom_rule.items.end(); ++it4) {
					Element *elem = *it4;
					/*
					prepare the request
					*/
					std::string request1 = utils::shell::isShellElementTag(elem->getClassTag()) ? "material" : "section";
					std::string request2 = "fiber";
					std::string request3 = "stress";
					int argc = 5;
					const char **argv = new const char*[argc];
					argv[0] = request1.c_str();
					argv[2] = request2.c_str();
					argv[4] = request3.c_str();
					/*
					fiber section data for each gauss point
					*/
					std::vector<int> elem_gauss_id;
					std::vector<int> elem_sec_id;
					std::vector<bool> elem_dummy_sec_flags;
					std::vector<mpco::element::FiberSectionData> elem_sections;
					std::vector<int> elem_fiber_base_index; // stores the index base for fibers, typically 0-based but some sections (LayeredShell for example) assume 1-based fiber indexing.
					/*
					for each section (gauss point)...
					*/
					int trial_num_sec = 0; /* note: in setResponse gp number is 1-based */
					while (true) {
						trial_num_sec++;
#ifdef MPCO_WRITE_SECTION_IS_VERBOSE
						std::cout << "trial sec: " << trial_num_sec << "\n";
#endif // MPCO_WRITE_SECTION_IS_VERBOSE
						if (trial_num_sec > MPCO_MAX_TRIAL_NSEC) {
							// we should never get here, or at least we hope, anyway we need a limit!
							//opserr << "MPCORecorder warning: iterative guess of ngp: reached maximum number of iteration, giving up...\n";
							break;
						}
						std::stringstream ss_trial_num_sec; ss_trial_num_sec << trial_num_sec;
						std::string s_trial_num_sec = ss_trial_num_sec.str();
						argv[1] = s_trial_num_sec.c_str();
						/*
						for each fiber
						*/
						int trial_num_fib = -1; /* note: in setResponse fiber index is 0-based */
						bool break_sec_loop = false;
						bool first_fiber_done = false;
						while (true) {
							trial_num_fib++;
							if (trial_num_fib > MPCO_MAX_TRIAL_NFIB) {
								// we should never get here, or at least we hope, anyway we need a limit!
								opserr << "MPCORecorder warning: iterative guess of num fibers: reached maximum number of iteration, giving up...\n";
								break;
							}
							std::stringstream ss_trial_num_fib; ss_trial_num_fib << trial_num_fib;
							std::string s_trial_num_fib = ss_trial_num_fib.str();
							argv[3] = s_trial_num_fib.c_str();
							/*
							get the element response for the ith section and the ith fiber
							*/
							mpco::element::OutputDescriptor eo_descriptor;
							mpco::element::OutputDescriptorStream eo_stream(&eo_descriptor);
							Response *eo_response = elem->setResponse(argv, argc, eo_stream);
							eo_stream.finalizeSetResponse();
							if (eo_response)
								delete eo_response; // we dont need it now
							/*
							post process the response descriptor
							*/
							std::vector<mpco::element::FiberData> trial_fiberdata;
							std::vector<int> trial_fiber_material_id;
							std::vector<int> trial_sec_id;
							std::vector<int> trial_gp_id;
							std::vector<bool> trial_dummy_flag;
							eo_descriptor.getFiberData(trial_fiberdata, 
								trial_fiber_material_id,
								trial_sec_id, trial_gp_id, trial_dummy_flag);
							/*
							first loop break
							*/
							if (trial_gp_id.size() == 0) {
								/*
								we reached the maximum number of gauss points for this element.
								exit this loop, and also the outer section loop
								*/
#ifdef MPCO_WRITE_SECTION_IS_VERBOSE
								std::cout << "MPCORecorder: exiting fiber loop (no gp). iter = " << trial_num_fib << "\n";
#endif // MPCO_WRITE_SECTION_IS_VERBOSE
								break_sec_loop = true;
								break;
							}
							if (trial_sec_id.size() == 0) {
								/*
								we reached the maximum number of fibers for this section, or this element does not have sections
								exit this loop.
								*/
#ifdef MPCO_WRITE_SECTION_IS_VERBOSE
								std::cout << "MPCORecorder: exiting fiber loop (no sec). iter = " << trial_num_fib << "\n";
#endif // MPCO_WRITE_SECTION_IS_VERBOSE
								break;
							}
							/*
							error checks
							*/
							if ((trial_fiberdata.size() > 0) && (trial_fiberdata.size() != trial_sec_id.size())) {
								// this should never happen!
								opserr << "MPCORecorder FATAL Error: trial_fiberdata.size() != trial_sec_id.size()\n";
								exit(-1);
							}
							if (trial_fiber_material_id.size() != trial_fiberdata.size()) {
								// this should never happen!
								opserr << "MPCORecorder FATAL Error: trial_fiber_material_id.size() != trial_fiberdata.size()\n";
								exit(-1);
							}
							if ((trial_sec_id.size() > 0) && (trial_gp_id.size() != trial_sec_id.size())) {
								// this should never happen!
								opserr << "MPCORecorder FATAL Error: trial_gp_id.size() != trial_sec_id.size()\n";
								exit(-1);
							}
							if (trial_sec_id.size() != trial_dummy_flag.size()) {
								// this should never happen!
								opserr << "MPCORecorder FATAL Error: trial_sec_id.size() != trial_dummy_flag.size()\n";
								exit(-1);
							}
							if ((trial_fiberdata.size() > 1) || (trial_sec_id.size() > 1) || (trial_gp_id.size() > 1) || (trial_dummy_flag.size() > 1)) {
								// this should never happen!
								opserr << "MPCORecorder FATAL Error: iterative guess of num fibers: expected 1 trial output, given = \n"
									<< "# fiber data: " << (int)trial_fiberdata.size() << "\n"
									<< "# sec id: " << (int)trial_sec_id.size() << "\n"
									<< "# gauss id: " << (int)trial_gp_id.size() << "\n"
									<< "# dummy flag: " << (int)trial_dummy_flag.size() << "\n";
								exit(-1);
							}
							/*
							now this can be assumed:
							1) # gauss = # sec_id = # dummy_flag = 1
							2) # fiber can be zero or 1
							*/
							/** workaournd for sections that assume 1-based fiber indexing (now only ShellLayeredFiberSection) 
							*/
							if (trial_fiberdata.size() == 0 && trial_num_fib == 0) {
								continue; // go to next fiber iteration with id = 1
							}
							/*
							set data on first fiber iteration
							*/
							if (!first_fiber_done) {
								elem_gauss_id.push_back(trial_gp_id[0]);
								elem_sec_id.push_back(trial_sec_id[0]);
								elem_dummy_sec_flags.push_back(trial_dummy_flag[0]);
								elem_sections.push_back(mpco::element::FiberSectionData());
								elem_fiber_base_index.push_back(trial_num_fib); // do it here after the previous check on 'trial_fiberdata.size() == 0 && trial_num_fib == 0'
								first_fiber_done = true;
							}
							/*
							second loop break
							*/
							if (trial_fiberdata.size() == 0) {
								/*
								we reached the maximum number of fibers for this section.
								exit this loop.
								*/
#ifdef MPCO_WRITE_SECTION_IS_VERBOSE
								std::cout << "MPCORecorder: exiting fiber loop. iter = " << trial_num_fib << "\n";
#endif // MPCO_WRITE_SECTION_IS_VERBOSE
								break;
							}
							/*
							some extra consistency checks... \todo: remove ? seems useless...
							*/
							int curr_gauss_id = trial_gp_id[0];
							int curr_sec_id = trial_sec_id[0];
							bool curr_dummy_flag = trial_dummy_flag[0];
							if (curr_gauss_id != trial_gp_id[0]) {
								opserr << "MPCORecorder Error: curr gauss id = " << curr_gauss_id << ", trial gauss id = " << trial_gp_id[0] << "\n";
								exit(-1);
							}
							if (curr_sec_id != trial_sec_id[0]) {
								opserr << "MPCORecorder Error: curr sec id = " << curr_sec_id << ", trial sec id = " << trial_sec_id[0] << "\n";
								exit(-1);
							}
							if (curr_dummy_flag != trial_dummy_flag[0]) {
								opserr << "MPCORecorder Error: curr dummy_flag = " << std::boolalpha << curr_dummy_flag << ", trial dummy flag = " << std::boolalpha << trial_dummy_flag[0] << "\n";
								exit(-1);
							}
							/*
							add this fiber and go on...
							*/
							if (trial_fiberdata.size() == 1) {
								mpco::element::FiberSectionData &curr_sec_fiberdata = elem_sections.back();
								curr_sec_fiberdata.fibers.push_back(trial_fiberdata[0]);
								curr_sec_fiberdata.materials.push_back(trial_fiber_material_id[0]);
							}
						} // end: fiber loop
						/*
						eventually break the section loop
						*/
						if (break_sec_loop) {
#ifdef MPCO_WRITE_SECTION_IS_VERBOSE
							std::cout << "MPCORecorder: exiting section loop. iter = " << trial_num_sec << "\n";
#endif // MPCO_WRITE_SECTION_IS_VERBOSE
							break;
						}
					} // end: section loop
					/*
					free argv
					*/
					delete[] argv;
					/*
					fiber section data for each gauss point has been obtained
					*/
#ifdef MPCO_WRITE_SECTION_IS_VERBOSE
					std::stringstream buff;
					buff << "\nElement: " << elem->getTag() << "\n";
					buff << "gauss:\n";
					for (size_t i = 0; i < elem_gauss_id.size(); i++)
						buff << elem_gauss_id[i] << "  ";
					buff << "\nsection tags:\n";
					for (size_t i = 0; i < elem_sec_id.size(); i++)
						buff << elem_sec_id[i] << "  ";
					buff << "\nsection dummy flags:\n";
					for (size_t i = 0; i < elem_dummy_sec_flags.size(); i++)
						buff << std::boolalpha << elem_dummy_sec_flags[i] << "  ";
					buff << "\nsection fiber data:\n";
					for (size_t i = 0; i < elem_sections.size(); i++) {
						mpco::element::FiberSectionData &idata = elem_sections[i];
						buff << "  ith section data:\n";
						buff << "    " << std::setw(10) << "Y" << std::setw(10) << "Z" << std::setw(10) << "A" << "\n";
						for (size_t j = 0; j < idata.fibers.size(); j++) {
							mpco::element::FiberData &jdata = idata.fibers[j];
							buff << "    " << std::setw(10) << jdata.y << std::setw(10) << jdata.z << std::setw(10) << jdata.a << "\n";
						}
					}
					std::cout << buff.str();
#endif // MPCO_WRITE_SECTION_IS_VERBOSE
					/*
					store info in element-gauss-fibers map
					*/
					std::vector<std::pair<int, int> > &nfibers_per_gauss_point = m_data->elem_ngauss_nfiber_info[elem->getTag()];
					nfibers_per_gauss_point.resize(elem_gauss_id.size());
					for (size_t igp = 0; igp < elem_gauss_id.size(); igp++) {
						mpco::element::FiberSectionData &igp_fibersec_data = elem_sections[igp];
						nfibers_per_gauss_point[igp].first = elem_fiber_base_index[igp]; // start witn standard 0-based
						nfibers_per_gauss_point[igp].second = static_cast<int>(igp_fibersec_data.fibers.size());
					}
					/*
					if the ith fiber section has a valid tag (i.e. not coming from a floating fiberoutput...)
					let's add it to a unique collection (set) based on the section tag.
					otherwise (floating fibers...) add this section to a map<section, new_tag>.
					then create the section assignment for this element.
					*/
					size_t ngauss = elem_gauss_id.size();
					if (ngauss != elem_sec_id.size() || ngauss != elem_dummy_sec_flags.size() || ngauss != elem_sections.size()) {
						opserr << "MPCORecorder FATAL Error: inconsistent dimensions in writeFiberSections():\n"
							<< "# gauss = " << (int)ngauss
							<< "# sec id = " << (int)elem_sec_id.size()
							<< "# sec dummy flags = " << (int)elem_dummy_sec_flags.size()
							<< "# sections = " << (int)elem_sections.size() << "\n";
						exit(-1);
					}
					for (size_t igp = 0; igp < ngauss; igp++) {
						/*
						this check should not be necessary. just to make sure we didn't mess things up
						*/
						if (igp != elem_gauss_id[igp]) {
							opserr << "MPCORecorder FATAL Error: igp != elem_gauss_id[igp]\n";
							exit(-1);
						}
						int igp_sec_id = elem_sec_id[igp];
						bool igp_dummy_sec_flag = elem_dummy_sec_flags[igp];
						mpco::element::FiberSectionData &igp_fibersec_data = elem_sections[igp];
						if (igp_dummy_sec_flag) {
							/*
							if this section is a dummy section (i.e. coming from a FiberOutput without SectionOutput parent),
							let's add it to the aux vector with a local indexing
							*/
							mpco::element::SectionAssignment &current_assignment = aux_sec_assignments[igp_fibersec_data];
							if (current_assignment.is_new) {
								current_assignment.fiber_section_data = igp_fibersec_data;
								current_assignment.is_new = false;
							}
							current_assignment.assignments.push_back(mpco::element::ElemGaussPair(elem->getTag(), (int)igp));
						}
						else {
							/*
							if this section has a valid tag,
							get or create the section assignment, and add this elem-gauss pair
							*/
							mpco::element::SectionAssignment &current_assignment = sec_assignments[igp_sec_id];
							if (current_assignment.is_new) {
								current_assignment.fiber_section_data = igp_fibersec_data;
								current_assignment.is_new = false;
							}
							current_assignment.assignments.push_back(mpco::element::ElemGaussPair(elem->getTag(), (int)igp));
						}
					}
				} // end loop for each element in same int rule group
			}
		}
	}
	/*
	here we call the OPS_GetSectionForceDeformation to get the class name. This is not necessary, just to
	give the user more informations about the written sections. warning: in this method (writeSections())
	we assume that all sections have been created by the user via the "section" command, so that they all have different tags.
	However, in some places in opensees, sections are manually created with hardcoded tags
	(see "beamWithHinges")
	*/
	for (std::map<int, mpco::element::SectionAssignment>::iterator
		it = sec_assignments.begin(); it != sec_assignments.end(); ++it) {
		SectionForceDeformation *sfd = OPS_getSectionForceDeformation(it->first);
		if (sfd) {
			it->second.name = sfd->getClassType();
			/*
			$MP(2017/04/20). some fiber sections change the sign of the yLoc in the setResponse method. let's nullify this bug.
			Note that we need to keep on checking whether this bug will be removed in future versions! otherwise WE will make the bug.
			FiberSection3d SEC_TAG_FiberSection3d
			FiberSectionGJ SEC_TAG_FiberSectionGJ
			*/
			switch (sfd->getClassTag()) {
			case SEC_TAG_FiberSection3d:
			case SEC_TAG_FiberSectionGJ: {
				mpco::element::SectionAssignment &i_sec_asgn = it->second;
				for (size_t fiber_id = 0; fiber_id < i_sec_asgn.fiber_section_data.fibers.size(); fiber_id++) {
					i_sec_asgn.fiber_section_data.fibers[fiber_id].y *= -1.0;
				}
			}
			}
		}
	}
	/*
	let's find the largest section tag from sec_assignments.
	note: we are using std::map, which is sorted with the default comparator (std::less),
	so the last element has the largest key (in this case the largest section tag).
	Then move the dummy section to the sec_assignments map, using tags starting from max_sec_tag+1
	*/
	int max_sec_tag = sec_assignments.size() > 0 ? (--sec_assignments.end())->first : 0;
	for (std::map<mpco::element::FiberSectionData, mpco::element::SectionAssignment>::iterator
		it = aux_sec_assignments.begin(); it != aux_sec_assignments.end(); ++it) {
		sec_assignments[++max_sec_tag] = it->second;
	}
	/*
	write each assignment
	*/
	for (std::map<int, mpco::element::SectionAssignment>::iterator
		it = sec_assignments.begin(); it != sec_assignments.end(); ++it) {
		/*
		create the name for this group with the following format: SECTION_<TAG>[<CLASS_TYPE>]
		*/
		mpco::element::SectionAssignment &sec_asn = it->second;
		if (sec_asn.assignments.size() > 0) {
			std::stringstream ss_gp_name_buffer;
			ss_gp_name_buffer << "SECTION_" << it->first << "[" << sec_asn.name << "]";
			std::string gp_name_buffer = ss_gp_name_buffer.str();
			/*
			create the group for this section assignment
			*/
			hid_t h_gp_isec_asn = h5::group::create(h_gp_section_assignments, gp_name_buffer.c_str(), H5P_DEFAULT, m_data->info.h_group_proplist, H5P_DEFAULT);
			/*
			write attributes
			*/
			status = h5::attribute::write(h_gp_isec_asn, "ID", it->first);
			status = h5::attribute::write(h_gp_isec_asn, "NAME", sec_asn.name);
			/*
			write assignment data
			*/
			{
				hid_t dset = h5::dataset::createAndWrite(h_gp_isec_asn, "ASSIGNMENT", sec_asn.assignments);
				status = h5::dataset::close(dset);
			}
			{
				hid_t dset = h5::dataset::createAndWrite(h_gp_isec_asn, "FIBER_DATA", sec_asn.fiber_section_data.fibers);
				status = h5::dataset::close(dset);
			}
			{
				hid_t dset = h5::dataset::createAndWrite(h_gp_isec_asn, "FIBER_MATERIALS", sec_asn.fiber_section_data.materials);
				status = h5::dataset::close(dset);
			}
			/*
			close the group
			*/
			status = h5::group::close(h_gp_isec_asn);
		}
	}
	/*
	close sections and section assignments groups
	*/
	status = h5::group::close(h_gp_section_assignments);
	/*
	return
	*/
#ifdef MPCO_TIMING
	timer.stop();
#endif // MPCO_TIMING

	return retval;
}

int MPCORecorder::writeSets()
{
#ifdef MPCO_TIMING
	mpco::Timer timer("writeSets"); timer.start();
#endif // MPCO_TIMING
	/*
	error flags
	*/
	int retval = 0;
	herr_t status = 0;
	/*
	get all regions.
	note if the user is saving results not on the whole model, but only on a subset,
	the of each reagion we write only the nodes/elements contained in the subset
	and if the region (so modified) is empty, we skip it.
	*/
	std::set<int> nset;
	std::set<int> eset;
	if (m_data->has_region) {
		// create ordered sets, search will be faster
		for (size_t i = 0; i < m_data->node_set.size(); i++)
			nset.insert(m_data->node_set[i]);
		for (size_t i = 0; i < m_data->elem_set.size(); i++)
			eset.insert(m_data->elem_set[i]);
	}
	ID tags;
	m_data->info.domain->getRegionTags(tags);
	if (tags.Size() > 0) {
		/*
		write the sets group
		*/
		std::stringstream ss_gp_sets;
		ss_gp_sets << "MODEL_STAGE[" << m_data->info.current_model_stage_id << "]/MODEL/SETS";
		std::string gp_sets = ss_gp_sets.str();
		hid_t h_gp_sets = h5::group::create(m_data->info.h_file_id, gp_sets.c_str(), H5P_DEFAULT, m_data->info.h_group_proplist, H5P_DEFAULT);
		/*
		for each region
		*/
		for (int i = 0; i < tags.Size(); i++) {
			/*
			get current region
			*/
			MeshRegion *region = m_data->info.domain->getRegion(tags(i));
			if (region == NULL) continue;
			/*
			write the current set group
			*/
			std::stringstream ss_gp_curr_set;
			ss_gp_curr_set << "SET_" << tags(i);
			std::string gp_curr_set = ss_gp_curr_set.str();
			hid_t h_gp_curr_set = h5::group::create(h_gp_sets, gp_curr_set.c_str(), H5P_DEFAULT, m_data->info.h_group_proplist, H5P_DEFAULT);
			/*
			write the current set NODES dataset
			*/
			{
				const ID &nodes = region->getNodes();
				if (nodes.Size() > 0) {
					std::vector<int> aux;
					if (m_data->has_region) {
						for (int i = 0; i < nodes.Size(); i++)
							if (nset.find(nodes(i)) != nset.end())
								aux.push_back(nodes(i));
					}
					else {
						aux.resize(static_cast<size_t>(nodes.Size()));
						for (int i = 0; i < nodes.Size(); i++)
							aux[static_cast<size_t>(i)] = nodes(i);
					}
					if (aux.size() > 0) {
						hid_t dset_id = h5::dataset::createAndWrite(h_gp_curr_set, "NODES", aux);
						status = h5::dataset::close(dset_id);
					}
				}
			}
			/*
			write the current set ELEMENTS dataset
			*/
			{
				const ID &elems = region->getElements();
				if (elems.Size() > 0) {
					std::vector<int> aux;
					if (m_data->has_region) {
						for (int i = 0; i < elems.Size(); i++)
							if (eset.find(elems(i)) != eset.end())
								aux.push_back(elems(i));
					}
					else {
						aux.resize(static_cast<size_t>(elems.Size()));
						for (int i = 0; i < elems.Size(); i++)
							aux[static_cast<size_t>(i)] = elems(i);
					}
					if (aux.size() > 0) {
						hid_t dset_id = h5::dataset::createAndWrite(h_gp_curr_set, "ELEMENTS", aux);
						status = h5::dataset::close(dset_id);
					}
				}
			}
			/*
			close the current set group
			*/
			status = h5::group::close(h_gp_curr_set);
		}
		/*
		close the sets group
		*/
		status = h5::group::close(h_gp_sets);
	}
	/*
	return
	*/
#ifdef MPCO_TIMING
	timer.stop();
#endif // MPCO_TIMING
	return retval;
}

int MPCORecorder::initNodeRecorders()
{
#ifdef MPCO_TIMING
	mpco::Timer timer("initNodeRecorders"); timer.start();
#endif // MPCO_TIMING
	/*
	clear previous recorders
	*/
	clearNodeRecorders();
	/*
	create all required requests for nodal results
	*/
	for (size_t i = 0; i < m_data->nodal_results_requests.size(); i++) {
		mpco::NodalResultType::Enum rtype = m_data->nodal_results_requests[i];
		switch (rtype)
		{
		case mpco::NodalResultType::Displacement:
			m_data->nodal_recorders[rtype] = new mpco::node::ResultRecorderDisplacement(m_data->info);
			break;
		case mpco::NodalResultType::Rotation:
			m_data->nodal_recorders[rtype] = new mpco::node::ResultRecorderRotation(m_data->info);
			break;
		case mpco::NodalResultType::Pressure:
			m_data->nodal_recorders[rtype] = new mpco::node::ResultRecorderPressure(m_data->info);
			break;
		case mpco::NodalResultType::ReactionForce:
			m_data->nodal_recorders[rtype] = new mpco::node::ResultRecorderReactionForce(m_data->info);
			break;
		case mpco::NodalResultType::ReactionMoment:
			m_data->nodal_recorders[rtype] = new mpco::node::ResultRecorderReactionMoment(m_data->info);
			break;
		case mpco::NodalResultType::ReactionForceIncludingInertia:
			m_data->nodal_recorders[rtype] = new mpco::node::ResultRecorderReactionForceIncIntertia(m_data->info);
			break;
		case mpco::NodalResultType::ReactionMomentIncludingInertia:
			m_data->nodal_recorders[rtype] = new mpco::node::ResultRecorderReactionMomentIncIntertia(m_data->info);
			break;
		case mpco::NodalResultType::RayleighForce:
			m_data->nodal_recorders[rtype] = new mpco::node::ResultRecorderReactionForceRayleigh(m_data->info);
			break;
		case mpco::NodalResultType::RayleighMoment:
			m_data->nodal_recorders[rtype] = new mpco::node::ResultRecorderReactionMomentRayleigh(m_data->info);
			break;
		case mpco::NodalResultType::Velocity:
			m_data->nodal_recorders[rtype] = new mpco::node::ResultRecorderVelocity(m_data->info);
			break;
		case mpco::NodalResultType::AngularVelocity:
			m_data->nodal_recorders[rtype] = new mpco::node::ResultRecorderAngularVelocity(m_data->info);
			break;
		case mpco::NodalResultType::Acceleration:
			m_data->nodal_recorders[rtype] = new mpco::node::ResultRecorderAcceleration(m_data->info);
			break;
		case mpco::NodalResultType::AngularAcceleration:
			m_data->nodal_recorders[rtype] = new mpco::node::ResultRecorderAngularAcceleration(m_data->info);
			break;
		case mpco::NodalResultType::ModesOfVibration:
			m_data->nodal_recorders[rtype] = new mpco::node::ResultRecorderModesOfVibration(m_data->info);
			break;
		case mpco::NodalResultType::ModesOfVibrationRotational:
			m_data->nodal_recorders[rtype] = new mpco::node::ResultRecorderModesOfVibrationRotational(m_data->info);
			break;
		case mpco::NodalResultType::DisplacementSensitivity:
			m_data->nodal_recorders[rtype] = new mpco::node::ResultRecorderDisplacementSensitivity(m_data->info, m_data->sens_grad_indices[i]);
			break;
		case mpco::NodalResultType::RotationSensitivity:
			m_data->nodal_recorders[rtype] = new mpco::node::ResultRecorderRotationSensitivity(m_data->info, m_data->sens_grad_indices[i]);
			break;
		case mpco::NodalResultType::VelocitySensitivity:
			m_data->nodal_recorders[rtype] = new mpco::node::ResultRecorderVelocitySensitivity(m_data->info, m_data->sens_grad_indices[i]);
			break;
		case mpco::NodalResultType::AngularVelocitySensitivity:
			m_data->nodal_recorders[rtype] = new mpco::node::ResultRecorderAngularVelocitySensitivity(m_data->info, m_data->sens_grad_indices[i]);
			break;
		case mpco::NodalResultType::AccelerationSensitivity:
			m_data->nodal_recorders[rtype] = new mpco::node::ResultRecorderAccelerationSensitivity(m_data->info, m_data->sens_grad_indices[i]);
			break;
		case mpco::NodalResultType::AngularAccelerationSensitivity:
			m_data->nodal_recorders[rtype] = new mpco::node::ResultRecorderAngularAccelerationSensitivity(m_data->info, m_data->sens_grad_indices[i]);
			break;
		default:
			break;
		}
	}
	/*
	return
	*/
#ifdef MPCO_TIMING
	timer.stop();
#endif // MPCO_TIMING
	return 0;
}

int MPCORecorder::clearNodeRecorders()
{
	for (mpco::node::ResultRecorderMap::iterator it = m_data->nodal_recorders.begin();
		it != m_data->nodal_recorders.end(); ++it) {
		mpco::node::ResultRecorder *recorder = it->second;
		if (recorder)
			delete recorder;
	}
	m_data->nodal_recorders.clear();
	return 0;
}

int MPCORecorder::initElementRecorders()
{
#ifdef MPCO_TIMING
	mpco::Timer timer("initElementRecorders"); timer.start();
#endif // MPCO_TIMING
	/*
	clear previous recorders
	*/
	clearElementRecorders();	
	/*
	generate elemental recorders, based on the element type, and the element output descriptor.
	The data size (num columns) of the same response may change even for element of the same class,
	for example (in some beamns) if the number of integration points is different from beam to beam...
	*/
	m_data->elemental_recorders.resize(m_data->elemental_results_requests.size());
	for (size_t cnt_request = 0; cnt_request < m_data->elemental_results_requests.size(); cnt_request++) {
		const std::vector<std::string> &request = m_data->elemental_results_requests[cnt_request];
		/*
		copy request to **argv
		*/
		std::vector<std::string> request_mod;
		bool do_all_materials = false;
		bool do_all_sections = false;
		bool do_all_fibers = false;
		size_t material_id_placeholder_index = 0;
		size_t section_id_placeholder_index = 0;
		size_t fiber_id_placeholder_index = 0;
		for (size_t i = 0; i < request.size(); i++) {
			const std::string &current_request = request[i];
			request_mod.push_back(current_request);
			if (i == 0) { // note: all sections or all materials available only as first options
				if (current_request == "section") {
					request_mod.push_back(""); // placeholder for section id
					section_id_placeholder_index = request_mod.size() - 1;
					do_all_sections = true;
				}
				else if (current_request == "material") {
					request_mod.push_back(""); // placeholder for material id
					material_id_placeholder_index = request_mod.size() - 1;
					do_all_materials = true;
				}
			}
			else {
				if (current_request == "fiber" && do_all_sections) {
					request_mod.push_back(""); // placeholder for fiber id
					fiber_id_placeholder_index = request_mod.size() - 1;
					do_all_fibers = true;
				}
			}
		}
		int argc = (int)request_mod.size();
		const char **argv = new const char*[argc];
		for (size_t i = 0; i < request_mod.size(); i++)
			argv[i] = request_mod[i].c_str();
		/**
		$WO:SHELL_SEC_KEYWORD
		*/
		std::string aux_section_keyword_standard = "section";
		std::string aux_section_keyword_for_shells = "material";
		/*
		populate the current elemental result
		*/
		mpco::element::ResultRecorder &recorder = m_data->elemental_recorders[cnt_request];
		recorder.result_request = request;
		/*
		loop over all element group-by-tag
		*/
		for (mpco::element::ElementCollection::submap_type::iterator it1 = m_data->elements.items.begin();
			it1 != m_data->elements.items.end(); ++it1) {
			mpco::element::ElementWithSameClassTagCollection &elem_by_tag = it1->second;
			mpco::element::OutputWithSameClassTagCollection &eo_by_tag = recorder.response_map[it1->first];
			/**
			$WO:SHELL_SEC_KEYWORD
			if tag is that of those shells expecting material keyword instead of section keyword...
			*/
			bool standard_section_keyword_modified(false);
			if (do_all_sections && utils::shell::isShellElementTag(it1->first)) {
				argv[0] = aux_section_keyword_for_shells.c_str();
				standard_section_keyword_modified = true;
			}
			/*
			loop over all element group-by-integration rule
			*/
			for (mpco::element::ElementWithSameClassTagCollection::submap_type::iterator it2 = elem_by_tag.items.begin();
				it2 != elem_by_tag.items.end(); ++it2) {
				mpco::element::ElementWithSameIntRuleCollection &elem_by_rule = it2->second;
				mpco::element::OutputWithSameIntRuleCollection &eo_by_rule = eo_by_tag.items[it2->first];
				/*
				loop over all element group-by-custom(or default) integration rule
				*/
				for (mpco::element::ElementWithSameIntRuleCollection::submap_type::iterator it3 = elem_by_rule.items.begin();
					it3 != elem_by_rule.items.end(); ++it3) {
					mpco::element::ElementWithSameCustomIntRuleCollection &elem_by_custom_rule = it3->second;
					mpco::element::OutputWithSameCustomIntRuleCollection &eo_by_custom_rule = eo_by_rule.items[it3->first];
					int header_local_index = 0;
					/*
					now loop over all mapped elements, call setResponse to 1) create the response and 2) create the descriptor
					*/
					for (mpco::element::ElementWithSameCustomIntRuleCollection::collection_type::iterator it4 = elem_by_custom_rule.items.begin();
						it4 != elem_by_custom_rule.items.end(); ++it4) {
						Element *elem = *it4;
						mpco::element::OutputDescriptor eo_descriptor;
						Response *eo_response = 0;
						mpco::element::OutputDescriptorStream eo_stream(&eo_descriptor);
						/*
						now loop over all integration points with cross section / material if a section/material response is required
						*/
						if (do_all_sections) {
							// find the vector of gauss points and number of fibers mapped to this element
							std::map<int, std::vector<std::pair<int, int> > >::const_iterator it_elem_ngauss_nfiber = m_data->elem_ngauss_nfiber_info.find(elem->getTag());
							if (it_elem_ngauss_nfiber != m_data->elem_ngauss_nfiber_info.end()) {
								CompositeResponse *sec_comp_response = new CompositeResponse();
								int num_sec_responses = 0;
								/*
								for each section (gauss point)...
								*/
								for(size_t section_id = 0; section_id < it_elem_ngauss_nfiber->second.size(); section_id++) { 
									// set section id in argv
									std::stringstream ss_section_id; ss_section_id << section_id + 1; /* note: in setResponse gp number is 1-based */
									std::string s_section_id = ss_section_id.str();
									argv[section_id_placeholder_index] = s_section_id.c_str();
									// set response
									Response *sec_response = 0;
									/*
									now loop over all fibers if a fiber response is required
									*/
									if (do_all_fibers) {
										CompositeResponse *fib_comp_response = new CompositeResponse();
										int num_fib_responses = 0;
										/*
										for each fiber...
										*/
										int fiber_base_id = it_elem_ngauss_nfiber->second[section_id].first;
										int num_fibers = it_elem_ngauss_nfiber->second[section_id].second;
										for(int fiber_id = 0; fiber_id < num_fibers; fiber_id++) { 
											// set fiber id in argv
											std::string s_fiber_id = utils::strings::to_string(fiber_id + fiber_base_id); /* note: in setResponse fiber index is 0-based */
											argv[fiber_id_placeholder_index] = s_fiber_id.c_str();
											// set response
											Response *fib_response = 0;
											fib_response = elem->setResponse(argv, argc, eo_stream);
											if (fib_response) {
												num_fib_responses = fib_comp_response->addResponse(fib_response);
											}
										} // end fiber while loop
										if (num_fib_responses == 0) { // no valid fiber responses found
											delete fib_comp_response; // here sec_response is left null -> exit also from section loop
										}
										else {
											sec_response = fib_comp_response;
										}
									}
									else {
										sec_response = elem->setResponse(argv, argc, eo_stream);
									}
									if (sec_response) {
										num_sec_responses = sec_comp_response->addResponse(sec_response);
									}
								} // end section loop
								if (num_sec_responses == 0) { // no valid section responses found
									delete sec_comp_response;
								}
								else {
									eo_response = sec_comp_response;
								}
							}
						}
						else if (do_all_materials) {
							/*
							for each material (gauss point)...
							*/
							CompositeResponse *mat_comp_response = new CompositeResponse();
							int num_mat_responses = 0;
							int material_id = 0; /* note: in setResponse gp number is 1-based */
							while (true) {
								material_id++;
								if (material_id > MPCO_MAX_TRIAL_NSEC) {
									break;
								}
								// set material id in argv
								std::stringstream ss_material_id; ss_material_id << material_id;
								std::string s_material_id = ss_material_id.str();
								argv[material_id_placeholder_index] = s_material_id.c_str();
								// set response
								Response *mat_response = elem->setResponse(argv, argc, eo_stream);
								if (mat_response) {
									num_mat_responses = mat_comp_response->addResponse(mat_response);
								}
								else {
									// not response at this material point, break the loop, probably we reached the element num of gauss...
									break;
								}
							} // end material loop
							if (num_mat_responses == 0) { // no valid material responses found
								delete mat_comp_response;
							}
							else {
								eo_response = mat_comp_response;
							}
						}
						else {
							eo_response = elem->setResponse(argv, argc, eo_stream);
						}
						eo_stream.finalizeSetResponse();
						if (do_all_fibers) {
							eo_descriptor.purge();
						}
						if (eo_response) {
							/*
							get (or create and get) the list of MPCORecorder_ElementResultRecorder mapped to this descriptor
							and add the new element-response pair
							*/
							mpco::element::OutputDescriptorHeader header(eo_descriptor.makeHeader());
							header.workaroundForSizeInconsistency(eo_response->getInformation().getData().Size());
							header.workaroundForDuplicatedComponents();
							mpco::element::OutputResponseCollection &eo_by_header = eo_by_custom_rule.items[header];
							if (eo_by_header.is_new) {
								/*
								create a name for this dataset using the following format <class_tag>-<class_name>[<integration_rule>:<custom_rule_index>:<header_index>]
								*/
								std::stringstream ss_buffer;
								ss_buffer << elem_by_tag.class_tag << "-" << elem_by_tag.class_name
									<< "[" << elem_by_rule.int_rule_type << ":"
									<< elem_by_custom_rule.custom_int_rule_index << ":"
									<< header_local_index++ << "]";
								eo_by_header.dir_name = ss_buffer.str();
								eo_by_header.is_new = false;
							}
							eo_by_header.items.push_back(mpco::element::OutputResponse(elem, eo_response));
							/*
							add this response pointer to the list of all responses. it will be used to release memory later on
							*/
							m_data->elemental_responses.push_back(eo_response);
						}
					}
				}
			}
			/**
			$WO:SHELL_SEC_KEYWORD
			reset standard keyword
			*/
			if (standard_section_keyword_modified) {
				argv[0] = aux_section_keyword_standard.c_str();
				standard_section_keyword_modified = false;
			}
		}
		/*
		delete argv, do not delete its contents! they are handled by the string object
		*/
		delete[] argv;
	}
	/*
	return
	*/
#ifdef MPCO_TIMING
	timer.stop();
#endif // MPCO_TIMING
	return 0;
}

int MPCORecorder::clearElementRecorders()
{
	/*
	manually destroy all responses, they are not handled by elemental recorders
	*/
	for (size_t i = 0; i < m_data->elemental_responses.size(); i++) {
		if (m_data->elemental_responses[i])
			delete m_data->elemental_responses[i];
	}
	m_data->elemental_responses.clear();
	/*
	clear recorders
	*/
	m_data->elemental_recorders.clear();
	return 0;
}

int MPCORecorder::recordResultsOnNodes()
{
#ifdef MPCO_TIMING
	mpco::Timer timer("recordResultsOnNodes"); timer.start();
#endif // MPCO_TIMING
	/*
	error flags
	*/
	int retval = 0;
	herr_t status = 0;
	/*
	quick return
	*/
	size_t num_nodes = m_data->nodes.size();
	if (num_nodes == 0) 
		return retval;
	/*
	some preprocessing before recording modes of vibration.
	*/
	m_data->info.record_eigen_on_this_step = false;
	int num_eigen = *OPS_GetNumEigen();
	if (num_eigen > 0) { 
		/*
		only if we have eigenvalues
		*/
		bool eigen_requested = false;
		if (m_data->nodal_recorders.find(mpco::NodalResultType::ModesOfVibration) != m_data->nodal_recorders.end())
			eigen_requested = true;
		else if (m_data->nodal_recorders.find(mpco::NodalResultType::ModesOfVibrationRotational) != m_data->nodal_recorders.end())
			eigen_requested = true;
		if (eigen_requested) {
			/*
			only if eigenvalues have been requested by the user
			*/
			double eigen_set_time = m_data->info.domain->getTimeEigenvaluesSet();
			const Vector &new_eigenvalues = m_data->info.domain->getEigenvalues();
			if (m_data->info.eigen_first_initialization_done) {
				/*
				check for a new eigenvalue analysis
				*/
				if (std::abs(eigen_set_time - m_data->info.eigen_last_time_set) > 1.0e-10) {
					m_data->info.record_eigen_on_this_step = true;
				}
				else if (num_eigen != m_data->info.eigen_last_values.Size()) {
					m_data->info.record_eigen_on_this_step = true;
				}
				else {
					/*
					check eigenvalues one-by-one
					*/
					for (int i = 0; i < num_eigen; i++) {
						if (std::abs(new_eigenvalues[i] - m_data->info.eigen_last_values[i]) > 1.0e-10) {
							m_data->info.record_eigen_on_this_step = true;
							break;
						}
					}
				}
			}
			else {
				/*
				write it anyway, this is the first time
				*/
				m_data->info.record_eigen_on_this_step = true;
				m_data->info.eigen_first_initialization_done = true;
			}
			if (m_data->info.record_eigen_on_this_step) {
				m_data->info.eigen_last_time_set = eigen_set_time;
				m_data->info.eigen_last_values = new_eigenvalues;
			}
		}
	}
	/*
	now we can record all nodal results
	*/
	int previous_reac_type(-1);
	for (mpco::node::ResultRecorderMap::iterator it = m_data->nodal_recorders.begin();
		it != m_data->nodal_recorders.end(); ++it) {
		mpco::node::ResultRecorder *nodal_recorder = it->second;
		if (nodal_recorder) {
			int curr_reac_type = nodal_recorder->getReactionFlag();
			if (curr_reac_type != previous_reac_type) {
				if (curr_reac_type > -1 && curr_reac_type < 3) {
					m_data->info.domain->calculateNodalReactions(curr_reac_type);
				}
				previous_reac_type = curr_reac_type;
			}
			nodal_recorder->record(m_data->info, m_data->nodes);
		}
	}
	/*
	return
	*/
#ifdef MPCO_TIMING
	timer.stop();
#endif // MPCO_TIMING
	return retval;
}

int MPCORecorder::recordResultsOnElements()
{
#ifdef MPCO_TIMING
	mpco::Timer timer("recordResultsOnElements"); timer.start();
#endif // MPCO_TIMING
	/*
	error flags
	*/
	int retval = 0;
	herr_t status = 0;
	/*
	loop over all requested elemental results
	*/
	for (mpco::element::ResultRecorderCollection::iterator it_recorder = m_data->elemental_recorders.begin();
		it_recorder != m_data->elemental_recorders.end(); ++it_recorder) {
		/*
		get the current recorder
		*/
		mpco::element::ResultRecorder &recorder = *it_recorder;
		/*
		create a representative name for this result
		*/
		std::string result_name;
		std::string result_display_name;
		{
			std::stringstream ss_buffer;
			for (int i = 0; i < recorder.result_request.size(); i++) {
				if (i > 0) ss_buffer << ".";
				ss_buffer << recorder.result_request[i];
			}
			result_display_name = ss_buffer.str();
		}
		{
			std::stringstream ss_buffer;
			ss_buffer << "MODEL_STAGE[" << m_data->info.current_model_stage_id << "]/RESULTS/ON_ELEMENTS/" << result_display_name;
			result_name = ss_buffer.str();
		}
		/*
		generate the result directory on the first record
		*/
		if (!recorder.initialized) {
			/*
			\todo all attributes are temporary, we need a user-defined metadata for elemental results...
			*/
			hid_t h_gp_result = h5::group::createResultGroup(m_data->info.h_file_id, m_data->info.h_group_proplist,
				result_name, result_display_name, "", 0, "", "",
				(int)mpco::ResultType::Generic, (int)mpco::ResultDataType::Scalar);
			status = h5::group::close(h_gp_result);
			recorder.initialized = true;
		}
		/*
		loop over all eo group-by-tag
		*/
		for (mpco::element::ResultRecorder::collection_type::iterator it1 = recorder.response_map.begin();
			it1 != recorder.response_map.end(); ++it1) {
			mpco::element::OutputWithSameClassTagCollection &eo_by_tag = it1->second;
			/*
			loop over all eo group-by-integration rule
			*/
			for (mpco::element::OutputWithSameClassTagCollection::collection_type::iterator it2 = eo_by_tag.items.begin();
				it2 != eo_by_tag.items.end(); ++it2) {
				mpco::element::OutputWithSameIntRuleCollection &eo_by_rule = it2->second;
				/*
				loop over all eo group-by-custom(or default) integration rule
				*/
				for (mpco::element::OutputWithSameIntRuleCollection::collection_type::iterator it3 = eo_by_rule.items.begin();
					it3 != eo_by_rule.items.end(); ++it3) {
					mpco::element::OutputWithSameCustomIntRuleCollection &eo_by_custom_rule = it3->second;
					/*
					loop over all eo group-by-header
					*/
					for (mpco::element::OutputWithSameCustomIntRuleCollection::collection_type::iterator it4 = eo_by_custom_rule.items.begin();
						it4 != eo_by_custom_rule.items.end(); ++it4) {
						const mpco::element::OutputDescriptorHeader &header = it4->first;
						mpco::element::OutputResponseCollection &eo_by_header = it4->second;
						size_t num_rows = eo_by_header.items.size();
						/*
						create the name for this group
						*/
						std::string header_dir_name;
						{
							std::stringstream ss_buffer;
							ss_buffer << result_name << "/" << eo_by_header.dir_name;
							header_dir_name = ss_buffer.str();
						}
						if (!eo_by_header.initialized) {
							/*
							create the header group
							*/
							hid_t h_gp_header = h5::group::create(m_data->info.h_file_id, header_dir_name.c_str(), H5P_DEFAULT, m_data->info.h_group_proplist, H5P_DEFAULT);
							status = h5::attribute::write(h_gp_header, "NUM_COLUMNS", header.num_columns);
							/*
							add attributes for metadata
							*/
							hid_t h_gp_meta = h5::group::create(h_gp_header, "META", H5P_DEFAULT, m_data->info.h_group_proplist, H5P_DEFAULT);
							hid_t h_dset_mult = h5::dataset::createAndWrite(h_gp_meta, "MULTIPLICITY", header.multiplicity, header.multiplicity.size(), 1);
							status = h5::dataset::close(h_dset_mult);
							hid_t h_dset_gauss = h5::dataset::createAndWrite(h_gp_meta, "GAUSS_IDS", header.gauss_id, header.gauss_id.size(), 1);
							status = h5::dataset::close(h_dset_gauss);
							hid_t h_dset_ncomp = h5::dataset::createAndWrite(h_gp_meta, "NUM_COMPONENTS", header.num_components, header.num_components.size(), 1);
							status = h5::dataset::close(h_dset_ncomp);
							std::stringstream ss_buffer_components;
							for (size_t i = 0; i < header.components_path.size(); i++) {
								if(i > 0) ss_buffer_components << ";";
								const std::vector<int> &ipath = header.components_path[i];
								const std::vector<std::string> &icmps = header.components[i];
								for (size_t j = 0; j < ipath.size(); j++) {
									ss_buffer_components << ipath[j] << ".";
								}
								for (size_t j = 0; j < icmps.size(); j++) {
									if (j > 0) ss_buffer_components << ",";
									ss_buffer_components << icmps[j];
								}
							}
							std::string buffer_components = ss_buffer_components.str();
							hid_t h_dset_comp = h5::dataset::createAndWrite(h_gp_meta, "COMPONENTS", buffer_components);
							status = h5::dataset::close(h_dset_comp);
							status = h5::group::close(h_gp_meta);
							/*
							create the id dataset
							*/
							std::vector<int> buffer_id(num_rows);
							for (size_t i = 0; i < num_rows; i++)
								buffer_id[i] = eo_by_header.items[i].element->getTag();
							hid_t h_dset_id = h5::dataset::createAndWrite(h_gp_header, "ID", buffer_id, num_rows, 1);
							/*
							create the data group
							*/
							hid_t h_gp_data = h5::group::create(h_gp_header, "DATA", H5P_DEFAULT, m_data->info.h_group_proplist, H5P_DEFAULT);
							/*
							done
							*/
							status = h5::group::close(h_gp_data);
							status = h5::dataset::close(h_dset_id);
							status = h5::group::close(h_gp_header);
							eo_by_header.initialized = true;
						}
						/*
						create the dataset for this timestep
						*/
						std::stringstream ss_dset_name;
						ss_dset_name << header_dir_name << "/DATA/STEP_" << m_data->info.current_time_step_id;
						std::string dset_name = ss_dset_name.str();
						std::vector<double> buffer_data(num_rows * header.num_columns);
						for (size_t i = 0; i < num_rows; i++) {
							mpco::element::OutputResponse &current_response = eo_by_header.items[i];
							current_response.response->getResponse();
							const Vector &current_data = current_response.response->getInformation().getData();
							if (current_data.Size() != header.num_columns) {
								opserr << "MPCORecorder Error: invalid response data size\n";
								continue;
							}
							size_t offset = i * header.num_columns;
							for (size_t j = 0; j < header.num_columns; j++)
								buffer_data[offset + j] = current_data[(int)j];
						}
						hid_t h_dset_data = h5::dataset::createAndWrite(m_data->info.h_file_id, dset_name.c_str(), buffer_data, num_rows, header.num_columns);
						status = h5::attribute::write(h_dset_data, "STEP", m_data->info.current_time_step_id);
						status = h5::attribute::write(h_dset_data, "TIME", m_data->info.current_time_step);
						status = h5::dataset::close(h_dset_data);
					}
				}
			}
		}		
	}
	/*
	return
	*/
#ifdef MPCO_TIMING
	timer.stop();
#endif // MPCO_TIMING
	return retval;
}

/*************************************************************************************

MPCORecorder generator

**************************************************************************************/

/* static class instance counter */
static int mpco_recorder_counter = 0;

void* OPS_MPCORecorder() 
{
	// on first call
	if (mpco_recorder_counter == 0) {
		mpco_recorder_counter++;
		opserr << 
			"MPCO recorder - Written by ASDEA Software Technology: M.Petracca, G.Camata\n"
			"ASDEA Software Technology: https://asdeasoft.net \n"
			"STKO (Scientific ToolKit for OpenSees): https://asdeasoft.net/stko/ \n"
			"If you use this tool, please cite us:\n"
			"Petracca, M., Candeloro, F., & Camata, G. (2017). \"STKO user manual\". ASDEA Software Technology, Pescara Italy.\n"
			;
	}

#ifdef MPCO_HDF5_LOADED_AT_RUNTIME
	if (!LibraryLoader::instance().loaded) {
		opserr << "MPCORecorder error loading HDF5\n";
		return 0;
	}
#endif

	// check the minimum number of arguments -> filename
	int numdata = OPS_GetNumRemainingInputArgs();
	if (numdata < 1) {
		opserr << "MPCORecorder error: insufficient number of arguments\n";
		return 0;
	}

	// filename
	const char *filename = OPS_GetString();
	numdata--;

	// get the current domain
	Domain *domain = OPS_GetDomain();
	if (domain == 0) {
		opserr << "MPCORecorder error: domain is not defined\n";
		return 0;
	}

	// parse for optional arguments
	utils::parsing::option_type curr_opt = utils::parsing::opt_none;
	std::vector<mpco::NodalResultType::Enum> nodal_results_requests;
	std::vector<int> sens_grad_indices;
	std::vector<std::vector<std::string> > elemental_results_requests;
	std::vector<std::string> tokens;
	mpco::OutputFrequency output_freq;
	bool has_region = false;
	std::set<int> node_set;
	std::set<int> elem_set;
	int one_item = 1;

	while (numdata > 0) {
		const char* data = OPS_GetString();
		numdata--;
		if (strcmp(data, "-N") == 0) {
			curr_opt = utils::parsing::opt_result_on_nodes;
		}
		else if (strcmp(data, "-NS") == 0) {
			curr_opt = utils::parsing::opt_result_on_nodes_sens;
		}
		else if (strcmp(data, "-E") == 0) {
			curr_opt = utils::parsing::opt_result_on_elements;
		}
		else if (strcmp(data, "-T") == 0) {
			curr_opt = utils::parsing::opt_time;
			output_freq.reset();
		}
		else if (strcmp(data, "-R") == 0) {
			curr_opt = utils::parsing::opt_region;
			if (numdata > 0) {
				int region_tag = 0;
				if (OPS_GetInt(&one_item, &region_tag) != 0) {
					opserr << "MPCORecorder error: option -R (region) requires an extra parameter (int) for the region tag. (cannot get int value)\n";
					return 0;
				}
				MeshRegion *region = domain->getRegion(region_tag);
				if (region == 0) {
					opserr << "MPCORecorder error: region " << region_tag << " is null\n";
					return 0;
				}
				const ID &node_ids = region->getNodes();
				const ID &elem_ids = region->getElements();
				for (int i = 0; i < node_ids.Size(); i++) 
					node_set.insert(node_ids(i));
				for (int i = 0; i < elem_ids.Size(); i++)
					elem_set.insert(elem_ids(i));
				has_region = true;
				numdata--;
			}
			else {
				opserr << "MPCORecorder error: option -R (region) requires an extra parameter (int) for the region tag\n";
				return 0;
			}
			has_region = true;
		}
		else {
			switch (curr_opt)
			{
			case utils::parsing::opt_result_on_nodes: {
				if (strcmp(data, "displacement") == 0)
					nodal_results_requests.push_back(mpco::NodalResultType::Displacement);
				else if (strcmp(data, "rotation") == 0)
					nodal_results_requests.push_back(mpco::NodalResultType::Rotation);
				else if (strcmp(data, "velocity") == 0)
					nodal_results_requests.push_back(mpco::NodalResultType::Velocity);
				else if (strcmp(data, "angularVelocity") == 0)
					nodal_results_requests.push_back(mpco::NodalResultType::AngularVelocity);
				else if (strcmp(data, "acceleration") == 0)
					nodal_results_requests.push_back(mpco::NodalResultType::Acceleration);
				else if (strcmp(data, "angularAcceleration") == 0)
					nodal_results_requests.push_back(mpco::NodalResultType::AngularAcceleration);
				else if (strcmp(data, "reactionForce") == 0)
					nodal_results_requests.push_back(mpco::NodalResultType::ReactionForce);
				else if (strcmp(data, "reactionMoment") == 0)
					nodal_results_requests.push_back(mpco::NodalResultType::ReactionMoment);
				else if (strcmp(data, "reactionForceIncludingInertia") == 0)
					nodal_results_requests.push_back(mpco::NodalResultType::ReactionForceIncludingInertia);
				else if (strcmp(data, "reactionMomentIncludingInertia") == 0)
					nodal_results_requests.push_back(mpco::NodalResultType::ReactionMomentIncludingInertia);
				else if (strcmp(data, "rayleighForce") == 0)
					nodal_results_requests.push_back(mpco::NodalResultType::RayleighForce);
				else if (strcmp(data, "rayleighMoment") == 0)
					nodal_results_requests.push_back(mpco::NodalResultType::RayleighMoment);
				else if (strcmp(data, "pressure") == 0)
					nodal_results_requests.push_back(mpco::NodalResultType::Pressure);
				else if (strcmp(data, "modesOfVibration") == 0)
					nodal_results_requests.push_back(mpco::NodalResultType::ModesOfVibration);
				else if (strcmp(data, "modesOfVibrationRotational") == 0)
					nodal_results_requests.push_back(mpco::NodalResultType::ModesOfVibrationRotational);
				else {
					opserr << "MPCORecorder error: option -N with unknown result type (" << data << ")\n";
					return 0;
				}
				sens_grad_indices.push_back(0); // just a place holder
				break;
			}
			case utils::parsing::opt_result_on_nodes_sens: {
				if (strcmp(data, "displacementSensitivity") == 0)
					nodal_results_requests.push_back(mpco::NodalResultType::DisplacementSensitivity);
				else if (strcmp(data, "rotationSensitivity") == 0)
					nodal_results_requests.push_back(mpco::NodalResultType::RotationSensitivity);
				else if (strcmp(data, "velocitySensitivity") == 0)
					nodal_results_requests.push_back(mpco::NodalResultType::VelocitySensitivity);
				else if (strcmp(data, "angularVelocitySensitivity") == 0)
					nodal_results_requests.push_back(mpco::NodalResultType::AngularVelocitySensitivity);
				else if (strcmp(data, "accelerationSensitivity") == 0)
					nodal_results_requests.push_back(mpco::NodalResultType::AccelerationSensitivity);
				else if (strcmp(data, "angularAccelerationSensitivity") == 0)
					nodal_results_requests.push_back(mpco::NodalResultType::AngularAccelerationSensitivity);
				else {
					opserr << "MPCORecorder error: option -N with unknown result type (" << data << ")\n";
					return 0;
				}
				if (numdata > 0) {
					int grad_index;
					if (OPS_GetInt(&one_item, &grad_index) != 0) {
						opserr << "MPCORecorder error: option -NS requires an extra parameter (int) for the sensitivity parameter index. (cannot get int value)\n";
						return 0;
					}
					numdata--;
					sens_grad_indices.push_back(grad_index);
				}
				else {
					opserr << "MPCORecorder error: option -NS requires an extra argument for the sensitivity parameter index\n";
					return 0;
				}
				break;
			}
			case utils::parsing::opt_result_on_elements: {
				std::string temp(data);
				utils::strings::split(temp, '.', tokens, true);
				if(tokens.size() > 0)
					elemental_results_requests.push_back(tokens);
				break;
			}
			case utils::parsing::opt_time: {
				if (strcmp(data, "dt") == 0) {
					output_freq.type = mpco::OutputFrequency::DeltaTime;
					output_freq.nsteps = 1;
					if (numdata > 0) {
						if (OPS_GetDouble(&one_item, &output_freq.dt) != 0) {
							opserr << "MPCORecorder error: invalid double argument for the delta time\n";
							return 0;
						}
						if (output_freq.dt < 0.0) output_freq.dt = 0.0;
						numdata--;
					}
					else {
						opserr << "MPCORecorder error: option -T with type dt requires an extra argument for the delta time\n";
						return 0;
					}
				}
				else if (strcmp(data, "nsteps") == 0) {
					output_freq.type = mpco::OutputFrequency::NumberOfSteps;
					output_freq.dt = 0.0;
					if (numdata > 0) {
						if (OPS_GetInt(&one_item, &output_freq.nsteps) != 0) {
							opserr << "MPCORecorder error: invalid int argument for the number of steps\n";
							return 0;
						}
						if (output_freq.nsteps < 1) output_freq.nsteps = 1; // make sure it's positive
						numdata--;
					}
					else {
						opserr << "MPCORecorder error: option -T with type dt requires an extra argument for the delta time\n";
						return 0;
					}
				}
				else {
					opserr << "MPCORecorder error: option -T with unknown frquency type (" << data << ")\n";
					return 0;
				}
				break;
			}
			default: {
				opserr << "MPCORecorder error: unknown arg with option none " << data << "\n";
				return 0;
			}
			}
		}
	}

	// create recorder
	MPCORecorder *new_recorder = new MPCORecorder();
	new_recorder->m_data->filename = filename;
	new_recorder->m_data->output_freq = output_freq;
	new_recorder->m_data->nodal_results_requests.swap(nodal_results_requests);
	new_recorder->m_data->sens_grad_indices.swap(sens_grad_indices);
	new_recorder->m_data->elemental_results_requests.swap(elemental_results_requests);
	new_recorder->m_data->has_region = has_region;
	for (std::set<int>::const_iterator it = node_set.begin(); it != node_set.end(); ++it)
		new_recorder->m_data->node_set.push_back(*it);
	for (std::set<int>::const_iterator it = elem_set.begin(); it != elem_set.end(); ++it)
		new_recorder->m_data->elem_set.push_back(*it);
	return new_recorder;
}

