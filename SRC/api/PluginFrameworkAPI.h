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

#ifndef PluginFrameworkAPI_h
#define PluginFrameworkAPI_h

// Written: Massimo Petracca 
// Created: 02/2020
// Revision: A
//
// Description: This file contains the API for the PluginFramework,
// that is, the interface between the PluginFramework and the 
// external libraries that will be used as plugins

#include <stdint.h>

/**
The jobs a plugin material should perform
*/
enum PluginMaterialJobType
{
	PF_MAT_GET_INIT_INFO = 0,
	PF_MAT_INITIALIZE,
	PF_MAT_FINALIZE,

	PF_MAT_COMMIT = 100,
	PF_MAT_REVERT,
	PF_MAT_REVERT_TO_START,

	PF_MAT_SERIALIZE = 200,
	PF_MAT_DESERIALIZE,

	PF_MAT_COMPUTE = 300,

	PF_MAT_GET_STRAIN = 400,
	PF_MAT_GET_STRAIN_RATE,
	PF_MAT_GET_STRESS,
	PF_MAT_GET_TANGENT,
	PF_MAT_GET_INITIAL_TANGENT,
	PF_MAT_GET_DAMP_TANGENT,
	PF_MAT_GET_RHO,
	PF_MAT_GET_ENERGY,
	PF_MAT_GET_IS_FAILED,
	PF_MAT_GET_RESPONSE,
	PF_MAT_GET_THERMAL_TANG_ELONG,
	PF_MAT_GET_TEMP_ELONG,

	PF_MAT_GET_VARIABLE = 500,
	PF_MAT_SET_VARIABLE,

	PF_MAT_UPDATE_PARAMETER = 600,
	PF_MAT_ACTIVATE_PARAMETER,

	PF_MAT_GET_RESPONSE_SENSITIVITY = 700,
	PF_MAT_GET_STRAIN_SENSITIVITY,
	PF_MAT_GET_STRESS_SENSITIVITY,
	PF_MAT_GET_TANGENT_SENSITIVITY,
	PF_MAT_GET_INITIAL_TANGENT_SENSITIVITY,
	PF_MAT_GET_DAMP_TANGENT_SENSITIVITY,
	PF_MAT_GET_RHO_SENSITIVITY,
	PF_MAT_COMMIT_SENSITIVITY

};

/**
Material types
*/

#define PF_MAT_TYPE_UNIAXIAL 1
#define PF_MAT_TYPE_ND_3D 2
#define PF_MAT_TYPE_ND_PLANE_STRAIN 3
#define PF_MAT_TYPE_ND_PLANE_STRESS 4
#define PF_MAT_TYPE_ND_AXISYMMETRIC 5
#define PF_MAT_TYPE_ND_PLATE_FIBER 6


typedef struct PluginMaterialData PluginMaterialData;

/**
The plugin material procedure.
Inputs:
data - contains all working data, input and output for material calculation and response requests
job - input tag of the current requested job the material should perform
rcode - the return code, 0 if successful, error code in any other case
*/
typedef void (*PluginMaterialProc)(PluginMaterialData* data, PluginMaterialJobType job, int* rcode);

/**
A simple structure containing exchange data for plugin unixial and nD materials.

This structure will be allocated/deallocated in OpenSees and will be passed as a pointer
to the material function. All dynamic data contained in this structure will
be set/unset in OpenSees as well. 

Only the "const char *message" 
member is a pointer to a string allocated in the plugin library and used by OpenSees
to get generic information about the plugin material. So the typical scenario would be
to store, in the plugin library, a static string, which is built and updated everytime
the plugin frameworks requests it, and the the plugin materials sets the "message"
member to point to the internal const data of that string. The plugin framework will only
read from it when requested, and will not use it in any other cases.

When job = PF_MAT_INITIALIZE:
	The format of "message" should be: 
	R1;R2;R2;...;RN
	where R1... RN are records, separated by ';'
	each record can be one of the following:
	------------------------------------------------------------------------------
	A|<name>|<default>
	R|<id>|<name>|<comp_1_name>|<comp_2_name>|...|<comp_N_name>
	P|<id>|<name>
	V|<id>|<name>|<n_comp>
	------------------------------------------------------------------------------
	notes:
	A = start of an input argument (material parameter) (double type assumed)
		It must have at least 1 field (the name), and at most 2 fields.
		The second field is optional and represents the default value 
		(use it in cases the input argument is optional)
	R = a response (result)
		It must have at 2+N fields, the id, the name, and N component,
	P = a parameter
		It must have 2 fields, the id and the name
	V = a variable
		It must have 3 fields, the id, the name, and the number of components
		
	fields in a record are separated by '|'
	------------------------------------------------------------------------------

When job = PF_MAT_SERIALIZE
	The "message" will be a custom string used to store in a text format
	the current status of the plugin material.

When job = PF_MAT_DESERIALIZE
	The "message" field will point to a string allocated in OpenSees, therefore
	the plugin can only use it to read from it for deserialziation.

\todo Add more documentation
*/
struct PluginMaterialData
{
	// Pointer to the external function.
	PluginMaterialProc proc;

	// Opaque pointer to custom data.
	// This is the only one that can must be allocated/released in the plugin library:
	//  - allocated in job = PF_MAT_INITIALIZE
	//  - deallocated in job = PF_MAT_FINALIZE
	// Typical usages: in a C++ plugin it can store the instance of a class, in C a plain structure...
	void* custom;

	// Description.
	// All data here should be set by the plugin function in job = PF_MAT_INITIALIZE
	// Note: the "message" field will store a pointer to a descriptive string in job = PF_MAT_INITIALIZE
	// while it will store a pointer to a 
	const char* message;
	int32_t tag;
	int32_t mat_type;
	int32_t n_param;
	
	// Output parameter/response/variable.
	// response_id and response set in the PluginFramework before the plugin
	// function call with job = PF_MAT_GET_<XXX>
	// note: response is not dynamically allocated but will point to existing
	// data of proper size
	int32_t response_id;
	uint64_t response_size;
	double* response;

	// for sensitivity
	double* sens_strain_gradient;
	int32_t sens_grad_index;
	int32_t sens_conditonal;
	int32_t sens_num_grad;
	int32_t sens_unused_placeholder; // for padding 8-bytes

	// Parameters.
	// Allocated only once by the PluginFramework after the plugin function
	// call with job = PF_MAT_INITIALIZE
	double* param;

	// Input.
	// Calculation input.
	// Not dynamically allocated, but set to point to existing data
	double* strain;
	double* strain_rate;
	double* temperature;
	double lch;
	double dT;
};

#endif
