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
                                                                        
// Original implementation: José Abell (UANDES), Massimo Petracca (ASDEA)
//
// ASDPlasticMaterial
//
// Fully general templated material class for plasticity modeling

#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include "AllASDPlasticMaterials.h"


void *OPS_AllASDPlasticMaterials(void)
{
	// some kudos
	static bool first_done = false;
	if (!first_done) {
		opserr << "Using ASDPlasticMaterial - Developed by: Jose Abell (UANDES), Massimo Petracca and Guido Camata (ASDEA Software Technology)\n";
		first_done = true;
	}

	// check arguments
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 1) {
		opserr << 
			"nDMaterial ASDPlasticMaterial Error: Few arguments (< 3).\n"
			"nDMaterial ASDPlasticMaterial $tag $E $v "
			"-Te $Te -Ts $Ts <-Td $Td> -Ce $Ce -Cs $Cs <-Cd $Cd> "
			"<-rho $rho> "
			"<-implex> <-implexControl $implexErrorTolerance $implexTimeReductionLimit> <-implexAlpha $alpha>"
			"<-crackPlanes $nct $ncc $smoothingAngle>"
			"<-eta $eta> <-tangent> <-autoRegularization>\n";
		return nullptr;
	}

	
	int tag_in = 1;
	double k0_in = 1;
	double H_alpha = 1;
	double H_k = 1;
	double E = 1;
	double nu = 0.0;
	double rho_ = 1000;

	NDMaterial* instance = new VonMisesLinearHardening(tag_in, k0_in, H_alpha, H_k, E, nu, rho_);

	return instance;

}

// 	// numData
// 	int numData = 1;

// 	// data
// 	int tag;
// 	double E;
// 	double v;
// 	double rho = 0.0;
// 	bool implex = false;
// 	bool implex_control = false;
// 	double implex_error_tolerance = 0.05;
// 	double implex_time_redution_limit = 0.01;
// 	double implex_alpha = 1.0;
// 	double eta = 0.0;
// 	bool tangent = false;
// 	bool auto_regularization = false;
// 	std::vector<double> Te, Ts, Td, Ce, Cs, Cd;
// 	int nct = 0;
// 	int ncc = 0;
// 	double smoothing_angle = 45.0;
	
// 	// get tag
// 	if (OPS_GetInt(&numData, &tag) != 0)  {
// 		opserr << "nDMaterial ASDConcrete3D Error: invalid 'tag'.\n";
// 		return nullptr;
// 	}

// 	// get Elasticity arguments
// 	if (OPS_GetDouble(&numData, &E) != 0) {
// 		opserr << "nDMaterial ASDConcrete3D Error: invalid 'E'.\n";
// 		return nullptr;
// 	}
// 	if (E <= 0.0) {
// 		opserr << "nDMaterial ASDConcrete3D Error: invalid value for 'E' (" << E << "). It should be strictly positive.\n";
// 		return nullptr;
// 	}
// 	if (OPS_GetDouble(&numData, &v) != 0) {
// 		opserr << "nDMaterial ASDConcrete3D Error: invalid 'v'.\n";
// 		return nullptr;
// 	}

// 	// utilities (code re-use)
// 	auto lam_optional_int = [&numData](const char* variable, int& value) -> bool {
// 		if (OPS_GetNumRemainingInputArgs() > 0) {
// 			if (OPS_GetInt(&numData, &value) < 0) {
// 				opserr << "nDMaterial ASDConcrete3D Error: failed to get '" << variable << "'.\n";
// 				return false;
// 			}
// 		}
// 		else {
// 			opserr << "nDMaterial ASDConcrete3D Error: '" << variable << "' requested but not provided.\n";
// 			return false;
// 		}
// 		return true;
// 	};
// 	auto lam_optional_double = [&numData](const char* variable, double& value) -> bool {
// 		if (OPS_GetNumRemainingInputArgs() > 0) {
// 			if (OPS_GetDouble(&numData, &value) < 0) {
// 				opserr << "nDMaterial ASDConcrete3D Error: failed to get '" << variable << "'.\n";
// 				return false;
// 			}
// 		}
// 		else {
// 			opserr << "nDMaterial ASDConcrete3D Error: '" << variable << "' requested but not provided.\n";
// 			return false;
// 		}
// 		return true;
// 	};
// 	auto lam_optional_list = [&numData](const char* variable, std::vector<double>& value) -> bool {
// 		// first try expanded list like {*}$the_list,
// 		// also used in python like *the_list
// 		value.clear();
// 		while (OPS_GetNumRemainingInputArgs() > 0) {
// 			double item;
// 			auto old_num_rem = OPS_GetNumRemainingInputArgs();
// 			if (OPS_GetDoubleInput(&numData, &item) < 0) {
// 				auto new_num_rem = OPS_GetNumRemainingInputArgs();
// 				if (new_num_rem < old_num_rem)
// 					OPS_ResetCurrentInputArg(-1);
// 				break;
// 			}
// 			value.push_back(item);
// 		}
// 		// try Tcl list (it's a string after all...)
// 		if (value.size() == 0 && OPS_GetNumRemainingInputArgs() > 0) {
// 			std::string list_string = OPS_GetString();
// 			if (!string_to_list_of_doubles(list_string, ' ', value)) {
// 				opserr << "nDMaterial ASDConcrete3D Error: cannot parse the '" << variable << "' list.\n";
// 				return false;
// 			}
// 		}
// 		return true;
// 	};

// 	// optional parameters
// 	while (OPS_GetNumRemainingInputArgs() > 0) {
// 		const char* value = OPS_GetString();
// 		if (strcmp(value, "-rho") == 0) {
// 			if (!lam_optional_double("rho", rho))
// 				return nullptr;
// 		}
// 		else if (strcmp(value, "-implex") == 0) {
// 			implex = true;
// 		}
// 		else if (strcmp(value, "-implexControl") == 0) {
// 			implex_control = true;
// 			if (OPS_GetNumRemainingInputArgs() < 2) {
// 				opserr << "nDMaterial ASDConcrete3D Error: '-implexControl' given without the next 2 arguments $implexErrorTolerance $implexTimeReductionLimit.\n";
// 				return nullptr;
// 			}
// 			if (!lam_optional_double("implexErrorTolerance", implex_error_tolerance))
// 				return nullptr;
// 			if (!lam_optional_double("implexTimeReductionLimit", implex_time_redution_limit))
// 				return nullptr;
// 		}
// 		else if (strcmp(value, "-implexAlpha") == 0) {
// 			if (!lam_optional_double("alpha", implex_alpha))
// 				return nullptr;
// 		}
// 		else if (strcmp(value, "-eta") == 0) {
// 			if (!lam_optional_double("eta", eta))
// 				return nullptr;
// 		}
// 		else if (strcmp(value, "-tangent") == 0) {
// 			tangent = true;
// 		}
// 		else if (strcmp(value, "-autoRegularization") == 0) {
// 			auto_regularization = true;
// 		}
// 		else if (strcmp(value, "-Te") == 0) {
// 			if (!lam_optional_list("Te", Te))
// 				return nullptr;
// 		}
// 		else if (strcmp(value, "-Ts") == 0) {
// 			if (!lam_optional_list("Ts", Ts))
// 				return nullptr;
// 		}
// 		else if (strcmp(value, "-Td") == 0) {
// 			if (!lam_optional_list("Td", Td))
// 				return nullptr;
// 		}
// 		else if (strcmp(value, "-Ce") == 0) {
// 			if (!lam_optional_list("Ce", Ce))
// 				return nullptr;
// 		}
// 		else if (strcmp(value, "-Cs") == 0) {
// 			if (!lam_optional_list("Cs", Cs))
// 				return nullptr;
// 		}
// 		else if (strcmp(value, "-Cd") == 0) {
// 			if (!lam_optional_list("Cd", Cd))
// 				return nullptr;
// 		}
// 		else if (strcmp(value, "-crackPlanes") == 0) {
// 			if (OPS_GetNumRemainingInputArgs() < 3) {
// 				opserr << "nDMaterial ASDConcrete3D Error: '-crackPlanes' given without the next 3 arguments $nct $ncc and $smoothingAngle.\n";
// 				return nullptr;
// 			}
// 			if (!lam_optional_int("nct", nct))
// 				return nullptr;
// 			if (!lam_optional_int("ncc", ncc))
// 				return nullptr;
// 			if (!lam_optional_double("$smoothingAngle", smoothing_angle))
// 				return nullptr;
// 		}
// 	}

// 	// check lists
// 	if (Te.size() < 1) {
// 		opserr << "nDMaterial ASDConcrete3D Error: 'Te' list is empty. At least 1 non-zero value should be provided.\n";
// 		return nullptr;
// 	}
// 	if (Ts.size() != Te.size()) {
// 		opserr << "nDMaterial ASDConcrete3D Error: 'Te' (size = " <<
// 			static_cast<int>(Te.size()) << ") and 'Ts' (size = " <<
// 			static_cast<int>(Ts.size()) << ") lists should have the same size.\n";
// 		return nullptr;
// 	}
// 	if (Td.size() == 0) {
// 		Td.resize(Te.size(), 0.0);
// 	}
// 	else if (Td.size() != Te.size()) {
// 		opserr << "nDMaterial ASDConcrete3D Error: 'Te' (size = " <<
// 			static_cast<int>(Te.size()) << ") and 'Td' (size = " <<
// 			static_cast<int>(Td.size()) << ") lists should have the same size.\n";
// 		return nullptr;
// 	}
// 	if (Ce.size() < 1) {
// 		opserr << "nDMaterial ASDConcrete3D Error: 'Tc' list is empty. At least 1 non-zero value should be provided.\n";
// 		return nullptr;
// 	}
// 	if (Cs.size() != Ce.size()) {
// 		opserr << "nDMaterial ASDConcrete3D Error: 'Ce' (size = " <<
// 			static_cast<int>(Ce.size()) << ") and 'Cs' (size = " <<
// 			static_cast<int>(Cs.size()) << ") lists should have the same size.\n";
// 		return nullptr;
// 	}
// 	if (Cd.size() == 0) {
// 		Cd.resize(Ce.size(), 0.0);
// 	}
// 	else if (Cd.size() != Ce.size()) {
// 		opserr << "nDMaterial ASDConcrete3D Error: 'Ce' (size = " <<
// 			static_cast<int>(Ce.size()) << ") and 'Cd' (size = " <<
// 			static_cast<int>(Cd.size()) << ") lists should have the same size.\n";
// 		return nullptr;
// 	}

// 	// build the hardening laws
// 	ASDConcrete3DMaterial::HardeningLaw HT(tag, ASDConcrete3DMaterial::HardeningLawType::Tension, E, Te, Ts, Td);
// 	if (!HT.isValid()) {
// 		opserr << "nDMaterial ASDConcrete3D Error: Tensile hardening law is not valid.\n";
// 		return nullptr;
// 	}
// 	ASDConcrete3DMaterial::HardeningLaw HC(tag, ASDConcrete3DMaterial::HardeningLawType::Compression, E, Ce, Cs, Cd);
// 	if (!HC.isValid()) {
// 		opserr << "nDMaterial ASDConcrete3D Error: Compressive hardening law is not valid.\n";
// 		return nullptr;
// 	}

// 	// create the material
// 	NDMaterial* instance = new ASDConcrete3DMaterial(
// 		tag, 
// 		E, v, rho, eta, 
// 		implex, implex_control, implex_error_tolerance, implex_time_redution_limit, implex_alpha,
// 		tangent, auto_regularization,
// 		HT, HC,
// 		nct, ncc, smoothing_angle);
// 	if (instance == nullptr) {
// 		opserr << "nDMaterial ASDConcrete3D Error: failed to allocate a new material.\n";
// 		return nullptr;
// 	}
// 	return instance;
// }