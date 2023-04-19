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

// $Revision: 1.0 $
// $Date: 2020-11-27 18:07:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/ResponseSpectrumAnalysis.h,v $

// Written: Massimo Petracca (ASDEA Software) 
// Created: Fri Nov 27 18:07:11: 2020
// Revision: A
//
// Description: This file contains the class definition for ResponseSpectrumAnalysis.
//
// What: "@(#) ResponseSpectrumAnalysis.h, revA"

#include <ResponseSpectrumAnalysis.h>
#include <AnalysisModel.h>
#include <Domain.h>
#include <DomainModalProperties.h>
#include <TimeSeries.h>
#include <elementAPI.h>
#include <Node.h>
#include <NodeIter.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#endif

namespace {

	bool string_to_list_of_doubles(const std::string& text, char sep, std::vector<double>& out) {
		auto to_double = [](const std::string& text, double& num) -> bool {
			num = 0.0;
			try {
				num = std::stod(text);
				return true;
			}
			catch (...) {
				return false;
			}
		};
		if (out.size() > 0) out.clear();
		std::size_t start = 0, end = 0;
		double value;
		while (true) {
			end = text.find(sep, start);
			if (end == std::string::npos) {
				if (start < text.size()) {
					if (!to_double(text.substr(start), value))
						return false;
					out.push_back(value);
				}
				break;
			}
			std::string subs = text.substr(start, end - start);
			if (subs.size() > 0) {
				if (!to_double(subs, value))
					return false;
				out.push_back(value);
			}
			start = end + 1;
		}
		return true;
	}

}

int
OPS_ResponseSpectrumAnalysis(void)
{
	// responseSpectrum $tsTag $dir <-scale $scale>

	// some kudos
	static bool first_done = false;
	if (!first_done) {
		opserr << "Using ResponseSpectrumAnalysis - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
		first_done = true;
	}

	// get analysis model
	AnalysisModel* theAnalysisModel = *OPS_GetAnalysisModel();
	if (theAnalysisModel == nullptr) {
		opserr << "ResponseSpectrumAnalysis Error: no AnalysisModel available.\n";
		return -1;
	}
	if (theAnalysisModel->getDomainPtr() == nullptr) {
		opserr << "ResponseSpectrumAnalysis Error: no Domain available.\n";
		return -1;
	}

	// init default arguments
	TimeSeries* ts = nullptr;
	int dir = 1;
	double scale = 1.0;
	std::vector<double> Tn;
	std::vector<double> Sa;
	int mode_id = 0;
	bool single_mode = false;

	// make sure eigenvalue and modal properties have been called before
	DomainModalProperties modal_props;
	if (theAnalysisModel->getDomainPtr()->getModalProperties(modal_props) < 0) {
		opserr << "ResponseSpectrumAnalysis - eigen and modalProperties have not been called" << endln;
		return -1;
	}
	int ndf = modal_props.totalMass().Size();

	// parse
	int nargs = OPS_GetNumRemainingInputArgs();
	if (nargs < 2) {
		opserr << "ResponseSpectrumAnalysis $tsTag $dir <-scale $scale> <-damp $damp>\n"
			<< "or\n"
			<< "ResponseSpectrumAnalysis $dir -Tn $TnValues -fn $fnValues -Sa $SaValues <-scale $scale> <-damp $damp>\n"
			"Error: at least 2 arguments should be provided.\n";
		return -1;
	}

	// search for -Tn -Sa, if both found we use them as lists
	// otherwise we fallback to the old implementation of the timeSeries
	bool found_Tn = false;
	bool found_Sa = false;
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* arg = OPS_GetString();
		if (!found_Tn && (strcmp(arg, "-Tn") == 0 || strcmp(arg, "-fn") == 0))
			found_Tn = true;
		if (!found_Sa && strcmp(arg, "-Sa") == 0)
			found_Sa = true;
	}
	bool use_lists = found_Tn && found_Sa;
	if (found_Tn && !found_Sa) {
		opserr << "ResponseSpectrumAnalysis Error: found -Tn or -fn but not -Sa, please specify both of them or use a timeSeries\n";
		return -1;
	}
	if (found_Sa && !found_Tn) {
		opserr << "ResponseSpectrumAnalysis Error: found -Sa but not -Tn or -fn, please specify both of them or use a timeSeries\n";
		return -1;
	}
	OPS_ResetCurrentInputArg(-nargs);

	// num data
	int numData = 1;

	// get time series
	if (!use_lists) {
		int tstag;
		if (OPS_GetInt(&numData, &tstag) < 0) {
			opserr << "ResponseSpectrumAnalysis Error: Failed to get timeSeries tag.\n";
			return -1;
		}
		ts = OPS_getTimeSeries(tstag);
		if (ts == nullptr) {
			opserr << "ResponseSpectrumAnalysis Error: Failed to get timeSeries with tag = " << tstag << ".\n";
			return -1;
		}
	}

	// get direction
	if (OPS_GetInt(&numData, &dir) < 0) {
		opserr << "ResponseSpectrumAnalysis Error: Failed to get direction.\n";
		return -1;
	}
	if (dir < 1 || dir > ndf) {
		opserr << "ResponseSpectrumAnalysis Error: provided direction (" << dir << ") should be in the range 1-" << ndf << ".\n";
		return -1;
	}

	// parse optional data
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* value = OPS_GetString();
		if (strcmp(value, "-scale") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetDouble(&numData, &scale) < 0) {
					opserr << "ResponseSpectrumAnalysis Error: Failed to get scale factor.\n";
					return -1;
				}
			}
			else {
				opserr << "ResponseSpectrumAnalysis Error: scale factor requested but not provided.\n";
				return -1;
			}
		}
		else if (strcmp(value, "-mode") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 0) {
				if (OPS_GetInt(&numData, &mode_id) < 0) {
					opserr << "ResponseSpectrumAnalysis Error: Failed to get the mode_id.\n";
					return -1;
				}
				--mode_id; // make it 0-based
				single_mode = true;
			}
			else {
				opserr << "ResponseSpectrumAnalysis Error: mode_id requested but not provided.\n";
				return -1;
			}
		}
		else if (strcmp(value, "-Tn") == 0 || strcmp(value, "-fn") == 0) {
			// first try expanded list like {*}$the_list,
			// also used in python like *the_list
			Tn.clear();
			while (OPS_GetNumRemainingInputArgs() > 0) {
				double item;
				auto old_num_rem = OPS_GetNumRemainingInputArgs();
				if (OPS_GetDoubleInput(&numData, &item) < 0) {
					auto new_num_rem = OPS_GetNumRemainingInputArgs();
					if (new_num_rem < old_num_rem)
						OPS_ResetCurrentInputArg(-1);
					break;
				}
				if (strcmp(value, "-fn") == 0 && item != 0) {
					item = 1.0 / item;
				}
				Tn.push_back(item);
			}
			// try Tcl list (it's a string after all...)
			if (Tn.size() == 0 && OPS_GetNumRemainingInputArgs() > 0) {
				std::string list_string = OPS_GetString();
				if (!string_to_list_of_doubles(list_string, ' ', Tn)) {
					opserr << "ResponseSpectrumAnalysis Error: cannot parse the Tn list.\n";
					return -1;
				}
			}
		}
		else if (strcmp(value, "-Sa") == 0) {
			// first try expanded list like {*}$the_list,
			// also used in python like *the_list
			Sa.clear();
			while (OPS_GetNumRemainingInputArgs() > 0) {
				double item;
				auto old_num_rem = OPS_GetNumRemainingInputArgs();
				if (OPS_GetDoubleInput(&numData, &item) < 0) {
					auto new_num_rem = OPS_GetNumRemainingInputArgs();
					if (new_num_rem < old_num_rem)
						OPS_ResetCurrentInputArg(-1);
					break;
				}
				Sa.push_back(item);
			}
			// try Tcl list (it's a string after all...)
			if (Sa.size() == 0 && OPS_GetNumRemainingInputArgs() > 0) {
				std::string list_string = OPS_GetString();
				if (!string_to_list_of_doubles(list_string, ' ', Sa)) {
					opserr << "ResponseSpectrumAnalysis Error: cannot parse the Sa list.\n";
					return -1;
				}
			}
		}
	}

	// check Tn and Sa vectors
	if (use_lists) {
		if (Tn.size() != Sa.size()) {
			opserr << "ResponseSpectrumAnalysis Error: Sa and Tn lists must have the same length\n";
			opserr << (int)Tn.size() << " != " << (int)Sa.size() << "\n";
			return -1;
		}
		if (Tn.size() == 0) {
			opserr << "ResponseSpectrumAnalysis Error: Sa and Tn lists cannot be empty\n";
			return -1;
		}
		for (std::size_t i = 0; i < Tn.size(); ++i) {
			if (Tn[i] < 0.0) {
				opserr << "ResponseSpectrumAnalysis Error: Tn values must be positive (found " << Tn[i] << ")\n";
				return -1;
			}
			if (i > 0 && Tn[i] <= Tn[i - 1]) {
				opserr << "ResponseSpectrumAnalysis Error: Tn values must be monotonically increasing (found " << Tn[i] << " after " << Tn[i - 1] << ")\n";
				return -1;
			}
		}
		for (double item : Sa) {
			if (item < 0.0) {
				opserr << "ResponseSpectrumAnalysis Error: Sa values must be positive (found " << item << ")\n";
				return -1;
			}
		}
	}

	// ok, create the response spectrum analysis and run it here... 
	// no need to store it
	ResponseSpectrumAnalysis rsa(theAnalysisModel, ts, Tn, Sa, dir, scale);
	int result;
	if (single_mode)
		result = rsa.analyze(mode_id);
	else
		result = rsa.analyze();

	return result;
}

ResponseSpectrumAnalysis::ResponseSpectrumAnalysis(
	AnalysisModel* theModel,
	TimeSeries* theFunction,
	const std::vector<double>& Tn,
	const std::vector<double>& Sa,
	int theDirection,
	double scale
)
	: m_model(theModel)
	, m_function(theFunction)
	, m_Tn(Tn)
	, m_Sa(Sa)
	, m_direction(theDirection)
	, m_scale(scale)
	, m_current_mode(0)
{

}

ResponseSpectrumAnalysis::~ResponseSpectrumAnalysis()
{
}

int ResponseSpectrumAnalysis::analyze()
{
	// get the domain
	Domain* domain = m_model->getDomainPtr();

	// get the modal properties
	DomainModalProperties mp;
	if (domain->getModalProperties(mp) < 0) {
		opserr << "ResponseSpectrumAnalysis::analyze() - failed to get modal properties" << endln;
		return -1;
	}

	// size info
	int num_eigen = domain->getEigenvalues().Size();

	// check consistency
	int error_code;
	error_code = check();
	if (error_code < 0) return error_code;

	// loop over all required eigen-modes, compute the modal displacement
	// and save the results.
	// we just compute modal displacements without doing any (SRSS, CQC, etc..)
	// modal combination otherwise derived results cannot be computed.
	// for each mode, this analysis produces a new analysis step.
	// modal combination of displacements (or any derived results)
	// it's up to the user.
	for (m_current_mode = 0; m_current_mode < num_eigen; ++m_current_mode)
	{
		// init the new step
		error_code = beginMode();
		if (error_code < 0) return error_code;

		// compute modal acceleration for this mode using the 
		// provided response spectrum function (time series)
		error_code = solveMode();
		if (error_code < 0) return error_code;

		// done with the current step.
		// here the modal displacements will be recorded (and all other results
		// if requested via recorders...)
		error_code = endMode();
		if (error_code < 0) return error_code;
	}

	return 0;
}

int ResponseSpectrumAnalysis::analyze(int mode_id)
{
	// get the domain
	Domain* domain = m_model->getDomainPtr();

	// get the modal properties
	//const DomainModalProperties& mp = domain->getModalProperties();
	DomainModalProperties mp;
	if (domain->getModalProperties(mp) < 0) {
		opserr << "ResponseSpectrumAnalysis::analyze(mode_id) - failed to get modal properties" << endln;
		return -1;
	}

	// size info
	int num_eigen = domain->getEigenvalues().Size();
	if (mode_id < 0 || mode_id >= num_eigen) {
		opserr << "ResponseSpectrumAnalysis::analyze(mode_id) - The provided mode_id (" << mode_id + 1 << ") is out of range (1, " << num_eigen << ")\n";
		return -1;
	}
	m_current_mode = mode_id;

	// check consistency
	int error_code;
	error_code = check();
	if (error_code < 0) return error_code;

	// init the new step
	error_code = beginMode();
	if (error_code < 0) return error_code;

	// compute modal acceleration for this mode using the 
	// provided response spectrum function (time series)
	error_code = solveMode();
	if (error_code < 0) return error_code;

	// done with the current step.
	// here the modal displacements will be recorded (and all other results
	// if requested via recorders...)
	error_code = endMode();
	if (error_code < 0) return error_code;

	return 0;
}

int ResponseSpectrumAnalysis::check()
{
	// get the domain
	Domain* domain = m_model->getDomainPtr();

	// get the modal properties
	//const DomainModalProperties& mp = domain->getModalProperties();
	DomainModalProperties mp;
	if (domain->getModalProperties(mp) < 0) {
		opserr << "ResponseSpectrumAnalysis::check() - failed to get modal properties" << endln;
		return -1;
	}

	// number of eigen-modes
	int num_eigen = domain->getEigenvalues().Size();
	if (num_eigen < 1) {
		opserr << "ResponseSpectrumAnalysis::check() - No Eigenvalue provided.\n";
		return -1;
	}

	// check consistency
	auto check_eigen = [&mp, domain]() -> bool {
		const Vector& ev = domain->getEigenvalues();
		if (ev.Size() != mp.eigenvalues().Size())
			return false;
		double tol = std::max(1.0e-15, 1.0e-12 * ev.Norm());
		for (int i = 0; i < ev.Size(); ++i) {
			double a = ev(i);
			double b = mp.eigenvalues()(i);
			if (std::abs(a - b) > tol)
				return false;
		}
		return true;
	};
	if (!check_eigen()) {
		opserr << "ResponseSpectrumAnalysis::check() - Eigenvalues stored in DomainModalProperties are not equal to the eigenvalues in the model.\n"
			"Make sure to call the 'modalProperties' command\n"
			"after the 'eigen' command, and right before the 'responseSpectrum' command.\n";
		return -1;
	}

	return 0;
}

int ResponseSpectrumAnalysis::beginMode()
{
	// new step... (do nothing)
	if (m_model->analysisStep() < 0) {
		opserr << "ResponseSpectrumAnalysis::analyze() - the AnalysisModel failed"
			" at mode " << m_current_mode << "\n";
		return -1;
	}

	return 0;
}

int ResponseSpectrumAnalysis::endMode()
{
	// update domain
	if (m_model->updateDomain() < 0) {
		opserr << "ResponseSpectrumAnalysis::analyze() - the AnalysisModel failed in updateDomain"
			" at mode " << m_current_mode << "\n";
		return -1;
	}

	// commit domain
	if (m_model->commitDomain() < 0) {
		opserr << "ResponseSpectrumAnalysis::analyze() - the AnalysisModel failed in commitDomain"
			" at mode " << m_current_mode << "\n";
		return -1;
	}

	return 0;
}

int ResponseSpectrumAnalysis::solveMode()
{
	// get the domain
	Domain* domain = m_model->getDomainPtr();

	// get the modal properties
	//const DomainModalProperties& mp = domain->getModalProperties();
	DomainModalProperties mp;
	if (domain->getModalProperties(mp) < 0) {
		opserr << "ResponseSpectrumAnalysis::solveMode() - failed to get modal properties" << endln;
		return -1;
	}

	// size info
	int ndf = mp.totalMass().Size();

	// excited DOF.
	// Note: now we assume that a RS acts along one of the global directions, so we need
	// to consider only the associated column of the modal participation factors.
	// in future versions we can implement a general direction vector.
	int exdof = m_direction - 1; // make it 0-based

	// compute modal acceleration for this mode using the 
	// provided response spectrum function (time series)
	double lambda = mp.eigenvalues()(m_current_mode);
	double omega = std::sqrt(lambda);
	double freq = omega / 2.0 / M_PI;
	double period = 1.0 / freq;
	double mga = getSa(period);
	double Vscale = mp.eigenVectorScaleFactors()(m_current_mode);
	double MPF = mp.modalParticipationFactors()(m_current_mode, exdof);

	// loop over all nodes and compute the modal displacements
	Node* node;
	NodeIter& theNodes = domain->getNodes();
	while ((node = theNodes()) != 0) {

		// get the nodal eigenvector, according to the ndf of modal properties
		// and its scaling
		const Matrix& node_evec = node->getEigenvectors();
		int node_ndf = node_evec.noRows();

		// for each DOF...
		for (int i = 0; i < std::min(node_ndf, ndf); ++i) {
			// exclude any dof >= node_ndf (if ndf > node_ndf... in case node is U-only)
			// exclude any dof >= ndf (if node_ndf > ndf... in case the node has more DOFs than U and R)

			// we also need to exclude pressure dofs... easy in 3D because node_ndf is 4,
			// but in 2D it is 3, as in U-R...
			if (ndf == 6 && node_ndf == 4 && i == 3)
				continue;

			// eigenvector
			double V = node_evec(i, m_current_mode) * Vscale;

			// compute modal displacements for the i-th DOF
			double u_modal = V * MPF * mga / lambda;

			// save this displacement at the i-th dof as new trial displacement
			node->setTrialDisp(u_modal, i);
		}

	}

	return 0;
}

double ResponseSpectrumAnalysis::getSa(double T) const
{
	// use the time series if provided
	if (m_function)
		return m_function->getFactor(T);
	// otherwise use the vectors
	std::size_t n = m_Tn.size();
	// check empty
	if (n < 1)
		return 0.0;
	// check constant
	if (n == 1)
		return m_Sa[0];
	// check bounds
	if (T <= m_Tn.front())
		return m_Sa.front();
	if (T >= m_Tn.back())
		return m_Sa.back();
	// interpolate
	for (std::size_t i = 1; i < n; ++i) {
		double t1 = m_Tn[i];
		if (T <= t1) {
			double t0 = m_Tn[i - 1];
			double s1 = m_Sa[i];
			double s0 = m_Sa[i - 1];
			double dT = t1 - t0;
			double dS = s1 - s0;
			double scale = dT > 0.0 ? (T - t0) / dT : 0.5;
			return s0 + scale * dS;
		}
	}
	return m_Sa.back();
}

