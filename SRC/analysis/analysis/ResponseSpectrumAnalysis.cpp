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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#endif

#define DMP_ERR_INFO "( function: " << __func__ << ", file: \"" << __FILE__ << "\", line: " << __LINE__ << " )\n"
#define DMP_ERR(X) opserr << "FATAL ERROR: " << X << DMP_ERR_INFO, exit(-1)

#define DMP_DBL_LARGE 1.0e200

void
OPS_ResponseSpectrumAnalysis(void)
{
	// responseSpectrum $tsTag $dir $mcType <-scale $scale> <-damp $damp>

	// some kudos
	static bool first_done = false;
	if (!first_done) {
		opserr << "Using ResponseSpectrumAnalysis - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
		first_done = true;
	}

	// get analysis model
	AnalysisModel* theAnalysisModel = *OPS_GetAnalysisModel();
	if (theAnalysisModel == nullptr) {
		opserr << "modalProperties Error: no AnalysisModel available.\n";
		exit(-1);
	}
	if (theAnalysisModel->getDomainPtr() == nullptr) {
		opserr << "modalProperties Error: no Domain available.\n";
		exit(-1);
	}

	// init default arguments
	TimeSeries* ts = nullptr;
	int dir = 1;
	ResponseSpectrumAnalysis::ModalCombinationType mctype = ResponseSpectrumAnalysis::SRSS;
	double scale = 1.0;
	double damp = 0.05;

	// make sure eigenvalue and modal properties have been called before
	const auto& modal_props = theAnalysisModel->getDomainPtr()->getModalProperties();
	int ndf = modal_props.totalMass().Size();

	// parse
	int nargs = OPS_GetNumRemainingInputArgs();
	if (nargs < 3) {
		opserr << "responseSpectrum $tsTag $dir $mcType <-scale $scale> <-damp $damp>\n"
			"Error: at least 3 arguments should be provided.\n";
		exit(-1);
	}
	int numData = 1;
	int tstag;
	if (OPS_GetInt(&numData, &tstag) < 0) {
		opserr << "responseSpectrum Error: Failed to get timeSeries tag.\n";
		exit(-1);
	}
	ts = OPS_getTimeSeries(tstag);
	if (ts == nullptr) {
		opserr << "responseSpectrum Error: Failed to get timeSeries with tag = " << tstag << ".\n";
		exit(-1);
	}
	if (OPS_GetInt(&numData, &dir) < 0) {
		opserr << "responseSpectrum Error: Failed to get direction.\n";
		exit(-1);
	}
	if (dir < 1 || dir > ndf) {
		opserr << "responseSpectrum Error: provided direction (" << dir << ") should be in the range 1-" << ndf << ".\n";
		exit(-1);
	}
	const char* mc = OPS_GetString();
	if ((strcmp(mc, "srss") == 0) || (strcmp(mc, "SRSS") == 0)) {
		mctype = ResponseSpectrumAnalysis::SRSS;
	}
	else if ((strcmp(mc, "cqc") == 0) || (strcmp(mc, "CQC") == 0)) {
		mctype = ResponseSpectrumAnalysis::CQC;
	}
	else {
		opserr << "responseSpectrum Error: Invalid modal combination type (" << mc << ")\n"
			"It should be either 'srss' (or 'SRSS') or 'cqc' (or 'CQC').\n";
		exit(-1);
	}

	// parse optional data
	nargs = OPS_GetNumRemainingInputArgs();
	int loc = 0;
	while (loc < nargs) {
		const char* value = OPS_GetString();
		if (strcmp(value, "-scale") == 0) {
			if (loc < nargs - 1) {
				if (OPS_GetDouble(&numData, &scale) < 0) {
					opserr << "responseSpectrum Error: Failed to get scale factor.\n";
					exit(-1);
				}
				++loc;
			}
			else {
				opserr << "responseSpectrum Error: scale factor requested but not provided.\n";
				exit(-1);
			}
		}
		else if (strcmp(value, "-damp") == 0) {
			if (loc < nargs - 1) {
				if (OPS_GetDouble(&numData, &damp) < 0) {
					opserr << "responseSpectrum Error: Failed to get damping factor.\n";
					exit(-1);
				}
				++loc;
			}
			else {
				opserr << "responseSpectrum Error: damping factor requested but not provided.\n";
				exit(-1);
			}
		}
		++loc;
	}

	// ok, create the response spectrum analysis and run it here... 
	// no need to store it
	ResponseSpectrumAnalysis rsa(theAnalysisModel, ts, dir, mctype, scale, damp);
	rsa.analyze();

}

namespace {

	// performs the CQC modal combination in the general case
	// of different damping ratios for each mode,
	// and eventually different scale factors for each mode
	double cqc(const Vector& mu, const Vector& eval, const Vector& dmp, const Vector& scalf) {
		double u = 0.0;
		int ne = mu.Size();
		for (int i = 0; i < ne; ++i) {
			for (int j = 0; j < ne; ++j) {
				double di = dmp(i);
				double dj = dmp(j);
				double bij = eval(i) / eval(j);
				double ro = (8.0 * std::sqrt(di * dj) * (di + bij * dj) * std::pow(bij, 3.0 / 2.0)) /
					(std::pow(1.0 - std::pow(bij, 2), 2) + 4.0 * di * dj * bij * (1.0 + std::pow(bij, 2)) +
						4.0 * (std::pow(di, 2) + std::pow(dj, 2)) * std::pow(bij, 2));
				u += ((scalf(i) * mu(i)) * (scalf(j) * mu(j)) * ro);
			}
		}
		u = std::sqrt(u);
		return u;
	}

	// performs the SRSS modal combination
	double srss(const Vector& mu, const Vector& scalf) {
		double u = 0.0;
		int ne = mu.Size();
		for (int i = 0; i < ne; ++i) {
			u += std::pow(scalf(i) * mu(i), 2);
		}
		u = std::sqrt(u);
		return u;
	}
}

ResponseSpectrumAnalysis::ResponseSpectrumAnalysis(
	AnalysisModel* theModel,
	TimeSeries* theFunction,
	int theDirection,
	ModalCombinationType theMCType,
	double scale,
	double damp
)
	: m_model(theModel)
	, m_function(theFunction)
	, m_direction(theDirection)
	, m_mc_type(theMCType)
	, m_scale(scale)
	, m_damp(damp)
{

}

ResponseSpectrumAnalysis::~ResponseSpectrumAnalysis()
{
}

void ResponseSpectrumAnalysis::analyze()
{
	// mimic what any other analysis would do
	// note: RS is a single-step analysis

	// get the domain
	Domain* domain = m_model->getDomainPtr();

	// result flag
	int result;

	// new step... (do nothing)
	result = m_model->analysisStep();
	if (result < 0) {
		DMP_ERR("ResponseSpectrumAnalysis::analyze() - the AnalysisModel failed"
			<< " with domain at load factor "
			<< domain->getCurrentTime() << "\n");
	}

	// compute displacement
	computeRSdisplacement();

	// update domain
	result = m_model->updateDomain();
	if (result < 0) {
		DMP_ERR("ResponseSpectrumAnalysis::analyze() - the AnalysisModel in updateDomain"
			<< " with domain at load factor "
			<< domain->getCurrentTime() << "\n");
	}

	// commit domain
	result = m_model->commitDomain();
	if (result < 0) {
		DMP_ERR("ResponseSpectrumAnalysis::analyze() - the AnalysisModel in commitDomain"
			<< " with domain at load factor "
			<< domain->getCurrentTime() << "\n");
	}
}

void ResponseSpectrumAnalysis::computeRSdisplacement()
{
	// get the domain
	Domain* domain = m_model->getDomainPtr();

	// get the modal properties
	const DomainModalProperties& mp = domain->getModalProperties();

	// number of eigen-modes
	int num_eigen = domain->getEigenvalues().Size();
	if (num_eigen < 1)
		DMP_ERR("No Eigenvalue provided.\n");

	// size info
	int ndm = mp.centerOfMass().Size();
	int ndf = mp.totalMass().Size();

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
	if (!check_eigen())
		DMP_ERR("Eigenvalues stored in DomainModalProperties are not equal to the eigenvalues in the model.\n"
			"Make sure to call the 'modalProperties' command\n"
			"after the 'eigen' command, and right before the 'responseSpectrum' command.\n");

	// compute accelerations for each mode using the 
	// provided response spectrum function (time series)
	Vector mga(num_eigen);
	for (int i = 0; i < num_eigen; ++i) {
		double lambda = mp.eigenvalues()(i);
		double omega = std::sqrt(lambda);
		double freq = omega / 2.0 / M_PI;
		double period = 1.0 / freq;
		mga(i) = m_function->getFactor(period);
	}

	// some auxiliary data:
	// -------------------
	// a vector that stores, for each node, and for each DOF, the displacement
	// value for each mode
	Vector u_modal(num_eigen);
	// the nodal eigenvector
	Matrix V(ndf, num_eigen);
	// now we have only 1 constant damping. use this vector for future implementation
	// of variable damping
	Vector damp(num_eigen);
	for (int i = 0; i < num_eigen; ++i)
		damp(i) = m_damp;
	// now we don't have scale factors for each mode. use this vector for future implementation
	Vector scalf(num_eigen);
	for (int i = 0; i < num_eigen; ++i)
		scalf(i) = 1.0;

	// excited DOF.
	// Note: now we assume that a RS acts along one of the global directions, so we need
	// to consider only the associated column of the modal participation factors.
	// in future versions we can implement a general direction vector.
	int exdof = m_direction - 1; // make it 0-based

	// for each node
	Node* node;
	NodeIter& theNodes = domain->getNodes();
	while ((node = theNodes()) != 0) {

		// store the nodal eigenvector, according to the ndf of modal properties
		// and its scaling
		V.Zero();
		const Matrix& node_evec = node->getEigenvectors();
		int node_ndf = node_evec.noRows();
		for (int i = 0; i < std::min(node_ndf, ndf); ++i) {
			// exclude any dof >= node_ndf (if ndf > node_ndf... in case node is U-only)
			// exclude any dof >= ndf (if node_ndf > ndf... in case the node has more DOFs than U and R)
			// we also need to exclude pressure dofs... easy in 3D because node_ndf is 4,
			// but in 2D it is 3, as in U-R...
			if (ndf == 6 && node_ndf == 4 && i == 3)
				continue;
			for (int j = 0; j < num_eigen; ++j) {
				double scale = mp.eigenVectorScaleFactors()(j);
				V(i, j) = node_evec(i, j) * scale;
			}
		}

		// for each DOF...
		for (int i = 0; i < ndf; ++i) {

			// compute modal displacements for the i-th DOF
			for (int j = 0; j < num_eigen; ++j) {
				u_modal(j) = V(i, j) * mp.modalParticipationFactors()(j, exdof) * mga(j) / mp.eigenvalues()(j);
			}

			// modal combination for the i-th DOF
			double u = 0.0;
			switch (m_mc_type)
			{
			case ResponseSpectrumAnalysis::SRSS:
				u = srss(u_modal, scalf);
				break;
			case ResponseSpectrumAnalysis::CQC:
				u = cqc(u_modal, mp.eigenvalues(), damp, scalf);
				break;
			default:
				DMP_ERR("Error: modal combination type not supported.\n");
				break;
			}

			// save this displacement at the i-th dof as new trial displacement
			if(i < node_ndf)
				node->setTrialDisp(u, i);
		}

	}
}
