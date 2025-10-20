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

#include <Domain.h>
#include "ResponseSpectrumAnalysis.h"
#include "DomainModalProperties.h"
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

ResponseSpectrumAnalysis::ResponseSpectrumAnalysis(
  DomainModalProperties& mp,
  TimeSeries* theFunction,
  const std::vector<double>& Tn,
  const std::vector<double>& Sa,
  int theDirection,
  double scale
)
  : mp(mp)
  , m_domain(mp.getDomain())
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

int
ResponseSpectrumAnalysis::analyze()
{
  // get the domain
  Domain* domain = m_domain;

  // size info
  int num_eigen = domain->getEigenvalues().Size();

  // check consistency
  int error_code;
  error_code = check();
  if (error_code < 0)
    return error_code;

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

int
ResponseSpectrumAnalysis::analyze(int mode_id)
{
  // get the domain
  Domain* domain = m_domain;

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

int
ResponseSpectrumAnalysis::check()
{
  // get the domain
  Domain* domain = m_domain;


  // number of eigen-modes
  int num_eigen = domain->getEigenvalues().Size();
  if (num_eigen < 1) {
    opserr << "ResponseSpectrumAnalysis::check() - No Eigenvalue provided.\n";
    return -1;
  }

  // check consistency
  {
    const Vector& ev = domain->getEigenvalues();
    if (ev.Size() != mp.eigenvalues().Size())
      return -2;
    double tol = std::max(1.0e-15, 1.0e-12 * ev.Norm());
    for (int i = 0; i < ev.Size(); ++i) {
      double a = ev(i);
      double b = mp.eigenvalues()(i);
      if (std::abs(a - b) > tol)
        return -2;
    }
    return 0;
  };

  return 0;
}

int
ResponseSpectrumAnalysis::beginMode()
{
  // new step... (do nothing)
  if (m_domain->analysisStep(0.0) < 0) {
    opserr << "ResponseSpectrumAnalysis::analyze() - the AnalysisModel failed"
      " at mode " << m_current_mode << "\n";
    return -1;
  }

  return 0;
}

int ResponseSpectrumAnalysis::endMode()
{
  // update domain
  if (m_domain->update() < 0) {
    opserr << "ResponseSpectrumAnalysis::analyze() - the AnalysisModel failed in updateDomain"
      " at mode " << m_current_mode << "\n";
    return -1;
  }
  // update Handler (cmp)

  // commit domain
  if (m_domain->commit() < 0) {
    opserr << "ResponseSpectrumAnalysis::analyze() - the AnalysisModel failed in commitDomain"
      " at mode " << m_current_mode << "\n";
    return -1;
  }

  return 0;
}

int
ResponseSpectrumAnalysis::solveMode()
{
  // get the domain
  Domain* domain = m_domain;


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
