/* *********************************************************************
**    Module:	QzSimple1.cpp 
**
**    Purpose:	Provide a simple Q-z material for OpenSees.
**
**    Developed by Ross W. Boulanger
** Copyright @ 2002 The Regents of the University of California (The Regents). All Rights Reserved.
**
** The Regents grants permission, without fee and without a written license agreement, for (a) use, 
** reproduction, modification, and distribution of this software and its documentation by educational, 
** research, and non-profit entities for noncommercial purposes only; and (b) use, reproduction and 
** modification of this software by other entities for internal purposes only. The above copyright 
** notice, this paragraph and the following three paragraphs must appear in all copies and modifications 
** of the software and/or documentation.
**
** Permission to incorporate this software into products for commercial distribution may be obtained 
** by contacting the University of California 
** Office of Technology Licensing 
** 2150 Shattuck Avenue #510, 
** Berkeley, CA 94720-1620, 
** (510) 643-7201.
**
** This software program and documentation are copyrighted by The Regents of the University of California. 
** The Regents does not warrant that the operation of the program will be uninterrupted or error-free. The 
** end-user understands that the program was developed for research purposes and is advised not to rely 
** exclusively on the program for any reason.
**
** IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR 
** CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS 
** DOCUMENTATION, EVEN IF REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. REGENTS GRANTS 
** NO EXPRESS OR IMPLIED LICENSE IN ANY PATENT RIGHTS OF REGENTS BUT HAS IMPLEMENTED AN INDIVIDUAL 
** CONTRIBUTOR LICENSE AGREEMENT FOR THE OPENSEES PROJECT AT THE UNIVERISTY OF CALIFORNIA, BERKELEY 
** TO BENEFIT THE END USER.
**
** REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES 
** OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION,
** IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, 
** SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
**
** ****************************************************************** */

// $Revision: 1.0
// $Date: 2001/1/22
// $Source: /OpenSees/SRC/material/uniaxial/QzSimple1.cpp

// Written: RWB
// Created: Jan 2002
// Revision: A
// tested and checked: Boris Jeremic (jeremic@ucdavis.edu) Spring 2002
//
// Description: This file contains the class implementation for QzSimple1

#include <stdlib.h>
#include <math.h>
#include "QzSimple1.h"
#include <Vector.h>
#include <Channel.h>

// Controls on internal iterations between spring components
const int QZmaxIterations = 20;
const double QZtolerance = 1.0e-12;

/////////////////////////////////////////////////////////////////////
//	Constructor with data

QzSimple1::QzSimple1(int tag, int qzChoice, double Q_ult, double z_50,
				 double suctionRatio, double dash_pot)
:UniaxialMaterial(tag,MAT_TAG_QzSimple1),
 QzType(qzChoice), Qult(Q_ult), z50(z_50), suction(suctionRatio), dashpot(dash_pot)
{
  // Initialize QzSimple variables and history variables
  //
  this->revertToStart();
  initialTangent = Ttangent;
}

/////////////////////////////////////////////////////////////////////
//	Default constructor

QzSimple1::QzSimple1()
:UniaxialMaterial(0,MAT_TAG_QzSimple1),
 QzType(0), Qult(0.0), z50(0.0), suction(0.0), dashpot(0.0)
{
  // Initialize variables .. WILL NOT WORK AS NOTHING SET
  // this->revertToStart();

  // need to set iterations and tolerance

  // BTW maxIterations and tolerance should not be private variables, they
  // should be static .. all PySimple1 materials share the same values & 
  // these values don't change

}

/////////////////////////////////////////////////////////////////////
//	Default destructor
QzSimple1::~QzSimple1()
{
    // Does nothing
}

/////////////////////////////////////////////////////////////////////
void QzSimple1::getGap(double zlast, double dz, double dz_old)
{
	// For stability in Closure spring, limit "dz" step size to avoid
	// overshooting on the "closing" or "opening" of the gap.
	//
	if(zlast > 0.0 && (zlast + dz) < -QZtolerance) dz = -QZtolerance - zlast;
	if(zlast < 0.0 && (zlast + dz) >  QZtolerance) dz =  QZtolerance - zlast;
	TGap_z = zlast + dz;

	// Combine the Suction and Closure elements in parallel
	//
	getClosure(zlast,dz);
	getSuction(zlast,dz);
	TGap_Q = TSuction_Q + TClose_Q;
	TGap_tang = TSuction_tang + TClose_tang;

	return;
}

/////////////////////////////////////////////////////////////////////
void QzSimple1::getFarField(double z)
{
	TFar_z   = z;
	TFar_tang= TFar_tang;
	TFar_Q   = TFar_tang * TFar_z;

	return;
}

/////////////////////////////////////////////////////////////////////
void QzSimple1::getClosure(double zlast, double dz)
{
	TClose_z = zlast + dz;
	
	// Loading on the stiff "closed gap"
	//
	if(TClose_z <= 0.0) 
	{
		TClose_tang = 1000.0*Qult/z50;
		TClose_Q    = TClose_z * TClose_tang;
	}

	// Loading on the soft "open gap"
	//
	if(TClose_z > 0.0) 
	{
		TClose_tang = 0.001*Qult/z50;
		TClose_Q    = TClose_z * TClose_tang;
	}

	return;
}

/////////////////////////////////////////////////////////////////////
void QzSimple1::getSuction(double zlast, double dz)
{
	TSuction_z = zlast + dz;
	double Qmax=suction*Qult;
	double dzTotal=TSuction_z - CSuction_z;

	// Treat as elastic if dzTotal is below QZtolerance
	//
	if(fabs(dzTotal*TSuction_tang/Qult) < 3.0*QZtolerance) 
	{
		TSuction_Q = TSuction_Q + dz*TSuction_tang;
		if(fabs(TSuction_Q) >= Qmax) 
			TSuction_Q =(TSuction_Q/fabs(TSuction_Q))*(1.0-1.0e-8)*Qmax;
		return;
	}

	// Reset the history terms to the last Committed values, and let them
	// reset if the reversal of loading persists in this step.
	//
	if(TSuction_Qin != CSuction_Qin)
	{
		TSuction_Qin = CSuction_Qin;
		TSuction_zin = CSuction_zin;
	}

	// Change from positive to negative direction
	//
	if(CSuction_z > CSuction_zin && dzTotal < 0.0)
	{
		TSuction_Qin = CSuction_Q;
		TSuction_zin = CSuction_z;
	}
	// Change from negative to positive direction
	//
	if(CSuction_z < CSuction_zin && dzTotal > 0.0)
	{
		TSuction_Qin = CSuction_Q;
		TSuction_zin = CSuction_z;
	}
	
	// Positive loading
	//
	if(dzTotal >= 0.0)
	{
		TSuction_Q=Qmax-(Qmax-TSuction_Qin)*pow(0.5*z50,nd)
					*pow(0.5*z50 + TSuction_z - TSuction_zin,-nd);
		TSuction_tang=nd*(Qmax-TSuction_Qin)*pow(0.5*z50,nd)
					*pow(0.5*z50 + TSuction_z - TSuction_zin,-nd-1.0);
	}

	// Negative loading
	//
	if(dzTotal < 0.0)
	{
		TSuction_Q=-Qmax+(Qmax+TSuction_Qin)*pow(0.5*z50,nd)
					*pow(0.5*z50 - TSuction_z + TSuction_zin,-nd);
		TSuction_tang=nd*(Qmax+TSuction_Qin)*pow(0.5*z50,nd)
					*pow(0.5*z50 - TSuction_z + TSuction_zin,-nd-1.0);
	}

	// Ensure that |Q|<Qmax and tangent not zero or negative.
	//
	if(fabs(TSuction_Q) >= (1.0-QZtolerance)*Qmax) {
		TSuction_Q =(TSuction_Q/fabs(TSuction_Q))*(1.0-QZtolerance)*Qmax;}
	if(TSuction_tang <=1.0e-4*Qult/z50) TSuction_tang = 1.0e-4*Qult/z50;

	return;
}

/////////////////////////////////////////////////////////////////////
void QzSimple1::getNearField(double zlast, double dz, double dz_old)
{
	// Limit "dz" step size if it is oscillating in sign and not shrinking
	//
	if(dz*dz_old < 0.0 && fabs(dz/dz_old) > 0.5) dz = -dz_old/2.0;

	// Set "dz" so "z" is at middle of elastic zone if oscillation is large.
	//
	if(dz*dz_old < -z50*z50) {
		dz = (TNF_zinr + TNF_zinl)/2.0 - zlast;
	}
	
	// Establish trial "z" and direction of loading (with NFdz) for entire step
	//
	TNF_z = zlast + dz;
	double NFdz = TNF_z - CNF_z;

	// Treat as elastic if NFdz is below QZtolerance
	//
	if(fabs(NFdz*TNF_tang/Qult) < 3.0*QZtolerance) 
	{
		TNF_Q = TNF_Q + dz*TNF_tang;
		if(fabs(TNF_Q) >=Qult) TNF_Q=(TNF_Q/fabs(TNF_Q))*(1.0-QZtolerance)*Qult;
		return;
	}

	// Reset the history terms to the last Committed values, and let them
	// reset if the reversal of loading persists in this step.
	//
	if(TNF_Qinr != CNF_Qinr || TNF_Qinl != CNF_Qinl)
	{
		TNF_Qinr = CNF_Qinr;
		TNF_Qinl = CNF_Qinl;
		TNF_zinr = CNF_zinr;
		TNF_zinl = CNF_zinl;
	}

	// For stability, may have to limit "dz" step size if direction changed.
	//
	bool changeDirection = false;
	
	// Direction change from a yield point triggers new Elastic range
	//
	if(CNF_Q > CNF_Qinr && NFdz <0.0){				// from pos to neg
		changeDirection = true;
		if((CNF_Q - CNF_Qinl) > 2.0*Qult*Elast) Elast=(CNF_Q - CNF_Qinl)/(2.0*Qult);
		if(2.0*Elast > maxElast) Elast=maxElast/2.0;
		TNF_Qinr = CNF_Q;
		TNF_Qinl = TNF_Qinr - 2.0*Qult*Elast;
		TNF_zinr = CNF_z;
		TNF_zinl = TNF_zinr - (TNF_Qinr-TNF_Qinl)/NFkrig; 
	}
	if(CNF_Q < CNF_Qinl && NFdz > 0.0){				// from neg to pos
		changeDirection = true;
		if((CNF_Qinr - CNF_Q) > 2.0*Qult*Elast) Elast=(CNF_Qinr - CNF_Q)/(2.0*Qult);
		if(2.0*Elast > maxElast) Elast=maxElast/2.0;
		TNF_Qinl = CNF_Q;
		TNF_Qinr = TNF_Qinl + 2.0*Qult*Elast;
		TNF_zinl = CNF_z;
		TNF_zinr = TNF_zinl + (TNF_Qinr-TNF_Qinl)/NFkrig; 
	}

	// Now if there was a change in direction, limit the step size "dz"
	//
	if(changeDirection == true) {
		double maxdz = Elast*Qult/NFkrig;
		if(fabs(dz) > maxdz) dz = (dz/fabs(dz))*maxdz;
	}

	// Now, establish the trial value of "z" for use in this function call.
	//
	TNF_z = zlast + dz;

	// Postive loading
	//
	if(NFdz >= 0.0){
		// Check if elastic using z < zinr
		if(TNF_z <= TNF_zinr){							// stays elastic
			TNF_tang = NFkrig;
			TNF_Q = TNF_Qinl + (TNF_z - TNF_zinl)*NFkrig;
		}
		else {
			TNF_tang = np * (Qult-TNF_Qinr) * pow(zref,np) 
				* pow(zref - TNF_zinr + TNF_z, -np-1.0);
			TNF_Q = Qult - (Qult-TNF_Qinr)* pow(zref/(zref-TNF_zinr+TNF_z),np);
		}
	}

	// Negative loading
	//
	if(NFdz < 0.0){
		// Check if elastic using z < zinl
		if(TNF_z >= TNF_zinl){							// stays elastic
			TNF_tang = NFkrig;
			TNF_Q = TNF_Qinr + (TNF_z - TNF_zinr)*NFkrig;
		}
		else {
			TNF_tang = np * (Qult+TNF_Qinl) * pow(zref,np) 
				* pow(zref + TNF_zinl - TNF_z, -np-1.0);
			TNF_Q = -Qult + (Qult+TNF_Qinl)* pow(zref/(zref+TNF_zinl-TNF_z),np);
		}
	}

	// Ensure that |Q|<Qult and tangent not zero or negative.
	//
	if(fabs(TNF_Q) >= (1.0-QZtolerance)*Qult) { 
		TNF_Q=(TNF_Q/fabs(TNF_Q))*(1.0-QZtolerance)*Qult;
		TNF_tang = 1.0e-4*Qult/z50;
	}
	if(TNF_tang <= 1.0e-4*Qult/z50) TNF_tang = 1.0e-4*Qult/z50;

    return;
}

/////////////////////////////////////////////////////////////////////
int 
QzSimple1::setTrialStrain (double newz, double zRate)
{
	// Set trial values for displacement and load in the material
	// based on the last Tangent modulus.
	//
	double dz = newz - Tz;
	double dQ = Ttangent * dz;
	TzRate    = zRate;

	// Limit the size of step (dz or dQ) that can be imposed. Prevents
	// numerical difficulties upon load reversal at high loads
	// where a soft loading modulus becomes a stiff unloading modulus.
	//
	int numSteps = 1;
	double stepSize = 1.0;
	if(fabs(dQ/Qult) > 0.5) numSteps = 1 + int(fabs(dQ/(0.5*Qult)));
	if(fabs(dz/z50)  > 1.0 ) numSteps = 1 + int(fabs(dz/(1.0*z50)));
	stepSize = 1.0/float(numSteps);
	if(numSteps > 100) numSteps = 100;

	dz = stepSize * dz;

	// Main loop over the required number of substeps
	//
	for(int istep=1; istep <= numSteps; istep++)
	{
		Tz = Tz + dz;
		dQ = Ttangent * dz;
		
	// May substep within Gap or NearField element if oscillating, which can happen
	// when they jump from soft to stiff. Initialize history terms here.
	//
		double dz_gap_old = ((TQ + dQ) - TGap_Q)/TGap_tang;
		double dz_nf_old  = ((TQ + dQ) - TNF_Q) /TNF_tang;

	// Iterate to distribute displacement among the series components.
	// Use the incremental iterative strain & iterate at this strain.
	//
	for (int j=1; j < QZmaxIterations; j++)
	{
		TQ = TQ + dQ;
		if(fabs(TQ) >(1.0-QZtolerance)*Qult) TQ=(1.0-QZtolerance)*Qult*(TQ/fabs(TQ));

		// Stress & strain update in Near Field element
		double dz_nf = (TQ - TNF_Q)/TNF_tang;
		getNearField(TNF_z,dz_nf,dz_nf_old);
		
		// Residuals in Near Field element
		double Q_unbalance = TQ - TNF_Q;
		double zres_nf = (TQ - TNF_Q)/TNF_tang;
		dz_nf_old = dz_nf;

		// Stress & strain update in Gap element
		double dz_gap = (TQ - TGap_Q)/TGap_tang;
		getGap(TGap_z,dz_gap,dz_gap_old);

		// Residuals in Gap element
		double Q_unbalance2 = TQ - TGap_Q;
		double zres_gap = (TQ - TGap_Q)/TGap_tang;
		dz_gap_old = dz_gap;

		// Stress & strain update in Far Field element
		double dz_far = (TQ - TFar_Q)/TFar_tang;
		TFar_z = TFar_z + dz_far;
		getFarField(TFar_z);

		// Residuals in Far Field element
		double Q_unbalance3 = TQ - TFar_Q;
		double zres_far = (TQ - TFar_Q)/TFar_tang;

		// Update the combined tangent modulus
		Ttangent = pow(1.0/TGap_tang + 1.0/TNF_tang + 1.0/TFar_tang, -1.0);

		// Residual deformation across combined element
		double dv = Tz - (TGap_z + zres_gap)
			- (TNF_z + zres_nf) - (TFar_z + zres_far);

		// Residual "Q" increment 
		dQ = Ttangent * dv;

		// Test for convergence
		double Qsum = (fabs(Q_unbalance) + fabs(Q_unbalance2) + fabs(Q_unbalance3))/3.0;
		if(Qsum/Qult < QZtolerance) break;
	}
	}

	return 0;
}
/////////////////////////////////////////////////////////////////////
double 
QzSimple1::getStress(void)
{
	// Dashpot force is only due to velocity in the far field.
	// If converged, proportion by Tangents.
	// If not converged, proportion by ratio of displacements in components.
	//
	double ratio_disp =(1.0/TFar_tang)/(1.0/TFar_tang + 1.0/TNF_tang + 1.0/TGap_tang);
	if(Tz != Cz) {
		ratio_disp = (TFar_z - CFar_z)/(Tz - Cz);
		if(ratio_disp > 1.0) ratio_disp = 1.0;
		if(ratio_disp < 0.0) ratio_disp = 0.0;
	}
	double dashForce = dashpot * TzRate * ratio_disp;

	// Limit the combined force to Qult.
	//
	if(fabs(TQ + dashForce) >= (1.0-QZtolerance)*Qult)
		return (1.0-QZtolerance)*Qult*(TQ+dashForce)/fabs(TQ+dashForce);
	else return TQ + dashForce;
}
/////////////////////////////////////////////////////////////////////
double 
QzSimple1::getTangent(void)
{
    return this->Ttangent;
}
/////////////////////////////////////////////////////////////////////
double 
QzSimple1::getInitialTangent(void)
{
    return this->initialTangent;
}
/////////////////////////////////////////////////////////////////////
double 
QzSimple1::getDampTangent(void)
{
	// Damping tangent is produced only by the far field component.
	// If converged, proportion by Tangents.
	// If not converged, proportion by ratio of displacements in components.
	//
	double ratio_disp =(1.0/TFar_tang)/(1.0/TFar_tang + 1.0/TNF_tang + 1.0/TGap_tang);
	if(Tz != Cz) {
		ratio_disp = (TFar_z - CFar_z)/(Tz - Cz);
		if(ratio_disp > 1.0) ratio_disp = 1.0;
		if(ratio_disp < 0.0) ratio_disp = 0.0;
	}

	double DampTangent = dashpot * ratio_disp;

	// Minimum damping tangent referenced against Farfield spring
	//
	if(DampTangent < TFar_tang * 1.0e-12) DampTangent = TFar_tang * 1.0e-12;

	return DampTangent;
}
/////////////////////////////////////////////////////////////////////
double 
QzSimple1::getStrain(void)
{
    return this->Tz;
}
/////////////////////////////////////////////////////////////////////
double 
QzSimple1::getStrainRate(void)
{
    return this->TzRate;
}
/////////////////////////////////////////////////////////////////////
int 
QzSimple1::commitState(void)
{
  // Commit trial history variable -- Combined element
    Cz       = Tz;
    CQ       = TQ;
    Ctangent = Ttangent;
    
    // Commit trial history variables for Near Field component
    CNF_Qinr   = TNF_Qinr;
    CNF_Qinl   = TNF_Qinl; 
    CNF_zinr   = TNF_zinr;
    CNF_zinl   = TNF_zinl;	
    CNF_Q      = TNF_Q;
    CNF_z      = TNF_z;
    CNF_tang   = TNF_tang;
    
    // Commit trial history variables for Suction component
    CSuction_Qin  = TSuction_Qin;
    CSuction_zin  = TSuction_zin;
    CSuction_Q    = TSuction_Q;
    CSuction_z    = TSuction_z;
    CSuction_tang = TSuction_tang;
    
    // Commit trial history variables for Closure component
    CClose_Q      = TClose_Q;
    CClose_z      = TClose_z;
    CClose_tang   = TClose_tang;
    
    // Commit trial history variables for the Gap
    CGap_z    = TGap_z;
    CGap_Q    = TGap_Q;
    CGap_tang = TGap_tang;
    
    // Commit trial history variables for the Far Field
    CFar_z    = TFar_z;
    CFar_Q    = TFar_Q;
    CFar_tang = TFar_tang;
    
    return 0;
}

/////////////////////////////////////////////////////////////////////
int 
QzSimple1::revertToLastCommit(void)
{
  // Nothing to do here -- WRONG -- have a look at setTrialStrain() .. everything
  // calculated based on trial values & trial values updated in method .. need to 
  // reset to committed values
  
  // for convenience i am just gonna do the reverse of commit 
  Tz       = Cz;
  TQ       = CQ;
  Ttangent = Ctangent;
  
  TNF_Qinr   = CNF_Qinr;
  TNF_Qinl   = CNF_Qinl; 
  TNF_zinr   = CNF_zinr;
  TNF_zinl   = CNF_zinl;	
  TNF_Q      = CNF_Q;
  TNF_z      = CNF_z;
  TNF_tang   = CNF_tang;
  
  TSuction_Qin  = CSuction_Qin;
  TSuction_zin  = CSuction_zin;
  TSuction_Q    = CSuction_Q;
  TSuction_z    = CSuction_z;
  TSuction_tang = CSuction_tang;
  
  TClose_Q      = CClose_Q;
  TClose_z      = CClose_z;
  TClose_tang   = CClose_tang;
  
  TGap_z    = CGap_z;
  TGap_Q    = CGap_Q;
  TGap_tang = CGap_tang;
  
  TFar_z    = CFar_z;
  TFar_Q    = CFar_Q;
  TFar_tang = CFar_tang;

  return 0;
}

/////////////////////////////////////////////////////////////////////
int 
QzSimple1::revertToStart(void)
{

	// Reset gap "suction" if zero (or negative) or exceeds max value of 0.1
	//
	if(suction <= QZtolerance) suction = QZtolerance;
	if(suction > 0.1){
	  suction = 0.1;
	  opserr << "QzSimple1::QzSimple1 -- setting suction to max value of 0.1\n";
	}

	// Only allow zero or positive dashpot values
	//
	if(dashpot < 0.0) dashpot = 0.0;

	// Do not allow zero or negative values for z50 or Qult.
	//
	if(Qult <= 0.0 || z50 <= 0.0) {
	  opserr << "QzSimple1::QzSimple1 -- only accepts positive nonzero Qult and z50\n";
	  exit(-1);
	}

	// Initialize variables for Near Field rigid-plastic spring
	//
	if(QzType ==1) {	// Approx Reese & O'Neill (1987) drilled shafts on clay
		zref	= 0.35*z50;
		np		= 1.2;
		Elast	= 0.2;
		maxElast= 0.7;
		nd		= 1.0;
		TFar_tang= 0.525*Qult/z50;
	}
	else if (QzType == 2){
		zref	= 12.3*z50;
		np		= 5.5;
		Elast	= 0.3;
		maxElast= 0.7;
		nd		= 1.0;
		TFar_tang= 1.39*Qult/z50;
	}
	else{
	  opserr << "QzSimple1::QzSimple1 -- only accepts QzType of 1 or 2\n";
	  exit(-1);
	}

	// Far Field components: TFar_tang was set under "soil type" statements.
	//
	TFar_Q  = 0.0;
	TFar_z  = 0.0;

	// Near Field components
	//
	NFkrig   = 10000.0 * Qult / z50;
    TNF_Qinr = Elast*Qult;
	TNF_Qinl = -TNF_Qinr;
	TNF_zinr = TNF_Qinr / NFkrig;
	TNF_zinl = -TNF_zinr;
	TNF_Q    = 0.0;
	TNF_z    = 0.0;
	TNF_tang = NFkrig;

	// Suction components
	//
	TSuction_Qin  = 0.0;
	TSuction_zin  = 0.0;
	TSuction_Q    = 0.0;
	TSuction_z    = 0.0;
	TSuction_tang = nd*(Qult*suction-TSuction_Q)*pow(z50/2.0,nd)
					*pow(z50/2.0 - TSuction_z + TSuction_zin,-nd-1.0);

	// Closure components
	//
	TClose_Q     = 0.0; 
	TClose_z     = 0.0;
	TClose_tang  = 100.0*Qult/z50;

	// Gap (Suction + Closure in parallel)
	//
	TGap_z   = 0.0;
	TGap_Q   = 0.0;
	TGap_tang= TClose_tang + TSuction_tang;

	// Entire element (Far field + Near field + Gap in series)
	//
	Tz       = 0.0;
	TQ       = 0.0;
	Ttangent = pow(1.0/TGap_tang + 1.0/TNF_tang + 1.0/TFar_tang, -1.0);
	TzRate   = 0.0;

	// Now get all the committed variables initiated
	//
	this->commitState();

    return 0;
}

/////////////////////////////////////////////////////////////////////
UniaxialMaterial *
QzSimple1::getCopy(void)
{
    QzSimple1 *theCopy =
	new QzSimple1(this->getTag(),QzType,Qult,z50,suction,dashpot);

	// Copy parameters
	theCopy->zref    = zref;
	theCopy->np      = np;
	theCopy->Elast   = Elast;
	theCopy->maxElast= maxElast;
	theCopy->nd      = nd;

	// Copy internal parameters or constants
	theCopy->NFkrig  = NFkrig;
    
	// Copy committed history variables for Near Field
    theCopy->CNF_Qinr = CNF_Qinr;
    theCopy->CNF_Qinl = CNF_Qinl;
    theCopy->CNF_zinr = CNF_zinr;
    theCopy->CNF_zinl = CNF_zinl;
    theCopy->CNF_Q    = CNF_Q;
    theCopy->CNF_z    = CNF_z;
    theCopy->CNF_tang = CNF_tang;

	// Copy trial history variables for Near Field
    theCopy->TNF_Qinr = TNF_Qinr;
    theCopy->TNF_Qinl = TNF_Qinl;
    theCopy->TNF_zinr = TNF_zinr;
    theCopy->TNF_zinl = TNF_zinl;
    theCopy->TNF_Q    = TNF_Q;
    theCopy->TNF_z    = TNF_z;
    theCopy->TNF_tang = TNF_tang;

	// Copy committed history variables for Suction component
    theCopy->CSuction_Qin  = CSuction_Qin;
    theCopy->CSuction_zin  = CSuction_zin;
	theCopy->CSuction_Q    = CSuction_Q;
	theCopy->CSuction_z    = CSuction_z;
	theCopy->CSuction_tang = CSuction_tang;

	// Copy trial history variables for Suction component
    theCopy->TSuction_Qin  = TSuction_Qin;
    theCopy->TSuction_zin  = TSuction_zin;
	theCopy->TSuction_Q    = TSuction_Q;
	theCopy->TSuction_z    = TSuction_z;
	theCopy->TSuction_tang = TSuction_tang;

	// Copy committed history variables for Closure component
	theCopy->CClose_Q     = CClose_Q; 
	theCopy->CClose_z     = CClose_z;
	theCopy->CClose_tang  = CClose_tang;
	
	// Copy trail history variables for Closure component	
	theCopy->TClose_Q     = TClose_Q; 
	theCopy->TClose_z     = TClose_z;
	theCopy->TClose_tang  = TClose_tang;

	// Copy committed history variables for Gap component
	theCopy->CGap_z    = CGap_z;
	theCopy->CGap_Q    = CGap_Q;
	theCopy->CGap_tang = CGap_tang;	

	// Copy trial history variables for Gap component
	theCopy->TGap_z    = TGap_z;
	theCopy->TGap_Q    = TGap_Q;
	theCopy->TGap_tang = TGap_tang;	

	// Copy committed history variables for Far Field component
	theCopy->CFar_z    = CFar_z;
	theCopy->CFar_Q    = CFar_Q;
	theCopy->CFar_tang = CFar_tang;	

	// Copy trial history variables for Far Field component
	theCopy->TFar_z    = TFar_z;
	theCopy->TFar_Q    = TFar_Q;
	theCopy->TFar_tang = TFar_tang;	

	// Copy committed history variables for Entire Material
	theCopy->Cz        = Cz;
	theCopy->CQ        = CQ;
	theCopy->Ctangent  = Ctangent;

	// Copy trial history variables for Entire Material
	theCopy->Tz        = Tz;
	theCopy->TQ        = TQ;
	theCopy->Ttangent  = Ttangent;
	theCopy->TzRate    = TzRate;

    return theCopy;
}

/////////////////////////////////////////////////////////////////////
int 
QzSimple1::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(38);
  
  data(0) = this->getTag();
  data(1) = QzType;
  data(2) = Qult;
  data(3) = z50;
  data(4) = suction;
  data(5) = dashpot;
  data(6) = zref;
  data(7) = np;
  data(8) = Elast;
  data(9) = maxElast;
  data(10)= nd;
  data(11)= NFkrig;

  data(12) = CNF_Qinr;
  data(13) = CNF_Qinl;
  data(14) = CNF_zinr;
  data(15) = CNF_zinl;
  data(16) = CNF_Q;
  data(17) = CNF_z;
  data(18) = CNF_tang;

  data(19) = CSuction_Qin;
  data(20) = CSuction_zin;
  data(21) = CSuction_Q;
  data(22) = CSuction_z;
  data(23) = CSuction_tang;

  data(24) = CClose_Q;
  data(25) = CClose_z;
  data(26) = CClose_tang;

  data(27) = CGap_z;
  data(28) = CGap_Q;
  data(29) = CGap_tang;

  data(30) = CFar_z;
  data(31) = CFar_Q;
  data(32) = CFar_tang;

  data(33) = Cz;
  data(34) = CQ;
  data(35) = Ctangent;
  data(36) = TzRate;

  data(37) = initialTangent;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "QzSimple1::sendSelf() - failed to send data\n";

  return res;
}

/////////////////////////////////////////////////////////////////////
int 
QzSimple1::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(38);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "QzSimple1::recvSelf() - failed to receive data\n";
      CNF_tang = 0; 
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
	QzType   = (int)data(1);
	Qult     = data(2);
	z50      = data(3);
	suction  = data(4);
	dashpot  = data(5);
	zref     = data(6);
	np       = data(7);
	Elast    = data(8);
	maxElast = data(9);
	nd       = data(10);
	NFkrig   = data(11);

	CNF_Qinr = data(12);
	CNF_Qinl = data(13);
	CNF_zinr = data(14);
	CNF_zinl = data(15);
	CNF_Q    = data(16);
	CNF_z    = data(17);
	CNF_tang = data(18);

	CSuction_Qin  = data(19);
	CSuction_zin  = data(20);
	CSuction_Q    = data(21);
	CSuction_z    = data(22);
	CSuction_tang = data(23);

	CClose_Q      = data(24);
	CClose_z      = data(25);
	CClose_tang   = data(26);

	CGap_z    = data(27);
	CGap_Q    = data(28);
	CGap_tang = data(29);

	CFar_z    = data(30);
	CFar_Q    = data(31);
	CFar_tang = data(32);

	Cz        = data(33);
	CQ        = data(34);
	Ctangent  = data(35);
	TzRate    = data(36);
	
	initialTangent = data(37);

	// set the trial quantities
	this->revertToLastCommit();
  }
    
  return res;
}

/////////////////////////////////////////////////////////////////////
void 
QzSimple1::Print(OPS_Stream &s, int flag)
{
    s << "QzSimple1, tag: " << this->getTag() << endln;
    s << "  QzType: " << QzType << endln;
    s << "  Qult: " << Qult << endln;
    s << "  z50: " << z50 << endln;
    s << "  suction: " << suction << endln;
	s << "  dashpot: " << dashpot << endln;
}

/////////////////////////////////////////////////////////////////////

