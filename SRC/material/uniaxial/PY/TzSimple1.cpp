/* *********************************************************************
**    Module:	TzSimple1.cpp 
**
**    Purpose:	Provide a simple t-z spring for OpenSees.
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
// $Date: 2002/1/19
// $Source: /OpenSees/SRC/material/uniaxial/TzSimple1.cpp

// Written: RWB
// Created: Jan 2002
// Revision: A
// tested and checked: Boris Jeremic (jeremic@ucdavis.edu) Spring 2002
//
// Description: This file contains the class implementation for TzSimple1

#include <stdlib.h>

#include "TzSimple1.h"
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <elementAPI.h>

// Controls on internal iteration between spring components
const int TZmaxIterations = 20;
const double TZtolerance = 1.0e-12;

void* OPS_TzSimple1()
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 4) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: uniaxialMaterial TzSimple1 tag? tzType? tult? z50? dashpot?\n";
	return 0;
    }
    
    int idata[2];
    numdata = 2;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
	opserr << "WARNING invalid int inputs\n";
	return 0;
    }
    
    double ddata[3] = {0,0,0};
    numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 3) numdata = 3;
    if (OPS_GetDoubleInput(&numdata, ddata) < 0) {
	opserr << "WARNING invalid double inputs\n";
	return 0;
    }
    
    UniaxialMaterial *theMaterial = 0;
    theMaterial = new TzSimple1(idata[0], MAT_TAG_PySimple1, idata[1], ddata[0], ddata[1],
    				ddata[2]);
    
    return theMaterial;
}

/////////////////////////////////////////////////////////////////////
//	Constructor with data

TzSimple1::TzSimple1(int tag,int classtag, int tz_type,double t_ult,double z_50,double dash_pot)
:UniaxialMaterial(tag,classtag),
 tzType(tz_type), tult(t_ult), z50(z_50), dashpot(dash_pot)
{
  // Initialize TzSimple variables and history variables
  //
  this->revertToStart();
  initialTangent = Ttangent;
}

/////////////////////////////////////////////////////////////////////
//	Default constructor

TzSimple1::TzSimple1()
:UniaxialMaterial(0,0),
 tzType(0), tult(0.0), z50(0.0), dashpot(0.0)
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
TzSimple1::~TzSimple1()
{
    // Does nothing
}

/////////////////////////////////////////////////////////////////////
void TzSimple1::getFarField(double z)
{
	TFar_z   = z;
	TFar_tang= TFar_tang;
	TFar_t   = TFar_tang * TFar_z;

	return;
}

/////////////////////////////////////////////////////////////////////
void TzSimple1::getNearField(double zlast, double dz, double dz_old)
{
	// Limit "dz" step size if it is osillating and not shrinking.
	//
	if(dz*dz_old < 0.0 && fabs(dz/dz_old) > 0.5) dz = -dz_old/2.0;

	// Establish trial "z" and direction of loading (with dzTotal) for entire step.
	//	
	TNF_z = zlast + dz;
	double dzTotal = TNF_z - CNF_z;

	// Treat as elastic if dzTotal is below TZtolerance
	//
	if(fabs(dzTotal*TNF_tang/tult) < 10.0*TZtolerance) 
	{
		TNF_t = TNF_t + dz*TNF_tang;
		if(fabs(TNF_t) >=(1.0-TZtolerance)*tult) 
			TNF_t =(TNF_t/fabs(TNF_t))*(1.0-TZtolerance)*tult;
		return;
	}

	// Reset the history terms to the last Committed values, and let them
	// reset if the reversal of loading persists in this step.
	//
	if(TNF_tin != CNF_tin)
	{
		TNF_tin = CNF_tin;
		TNF_zin = CNF_zin;
	}

	// Change from positive to negative direction
	//
	if(CNF_z > CNF_zin && dzTotal < 0.0)
	{
		TNF_tin = CNF_t;
		TNF_zin = CNF_z;
	}

	// Change from negative to positive direction
	//
	if(CNF_z < CNF_zin && dzTotal > 0.0)
	{
		TNF_tin = CNF_t;
		TNF_zin = CNF_z;
	}
	
	// Positive loading
	//
	if(dzTotal > 0.0)
	{
		TNF_t=tult-(tult-TNF_tin)*pow(zref,np)
					*pow(zref + TNF_z - TNF_zin,-np);
		TNF_tang=np*(tult-TNF_tin)*pow(zref,np)
					*pow(zref + TNF_z - TNF_zin,-np-1.0);
	}
	// Negative loading
	//
	if(dzTotal < 0.0)
	{
		TNF_t=-tult+(tult+TNF_tin)*pow(zref,np)
					*pow(zref - TNF_z + TNF_zin,-np);
		TNF_tang=np*(tult+TNF_tin)*pow(zref,np)
					*pow(zref - TNF_z + TNF_zin,-np-1.0);
	}

	// Ensure that |t|<tult and tangent not zero or negative.
	//
	if(fabs(TNF_t) >=tult) {
		TNF_t =(TNF_t/fabs(TNF_t))*(1.0-TZtolerance)*tult;}
	if(TNF_tang <=1.0e-4*tult/z50) TNF_tang = 1.0e-4*tult/z50;

	return;
}

/////////////////////////////////////////////////////////////////////
int 
TzSimple1::setTrialStrain (double newz, double zRate)
{
	// Set trial values for displacement and load in the material
	// based on the last Tangent modulus.
	//
	double dz = newz - Tz;
	double dt = Ttangent * dz;
	TzRate    = zRate;

	// Limit the size of step (dz or dt) that can be imposed. Prevents
	// oscillation under large load reversal steps
	//
	int numSteps = 1;
	double stepSize = 1.0;
	if(fabs(dt/tult) > 0.5)  numSteps = 1 + int(fabs(dt/(0.5*tult)));
	if(fabs(dz/z50)  > 1.0 ) numSteps = 1 + int(fabs(dz/(1.0*z50)));
	stepSize = 1.0/float(numSteps);
	if(numSteps > 100) numSteps = 100;

	dz = stepSize * dz;

	// Main loop over the required number of substeps
	//
	for(int istep=1; istep <= numSteps; istep++)
	{
		Tz = Tz + dz;
		dt = Ttangent * dz;
		
		// May substep in NearField component if not making progress due to oscillation
		// The following history term is initialized here.
		//
		double dz_nf_old = ((Tt+dt) - TNF_t)/TNF_tang;
		
	// Iterate to distribute displacement between elastic & plastic components.
	// Use the incremental iterative strain & iterate at this strain.
	//
	for (int j=1; j < TZmaxIterations; j++)
	{
		Tt = Tt + dt;
		if(fabs(Tt) >(1.0-TZtolerance)*tult) Tt=(1.0-TZtolerance)*tult*(Tt/fabs(Tt));

		// Stress & strain update in Near Field element
		double dz_nf = (Tt - TNF_t)/TNF_tang;
		getNearField(TNF_z,dz_nf,dz_nf_old);
		
		// Residuals in Near Field element
		double t_unbalance = Tt - TNF_t;
		double zres_nf = (Tt - TNF_t)/TNF_tang;
		dz_nf_old = dz_nf;

		// Stress & strain update in Far Field element
		double dz_far = (Tt - TFar_t)/TFar_tang;
		TFar_z = TFar_z + dz_far;
		getFarField(TFar_z);

		// Residuals in Far Field element
		double t_unbalance2 = Tt - TFar_t;
		double zres_far = (Tt - TFar_t)/TFar_tang;

		// Update the combined tangent modulus
		Ttangent = pow(1.0/TNF_tang + 1.0/TFar_tang, -1.0);

		// Residual deformation across combined element
		double dv = Tz - (TNF_z + zres_nf) - (TFar_z + zres_far);

		// Residual "t" increment 
		dt = Ttangent * dv;

		// Test for convergence
		double tsum = fabs(t_unbalance) + fabs(t_unbalance2);
		if(tsum/tult < TZtolerance) break;
	}
	}

	return 0;
}
/////////////////////////////////////////////////////////////////////
double 
TzSimple1::getStress(void)
{
	// Dashpot force is only due to velocity in the far field.
	// If converged, proportion by Tangents.
	// If not converged, proportion by ratio of displacements in components.
	//
	double ratio_disp =(1.0/TFar_tang)/(1.0/TFar_tang + 1.0/TNF_tang);
	if(Tz != Cz) {
		ratio_disp = (TFar_z - CFar_z)/(Tz - Cz);
		if(ratio_disp > 1.0) ratio_disp = 1.0;
		if(ratio_disp < 0.0) ratio_disp = 0.0;
	}
	double dashForce = dashpot * TzRate * ratio_disp;

	// Limit the combined force to tult.
	//
	if(fabs(Tt + dashForce) >= (1.0-TZtolerance)*tult)
		return (1.0-TZtolerance)*tult*(Tt+dashForce)/fabs(Tt+dashForce);
	else return Tt + dashForce;
}
/////////////////////////////////////////////////////////////////////
double 
TzSimple1::getTangent(void)
{
    return this->Ttangent;
}
/////////////////////////////////////////////////////////////////////
double 
TzSimple1::getInitialTangent(void)
{
    return this->initialTangent;
}
/////////////////////////////////////////////////////////////////////
double 
TzSimple1::getDampTangent(void)
{
	// Damping tangent is produced only by the far field component.
	// If converged, proportion by Tangents.
	// If not converged, proportion by ratio of displacements in components.
	//
	double ratio_disp =(1.0/TFar_tang)/(1.0/TFar_tang + 1.0/TNF_tang);
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
TzSimple1::getStrain(void)
{
    return this->Tz;
}
/////////////////////////////////////////////////////////////////////
double 
TzSimple1::getStrainRate(void)
{
    return this->TzRate;
}
/////////////////////////////////////////////////////////////////////
int
TzSimple1::commitState(void)
{
  // Commit trial history variable -- Combined element
  Cz       = Tz;
  Ct       = Tt;
  Ctangent = Ttangent;
  
  // Commit trial history variables for Near Field component
  CNF_tin = TNF_tin;
  CNF_zin = TNF_zin;
  CNF_t   = TNF_t;
  CNF_z   = TNF_z;
  CNF_tang= TNF_tang;
  
  // Commit trial history variables for the Far Field
  CFar_z    = TFar_z;
  CFar_t    = TFar_t;
  CFar_tang = TFar_tang;
    
  return 0;
}

/////////////////////////////////////////////////////////////////////
int 
TzSimple1::revertToLastCommit(void)
{
  // Nothing to do here -- WRONG -- have a look at setTrialStrain() .. everything
  // calculated based on trial values & trial values updated in method .. need to 
  // reset to committed values
  
  // for convenience i am just gonna do the reverse of commit
  Tz       = Cz;
  Tt       = Ct;
  Ttangent = Ctangent;
  
  TNF_tin = CNF_tin;
  TNF_zin = CNF_zin;
  TNF_t   = CNF_t;
  TNF_z   = CNF_z;
  TNF_tang= CNF_tang;
  
  TFar_z    = CFar_z;
  TFar_t    = CFar_t;
  TFar_tang = CFar_tang;

  return 0;
}

/////////////////////////////////////////////////////////////////////
int 
TzSimple1::revertToStart(void)
{

	// If tzType = 0, then it is entering with the default constructor.
	// To avoid division by zero, set small nonzero values for terms.
	//
	if(tzType == 0){
		tult = 1.0e-12;
		z50  = 1.0e12;
	}

	// Only allow zero or positive dashpot values
	//
	if(dashpot < 0.0) dashpot = 0.0;

	// Do not allow zero or negative values for z50 or tult.
	//
	if(tult <= 0.0 || z50 <= 0.0){
		opserr << "WARNING -- only accepts positive nonzero tult and z50" << endln;
		opserr << "TzLiq1: " << endln;
		opserr << "tzType: " << tzType << endln;
		exit(-1);
	}
		
	// Initialize variables for Near Field plastic component
	//
	if(tzType ==0) {			// This will happen with default constructor
		zref  = 0.5*z50;
		np    = 1.5;
		TFar_tang   = 0.70791*tult/(z50);
	}
	else if(tzType ==1) {		// Backbone approximates Reese & O'Neill 1987
		zref  = 0.5*z50;
		np    = 1.5;
		TFar_tang	= 0.70791*tult/(z50);
	}
	else if (tzType == 2){		// Backbone approximates Mosher 1984
		zref  = 0.6*z50;
		np    = 0.85;
		TFar_tang   = 2.0504*tult/z50;
	}
	else{
		opserr << "WARNING -- only accepts tzType of 1 or 2" << endln;
		opserr << "TzLiq1: " << endln;
		opserr << "tzType: " << tzType << endln;
		exit(-1);
	}

	// Far Field components: TFar_tang was set under "tzType" statements.
	//
	TFar_t  = 0.0;
	TFar_z  = 0.0;

	// Near Field components
	//
	TNF_tin = 0.0;
	TNF_zin = 0.0;
	TNF_t   = 0.0;
	TNF_z   = 0.0;
	TNF_tang= np*tult*pow(zref,np)*pow(zref,-np-1.0);

	// Entire element (Far field + Near field + Gap in series)
	//
	Tz       = 0.0;
	Tt       = 0.0;
	Ttangent = pow(1.0/TNF_tang + 1.0/TFar_tang, -1.0);
	TzRate   = 0.0;

	// Now get all the committed variables initiated
	//
	this->commitState();

    return 0;
}

/////////////////////////////////////////////////////////////////////
UniaxialMaterial *
TzSimple1::getCopy(void)
{
    TzSimple1 *theCopy;			// pointer to a TzSimple1 class
	theCopy = new TzSimple1();	// new instance of this class
	*theCopy= *this;			// theCopy (dereferenced) = this (dereferenced pointer)
	return theCopy;
}

/////////////////////////////////////////////////////////////////////
int 
TzSimple1::sendSelf(int cTag, Channel &theChannel)
{
	int res = 0;
  
	static Vector data(20);
  
	data(0) = this->getTag();
	data(1) = tzType;
	data(2) = tult;
	data(3) = z50;
	data(4) = dashpot;
	data(5) = zref;
	data(6) = np;

	data(7)  = CNF_tin;
	data(8)  = CNF_zin;
	data(9)  = CNF_t;
	data(10) = CNF_z;
	data(11) = CNF_tang;

	data(12) = CFar_z;
	data(13) = CFar_t;
	data(14) = CFar_tang;

	data(15) = Cz;
	data(16) = Ct;
	data(17) = Ctangent;
	data(18) = TzRate;

	data(19) = initialTangent;

	res = theChannel.sendVector(this->getDbTag(), cTag, data);
	if (res < 0) 
		opserr << "TzSimple1::sendSelf() - failed to send data\n";

	return res;
}

/////////////////////////////////////////////////////////////////////
int 
TzSimple1::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(20);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "TzSimple1::recvSelf() - failed to receive data\n";
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
	tzType = (int)data(1);
	tult     = data(2);
	z50      = data(3);
	dashpot  = data(4);
	zref     = data(5);
	np       = data(6);
	
	CNF_tin  = data(7);
    CNF_zin	 = data(8);
	CNF_t	 = data(9);
	CNF_z	 = data(10);
	CNF_tang = data(11);

	CFar_z    = data(12);
	CFar_t    = data(13);
	CFar_tang = data(14);

	Cz        = data(15);
	Ct        = data(16);
	Ctangent  = data(17);
	TzRate    = data(18);
	
	initialTangent = data(19);
	
	this->revertToLastCommit();
  }
    
  return res;
}

/////////////////////////////////////////////////////////////////////
void 
TzSimple1::Print(OPS_Stream &s, int flag)
{
    s << "TzSimple1, tag: " << this->getTag() << endln;
    s << "  tzType: " << tzType << endln;
    s << "  tult: " << tult << endln;
    s << "  z50: " << z50 << endln;
    s << "  dashpot: " << dashpot << endln;
}

/////////////////////////////////////////////////////////////////////

