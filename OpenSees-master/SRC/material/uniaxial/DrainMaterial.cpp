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
                                                                        
// $Revision: 1.5 $
// $Date: 2003-04-02 22:02:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/DrainMaterial.cpp,v $
                                                                        
// Written: MHS
// Created: June 2001
//
// Description: This file contains the class definition for 
// DrainMaterial. DrainMaterial wraps a Drain spring element subroutine
// and converts it to the UniaxialMaterial interface for use in
// zero length elements, beam sections, or anywhere else
// UniaxialMaterials may be used.
//
// For more information see the Drain-2DX user guide:
//    Allahabadi, R.; Powell, G. H.
//    UCB/EERC-88/06, Berkeley: Earthquake Engineering Research Center,
//    University of California, Mar. 1988, 1 vol.

#include <OPS_Globals.h>
#include <DrainMaterial.h>
#include <Vector.h>
#include <ID.h>
#include <Channel.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>

DrainMaterial::DrainMaterial(int tag, int classTag, int nhv, int ndata, double b)
:UniaxialMaterial(tag,classTag), data(0), hstv(0),
numData(ndata), numHstv(nhv), epsilonP(0.0), sigmaP(0.0), tangentP(0.0), beto(b),
epsilon(0.0), epsilonDot(0.0), sigma(0.0), tangent(0.0)
{
  if (numHstv < 0)
    numHstv = 0;

  if (numHstv > 0) {
    // Allocate history array
    hstv = new double[2*numHstv];
    if (hstv == 0) {
      opserr << "DrainMaterial::DrainMaterial -- failed to allocate history array -- type : " <<
	this->getClassTag() << endln;
      exit(-1);
    }
			    
    
    // Initialize to zero
    for (int i = 0; i < 2*numHstv; i++)
      hstv[i] = 0.0;
  }
  
  if (numData < 0)
    numData = 0;
  
  if (numData > 0) {
    // Allocate material parameter array
    data = new double[numData];
    if (data == 0) {
      opserr << "DrainMaterial::DrainMaterial -- failed to allocate data array -- type: " <<
	this->getClassTag() << endln;
      exit(-1);
    }
    
    // Initialize to zero
    for (int i = 0; i < numData; i++)
      data[i] = 0.0;
  }

  // determine initial tangent
  this->invokeSubroutine();
  initialTangent = tangent;
}

DrainMaterial::~DrainMaterial()
{
  if (hstv != 0)
    delete [] hstv;
  
  if (data != 0)
    delete [] data;
}

int
DrainMaterial::setTrialStrain(double strain, double strainRate)
{
	// Store the strain
	epsilon = strain;
	epsilonDot = strainRate;

	// Invoke Drain subroutine
	return this->invokeSubroutine();
}

int
DrainMaterial::setTrial(double strain, double &stress, double &stiff, double strainRate)
{
	// Store the strain
	epsilon = strain;
	epsilonDot = strainRate;

	// Invoke Drain subroutine
	int res = this->invokeSubroutine();

	stress = sigma;
	stiff = tangent;

	return res;
}

double
DrainMaterial::getStrain(void)
{
	return epsilon;
}

double
DrainMaterial::getStrainRate(void)
{
	return epsilonDot;
}

double
DrainMaterial::getStress(void)
{
	return sigma;
}

double
DrainMaterial::getTangent(void)
{
	return tangent;
}

double
DrainMaterial::getInitialTangent(void)
{
	return initialTangent;
}

double
DrainMaterial::getDampTangent(void)
{
	// Damping computed here using the last committed tangent
	// rather than the initial tangent!
	return beto*tangentP;
}

int
DrainMaterial::commitState(void)
{
	// Set committed values equal to corresponding trial values
	for (int i = 0; i < numHstv; i++)
		hstv[i] = hstv[i+numHstv];

	epsilonP = epsilon;
	sigmaP   = sigma;
	tangentP = tangent;

	return 0;
}

int
DrainMaterial::revertToLastCommit(void)
{
        // Set trial values equal to corresponding committed values
	for (int i = 0; i < numHstv; i++)
		hstv[i+numHstv] = hstv[i];

	epsilon = epsilonP;
	sigma   = sigmaP;
	tangent = tangentP;

	return 0;
}

int
DrainMaterial::revertToStart(void)
{
	// Set all trial and committed values to zero
	for (int i = 0; i < 2*numHstv; i++)
		hstv[i] = 0.0;

	epsilonP = 0.0;
	sigmaP   = 0.0;
	tangentP = 0.0;

	return 0;
}

// WARNING -- if you wish to override any method in this base class, you must
// also override the getCopy method to return a pointer to the derived class!!!
UniaxialMaterial*
DrainMaterial::getCopy(void)
{
	DrainMaterial *theCopy = 
		new DrainMaterial(this->getTag(), this->getClassTag(), numHstv, numData, beto);

	int i;
	for (i = 0; i < 2*numHstv; i++)
		theCopy->hstv[i] = hstv[i];

	for (i = 0; i < numData; i++)
		theCopy->data[i] = data[i];

	theCopy->epsilonP = epsilonP;
	theCopy->sigmaP   = sigmaP;
	theCopy->tangentP = tangentP;

	return theCopy;
}

int 
DrainMaterial::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	Vector vecData(numHstv+numData+5);

	int i, j;
	// Copy history variables into vector
	for (i = 0; i < numHstv; i++)
		vecData(i) = hstv[i];

	// Copy material properties into vector
	for (i = 0, j = numHstv; i < numData; i++, j++)
		vecData(j) = data[i];

	vecData(j++) = epsilonP;
	vecData(j++) = sigmaP;
	vecData(j++) = tangentP;
	vecData(j++) = beto;
	vecData(j++) = this->getTag();

	res += theChannel.sendVector(this->getDbTag(), commitTag, vecData);
	if (res < 0) 
		opserr << "DrainMaterial::sendSelf() - failed to send Vector data\n";

	return res;
}

int
DrainMaterial::recvSelf(int commitTag, Channel &theChannel,
							FEM_ObjectBroker &theBroker)
{
	int res = 0;

	Vector vecData(numHstv+numData+5);

	res += theChannel.recvVector(this->getDbTag(), commitTag, vecData);
	if (res < 0) {
		opserr << "DrainMaterial::recvSelf() - failed to receive Vector data\n";
		return res;
	}

	int i, j;
	// Copy history variables from vector
	for (i = 0; i < numHstv; i++) {
	  hstv[i] = vecData(i);
	  hstv[i+numHstv] = vecData(i);
	}

	// Copy material properties from vector
	for (i = 0, j = numHstv; i < numData; i++, j++)
	  data[i] = vecData(j);

	epsilonP = vecData(j++);
	sigmaP   = vecData(j++);
	tangentP = vecData(j++);
	beto     = vecData(j++);

	this->setTag((int)vecData(j));

	epsilon = epsilonP;
	sigma   = sigmaP;
	tangent = tangentP;

	return res;
}
    
void
DrainMaterial::Print(OPS_Stream &s, int flag)
{
	s << "DrainMaterial, type: ";
	
	switch (this->getClassTag()) {
	case MAT_TAG_DrainHardening:
		s << "Hardening" << endln;
		break;
	case MAT_TAG_DrainBilinear:
		s << "Bilinear" << endln;
		break;
	case MAT_TAG_DrainClough1:
		s << "Clough1" << endln;
		break;
	case MAT_TAG_DrainClough2:
		s << "Clough2" << endln;
		break;
	case MAT_TAG_DrainPinch1:
		s << "Pinch1" << endln;
		break;
	// Add more cases as needed

	default:
		s << "Material identifier = " << this->getClassTag() << endln;
		break;
	}
}

#ifdef _WIN32

// Declarations for the Hardening subroutines
extern "C" int FILL00(double *data, double *hstv, double *stateP);
extern "C" int RESP00(int *kresis, int *ksave, int *kgem, int *kstep,
							   int *ndof, int *kst, int *kenr,
							   double *ener, double *ened, double *enso, double *beto,
							   double *relas, double *rdamp, double *rinit,
							   double *ddise, double *dise, double *vele);
extern "C" int STIF00(int *kstt, int *ktype, int *ndof, double *fk);
extern "C" int GET00(double *hstv);

#define fill00_		FILL00
#define resp00_		RESP00
#define stif00_		STIF00
#define get00_		GET00


// I don't know which subroutines to call, so fill in the XX for Bilinear later -- MHS
// Declarations for the Bilinear subroutines
//extern "C" int _stdcall FILLXX(double *data, double *hstv, double *stateP);
//extern "C" int _stdcall RESPXX(int *kresis, int *ksave, int *kgem, int *kstep,
//							   int *ndof, int *kst, int *kenr,
//							   double *ener, double *ened, double *enso, double *beto,
//							   double *relas, double *rdamp, double *rinit,
//							   double *ddise, double *dise, double *vele);
//extern "C" int _stdcall STIFXX(int *kstt, int *ktype, int *ndof, double *fk);
//extern "C" int _stdcall GETXX(double *hstv);

//#define fillXX_		FILLXX
//#define respXX_		RESPXX
//#define stifXX_		STIFXX
//#define getXX_		GETXX


// I don't know which subroutines to call, so fill in the XX for Clough1 later -- MHS
// Declarations for the Clough1 subroutines
//extern "C" int _stdcall FILLXX(double *data, double *hstv, double *stateP);
//extern "C" int _stdcall RESPXX(int *kresis, int *ksave, int *kgem, int *kstep,
//							   int *ndof, int *kst, int *kenr,
//							   double *ener, double *ened, double *enso, double *beto,
//							   double *relas, double *rdamp, double *rinit,
//							   double *ddise, double *dise, double *vele);
//extern "C" int _stdcall STIFXX(int *kstt, int *ktype, int *ndof, double *fk);
//extern "C" int _stdcall GETXX(double *hstv);

//#define fillXX_		FILLXX
//#define respXX_		RESPXX
//#define stifXX_		STIFXX
//#define getXX_		GETXX


// I don't know which subroutines to call, so fill in the XX for Clough2 later -- MHS
// Declarations for the Clough2 subroutines
//extern "C" int _stdcall FILLXX(double *data, double *hstv, double *stateP);
//extern "C" int _stdcall RESPXX(int *kresis, int *ksave, int *kgem, int *kstep,
//							   int *ndof, int *kst, int *kenr,
//							   double *ener, double *ened, double *enso, double *beto,
//							   double *relas, double *rdamp, double *rinit,
//							   double *ddise, double *dise, double *vele);
//extern "C" int _stdcall STIFXX(int *kstt, int *ktype, int *ndof, double *fk);
//extern "C" int _stdcall GETXX(double *hstv);

//#define fillXX_		FILLXX
//#define respXX_		RESPXX
//#define stifXX_		STIFXX
//#define getXX_		GETXX


// I don't know which subroutines to call, so fill in the XX for Pinch1 later -- MHS
// Declarations for the Pinch1 subroutines
//extern "C" int _stdcall FILLXX(double *data, double *hstv, double *stateP);
//extern "C" int _stdcall RESPXX(int *kresis, int *ksave, int *kgem, int *kstep,
//							   int *ndof, int *kst, int *kenr,
//							   double *ener, double *ened, double *enso, double *beto,
//							   double *relas, double *rdamp, double *rinit,
//							   double *ddise, double *dise, double *vele);
//extern "C" int _stdcall STIFXX(int *kstt, int *ktype, int *ndof, double *fk);
//extern "C" int _stdcall GETXX(double *hstv);

//#define fillXX_		FILLXX
//#define respXX_		RESPXX
//#define stifXX_		STIFXX
//#define getXX_		GETXX


// Add more declarations as needed

#else

// Declarations for the Hardening subroutines
extern "C" int fill00_(double *data, double *hstv, double *stateP);
extern "C" int resp00_(int *kresis, int *ksave, int *kgem, int *kstep,
								int *ndof, int *kst, int *kenr,
								double *ener, double *ened, double *enso, double *beto,
								double *relas, double *rdamp, double *rinit,
								double *ddise, double *dise, double *vele);
extern "C" int stif00_(int *kstt, int *ktype, int *ndof, double *fk);
extern "C" int get00_(double *hstv);


// I don't know which subroutines to call, so fill in the XX for Bilinear later -- MHS
// Declarations for the Bilinear subroutines
//extern "C" int fillXX_(double *data, double *hstv, double *stateP);
//extern "C" int respXX_(int *kresis, int *ksave, int *kgem, int *kstep,
//								int *ndof, int *kst, int *kenr,
//								double *ener, double *ened, double *enso, double *beto,
//								double *relas, double *rdamp, double *rinit,
//								double *ddise, double *dise, double *vele);
//extern "C" int stifXX_(int *kstt, int *ktype, int *ndof, double *fk);
//extern "C" int getXX_(double *hstv);


// I don't know which subroutines to call, so fill in the XX for Clough1 later -- MHS
// Declarations for the Clough1 subroutines
//extern "C" int fillXX_(double *data, double *hstv, double *stateP);
//extern "C" int respXX_(int *kresis, int *ksave, int *kgem, int *kstep,
//								int *ndof, int *kst, int *kenr,
//								double *ener, double *ened, double *enso, double *beto,
//								double *relas, double *rdamp, double *rinit,
//								double *ddise, double *dise, double *vele);
//extern "C" int stifXX_(int *kstt, int *ktype, int *ndof, double *fk);
//extern "C" int getXX_(double *hstv);


// I don't know which subroutines to call, so fill in the XX for Clough2 later -- MHS
// Declarations for the Clough2 subroutines
//extern "C" int fillXX_(double *data, double *hstv, double *stateP);
//extern "C" int respXX_(int *kresis, int *ksave, int *kgem, int *kstep,
//								int *ndof, int *kst, int *kenr,
//								double *ener, double *ened, double *enso, double *beto,
//								double *relas, double *rdamp, double *rinit,
//								double *ddise, double *dise, double *vele);
//extern "C" int stifXX_(int *kstt, int *ktype, int *ndof, double *fk);
//extern "C" int getXX_(double *hstv);


// I don't know which subroutines to call, so fill in the XX for Pinch1 later -- MHS
// Declarations for the Pinch1 subroutines
//extern "C" int fillXX_(double *data, double *hstv, double *stateP);
//extern "C" int respXX_(int *kresis, int *ksave, int *kgem, int *kstep,
//								int *ndof, int *kst, int *kenr,
//								double *ener, double *ened, double *enso, double *beto,
//								double *relas, double *rdamp, double *rinit,
//								double *ddise, double *dise, double *vele);
//extern "C" int stifXX_(int *kstt, int *ktype, int *ndof, double *fk);
//extern "C" int getXX_(double *hstv);


// Add more declarations as needed

#endif

int
DrainMaterial::invokeSubroutine(void)
{
	// Number of degrees of freedom
	static const int NDOF = 2;

	// Flags sent into RESPXX subroutine
	int kresis = 2;		// Compute static and damping forces
	int ksave  = 0;		// Do not save results
	int kgem   = 0;		// Geometric effects (not used)
	int kstep  = 1;		// Step number (set by subroutine)
	int ndof   = NDOF;	// Number of degrees of freedom
	int kst    = 1;		// Stiffness formation code
	int kenr   = 2;		// Calculate static and dynamic energy

	// Energy terms computed in RESPXX subroutine
	double ener = 0.0;	// Change in elasto-plastic energy
	double ened = 0.0;	// Change in damping energy
	double enso = 0.0;	// Change in second-order energy (not used)

	// Force terms computed in RESPXX subroutine
	static double relas[NDOF];	// Resisting force vector
	static double rdamp[NDOF];	// Damping force vector
	static double rinit[NDOF];	// Initial force vector (not used)

	// Total displacement vector
	static double dise[NDOF];
	dise[0] = 0.0;
	dise[1] = epsilon;

	// Incremental displacement vector
	static double ddise[NDOF];
	ddise[0] = 0.0;
	ddise[1] = epsilon-epsilonP;

	// Velocity vector
	static double vele[NDOF];
	vele[0] = 0.0;
	vele[1] = epsilonDot;

	// Fill in committed state array
	static double stateP[3];
	stateP[0] = epsilonP;
	stateP[1] = sigmaP;
	stateP[2] = tangentP;

	// Flags sent into STIFXX subroutine
	int kstt  = 1;			// Total stiffness
	int ktype = 1;			// Elastic stiffness only

	// Stiffness computed in STIFXX subroutine
	static double fk[NDOF*NDOF];

	switch (this->getClassTag()) {
	case MAT_TAG_DrainHardening:
		// Fill the common block with parameters and history variables
		fill00_(data, hstv, stateP);

		// Call the response subroutine
		resp00_(&kresis, &ksave, &kgem, &kstep, &ndof, &kst, &kenr,
			&ener, &ened, &enso, &beto, relas, rdamp, rinit, ddise, dise, vele);
		
		// Call the stiffness subroutine
		stif00_(&kstt, &ktype, &ndof, fk);
		
		// Get the trial history variables
		get00_(&hstv[numHstv]);
		break;

	case MAT_TAG_DrainBilinear:
		// I don't know which subroutines to call, so fill in the XX for Bilinear later -- MHS
		opserr << "DrainMaterial::invokeSubroutine -- Bilinear subroutine not yet linked\n"; exit(-1);


		//fillXX_(data, hstv, stateP);
		//respXX_(&kresis, &ksave, &kgem, &kstep, &ndof, &kst, &kenr,
		//	&ener, &ened, &enso, &beto, relas, rdamp, rinit, ddise, dise, vele);
		//stifXX_(&kstt, &ktype, &ndof, fk);
		//getXX_(&hstv[numHstv]);
		break;

	case MAT_TAG_DrainClough1:
		// I don't know which subroutines to call, so fill in the XX for Clough1 later -- MHS
		opserr << "DrainMaterial::invokeSubroutine -- Clough1 subroutine not yet linked\n"; exit(-1);

		//fillXX_(data, hstv, stateP);
		//respXX_(&kresis, &ksave, &kgem, &kstep, &ndof, &kst, &kenr,
		//	&ener, &ened, &enso, &beto, relas, rdamp, rinit, ddise, dise, vele);
		//stifXX_(&kstt, &ktype, &ndof, fk);
		//getXX_(&hstv[numHstv]);
		break;

	case MAT_TAG_DrainClough2:
		// I don't know which subroutines to call, so fill in the XX for Clough2 later -- MHS
		opserr << "DrainMaterial::invokeSubroutine -- Clough2 subroutine not yet linked\n"; exit(-1);
			  
		//fillXX_(data, hstv, stateP);
		//respXX_(&kresis, &ksave, &kgem, &kstep, &ndof, &kst, &kenr,
		//	&ener, &ened, &enso, &beto, relas, rdamp, rinit, ddise, dise, vele);
		//stifXX_(&kstt, &ktype, &ndof, fk);
		//getXX_(&hstv[numHstv]);
		break;

	case MAT_TAG_DrainPinch1:
		// I don't know which subroutines to call, so fill in the XX for Pinch1 later -- MHS
		opserr << "DrainMaterial::invokeSubroutine -- Pinch1 subroutine not yet linked\n"; exit(-1);
		
		//fillXX_(data, hstv, stateP);
		//respXX_(&kresis, &ksave, &kgem, &kstep, &ndof, &kst, &kenr,
		//	&ener, &ened, &enso, &beto, relas, rdamp, rinit, ddise, dise, vele);
		//stifXX_(&kstt, &ktype, &ndof, fk);
		//getXX_(&hstv[numHstv]);
		break;

	// Add more cases as needed

	default:
		opserr << "DrainMaterial::invokeSubroutine -- unknown material type\n"; exit(-1);
		return -1;
	}

	// Total stress is elastic plus damping force
	sigma = relas[1] + rdamp[1];

	// Get tangent stiffness
	tangent = fk[0];

	return 0;
}
