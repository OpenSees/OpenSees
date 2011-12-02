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
                                                                        
// $Revision: 1.18 $
// $Date: 2005-04-04 23:53:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/FedeasMaterial.cpp,v $
                                                                        
// Written: MHS
// Created: Jan 2001
//
// Description: This file contains the class definition for 
// FedeasMaterial. FedeasMaterial provides a FORTRAN interface
// for programming uniaxial material models, using the subroutine
// interface from the FEDEAS ML1D library, developed by F.C. Filippou.
//
// For more information visit the FEDEAS web page:
//    http://www.ce.berkeley.edu/~filippou/Research/fedeas.htm

#include <OPS_Globals.h>
#include <FedeasMaterial.h>
#include <Vector.h>
#include <ID.h>
#include <Channel.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

FedeasMaterial::FedeasMaterial(int tag, int classTag, int nhv, int ndata)
  :UniaxialMaterial(tag,classTag),
   data(0), hstv(0), numData(ndata), numHstv(nhv),
   epsilonP(0.0), sigmaP(0.0), tangentP(0),
   epsilon(0.0), sigma(0.0), tangent(0)
{
  if (numHstv < 0)
    numHstv = 0;
  
  if (numHstv > 0) {
    // Allocate history array
    hstv = new double[2*numHstv];
    if (hstv == 0) {
      opserr << "FedeasMaterial::FedeasMaterial -- failed to allocate history array -- type " << 
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
      opserr << "FedeasMaterial::FedeasMaterial -- failed to allocate data array -- type : " <<
	this->getClassTag() << endln;
      exit(-1);
    }
			    
    // Initialize to zero
    for (int i = 0; i < numData; i++)
      data[i] = 0.0;
  }
}

FedeasMaterial::~FedeasMaterial()
{
  if (hstv != 0)
    delete [] hstv;
  
  if (data != 0)
    delete [] data;
}

int
FedeasMaterial::setTrialStrain(double strain, double strainRate)
{

  if (fabs(strain-epsilon) > DBL_EPSILON) {
    // Store the strain
    epsilon = strain;
    
    // Tells subroutine to do normal operations for stress and tangent
    int ist = 1;
    
    // Call the subroutine
    return this->invokeSubroutine(ist);
  }
  return 0;
}

int
FedeasMaterial::setTrial(double strain, double &stress, double &stiff, double strainRate)
{
  int res = 0;
  if (fabs(strain-epsilon) > DBL_EPSILON) {

    // Store the strain
    epsilon = strain;
    
    // Tells subroutine to do normal operations for stress and tangent
    int ist = 1;
    
    // Call the subroutine
    res = this->invokeSubroutine(ist);
  }
    
  stress = sigma;
  stiff = tangent;
  
  return res;
}

double
FedeasMaterial::getStrain(void)
{
  return epsilon;
}

double
FedeasMaterial::getStress(void)
{
  return sigma;
}

double
FedeasMaterial::getTangent(void)
{
  return tangent;
}

int
FedeasMaterial::commitState(void)
{
  // Set committed values equal to corresponding trial values
  for (int i = 0; i < numHstv; i++)
    hstv[i] = hstv[i+numHstv];
  
  epsilonP = epsilon;
  sigmaP = sigma;
  tangentP = tangent;

  return 0;
}

int
FedeasMaterial::revertToLastCommit(void)
{
  // Set trial values equal to corresponding committed values
  for (int i = 0; i < numHstv; i++)
    hstv[i+numHstv] = hstv[i];
  
  epsilon = epsilonP;
  sigma = sigmaP;
  tangent = tangentP;

  return 0;
}

int
FedeasMaterial::revertToStart(void)
{
  // Set all trial and committed values to zero
  for (int i = 0; i < 2*numHstv; i++)
    hstv[i] = 0.0;
  
  epsilonP = 0.0;
  sigmaP = 0.0;
  tangentP = this->getInitialTangent();

  epsilon = 0.0;
  sigma = 0.0;
  tangent = tangentP;

  return 0;
}

int 
FedeasMaterial::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  Vector vecData(numHstv+numData+4);
  
  int i, j;
  // Copy only the committed history variables into vector
  for (i = 0; i < numHstv; i++)
    vecData(i) = hstv[i];
  
  // Copy material properties into vector
  for (i = 0, j = numHstv; i < numData; i++, j++)
    vecData(j) = data[i];
  
  vecData(j++) = epsilonP;
  vecData(j++) = sigmaP;
  vecData(j++) = tangentP;
  vecData(j++) = this->getTag();
  
  res += theChannel.sendVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) 
    opserr << "FedeasMaterial::sendSelf() - failed to send Vector data\n";
  
  return res;
}

int
FedeasMaterial::recvSelf(int commitTag, Channel &theChannel,
			 FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  Vector vecData(numHstv+numData+4);
  
  res += theChannel.recvVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) {
    opserr << "FedeasMaterial::recvSelf() - failed to receive Vector data\n";
    return res;
  }
  
  int i, j;
  // Copy committed history variables from vector
  for (i = 0; i < numHstv; i++)
    hstv[i] = vecData(i);
  
  // Copy material properties from vector
  for (i = 0, j = numHstv; i < numData; i++, j++)
    data[i] = vecData(j);
  
  epsilonP = vecData(j++);
  sigmaP   = vecData(j++);
  tangentP = vecData(j++);
  this->setTag((int)vecData(j++));

  tangent = tangentP;
  sigma = sigmaP;
  epsilon = epsilonP;

  return res;
}

void
FedeasMaterial::Print(OPS_Stream &s, int flag)
{
  s << "FedeasMaterial, type: ";
	
  switch (this->getClassTag()) {
  case MAT_TAG_FedeasHardening:
    s << "Hardening" << endln;
    break;
  case MAT_TAG_FedeasBond1:
    s << "Bond1" << endln;
    break;
  case MAT_TAG_FedeasBond2:
    s << "Bond2" << endln;
    break;
  case MAT_TAG_FedeasConcrete1:
    s << "Concrete1" << endln;
    break;
  case MAT_TAG_FedeasConcrete2:
    s << "Concrete2" << endln;
    break;
  case MAT_TAG_FedeasConcrete3:
    s << "Concrete3" << endln;
    break;
  case MAT_TAG_FedeasHysteretic1:
    s << "Hysteretic1" << endln;
    break;
  case MAT_TAG_FedeasHysteretic2:
    s << "Hysteretic2" << endln;
    break;
  case MAT_TAG_FedeasSteel1:
    s << "Steel1" << endln;
    break;
  case MAT_TAG_FedeasSteel2:
    s << "Steel2" << endln;
    break;
    // Add more cases as needed
    
  default:
    s << "Material identifier = " << this->getClassTag() << endln;
    break;
  }
}


#ifdef _pgCC
#define hard_1__ hard1_1_
#endif



#ifdef _WIN32

extern "C" int BOND_1(double *matpar, double *hstvP, double *hstv,
			       double *strainP, double *stressP, double *dStrain,
			       double *tangent, double *stress, int *ist);

extern "C" int BOND_2(double *matpar, double *hstvP, double *hstv,
			       double *strainP, double *stressP, double *dStrain,
			       double *tangent, double *stress, int *ist);

extern "C" int CONCRETE_1(double *matpar, double *hstvP, double *hstv,
				   double *strainP, double *stressP, double *dStrain,
				   double *tangent, double *stress, int *ist);

extern "C" int CONCRETE_2(double *matpar, double *hstvP, double *hstv,
				   double *strainP, double *stressP, double *dStrain,
				   double *tangent, double *stress, int *ist);

extern "C" int CONCRETE_3(double *matpar, double *hstvP, double *hstv,
				   double *strainP, double *stressP, double *dStrain,
				   double *tangent, double *stress, int *ist);

extern "C" int HARD_1(double *matpar, double *hstvP, double *hstv,
			       double *strainP, double *stressP, double *dStrain,
			       double *tangent, double *stress, int *ist);

extern "C" int HYSTER_1(double *matpar, double *hstvP, double *hstv,
				 double *strainP, double *stressP, double *dStrain,
				 double *tangent, double *stress, int *ist);

extern "C" int HYSTER_2(double *matpar, double *hstvP, double *hstv,
				 double *strainP, double *stressP, double *dStrain,
				 double *tangent, double *stress, int *ist);

extern "C" int STEEL_1(double *matpar, double *hstvP, double *hstv,
				double *strainP, double *stressP, double *dStrain,
				double *tangent, double *stress, int *ist);

extern "C" int  STEEL_2(double *matpar, double *hstvP, double *hstv,
				double *strainP, double *stressP, double *dStrain,
				double *tangent, double *stress, int *ist);

// Add more declarations as needed

#define bond_1__	BOND_1
#define bond_2__	BOND_2
#define concrete_1__	CONCRETE_1
#define concrete_2__	CONCRETE_2
#define concrete_3__	CONCRETE_3
#define hard_1__ 	HARD_1
#define hyster_1__	HYSTER_1
#define hyster_2__	HYSTER_2
#define steel_1__	STEEL_1
#define steel_2__	STEEL_2

#else

extern "C" int bond_1__(double *matpar, double *hstvP, double *hstv,
			double *strainP, double *stressP, double *dStrain,
			double *tangent, double *stress, int *ist);

extern "C" int bond_2__(double *matpar, double *hstvP, double *hstv,
			double *strainP, double *stressP, double *dStrain,
			double *tangent, double *stress, int *ist);

extern "C" int concrete_1__(double *matpar, double *hstvP, double *hstv,
			    double *strainP, double *stressP, double *dStrain,
			    double *tangent, double *stress, int *ist);

extern "C" int concrete_2__(double *matpar, double *hstvP, double *hstv,
			    double *strainP, double *stressP, double *dStrain,
			    double *tangent, double *stress, int *ist);

extern "C" int concrete_3__(double *matpar, double *hstvP, double *hstv,
			    double *strainP, double *stressP, double *dStrain,
			    double *tangent, double *stress, int *ist);

extern "C" int hard_1__(double *matpar, double *hstvP, double *hstv,
			double *strainP, double *stressP, double *dStrain,
			double *tangent, double *stress, int *ist);

extern "C" int hyster_1__(double *matpar, double *hstvP, double *hstv,
			  double *strainP, double *stressP, double *dStrain,
			  double *tangent, double *stress, int *ist);

extern "C" int hyster_2__(double *matpar, double *hstvP, double *hstv,
			  double *strainP, double *stressP, double *dStrain,
			  double *tangent, double *stress, int *ist);

extern "C" int steel_1__(double *matpar, double *hstvP, double *hstv,
			 double *strainP, double *stressP, double *dStrain,
			 double *tangent, double *stress, int *ist);

extern "C" int steel_2__(double *matpar, double *hstvP, double *hstv,
			 double *strainP, double *stressP, double *dStrain,
			 double *tangent, double *stress, int *ist);

// Add more declarations as needed

#endif


int
FedeasMaterial::invokeSubroutine(int ist)
{
  // Compute strain increment
  double dEpsilon = epsilon-epsilonP;
  
  switch (this->getClassTag()) {
  case MAT_TAG_FedeasHardening:
    hard_1__(data, hstv, &hstv[numHstv], &epsilonP, &sigmaP, &dEpsilon, 
	     &sigma, &tangent, &ist);
    break;

  case MAT_TAG_FedeasBond1:
#ifdef _WIN32
    bond_1__(data, hstv, &hstv[numHstv], &epsilonP, &sigmaP, &dEpsilon,
	     &sigma, &tangent, &ist);
#else
	opserr << "FedeasMaterial::invokeSubroutine -- Bond1 subroutine not yet linked\n";
#endif
    break;
    
  case MAT_TAG_FedeasBond2:
#ifdef _WIN32
    bond_2__(data, hstv, &hstv[numHstv], &epsilonP, &sigmaP, &dEpsilon,
	     &sigma, &tangent, &ist);
#else
	opserr << "FedeasMaterial::invokeSubroutine -- Bond2 subroutine not yet linked\n";
#endif
    break;
    
  case MAT_TAG_FedeasConcrete1:
#ifdef _WIN32
    concrete_1__(data, hstv, &hstv[numHstv], &epsilonP, &sigmaP, &dEpsilon, 
		 &sigma, &tangent, &ist);
#else
	opserr << "FedeasMaterial::invokeSubroutine -- Concrete1 subroutine not yet linked\n";
#endif
    break;
    
  case MAT_TAG_FedeasConcrete2:
#ifdef _WIN32
    concrete_2__(data, hstv, &hstv[numHstv], &epsilonP, &sigmaP, &dEpsilon, 
		 &sigma, &tangent, &ist);
#else
	opserr << "FedeasMaterial::invokeSubroutine -- Concrete2 subroutine not yet linked\n";
#endif
    break;
    
  case MAT_TAG_FedeasConcrete3:
#ifdef _WIN32
    concrete_3__(data, hstv, &hstv[numHstv], &epsilonP, &sigmaP, &dEpsilon, 
		 &sigma, &tangent, &ist);
#else
	opserr << "FedeasMaterial::invokeSubroutine -- Concrete3 subroutine not yet linked\n";
#endif
    break;
        
  case MAT_TAG_FedeasHysteretic1:
#ifdef _WIN32
    hyster_1__(data, hstv, &hstv[numHstv], &epsilonP, &sigmaP, &dEpsilon, 
	       &sigma, &tangent, &ist);
#else
	opserr << "FedeasMaterial::invokeSubroutine -- Hysteretic1 subroutine not yet linked\n";
#endif
    break;
    
  case MAT_TAG_FedeasHysteretic2:
#ifdef _WIN32
    hyster_2__(data, hstv, &hstv[numHstv], &epsilonP, &sigmaP, &dEpsilon, 
	       &sigma, &tangent, &ist);
#else
	opserr << "FedeasMaterial::invokeSubroutine -- Hysteretic2 subroutine not yet linked\n";
#endif
    break;
    
  case MAT_TAG_FedeasSteel1:
#ifdef _WIN32
    steel_1__(data, hstv, &hstv[numHstv], &epsilonP, &sigmaP, &dEpsilon, 
	      &sigma, &tangent, &ist);
#else
	opserr << "FedeasMaterial::invokeSubroutine -- Steel1 subroutine not yet linked\n";
#endif
    break;
    
  case MAT_TAG_FedeasSteel2:
#ifdef _WIN32
    steel_2__(data, hstv, &hstv[numHstv], &epsilonP, &sigmaP, &dEpsilon, 
	      &sigma, &tangent, &ist);
#else
	opserr << "FedeasMaterial::invokeSubroutine -- Steel2 subroutine not yet linked\n";
#endif
    break;
    
    // Add more cases as needed
  default:
    opserr << "FedeasMaterial::invokeSubroutine -- unknown material type\n";
    return -1;
  }
  
  return 0;
}

