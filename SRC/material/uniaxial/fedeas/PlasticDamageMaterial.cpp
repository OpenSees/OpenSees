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
                                                                        
// $Revision: 1.6 $
// $Date: 2007-02-23 01:04:40 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/fedeas/PlasticDamageMaterial.cpp,v $
                                                                      
// Written: Jeeho Lee
// Created: Feb 2007
//
// Description: This file contains the class/methods definition for 
// PlasticDamageMaterial. PlasticDamageMaterial wraps the FEDEAS
// 1d Material subroutine: PD_1

#include <stdlib.h>
#include <PlasticDamageMaterial.h>

PlasticDamageMaterial::PlasticDamageMaterial(int tag, double E, double Ft, double Fc,
					     double ft_max, double fcy, double fc_max, double kt_crit, double Relax):
  // 11 history variables and 8 material parameters
    FedeasMaterial(tag, MAT_TAG_PlasticDamage, 11, 8)
{
  // Fill in material parameters
  data[0] = E;
  data[1] = Ft;
  data[2] = Fc;
  data[3] = ft_max;
  data[4] = fcy;
  data[5] = fc_max;
  data[6] = kt_crit;
  data[7] = Relax;
  
  tangentP = E;
  tangent = tangentP;
}

PlasticDamageMaterial::PlasticDamageMaterial(int tag, const Vector &d):
  // 11 history variables and 8 material parameters
  FedeasMaterial(tag, MAT_TAG_PlasticDamage, 11, 8)
{
  if (d.Size() != numData) {
    opserr << "PlasticDamageMaterial::PlasticDamageMaterial -- not enough input arguments\n";
    exit(-1);	
  }
  
  for (int i = 0; i < numData; i++)
    data[i] = d(i);
  
  tangentP = data[0];
  tangent = tangentP;
}

PlasticDamageMaterial::PlasticDamageMaterial(void):
FedeasMaterial(0, MAT_TAG_PlasticDamage, 11, 8)
{
  // Does nothing
}

PlasticDamageMaterial::~PlasticDamageMaterial(void)
{
  // Does nothing
}

UniaxialMaterial*
PlasticDamageMaterial::getCopy(void)
{
  Vector d(data, numData);
  
  PlasticDamageMaterial *theCopy = new PlasticDamageMaterial(this->getTag(), d);
  
  // Copy history variables
  for (int i = 0; i < 2*numHstv; i++)
    theCopy->hstv[i] = hstv[i];
  
  theCopy->epsilonP = epsilonP;
  theCopy->sigmaP   = sigmaP;
  theCopy->tangentP = tangentP;
  
  theCopy->epsilon = epsilonP;
  theCopy->sigma = sigmaP;
  theCopy->tangent = tangentP;  
  
  return theCopy;
}   

double
PlasticDamageMaterial::getInitialTangent(void)
{
  //return E;
  return data[0];
}


#ifdef _WIN32

extern "C" int pd_(double *matpar, double *hstvP, double *hstv,
		  double *strainP, double *stressP, double *dStrain,
		  double *tangent, double *stress, int *ist);

// Add more declarations as needed

//#define pd_	PD

#else

extern "C" int pd_(double *matpar, double *hstvP, double *hstv,
		   double *strainP, double *stressP, double *dStrain,
		   double *tangent, double *stress, int *ist);

// Add more declarations as needed

#endif


int
PlasticDamageMaterial::invokeSubroutine(int ist)
{
  // Compute strain increment
  double dEpsilon = epsilon-epsilonP;
  
  pd_(data, hstv, &hstv[numHstv], &epsilonP, &sigmaP, &dEpsilon, &sigma, &tangent, &ist);
	   
  
  return 0;
}

