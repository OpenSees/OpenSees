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
                                                                        
// $Revision: 1.19 $                                                              
// $Date: 2004-02-24 20:50:58 $                                                                  
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropicMaterial.cpp,v $                                                                
                                                                        
                                                                        
// File: ~/material/ElasticIsotropicMaterial.C
//
// Written: MHS 
// Created: Feb 2000
// Revision: A
// Boris Jeremic (@ucdavis.edu) 19June2002 added getE, getnu
//
// Description: This file contains the class implementation for ElasticIsotropicMaterial.
//
// What: "@(#) ElasticIsotropicMaterial.C, revA"

#include <string.h>

#include <ElasticIsotropicMaterial.h>
#include <ElasticIsotropicPlaneStress2D.h>
#include <ElasticIsotropicPlaneStrain2D.h>
#include <ElasticIsotropicAxiSymm.h>
#include <ElasticIsotropic3D.h>
#include <PressureDependentElastic3D.h>
#include <ElasticIsotropicPlateFiber.h>
#include <ElasticIsotropicBeamFiber.h>

#include <Tensor.h>
#include <Channel.h>

#include <OPS_Globals.h>

ElasticIsotropicMaterial::ElasticIsotropicMaterial
(int tag, int classTag, double e, double nu, double r)
  :NDMaterial(tag, classTag), E(e), v(nu), rho(r)
{

}

ElasticIsotropicMaterial::ElasticIsotropicMaterial
(int tag, double e, double nu, double r)
  :NDMaterial(tag, ND_TAG_ElasticIsotropic), E(e), v(nu), rho(r)
{

}

ElasticIsotropicMaterial::~ElasticIsotropicMaterial()
{
	
}

double
ElasticIsotropicMaterial::getRho() 
{ 
  return rho ;
}


// Boris Jeremic (@ucdavis.edu) 19June2002
double ElasticIsotropicMaterial::getE() 
{ 
  return E;
}

// Boris Jeremic (@ucdavis.edu) 19June2002
double ElasticIsotropicMaterial::getnu() 
{ 
  return v;
}



NDMaterial*
ElasticIsotropicMaterial::getCopy (const char *type)
{
    if (strcmp(type,"PlaneStress2D") == 0 || strcmp(type,"PlaneStress") == 0)
    {
	ElasticIsotropicPlaneStress2D *theModel;
	theModel = new ElasticIsotropicPlaneStress2D (this->getTag(), E, v, rho);
		// DOES NOT COPY sigma, D, and epsilon ...
		// This function should only be called during element instantiation, so
		// no state determination is performed on the material model object
		// prior to copying the material model (calling this function)
	return theModel;
    }

    else if (strcmp(type,"PlaneStrain2D") == 0 || strcmp(type,"PlaneStrain") == 0)
    {
	ElasticIsotropicPlaneStrain2D *theModel;
	theModel = new ElasticIsotropicPlaneStrain2D (this->getTag(), E, v, rho);
		// DOES NOT COPY sigma, D, and epsilon ...
		// This function should only be called during element instantiation, so
		// no state determination is performed on the material model object
		// prior to copying the material model (calling this function)
	return theModel;
    }
    else if (strcmp(type,"AxiSymmetric2D") == 0 || strcmp(type,"AxiSymmetric") == 0)
    {
	ElasticIsotropicAxiSymm *theModel;
	theModel = new ElasticIsotropicAxiSymm(this->getTag(), E, v, rho);
		// DOES NOT COPY sigma, D, and epsilon ...
		// This function should only be called during element instantiation, so
		// no state determination is performed on the material model object
		// prior to copying the material model (calling this function)
	return theModel;
    }
///////////////////////////////
    else if (strcmp(type,"ThreeDimensional") == 0 || 
	     strcmp(type,"3D") == 0)
    {
	ElasticIsotropic3D *theModel;
	theModel = new ElasticIsotropic3D (this->getTag(), E, v, rho);
		// DOES NOT COPY sigma, D, and epsilon ...
		// This function should only be called during element instantiation, so
		// no state determination is performed on the material model object
		// prior to copying the material model (calling this function)
	return theModel;
    }
///////////////////////////////
    else if (strcmp(type,"PlateFiber") == 0)
    {
	ElasticIsotropicPlateFiber *theModel;
	theModel = new ElasticIsotropicPlateFiber(this->getTag(), E, v, rho);
		// DOES NOT COPY sigma, D, and epsilon ...
		// This function should only be called during element instantiation, so
		// no state determination is performed on the material model object
		// prior to copying the material model (calling this function)
	return theModel;
    }
    else if (strcmp(type,"BeamFiber") == 0)
    {
	ElasticIsotropicBeamFiber *theModel;
	theModel = new ElasticIsotropicBeamFiber(this->getTag(), E, v, rho);
		// DOES NOT COPY sigma, D, and epsilon ...
		// This function should only be called during element instantiation, so
		// no state determination is performed on the material model object
		// prior to copying the material model (calling this function)
	return theModel;
    }


    // Handle other cases
    else
    {
      opserr << "ElasticIsotropicMaterial::getModel failed to get model: " << type << endln;
      exit(-1);
    }
    
    return 0;
}

int
ElasticIsotropicMaterial::setTrialStrain (const Vector &v)
{
    opserr << "ElasticIsotropicMaterial::setTrialStrain -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
ElasticIsotropicMaterial::setTrialStrain (const Vector &v, const Vector &rate)
{
    opserr << "ElasticIsotropicMaterial::setTrialStrain -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
ElasticIsotropicMaterial::setTrialStrainIncr (const Vector &v)
{
    opserr << "ElasticIsotropicMaterial::setTrialStrainIncr -- subclass responsibility\n";
    exit(-1);
    return -1;
}

int
ElasticIsotropicMaterial::setTrialStrainIncr (const Vector &v, const Vector &rate)
{
    opserr << "ElasticIsotropicMaterial::setTrialStrainIncr -- subclass responsibility\n";
    exit(-1);
    return -1;
}

const Matrix&
ElasticIsotropicMaterial::getTangent (void)
{
  opserr << "ElasticIsotropicMaterial::getTangent -- subclass responsibility\n";
  exit(-1);

  // Just to make it compile
  Matrix *ret = new Matrix();
  return *ret;
}

const Matrix&
ElasticIsotropicMaterial::getInitialTangent (void)
{
  opserr << "ElasticIsotropicMaterial::getInitialTangent -- subclass responsibility\n";
  exit(-1);

  // Just to make it compile
  Matrix *ret = new Matrix();
  return *ret;
}

const Vector&
ElasticIsotropicMaterial::getStress (void)
{
  opserr << "ElasticIsotropicMaterial::getStress -- subclass responsibility\n";
  exit(-1);
    
  // Just to make it compile
  Vector *ret = new Vector();
  return *ret;
}

const Vector&
ElasticIsotropicMaterial::getStrain (void)
{
  opserr << "ElasticIsotropicMaterial::getStrain -- subclass responsibility\n";
  exit(-1);

  // Just to make it compile
  Vector *ret = new Vector();
  return *ret;
}

int
ElasticIsotropicMaterial::setTrialStrain (const Tensor &v)
{
    opserr << "ElasticIsotropicMaterial::setTrialStrain -- subclass responsibility\n";
    exit(-1);

    return -1;
}

int
ElasticIsotropicMaterial::setTrialStrain (const Tensor &v, const Tensor &r)
{
    opserr << "ElasticIsotropicMaterial::setTrialStrain -- subclass responsibility\n";
    exit(-1);

    return -1;
}

int
ElasticIsotropicMaterial::setTrialStrainIncr (const Tensor &v)
{
    opserr << "ElasticIsotropicMaterial::setTrialStrainIncr -- subclass responsibility\n";
    exit(-1);

    return -1;
}

int
ElasticIsotropicMaterial::setTrialStrainIncr (const Tensor &v, const Tensor &r)
{
    opserr << "ElasticIsotropicMaterial::setTrialStrainIncr -- subclass responsibility\n";

    return -1;
}

const Tensor&
ElasticIsotropicMaterial::getTangentTensor (void)
{
  opserr << "ElasticIsotropicMaterial::getTangentTensor -- subclass responsibility\n";
  exit(-1);
  
  // Just to make it compile
  Tensor *t = new Tensor;
  return *t;
}

const stresstensor ElasticIsotropicMaterial::getStressTensor (void)
{
  opserr << "ElasticIsotropicMaterial::getStressTensor -- subclass responsibility\n";
  exit(-1);

  // Just to make it compile
  stresstensor t;
  return t;
}

const straintensor ElasticIsotropicMaterial::getStrainTensor (void)
{
  opserr << "ElasticIsotropicMaterial::getStrainTensor -- subclass responsibility\n";
  exit(-1);

  // Just to make it compile
  straintensor t;
  return t;
}

const straintensor ElasticIsotropicMaterial::getPlasticStrainTensor (void)
{
  opserr << "ElasticIsotropicMaterial::getPlasticStrainTensor -- subclass responsibility\n";
  exit(-1);
	
  // Just to make it compile
  straintensor t;
  return t;
}

int
ElasticIsotropicMaterial::commitState (void)
{
  opserr << "ElasticIsotropicMaterial::commitState -- subclass responsibility\n";
  exit(-1);
  return -1;
}

int
ElasticIsotropicMaterial::revertToLastCommit (void)
{
  opserr << "ElasticIsotropicMaterial::revertToLastCommit -- subclass responsibility\n";
  exit(-1);
    
  return -1;
}

int
ElasticIsotropicMaterial::revertToStart (void)
{
  opserr << "ElasticIsotropicMaterial::revertToStart -- subclass responsibility\n";
  exit(-1);
  return -1;
}

NDMaterial*
ElasticIsotropicMaterial::getCopy (void)
{
  opserr << "ElasticIsotropicMaterial::getCopy -- subclass responsibility\n";
  exit(-1);
  return 0;
}

const char*
ElasticIsotropicMaterial::getType (void) const
{
  opserr << "ElasticIsotropicMaterial::getType -- subclass responsibility\n";
  exit(-1);	

  return 0;
}

int
ElasticIsotropicMaterial::getOrder (void) const
{
  opserr << "ElasticIsotropicMaterial::getOrder -- subclass responsibility\n";
  exit(-1);
  return -1;
}

int
ElasticIsotropicMaterial::sendSelf (int commitTag, Channel &theChannel)
{
  int res = 0;

  static Vector data(4);
  
  data(0) = this->getTag();
  data(1) = E;
  data(2) = v;
  data(3) = rho;
  
 res += theChannel.sendVector(this->getDbTag(), commitTag, data);
 if (res < 0) {
   opserr << "ElasticIsotropicMaterial::sendSelf -- could not send Vector\n";
   return res;
 }

 return res;
}

int
ElasticIsotropicMaterial::recvSelf (int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(4);
  
  res += theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
   opserr << "ElasticIsotropicMaterial::recvSelf -- could not recv Vector\n";
   return res;
  }
    
  this->setTag((int)data(0));
  E = data(1);
  v = data(2);
  rho = data(3);
  
  return res;
}

void
ElasticIsotropicMaterial::Print (OPS_Stream &s, int flag)
{
	s << "Elastic Isotropic Material Model" << endln;
	s << "\tE:  " << E << endln;
	s << "\tv:  " << v << endln;
	s << "\trho:  " << rho << endln;

	return;
}

int 
ElasticIsotropicMaterial::setParameter(char **argv, int argc, Information &info)
{
  return -1;
}

int 
ElasticIsotropicMaterial::updateParameter(int parameterID, Information &info)
{ 
  return -1;
}
